#!/usr/bin/env nextflow

/*
 * SLE Methylation Pipeline — DSL2
 * FASTQ → alignment → methylation calling → batch correction →
 * cell deconvolution → DMR detection → patient stratification
 */

nextflow.enable.dsl = 2

// --- Module imports ---
include { FETCH_SRA }        from './modules/fetch_sra'
include { FASTP }             from './modules/fastp'
include { BWAMETH_INDEX }    from './modules/bwameth_index'
include { BWAMETH_ALIGN }    from './modules/bwameth_align'
include { MARK_DUPLICATES }  from './modules/mark_duplicates'
include { BAM_TO_CRAM }     from './modules/bam_to_cram'
include { METHYLDACKEL }     from './modules/methyldackel'
include { COMBAT_METH }      from './modules/combat_meth'
include { PCA_PLOT }          from './modules/pca_plot'
include { HOUSEMAN_DECONV }  from './modules/houseman_deconv'
include { REGION_DETECT }    from './modules/region_detect'
include { NMF_STRATIFY }     from './modules/nmf_stratify'

// --- Help message ---
if (params.help) {
    log.info """
    ===========================================
     SLE Methylation Pipeline v${workflow.manifest.version}
    ===========================================
    Usage:
      nextflow run main.nf -profile conda \\
        --sample_sheet samples.csv \\
        --genome /path/to/hg38.fa \\
        --outdir results

    Required:
      --sample_sheet   CSV with columns: sample_id, fastq_1, fastq_2, condition, batch
                       For SRA downloads, set fastq_1 to an SRR accession and leave fastq_2 empty
      --genome         Path to reference genome FASTA

    Options:
      --outdir         Output directory [default: ${params.outdir}]
      --min_depth      MethylDackel minimum depth [default: ${params.min_depth}]
      --nmf_kmin       NMF minimum rank [default: ${params.nmf_kmin}]
      --nmf_kmax       NMF maximum rank [default: ${params.nmf_kmax}]
      --nmf_nruns      NMF runs per rank [default: ${params.nmf_nruns}]
      --window_size    DMR window size [default: ${params.window_size}]
      --step_size      DMR step size [default: ${params.step_size}]
      --min_cpgs       DMR min CpGs per window [default: ${params.min_cpgs}]
      --alignment_only      Stop after methylation extraction (skip cohort steps) [default: false]
      --clinical_metadata   TSV with numeric clinical variables for NMF correlation [default: none]
    """.stripIndent()
    exit 0
}

// --- Parameter validation ---
if (!params.sample_sheet) {
    error "Please provide --sample_sheet"
}
if (!params.genome) {
    error "Please provide --genome"
}

// Validate sample sheet exists and has required columns
def ss_file = file(params.sample_sheet)
if (!ss_file.exists()) {
    error "Sample sheet not found: ${params.sample_sheet}"
}
def ss_header = ss_file.readLines()[0].split(',') as List
def required_cols = ['sample_id', 'fastq_1', 'fastq_2', 'condition', 'batch']
def missing_cols = required_cols.findAll { !ss_header.contains(it) }
if (missing_cols) {
    error "Sample sheet missing required columns: ${missing_cols.join(', ')}. Found: ${ss_header.join(', ')}"
}

// Validate genome exists (skip for S3/GCS paths)
if (!params.genome.startsWith('s3://') && !params.genome.startsWith('gs://')) {
    def genome_file = file(params.genome)
    if (!genome_file.exists()) {
        error "Genome file not found: ${params.genome}"
    }
}

// --- Main workflow ---
workflow {

    // Sample sheet as file for cohort-level processes
    sample_sheet_ch = Channel.fromPath(params.sample_sheet)

    if (params.skip_alignment) {
        // ---- Skip alignment mode ----
        // Use pre-made bedGraph files (e.g., from test data generator)
        log.info "Skipping alignment — using pre-made bedGraphs from: ${params.bedgraph_dir}"

        if (!params.bedgraph_dir) {
            error "skip_alignment=true requires --bedgraph_dir"
        }

        bedgraphs_ch = Channel
            .fromPath("${params.bedgraph_dir}/*_CpG.bedGraph")
            .collect()

    } else {
        // ---- Full pipeline: FASTQ → alignment → methylation calling ----

        // Parse sample sheet — split into SRA accessions vs local FASTQs
        Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .branch { row ->
                sra:   row.fastq_1 =~ /^SRR/ && (!row.fastq_2 || row.fastq_2 == '')
                local: true
            }
            .set { sample_rows }

        // ---- Step 0: Fetch SRA data (if any) ----
        sra_ch = sample_rows.sra
            .map { row -> tuple(row.sample_id, row.fastq_1) }

        // FETCH_SRA disabled — use local FASTQs only (samples_downsampled.csv has no SRR rows)
        // FETCH_SRA(sra_ch)

        // Local FASTQs
        local_reads_ch = sample_rows.local
            .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)]) }

        reads_ch = local_reads_ch

        // Reference genome
        genome_ch = Channel.fromPath(params.genome)

        // Genome index (.fai) — derive from genome path if not provided
        if (params.genome_fai) {
            fai_ch = Channel.fromPath(params.genome_fai)
        } else {
            fai_ch = Channel.fromPath("${params.genome}.fai")
        }

        // ---- Step 1: QC + Trim (fastp) ----
        FASTP(reads_ch)

        // ---- Step 2: Index genome ----
        BWAMETH_INDEX(genome_ch)

        // ---- Step 3: Align ----
        BWAMETH_ALIGN(
            FASTP.out.trimmed_reads,
            genome_ch.first(),
            BWAMETH_INDEX.out.index.first()
        )

        // ---- Step 4: Mark duplicates ----
        MARK_DUPLICATES(BWAMETH_ALIGN.out.bam)

        // ---- Step 4b: Convert BAM → CRAM (40-60% smaller) ----
        dedup_bam_bai = MARK_DUPLICATES.out.bam
            .join(MARK_DUPLICATES.out.bai)
            .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }

        BAM_TO_CRAM(
            dedup_bam_bai,
            genome_ch.first(),
            fai_ch.first()
        )

        // ---- Step 5: Methylation extraction (from CRAMs) ----
        cram_crai = BAM_TO_CRAM.out.cram
            .join(BAM_TO_CRAM.out.crai)
            .map { sample_id, cram, crai -> tuple(sample_id, cram, crai) }

        METHYLDACKEL(
            cram_crai,
            genome_ch.first(),
            fai_ch.first()
        )

        // Collect all bedGraphs (strip sample_id, just get the files)
        bedgraphs_ch = METHYLDACKEL.out.bedgraph
            .map { sample_id, bg -> bg }
            .collect()
    }

    // ---- Cohort-level analysis (requires multiple samples) ----
    if (!params.alignment_only) {

        // ---- Step 7: ComBat-meth batch correction ----
        COMBAT_METH(
            bedgraphs_ch,
            sample_sheet_ch.first()
        )

        // ---- Step 8: PCA plots ----
        PCA_PLOT(
            COMBAT_METH.out.mvalues_raw,
            COMBAT_METH.out.mvalues_corrected,
            sample_sheet_ch.first()
        )

        // ---- Step 9: Houseman deconvolution ----
        HOUSEMAN_DECONV(
            COMBAT_METH.out.beta_matrix,
            sample_sheet_ch.first()
        )

        // ---- Step 10: Region detection (dmrseq on raw counts) ----
        REGION_DETECT(
            bedgraphs_ch,
            sample_sheet_ch.first()
        )

        // ---- Step 11: NMF patient stratification ----
        // Now receives cell fractions (for regression) and DMRs (for feature selection)
        NMF_STRATIFY(
            COMBAT_METH.out.mvalues_corrected,
            sample_sheet_ch.first(),
            HOUSEMAN_DECONV.out.fractions,
            REGION_DETECT.out.dmrs
        )
    }
}

// --- Completion message ---
workflow.onComplete {
    log.info """
    ===========================================
     Pipeline Complete
    ===========================================
     Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
     Duration : ${workflow.duration}
     Output   : ${params.outdir}
    ===========================================
    """.stripIndent()
}
