#!/usr/bin/env nextflow

/*
 * SLE Methylation Pipeline — DSL2
 * FASTQ → alignment → methylation calling → batch correction →
 * cell deconvolution → DMR detection → patient stratification
 */

nextflow.enable.dsl = 2

// --- Module imports ---
include { FASTQC }           from './modules/fastqc'
include { TRIM_GALORE }      from './modules/trim_galore'
include { BWAMETH_INDEX }    from './modules/bwameth_index'
include { BWAMETH_ALIGN }    from './modules/bwameth_align'
include { MARK_DUPLICATES }  from './modules/mark_duplicates'
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

// --- Main workflow ---
workflow {

    // Parse sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)])
        }
        .set { reads_ch }

    // Sample sheet as file for cohort-level processes
    sample_sheet_ch = Channel.fromPath(params.sample_sheet)

    // Reference genome
    genome_ch = Channel.fromPath(params.genome)

    // Genome index (.fai) — derive from genome path if not provided
    if (params.genome_fai) {
        fai_ch = Channel.fromPath(params.genome_fai)
    } else {
        fai_ch = Channel.fromPath("${params.genome}.fai")
    }

    // ---- Step 1: QC ----
    FASTQC(reads_ch)

    // ---- Step 2: Trim ----
    TRIM_GALORE(reads_ch)

    // ---- Step 3: Index genome ----
    BWAMETH_INDEX(genome_ch)

    // ---- Step 4: Align ----
    BWAMETH_ALIGN(
        TRIM_GALORE.out.trimmed_reads,
        genome_ch.first(),
        BWAMETH_INDEX.out.index.first()
    )

    // ---- Step 5: Mark duplicates ----
    MARK_DUPLICATES(BWAMETH_ALIGN.out.bam)

    // ---- Step 6: Methylation extraction ----
    // Join dedup BAM with its BAI by sample_id
    dedup_bam_bai = MARK_DUPLICATES.out.bam
        .join(MARK_DUPLICATES.out.bai)
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }

    METHYLDACKEL(
        dedup_bam_bai,
        genome_ch.first(),
        fai_ch.first()
    )

    // ---- Step 7: ComBat-meth batch correction ----
    // Collect all bedGraphs (strip sample_id, just get the files)
    bedgraphs_ch = METHYLDACKEL.out.bedgraph
        .map { sample_id, bg -> bg }
        .collect()

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

    // ---- Step 10: Region detection (DMRs) ----
    REGION_DETECT(
        COMBAT_METH.out.mvalues_corrected,
        sample_sheet_ch.first()
    )

    // ---- Step 11: NMF patient stratification ----
    NMF_STRATIFY(
        COMBAT_METH.out.mvalues_corrected,
        sample_sheet_ch.first()
    )
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
