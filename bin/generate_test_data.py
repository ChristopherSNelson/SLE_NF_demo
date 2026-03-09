#!/usr/bin/env python3
"""
Generate simulated chr22 bisulfite sequencing test data.
Creates a small reference FASTA with planted CpG islands,
paired-end FASTQs (4 samples), sample sheet, and pre-made bedGraph files.

The bedGraph files allow testing downstream steps (ComBat → NMF) without
requiring successful bisulfite alignment, which is unreliable with tiny
simulated data.

The reference is designed to produce actual CpG coverage:
- 10kb sequence with 3 CpG-dense regions (~20 CpGs each)
- 5000 reads per sample concentrated on CpG regions
"""

import os
import random
import gzip
import argparse
import math


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='test_data', help='Output directory')
    parser.add_argument('--n_reads', type=int, default=5000, help='Reads per sample')
    parser.add_argument('--read_len', type=int, default=100, help='Read length')
    parser.add_argument('--skip_fastqs', action='store_true',
                        help='Skip FASTQ generation (only produce bedGraphs and sample sheet)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    random.seed(42)

    # --- Sample definitions ---
    # 6 samples (3 per batch) — ComBatMet needs >=3 per batch for variance estimation
    samples = [
        {'sample_id': 'SLE_01',     'condition': 'SLE',     'batch': 'batch1'},
        {'sample_id': 'SLE_02',     'condition': 'SLE',     'batch': 'batch2'},
        {'sample_id': 'SLE_03',     'condition': 'SLE',     'batch': 'batch1'},
        {'sample_id': 'Control_01', 'condition': 'Control', 'batch': 'batch1'},
        {'sample_id': 'Control_02', 'condition': 'Control', 'batch': 'batch2'},
        {'sample_id': 'Control_03', 'condition': 'Control', 'batch': 'batch2'},
    ]

    # --- CpG site positions (shared across all outputs) ---
    # 200 CpG sites spread across chr22 in 3 regions, mimicking CpG islands
    cpg_regions = [
        (1000, 1500),   # Region 1: positions 1000-1500
        (4000, 4500),   # Region 2: positions 4000-4500
        (7000, 7500),   # Region 3: positions 7000-7500
    ]

    # Generate CpG positions within regions (every ~10-20 bp)
    cpg_positions = []
    for start, end in cpg_regions:
        pos = start
        while pos < end - 1:
            cpg_positions.append(pos)
            pos += random.randint(8, 20)

    # Add some scattered CpGs outside islands
    for _ in range(40):
        pos = random.randint(0, 9998)
        # Avoid duplicates and region overlaps
        if pos not in cpg_positions:
            cpg_positions.append(pos)
    cpg_positions.sort()

    print(f"Generated {len(cpg_positions)} CpG site positions")

    # --- Generate reference FASTA with planted CpGs ---
    ref_len = 10000
    ref_list = [random.choice('AT') for _ in range(ref_len)]
    for pos in cpg_positions:
        if pos + 1 < ref_len:
            ref_list[pos] = 'C'
            ref_list[pos + 1] = 'G'
    ref_seq = ''.join(ref_list)

    ref_path = os.path.join(args.outdir, 'chr22.fa')
    with open(ref_path, 'w') as f:
        f.write('>chr22\n')
        for i in range(0, len(ref_seq), 80):
            f.write(ref_seq[i:i+80] + '\n')

    fai_path = ref_path + '.fai'
    with open(fai_path, 'w') as f:
        f.write(f'chr22\t{ref_len}\t6\t80\t81\n')

    print(f"Reference: {ref_len} bp with {len(cpg_positions)} CpG sites")

    # --- Generate synthetic bedGraph files ---
    # These allow testing ComBat → NMF without requiring successful alignment
    bedgraph_dir = os.path.join(args.outdir, 'bedgraphs')
    os.makedirs(bedgraph_dir, exist_ok=True)

    for sample in samples:
        sid = sample['sample_id']
        bg_path = os.path.join(bedgraph_dir, f'{sid}_CpG.bedGraph')

        # Base methylation differs by condition (SLE = hypomethylated)
        if sample['condition'] == 'SLE':
            base_beta = 0.30
        else:
            base_beta = 0.70

        # Batch effect: shift beta by ±0.05 depending on batch
        if sample['batch'] == 'batch1':
            batch_shift = 0.05
        else:
            batch_shift = -0.05

        with open(bg_path, 'w') as f:
            # MethylDackel track header line (combat_meth.R skips this with skip=1)
            f.write(f'track type=bedGraph description="{sid} CpG methylation"\n')
            for pos in cpg_positions:
                # Per-site beta with biological + batch + noise variation
                # Region-specific effects: region 2 has stronger SLE signal
                region_boost = 0.0
                if 4000 <= pos <= 4500 and sample['condition'] == 'SLE':
                    region_boost = -0.15  # extra hypomethylation in region 2 for SLE

                beta = base_beta + batch_shift + region_boost + random.gauss(0, 0.08)
                beta = max(0.01, min(0.99, beta))  # clamp to (0,1)

                # Simulate coverage (20-80x) and compute count_m / count_u
                total_reads = random.randint(20, 80)
                count_m = round(beta * total_reads)
                count_u = total_reads - count_m

                # MethylDackel format: chr start end beta% count_m count_u
                beta_pct = (count_m / total_reads) * 100
                f.write(f'chr22\t{pos}\t{pos + 1}\t{beta_pct:.1f}\t{count_m}\t{count_u}\n')

    print(f"  bedGraphs: {bedgraph_dir}/ ({len(cpg_positions)} CpGs x {len(samples)} samples)")

    # --- Generate FASTQs (optional, for full pipeline testing) ---
    if not args.skip_fastqs:
        _generate_fastqs(args, samples, ref_seq, cpg_regions, ref_len)

    # --- Write sample sheet ---
    csv_path = os.path.join(args.outdir, 'samples.csv')
    with open(csv_path, 'w') as f:
        f.write('sample_id,fastq_1,fastq_2,condition,batch\n')
        for s in samples:
            fq1 = s.get('fastq_1', '')
            fq2 = s.get('fastq_2', '')
            f.write(f"{s['sample_id']},{fq1},{fq2},{s['condition']},{s['batch']}\n")

    print(f"Test data generated in {args.outdir}/")
    print(f"  Reference: {ref_path}")
    print(f"  Sample sheet: {csv_path}")
    print(f"  {len(samples)} samples, {len(cpg_positions)} CpG sites")


def _generate_fastqs(args, samples, ref_seq, cpg_regions, ref_len):
    """Generate simulated bisulfite-converted paired-end FASTQs."""
    fastq_dir = os.path.join(args.outdir, 'fastqs')
    os.makedirs(fastq_dir, exist_ok=True)

    def bisulfite_convert(seq, methylation_rate=0.5):
        result = list(seq)
        for i in range(len(result)):
            if result[i] == 'C':
                if i + 1 < len(result) and result[i + 1] == 'G':
                    if random.random() > methylation_rate:
                        result[i] = 'T'
                else:
                    result[i] = 'T'
        return ''.join(result)

    def reverse_complement(seq):
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(comp.get(b, 'N') for b in reversed(seq))

    def random_quality(length, min_q=30, max_q=40):
        return ''.join(chr(random.randint(min_q, max_q) + 33) for _ in range(length))

    for sample in samples:
        sid = sample['sample_id']
        fq1_path = os.path.join(fastq_dir, f'{sid}_1.fastq.gz')
        fq2_path = os.path.join(fastq_dir, f'{sid}_2.fastq.gz')

        if sample['condition'] == 'SLE':
            meth_rate = 0.3
        else:
            meth_rate = 0.7

        with gzip.open(fq1_path, 'wt') as fq1, gzip.open(fq2_path, 'wt') as fq2:
            for i in range(args.n_reads):
                if random.random() < 0.8:
                    region_start, region_end = random.choice(cpg_regions)
                    frag_start = random.randint(
                        max(0, region_start - args.read_len),
                        min(region_end, ref_len - 200)
                    )
                else:
                    frag_start = random.randint(0, ref_len - 200)

                frag_len = random.randint(150, 200)
                frag_end = min(frag_start + frag_len, ref_len)
                fragment = ref_seq[frag_start:frag_end]

                if len(fragment) < args.read_len:
                    continue

                r1_seq = bisulfite_convert(fragment[:args.read_len], meth_rate)
                r2_raw = reverse_complement(fragment[-args.read_len:])
                r2_seq = bisulfite_convert(r2_raw, meth_rate)

                qual1 = random_quality(args.read_len)
                qual2 = random_quality(args.read_len)

                fq1.write(f'@{sid}_{i}/1\n{r1_seq}\n+\n{qual1}\n')
                fq2.write(f'@{sid}_{i}/2\n{r2_seq}\n+\n{qual2}\n')

        sample['fastq_1'] = os.path.abspath(fq1_path)
        sample['fastq_2'] = os.path.abspath(fq2_path)


if __name__ == '__main__':
    main()
