#!/usr/bin/env python3
"""
Generate simulated chr22 bisulfite sequencing test data.
Creates a small reference FASTA with planted CpG islands,
paired-end FASTQs (4 samples), and sample sheet.

The reference is designed to produce actual CpG coverage:
- 10kb sequence with 3 CpG-dense regions (~20 CpGs each)
- 5000 reads per sample concentrated on CpG regions
- This ensures MethylDackel finds CpG sites with sufficient depth
"""

import os
import random
import gzip
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='test_data', help='Output directory')
    parser.add_argument('--n_reads', type=int, default=5000, help='Reads per sample')
    parser.add_argument('--read_len', type=int, default=100, help='Read length')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    random.seed(42)

    # --- Generate reference with planted CpG islands ---
    # Build a 10kb sequence with 3 CpG-rich regions
    ref_len = 10000

    # Start with random non-CpG bases (A, T only to avoid accidental CpGs)
    ref_list = [random.choice('AT') for _ in range(ref_len)]

    # Plant CpG islands at known positions
    cpg_regions = [
        (1000, 1500),   # Region 1: positions 1000-1500
        (4000, 4500),   # Region 2: positions 4000-4500
        (7000, 7500),   # Region 3: positions 7000-7500
    ]

    for start, end in cpg_regions:
        pos = start
        while pos < end - 1:
            # Place CG dinucleotides every ~10-20 bp with some flanking sequence
            ref_list[pos] = 'C'
            ref_list[pos + 1] = 'G'
            # Add some non-CpG context around it
            gap = random.randint(8, 20)
            pos += gap

    ref_seq = ''.join(ref_list)

    # Count actual CpGs for verification
    n_cpgs = sum(1 for i in range(len(ref_seq) - 1) if ref_seq[i] == 'C' and ref_seq[i+1] == 'G')
    print(f"Reference: {ref_len} bp with {n_cpgs} CpG sites")

    ref_path = os.path.join(args.outdir, 'chr22.fa')
    with open(ref_path, 'w') as f:
        f.write('>chr22\n')
        for i in range(0, len(ref_seq), 80):
            f.write(ref_seq[i:i+80] + '\n')

    # Create .fai index
    fai_path = ref_path + '.fai'
    with open(fai_path, 'w') as f:
        f.write(f'chr22\t{ref_len}\t6\t80\t81\n')

    # --- Sample definitions ---
    samples = [
        {'sample_id': 'SLE_01',     'condition': 'SLE',     'batch': 'batch1'},
        {'sample_id': 'SLE_02',     'condition': 'SLE',     'batch': 'batch2'},
        {'sample_id': 'Control_01', 'condition': 'Control', 'batch': 'batch1'},
        {'sample_id': 'Control_02', 'condition': 'Control', 'batch': 'batch2'},
    ]

    def bisulfite_convert(seq, methylation_rate=0.5):
        """Simulate bisulfite conversion: C->T except at methylated CpG sites."""
        result = list(seq)
        for i in range(len(result)):
            if result[i] == 'C':
                if i + 1 < len(result) and result[i + 1] == 'G':
                    # CpG context: convert only if unmethylated
                    if random.random() > methylation_rate:
                        result[i] = 'T'
                else:
                    # Non-CpG C: always convert (bisulfite treatment)
                    result[i] = 'T'
        return ''.join(result)

    def reverse_complement(seq):
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(comp.get(b, 'N') for b in reversed(seq))

    def random_quality(length, min_q=30, max_q=40):
        return ''.join(chr(random.randint(min_q, max_q) + 33) for _ in range(length))

    # --- Generate FASTQs ---
    fastq_dir = os.path.join(args.outdir, 'fastqs')
    os.makedirs(fastq_dir, exist_ok=True)

    for sample in samples:
        sid = sample['sample_id']
        fq1_path = os.path.join(fastq_dir, f'{sid}_1.fastq.gz')
        fq2_path = os.path.join(fastq_dir, f'{sid}_2.fastq.gz')

        # SLE samples have different methylation rates than controls
        if sample['condition'] == 'SLE':
            meth_rate = 0.3  # hypomethylated in SLE
        else:
            meth_rate = 0.7  # normal methylation in controls

        with gzip.open(fq1_path, 'wt') as fq1, gzip.open(fq2_path, 'wt') as fq2:
            for i in range(args.n_reads):
                # 80% of reads target CpG regions, 20% random
                if random.random() < 0.8:
                    # Pick a CpG region
                    region_start, region_end = random.choice(cpg_regions)
                    # Start read near the region
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

                # Read 1: forward, bisulfite converted
                r1_seq = bisulfite_convert(fragment[:args.read_len], meth_rate)
                # Read 2: reverse complement of fragment end, bisulfite converted
                r2_raw = reverse_complement(fragment[-args.read_len:])
                r2_seq = bisulfite_convert(r2_raw, meth_rate)

                qual1 = random_quality(args.read_len)
                qual2 = random_quality(args.read_len)

                fq1.write(f'@{sid}_{i}/1\n{r1_seq}\n+\n{qual1}\n')
                fq2.write(f'@{sid}_{i}/2\n{r2_seq}\n+\n{qual2}\n')

        sample['fastq_1'] = os.path.abspath(fq1_path)
        sample['fastq_2'] = os.path.abspath(fq2_path)

    # --- Write sample sheet ---
    csv_path = os.path.join(args.outdir, 'samples.csv')
    with open(csv_path, 'w') as f:
        f.write('sample_id,fastq_1,fastq_2,condition,batch\n')
        for s in samples:
            f.write(f"{s['sample_id']},{s['fastq_1']},{s['fastq_2']},{s['condition']},{s['batch']}\n")

    print(f"Test data generated in {args.outdir}/")
    print(f"  Reference: {ref_path}")
    print(f"  Sample sheet: {csv_path}")
    print(f"  {len(samples)} samples x {args.n_reads} reads each")
    print(f"  CpG regions: {cpg_regions}")


if __name__ == '__main__':
    main()
