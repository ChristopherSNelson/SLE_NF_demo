#!/usr/bin/env python3
"""
Generate simulated chr22 bisulfite sequencing test data.
Creates a small reference FASTA, paired-end FASTQs (4 samples), and sample sheet.
"""

import os
import random
import gzip
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='test_data', help='Output directory')
    parser.add_argument('--n_reads', type=int, default=1000, help='Reads per sample')
    parser.add_argument('--read_len', type=int, default=100, help='Read length')
    parser.add_argument('--ref_len', type=int, default=50000, help='Reference length')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    random.seed(42)

    # --- Generate chr22 reference ---
    ref_seq = ''.join(random.choices('ACGT', k=args.ref_len))
    ref_path = os.path.join(args.outdir, 'chr22.fa')
    with open(ref_path, 'w') as f:
        f.write('>chr22\n')
        for i in range(0, len(ref_seq), 80):
            f.write(ref_seq[i:i+80] + '\n')

    # Create .fai index
    fai_path = ref_path + '.fai'
    with open(fai_path, 'w') as f:
        # chr22\tlength\toffset\tlinebases\tlinewidth
        n_lines = (args.ref_len + 79) // 80
        last_line = args.ref_len % 80 if args.ref_len % 80 != 0 else 80
        f.write(f'chr22\t{args.ref_len}\t6\t80\t81\n')

    # --- Sample definitions ---
    samples = [
        {'sample_id': 'SLE_01',     'condition': 'SLE',     'batch': 'batch1'},
        {'sample_id': 'SLE_02',     'condition': 'SLE',     'batch': 'batch2'},
        {'sample_id': 'Control_01', 'condition': 'Control', 'batch': 'batch1'},
        {'sample_id': 'Control_02', 'condition': 'Control', 'batch': 'batch2'},
    ]

    def bisulfite_convert(seq, conversion_rate=0.99):
        """Simulate bisulfite conversion: C->T except at CpG sites (partially)."""
        result = list(seq)
        for i in range(len(result)):
            if result[i] == 'C':
                # If CpG context, partially convert based on methylation
                if i + 1 < len(result) and result[i + 1] == 'G':
                    if random.random() > 0.5:  # 50% methylation at CpGs
                        result[i] = 'T'
                else:
                    if random.random() < conversion_rate:
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

        with gzip.open(fq1_path, 'wt') as fq1, gzip.open(fq2_path, 'wt') as fq2:
            for i in range(args.n_reads):
                # Random fragment from reference
                frag_len = random.randint(150, 300)
                max_start = max(0, args.ref_len - frag_len)
                start = random.randint(0, max_start)
                fragment = ref_seq[start:start + frag_len]

                # Read 1: forward, bisulfite converted
                r1_seq = bisulfite_convert(fragment[:args.read_len])
                # Read 2: reverse complement of fragment end, bisulfite converted
                r2_raw = reverse_complement(fragment[-args.read_len:])
                r2_seq = bisulfite_convert(r2_raw)

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


if __name__ == '__main__':
    main()
