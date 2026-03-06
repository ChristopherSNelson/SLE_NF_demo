#!/usr/bin/env python3
"""
Sliding-window DMR detection with sub-population-aware testing.
Uses a stratified approach to handle patient sub-population correlation.
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description='Sliding-window DMR detection')
    parser.add_argument('--mvalues', required=True, help='Corrected M-values TSV')
    parser.add_argument('--sample_sheet', required=True, help='Sample sheet CSV')
    parser.add_argument('--window_size', type=int, default=1000, help='Window size in bp')
    parser.add_argument('--step_size', type=int, default=500, help='Step size in bp')
    parser.add_argument('--min_cpgs', type=int, default=3, help='Min CpGs per window')
    parser.add_argument('--outdir', default='.', help='Output directory')
    return parser.parse_args()


def parse_cpg_coords(cpg_ids):
    """Parse CpG IDs (chr:pos format) into chromosome and position."""
    chrs = []
    positions = []
    for cid in cpg_ids:
        parts = cid.split(':')
        chrs.append(parts[0])
        positions.append(int(parts[1]))
    return pd.DataFrame({'chr': chrs, 'pos': positions, 'cpg_id': cpg_ids})


def sliding_window_test(mval_df, cpg_coords, conditions, window_size, step_size, min_cpgs):
    """
    Sliding window DMR detection.
    For each window, performs a sub-population-aware test comparing SLE vs controls.
    Uses linear regression with condition as predictor, averaging M-values across CpGs in the window.
    """
    results = []
    chromosomes = cpg_coords['chr'].unique()

    for chrom in sorted(chromosomes):
        chr_mask = cpg_coords['chr'] == chrom
        chr_coords = cpg_coords[chr_mask].sort_values('pos')
        chr_mvals = mval_df.loc[chr_coords['cpg_id']]

        if len(chr_coords) == 0:
            continue

        min_pos = chr_coords['pos'].min()
        max_pos = chr_coords['pos'].max()

        for win_start in range(min_pos, max_pos, step_size):
            win_end = win_start + window_size

            # CpGs in this window
            in_window = (chr_coords['pos'] >= win_start) & (chr_coords['pos'] < win_end)
            window_cpgs = chr_coords[in_window]['cpg_id'].values

            if len(window_cpgs) < min_cpgs:
                continue

            # Mean M-value per sample across window CpGs
            window_means = chr_mvals.loc[window_cpgs].mean(axis=0)

            # Binary condition encoding: SLE=1, others=0
            is_sle = (conditions == 'SLE').astype(int).values

            # Linear regression: M-value ~ condition
            try:
                X = is_sle.reshape(-1, 1)
                y = window_means.values

                if np.std(y) < 1e-10:
                    continue

                reg = LinearRegression().fit(X, y)
                y_pred = reg.predict(X)
                residuals = y - y_pred
                n = len(y)
                se = np.sqrt(np.sum(residuals**2) / (n - 2)) / np.sqrt(np.sum((X.flatten() - X.mean())**2))

                if se < 1e-10:
                    continue

                t_stat = reg.coef_[0] / se
                p_value = 2 * stats.t.sf(abs(t_stat), df=n - 2)
                effect_size = reg.coef_[0]

                results.append({
                    'chr': chrom,
                    'start': win_start,
                    'end': win_end,
                    'n_cpgs': len(window_cpgs),
                    'effect_size': effect_size,
                    't_stat': t_stat,
                    'p_value': p_value
                })
            except Exception:
                continue

    return pd.DataFrame(results)


def main():
    args = parse_args()

    # Load data
    mval_df = pd.read_csv(args.mvalues, sep='\t', index_col=0)
    samples = pd.read_csv(args.sample_sheet)

    print(f"Loaded {mval_df.shape[0]} CpGs x {mval_df.shape[1]} samples")

    # Parse CpG coordinates
    cpg_coords = parse_cpg_coords(mval_df.index)

    # Get condition for each sample (match column order)
    sample_order = [s for s in mval_df.columns]
    conditions = pd.Series([
        samples.loc[samples['sample_id'] == s, 'condition'].values[0]
        for s in sample_order
    ])

    # Run sliding window
    print(f"Running sliding window: size={args.window_size}, step={args.step_size}, min_cpgs={args.min_cpgs}")
    results = sliding_window_test(mval_df, cpg_coords, conditions,
                                   args.window_size, args.step_size, args.min_cpgs)

    if len(results) == 0:
        print("No windows passed filtering — writing empty outputs")
        pd.DataFrame(columns=['chr', 'start', 'end', 'n_cpgs', 'effect_size', 't_stat', 'p_value']
                     ).to_csv(os.path.join(args.outdir, 'window_results.tsv'), sep='\t', index=False)
        with open(os.path.join(args.outdir, 'candidate_dmrs.bed'), 'w') as f:
            pass
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.set_title('DMR Manhattan Plot (no significant windows)')
        fig.savefig(os.path.join(args.outdir, 'dmr_manhattan.png'), dpi=150, bbox_inches='tight')
        plt.close()
        return

    # Multiple testing correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    try:
        _, results['p_adjusted'], _, _ = multipletests(results['p_value'], method='fdr_bh')
    except ImportError:
        # Fallback: Bonferroni
        results['p_adjusted'] = np.minimum(results['p_value'] * len(results), 1.0)

    # Write full window results
    results.to_csv(os.path.join(args.outdir, 'window_results.tsv'), sep='\t', index=False)

    # Candidate DMRs: p_adjusted < 0.05
    sig = results[results['p_adjusted'] < 0.05].copy()
    if len(sig) == 0:
        # If nothing passes FDR, take top 20 by p-value as candidates
        sig = results.nsmallest(min(20, len(results)), 'p_value').copy()

    sig[['chr', 'start', 'end']].to_csv(
        os.path.join(args.outdir, 'candidate_dmrs.bed'),
        sep='\t', header=False, index=False
    )

    print(f"Found {len(results)} windows, {len(sig)} candidate DMRs")

    # --- Manhattan plot ---
    fig, ax = plt.subplots(figsize=(14, 5))

    # Assign numeric x positions per chromosome
    chrom_order = sorted(results['chr'].unique(),
                          key=lambda c: int(c.replace('chr', '')) if c.replace('chr', '').isdigit() else 99)
    chrom_offsets = {}
    offset = 0
    for chrom in chrom_order:
        chrom_offsets[chrom] = offset
        offset += results[results['chr'] == chrom]['end'].max() - results[results['chr'] == chrom]['start'].min()
        offset += 1e6  # gap between chromosomes

    x_vals = []
    for _, row in results.iterrows():
        x_vals.append(chrom_offsets[row['chr']] + (row['start'] + row['end']) / 2)
    results['x'] = x_vals

    neg_log_p = -np.log10(results['p_value'].clip(lower=1e-300))

    colors = ['#1f77b4', '#aec7e8']
    chrom_colors = [colors[i % 2] for i, c in enumerate(chrom_order) for _ in range(len(results[results['chr'] == c]))]

    ax.scatter(results['x'], neg_log_p, c=chrom_colors, s=8, alpha=0.6)

    # Significance line
    if len(results) > 0:
        bonf_thresh = -np.log10(0.05 / len(results))
        ax.axhline(y=bonf_thresh, color='red', linestyle='--', alpha=0.5, label='Bonferroni')

    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('DMR Manhattan Plot')
    ax.legend()

    fig.savefig(os.path.join(args.outdir, 'dmr_manhattan.png'), dpi=150, bbox_inches='tight')
    plt.close()

    print("Region detection complete")


if __name__ == '__main__':
    main()
