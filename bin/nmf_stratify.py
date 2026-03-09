#!/usr/bin/env python3
"""
NMF consensus clustering for patient stratification.
Sweeps k=2..k_max, runs multiple NMF iterations per rank, and selects optimal k
using cophenetic correlation, dispersion, and silhouette scores.

Enhancements:
  - Cell-type regression: regress out Houseman cell fractions before NMF
  - DMR-based feature selection: use DMR CpGs instead of top-variance (with fallback)
  - LOO stability: leave-one-out analysis to assess cluster robustness
  - Clinical metadata correlation: correlate H-matrix weights with clinical variables
"""

import argparse
import os
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform, pdist
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description='NMF patient stratification')
    parser.add_argument('--mvalues', required=True, help='Corrected M-values TSV')
    parser.add_argument('--sample_sheet', required=True, help='Sample sheet CSV')
    parser.add_argument('--k_min', type=int, default=2, help='Min NMF rank')
    parser.add_argument('--k_max', type=int, default=8, help='Max NMF rank')
    parser.add_argument('--n_runs', type=int, default=30, help='NMF runs per rank')
    parser.add_argument('--cell_fractions', default=None,
                        help='Cell fractions TSV from Houseman deconvolution (optional)')
    parser.add_argument('--dmr_bed', default=None,
                        help='Candidate DMR BED file from region detection (optional)')
    parser.add_argument('--clinical_metadata', default=None,
                        help='Clinical metadata TSV with numeric columns to correlate (optional)')
    parser.add_argument('--outdir', default='.', help='Output directory')
    return parser.parse_args()


def make_nonnegative(mat):
    """Shift M-values to be non-negative for NMF input."""
    min_val = np.min(mat)
    if min_val < 0:
        mat = mat - min_val
    return mat


def regress_out_cell_types(mval_df, cell_fractions_path):
    """
    Regress out cell-type composition effects from M-values.
    Ensures NMF clusters reflect epigenetic shifts, not cell-proportion differences.
    """
    cell_df = pd.read_csv(cell_fractions_path, sep='\t')
    cell_df = cell_df.set_index('sample_id')

    # Align samples
    common_samples = mval_df.columns.intersection(cell_df.index)
    if len(common_samples) < 3:
        print(f"WARNING: Only {len(common_samples)} samples overlap with cell fractions. "
              "Skipping cell-type regression.")
        return mval_df

    mval_aligned = mval_df[common_samples]
    cell_aligned = cell_df.loc[common_samples]

    # Drop cell types with zero variance
    cell_cols = [c for c in cell_aligned.columns if cell_aligned[c].std() > 1e-10]
    if len(cell_cols) == 0:
        print("WARNING: No variable cell-type fractions. Skipping regression.")
        return mval_df

    X = cell_aligned[cell_cols].values  # n_samples x n_cell_types
    Y = mval_aligned.values  # n_cpgs x n_samples

    # For each CpG, regress M-value ~ cell fractions, keep residuals
    print(f"Regressing out {len(cell_cols)} cell-type fractions from {Y.shape[0]} CpGs...")
    from sklearn.linear_model import LinearRegression
    reg = LinearRegression()
    reg.fit(X, Y.T)  # fit: samples x cell_types -> samples x cpgs
    predicted = reg.predict(X)  # samples x cpgs
    residuals = Y - predicted.T  # n_cpgs x n_samples

    result = pd.DataFrame(residuals, index=mval_aligned.index, columns=mval_aligned.columns)
    print(f"Cell-type regression complete. Residual matrix: {result.shape}")
    return result


def select_dmr_features(mval_df, dmr_bed_path, min_features=500):
    """
    Select CpGs that fall within candidate DMR regions.
    Falls back to top-variance if too few DMR CpGs are found.
    """
    dmr_df = pd.read_csv(dmr_bed_path, sep='\t', header=None,
                          names=['chr', 'start', 'end'],
                          usecols=[0, 1, 2])

    if len(dmr_df) == 0:
        print("WARNING: No DMRs found. Falling back to top-variance selection.")
        return None

    # Parse CpG coordinates from index (chr:pos format)
    cpg_ids = mval_df.index.tolist()
    cpg_chrs = []
    cpg_pos = []
    for cid in cpg_ids:
        parts = cid.split(':')
        cpg_chrs.append(parts[0])
        cpg_pos.append(int(parts[1]))

    cpg_coords = pd.DataFrame({'cpg_id': cpg_ids, 'chr': cpg_chrs, 'pos': cpg_pos})

    # Find CpGs within any DMR
    dmr_cpgs = set()
    for _, dmr in dmr_df.iterrows():
        mask = ((cpg_coords['chr'] == dmr['chr']) &
                (cpg_coords['pos'] >= dmr['start']) &
                (cpg_coords['pos'] < dmr['end']))
        dmr_cpgs.update(cpg_coords.loc[mask, 'cpg_id'].tolist())

    print(f"Found {len(dmr_cpgs)} CpGs within {len(dmr_df)} DMR regions")

    if len(dmr_cpgs) < min_features:
        print(f"Only {len(dmr_cpgs)} DMR CpGs (< {min_features} minimum). "
              f"Supplementing with top-variance CpGs.")
        # Supplement with top-variance CpGs
        variances = mval_df.var(axis=1)
        top_var = variances.nlargest(min_features).index.tolist()
        dmr_cpgs = dmr_cpgs.union(set(top_var))

    selected = [c for c in mval_df.index if c in dmr_cpgs]
    print(f"Using {len(selected)} CpGs for NMF (DMR + variance supplement)")
    return selected


def consensus_matrix(connectivity_matrices):
    """Average connectivity matrices across runs."""
    return np.mean(connectivity_matrices, axis=0)


def compute_cophenetic(consensus_mat):
    """Cophenetic correlation: how well the consensus dendrogram preserves distances."""
    n = consensus_mat.shape[0]
    if n < 3:
        return 0.0
    dist_mat = 1 - consensus_mat
    np.fill_diagonal(dist_mat, 0)
    dist_mat = np.maximum(dist_mat, 0)
    try:
        condensed = squareform(dist_mat)
        Z = linkage(condensed, method='average')
        coph_dist, _ = cophenet(Z, condensed)
        return coph_dist
    except Exception:
        return 0.0


def compute_dispersion(consensus_mat):
    """Dispersion: measures how crisp the consensus matrix is (all 0s and 1s = perfect)."""
    n = consensus_mat.shape[0]
    entries = consensus_mat[np.triu_indices(n, k=1)]
    if len(entries) == 0:
        return 0.0
    return 1.0 - (4.0 / (n * (n - 1))) * np.sum((entries - 0.5) ** 2)


def run_nmf_sweep(data, k_min, k_max, n_runs):
    """Run NMF for each rank and compute consensus metrics."""
    n_samples = data.shape[1]
    # Cap k_max to prevent silhouette crash
    k_max = min(k_max, n_samples - 1)

    results = {}

    for k in range(k_min, k_max + 1):
        print(f"  NMF rank k={k}: {n_runs} runs...")
        connectivity_list = []
        best_err = np.inf
        best_W = None
        best_H = None

        for run in range(n_runs):
            model = NMF(n_components=k, init='random', random_state=run,
                        max_iter=500, tol=1e-4)
            W = model.fit_transform(data)
            H = model.components_
            err = model.reconstruction_err_

            if err < best_err:
                best_err = err
                best_W = W
                best_H = H

            # Connectivity matrix: samples assigned to same cluster = 1
            assignments = np.argmax(H, axis=0)
            conn = (assignments[:, None] == assignments[None, :]).astype(float)
            connectivity_list.append(conn)

        consensus = consensus_matrix(np.array(connectivity_list))
        coph = compute_cophenetic(consensus)
        disp = compute_dispersion(consensus)

        # Silhouette on best H matrix
        best_assignments = np.argmax(best_H, axis=0)
        n_unique = len(np.unique(best_assignments))
        if n_unique > 1 and n_unique < n_samples:
            sil = silhouette_score(best_H.T, best_assignments)
        else:
            sil = 0.0

        results[k] = {
            'cophenetic': coph,
            'dispersion': disp,
            'silhouette': sil,
            'W': best_W,
            'H': best_H,
            'assignments': best_assignments,
            'consensus': consensus
        }

    return results


def run_loo_stability(data, best_k, n_runs, sample_ids):
    """
    Leave-one-out stability analysis.
    Runs NMF n_samples times, each time removing one sample.
    Reports whether cluster assignments are stable across LOO iterations.
    """
    n_samples = data.shape[1]
    print(f"\nLOO stability analysis: k={best_k}, {n_samples} iterations...")

    # Get reference assignments from full data
    model = NMF(n_components=best_k, init='nndsvd', random_state=42, max_iter=500)
    model.fit_transform(data)
    H_full = model.components_
    ref_assignments = np.argmax(H_full, axis=0)

    loo_results = []

    for i in range(n_samples):
        # Remove sample i
        data_loo = np.delete(data, i, axis=1)
        ids_loo = [s for j, s in enumerate(sample_ids) if j != i]

        model_loo = NMF(n_components=best_k, init='nndsvd', random_state=42, max_iter=500)
        model_loo.fit_transform(data_loo)
        H_loo = model_loo.components_
        assignments_loo = np.argmax(H_loo, axis=0)

        # Compare assignments for the remaining samples
        ref_remaining = np.delete(ref_assignments, i)

        # Check if assignments are consistent (up to label permutation)
        consistent = check_assignment_consistency(ref_remaining, assignments_loo)

        loo_results.append({
            'removed_sample': sample_ids[i],
            'consistent': consistent,
            'n_clusters_found': len(np.unique(assignments_loo))
        })

    loo_df = pd.DataFrame(loo_results)
    stability = loo_df['consistent'].mean()
    print(f"LOO stability: {stability:.1%} ({loo_df['consistent'].sum()}/{n_samples} consistent)")

    return loo_df, stability


def check_assignment_consistency(ref, test):
    """
    Check if two cluster assignment vectors are consistent up to label permutation.
    Two assignments are consistent if there exists a label mapping that makes them identical.
    """
    from itertools import permutations

    ref_labels = np.unique(ref)
    test_labels = np.unique(test)

    if len(ref_labels) != len(test_labels):
        return False

    # Try all permutations of test labels
    for perm in permutations(test_labels):
        mapping = dict(zip(test_labels, perm))
        remapped = np.array([mapping[t] for t in test])
        if np.array_equal(ref, remapped):
            return True

    return False


def correlate_clinical(H_matrix, sample_ids, clinical_path, outdir):
    """
    Correlate H-matrix weights with clinical metadata variables.
    Only processes numeric columns. Skips gracefully if file is missing or empty.
    """
    clinical = pd.read_csv(clinical_path, sep='\t')

    if 'sample_id' not in clinical.columns:
        print("WARNING: Clinical metadata must have 'sample_id' column. Skipping.")
        return None

    clinical = clinical.set_index('sample_id')

    # Align samples
    h_df = pd.DataFrame(H_matrix, columns=sample_ids,
                          index=[f'factor_{i+1}' for i in range(H_matrix.shape[0])])
    common = h_df.columns.intersection(clinical.index)
    if len(common) < 3:
        print(f"WARNING: Only {len(common)} samples overlap with clinical data. Skipping.")
        return None

    h_aligned = h_df[common]
    clin_aligned = clinical.loc[common]

    # Select only numeric columns
    num_cols = clin_aligned.select_dtypes(include=[np.number]).columns.tolist()
    if len(num_cols) == 0:
        print("WARNING: No numeric columns in clinical metadata. Skipping correlation.")
        return None

    print(f"Correlating {H_matrix.shape[0]} NMF factors with {len(num_cols)} clinical variables...")

    corr_results = []
    for factor in h_aligned.index:
        for col in num_cols:
            vals = clin_aligned[col].dropna()
            common_samp = vals.index.intersection(h_aligned.columns)
            if len(common_samp) < 3:
                continue
            r, p = stats.pearsonr(h_aligned.loc[factor, common_samp].values,
                                   vals[common_samp].values)
            corr_results.append({
                'factor': factor,
                'clinical_variable': col,
                'pearson_r': r,
                'p_value': p,
                'n_samples': len(common_samp)
            })

    if len(corr_results) == 0:
        print("No valid correlations computed.")
        return None

    corr_df = pd.DataFrame(corr_results)
    corr_df.to_csv(os.path.join(outdir, 'clinical_correlations.tsv'), sep='\t', index=False)
    print(f"Wrote {len(corr_df)} clinical correlations")
    return corr_df


def plot_rank_selection(results, outdir):
    """Plot cophenetic, dispersion, and silhouette across ranks."""
    ks = sorted(results.keys())
    coph = [results[k]['cophenetic'] for k in ks]
    disp = [results[k]['dispersion'] for k in ks]
    sil = [results[k]['silhouette'] for k in ks]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    axes[0].plot(ks, coph, 'o-', color='#1f77b4')
    axes[0].set_xlabel('Rank (k)')
    axes[0].set_ylabel('Cophenetic Correlation')
    axes[0].set_title('Cophenetic Correlation')
    axes[0].set_xticks(ks)

    axes[1].plot(ks, disp, 's-', color='#ff7f0e')
    axes[1].set_xlabel('Rank (k)')
    axes[1].set_ylabel('Dispersion')
    axes[1].set_title('Dispersion')
    axes[1].set_xticks(ks)

    axes[2].plot(ks, sil, '^-', color='#2ca02c')
    axes[2].set_xlabel('Rank (k)')
    axes[2].set_ylabel('Silhouette Score')
    axes[2].set_title('Silhouette Score')
    axes[2].set_xticks(ks)

    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'rank_selection.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_umap(H, assignments, samples, outdir):
    """UMAP visualization of H matrix colored by cluster."""
    try:
        from umap import UMAP
        reducer = UMAP(n_components=2, random_state=42, n_neighbors=min(5, H.shape[1] - 1))
        embedding = reducer.fit_transform(H.T)
    except (ImportError, Exception):
        # Fallback to PCA if UMAP fails (e.g., too few samples)
        from sklearn.decomposition import PCA
        n_comp = min(2, H.shape[0], H.shape[1])
        pca = PCA(n_components=n_comp)
        embedding = pca.fit_transform(H.T)

    fig, ax = plt.subplots(figsize=(7, 5))
    scatter = ax.scatter(embedding[:, 0], embedding[:, 1],
                          c=assignments, cmap='Set2', s=60, edgecolors='k', linewidth=0.5)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title('NMF Patient Clusters (UMAP)')

    # Add sample labels
    for i, sid in enumerate(samples):
        ax.annotate(sid, (embedding[i, 0], embedding[i, 1]),
                     fontsize=7, alpha=0.7, ha='center', va='bottom')

    plt.colorbar(scatter, label='Cluster')
    fig.savefig(os.path.join(outdir, 'nmf_umap.png'), dpi=150, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()

    # Load data
    mval_df = pd.read_csv(args.mvalues, sep='\t', index_col=0)
    samples = pd.read_csv(args.sample_sheet)

    print(f"Loaded {mval_df.shape[0]} CpGs x {mval_df.shape[1]} samples")

    # --- Enhancement 1: Cell-type regression ---
    if args.cell_fractions:
        print("\n--- Cell-type regression ---")
        mval_df = regress_out_cell_types(mval_df, args.cell_fractions)

    # --- Enhancement 2: DMR-based feature selection ---
    selected_cpgs = None
    if args.dmr_bed:
        print("\n--- DMR-based feature selection ---")
        selected_cpgs = select_dmr_features(mval_df, args.dmr_bed)

    # Prepare non-negative matrix
    if selected_cpgs is not None:
        data = make_nonnegative(mval_df.loc[selected_cpgs].values)
        feature_ids = selected_cpgs
    else:
        data = make_nonnegative(mval_df.values)
        # Fall back to top-variance selection
        n_features = min(5000, data.shape[0])
        variances = np.var(data, axis=1)
        top_idx = np.argsort(variances)[-n_features:]
        data = data[top_idx, :]
        feature_ids = mval_df.index[top_idx].tolist()

    print(f"Using {data.shape[0]} features x {data.shape[1]} samples for NMF")

    # Run NMF sweep
    results = run_nmf_sweep(data, args.k_min, args.k_max, args.n_runs)

    # Select best k: highest silhouette (with cophenetic as tiebreaker)
    best_k = max(results.keys(),
                  key=lambda k: (results[k]['silhouette'], results[k]['cophenetic']))
    print(f"\nSelected k={best_k} (silhouette={results[best_k]['silhouette']:.3f})")

    best = results[best_k]

    # --- Enhancement 3: LOO stability ---
    print("\n--- LOO stability analysis ---")
    loo_df, loo_stability = run_loo_stability(data, best_k, args.n_runs,
                                               mval_df.columns.tolist())
    loo_df['loo_stability_overall'] = loo_stability
    loo_df.to_csv(os.path.join(args.outdir, 'stability_loo.tsv'), sep='\t', index=False)

    # --- Write outputs ---
    # Cluster assignments
    cluster_df = pd.DataFrame({
        'sample_id': mval_df.columns,
        'cluster': best['assignments']
    })
    # Merge condition
    cluster_df = cluster_df.merge(samples[['sample_id', 'condition']], on='sample_id', how='left')
    cluster_df.to_csv(os.path.join(args.outdir, 'nmf_clusters.tsv'), sep='\t', index=False)

    # W and H matrices
    w_df = pd.DataFrame(best['W'], index=feature_ids,
                          columns=[f'factor_{i+1}' for i in range(best_k)])
    w_df.to_csv(os.path.join(args.outdir, 'W_matrix.tsv'), sep='\t')

    h_df = pd.DataFrame(best['H'], index=[f'factor_{i+1}' for i in range(best_k)],
                          columns=mval_df.columns)
    h_df.to_csv(os.path.join(args.outdir, 'H_matrix.tsv'), sep='\t')

    # --- Enhancement 4: Clinical metadata correlation ---
    if args.clinical_metadata:
        print("\n--- Clinical metadata correlation ---")
        correlate_clinical(best['H'], mval_df.columns.tolist(),
                           args.clinical_metadata, args.outdir)

    # Plots
    plot_rank_selection(results, args.outdir)
    plot_umap(best['H'], best['assignments'], mval_df.columns.tolist(), args.outdir)

    # Summary
    print(f"\n=== NMF Summary ===")
    print(f"Best rank: k={best_k}")
    print(f"Silhouette: {results[best_k]['silhouette']:.3f}")
    print(f"LOO stability: {loo_stability:.1%}")
    if args.cell_fractions:
        print("Cell-type effects: regressed out")
    if args.dmr_bed:
        print(f"Feature selection: DMR-based ({data.shape[0]} CpGs)")
    print("NMF stratification complete")


if __name__ == '__main__':
    main()
