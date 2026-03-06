#!/usr/bin/env python3
"""
NMF consensus clustering for patient stratification.
Sweeps k=2..k_max, runs multiple NMF iterations per rank, and selects optimal k
using cophenetic correlation, dispersion, and silhouette scores.
"""

import argparse
import os
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform, pdist
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
    parser.add_argument('--outdir', default='.', help='Output directory')
    return parser.parse_args()


def make_nonnegative(mat):
    """Shift M-values to be non-negative for NMF input."""
    min_val = np.min(mat)
    if min_val < 0:
        mat = mat - min_val
    return mat


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

    # Prepare non-negative matrix
    data = make_nonnegative(mval_df.values)

    # Select top variable features for speed
    n_features = min(5000, data.shape[0])
    variances = np.var(data, axis=1)
    top_idx = np.argsort(variances)[-n_features:]
    data_subset = data[top_idx, :]

    print(f"Using top {n_features} variable CpGs for NMF")

    # Run NMF sweep
    results = run_nmf_sweep(data_subset, args.k_min, args.k_max, args.n_runs)

    # Select best k: highest silhouette (with cophenetic as tiebreaker)
    best_k = max(results.keys(),
                  key=lambda k: (results[k]['silhouette'], results[k]['cophenetic']))
    print(f"Selected k={best_k} (silhouette={results[best_k]['silhouette']:.3f})")

    best = results[best_k]

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
    w_df = pd.DataFrame(best['W'], index=mval_df.index[top_idx],
                          columns=[f'factor_{i+1}' for i in range(best_k)])
    w_df.to_csv(os.path.join(args.outdir, 'W_matrix.tsv'), sep='\t')

    h_df = pd.DataFrame(best['H'], index=[f'factor_{i+1}' for i in range(best_k)],
                          columns=mval_df.columns)
    h_df.to_csv(os.path.join(args.outdir, 'H_matrix.tsv'), sep='\t')

    # Plots
    plot_rank_selection(results, args.outdir)
    plot_umap(best['H'], best['assignments'], mval_df.columns.tolist(), args.outdir)

    print("NMF stratification complete")


if __name__ == '__main__':
    main()
