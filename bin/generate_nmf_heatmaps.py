#!/usr/bin/env python3
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Patch

def main():
    parser = argparse.ArgumentParser(description='Generate NMF heatmaps for a single result set')
    parser.add_argument('--h_matrix', required=True, help='Path to H_matrix.tsv')
    parser.add_argument('--clusters', required=True, help='Path to nmf_clusters.tsv')
    parser.add_argument('--outdir', default='.', help='Output directory for plots')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    try:
        # Load the H-matrix
        h_df = pd.read_csv(args.h_matrix, sep='\t', index_col=0)
        clusters = pd.read_csv(args.clusters, sep='\t')
        
        # Align samples
        common_samples = h_df.columns.intersection(clusters['sample_id'])
        h_aligned = h_df[common_samples]
        clusters_aligned = clusters.set_index('sample_id').loc[common_samples]

        # --- 1. NMF-Ordered Heatmap ---
        sort_order = clusters_aligned.sort_values('cluster').index
        h_sorted = h_aligned[sort_order]
        clusters_sorted = clusters_aligned.loc[sort_order]

        cluster_ids = clusters_sorted['cluster'].unique()
        palette = sns.color_palette("Set2", len(cluster_ids))
        lut = dict(zip(cluster_ids, palette))
        col_colors = clusters_sorted['cluster'].map(lut)

        g = sns.clustermap(h_sorted, row_cluster=True, col_cluster=False, 
                           col_colors=col_colors, cmap='magma', figsize=(10, 8),
                           annot=True, fmt=".2e", cbar_kws={'label': 'NMF Factor Weight'})
        
        handles = [Patch(facecolor=lut[cid]) for cid in cluster_ids]
        plt.legend(handles, [f"Cluster {cid}" for cid in cluster_ids], 
                   title='NMF Cluster', bbox_to_anchor=(1, 1), 
                   bbox_transform=plt.gcf().transFigure, loc='upper right')

        g.fig.suptitle('NMF Factor Weights (Ordered by Cluster)', fontsize=16, y=1.05)
        g.savefig(os.path.join(args.outdir, "H_matrix_nmf_ordered.png"), dpi=300, bbox_inches='tight')
        plt.close('all')

        # --- 2. Clustered Heatmap ---
        if h_df.shape[0] >= 2 and h_df.shape[1] >= 2:
            g = sns.clustermap(h_df, cmap='magma', method='ward', metric='euclidean',
                               figsize=(10, 8), annot=True, fmt=".2e")
            g.fig.suptitle('NMF Factor Weights (Hierarchical Clustering)', fontsize=16, y=1.05)
            g.savefig(os.path.join(args.outdir, "H_matrix_clustered.png"), dpi=300, bbox_inches='tight')
            plt.close('all')

        # --- 3. Unclustered Heatmap ---
        plt.figure(figsize=(10, 6))
        sns.heatmap(h_df, cmap='magma', annot=True, fmt=".2e")
        plt.title('NMF Factor Weights (Raw Order)', fontsize=16)
        plt.savefig(os.path.join(args.outdir, "H_matrix_unclustered.png"), dpi=300, bbox_inches='tight')
        plt.close('all')
        
        print(f"Heatmaps successfully generated in {args.outdir}")

    except Exception as e:
        print(f"Error generating heatmaps: {e}")
        exit(1)

if __name__ == '__main__':
    main()
