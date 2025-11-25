#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regenerate volcano plots with improved formatting"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def create_volcano(de_result, title, filename):
    fig, ax = plt.subplots(figsize=(14, 10))
    
    x_limit = 5  # Limit x-axis to [-5, 5]
    
    # Handle pvals_adj = 0 by replacing with a small value
    de_result = de_result.copy()
    de_result['pvals_adj'] = de_result['pvals_adj'].replace(0, 1e-310)
    
    # Calculate -log10 p-value
    neg_log10_pval = -np.log10(de_result['pvals_adj'])
    y_max = neg_log10_pval[np.isfinite(neg_log10_pval)].max() * 1.1
    
    colors = np.where(
        (de_result['pvals_adj'] < 0.05) & (de_result['log2fc'] > 0.5), 'red',
        np.where((de_result['pvals_adj'] < 0.05) & (de_result['log2fc'] < -0.5), 'blue', 'lightgray')
    )
    
    ax.scatter(de_result['log2fc'], neg_log10_pval,
              c=colors, alpha=0.5, s=15)
    
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(-0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # Set axis limits
    ax.set_xlim(-x_limit, x_limit)
    ax.set_ylim(0, y_max)
    
    # Get significant genes
    sig_genes = de_result[(de_result['pvals_adj'] < 0.05) & (abs(de_result['log2fc']) > 0.5)].copy()
    sig_genes['neg_log10_pval'] = -np.log10(sig_genes['pvals_adj'])
    
    # Top 10 upregulated genes (most significant)
    top_up = sig_genes[sig_genes['log2fc'] > 0.5].nsmallest(10, 'pvals_adj')
    
    # Top 10 downregulated genes (most significant)
    top_down = sig_genes[sig_genes['log2fc'] < -0.5].nsmallest(10, 'pvals_adj')
    
    # Calculate y positions for labels to avoid overlap
    y_positions_up = np.linspace(y_max * 0.95, y_max * 0.3, len(top_up))
    y_positions_down = np.linspace(y_max * 0.95, y_max * 0.3, len(top_down))
    
    # Label upregulated genes (right side)
    for i, (_, row) in enumerate(top_up.iterrows()):
        x_pos = min(row['log2fc'], x_limit - 0.3)
        y_pos = -np.log10(row['pvals_adj'])
        y_text = y_positions_up[i]
        ax.annotate(row['names'],
                   xy=(x_pos, y_pos),
                   xytext=(x_limit + 0.5, y_text),
                   fontsize=10, color='darkred', fontweight='bold',
                   arrowprops=dict(arrowstyle='->', color='darkred', alpha=0.6, lw=1),
                   va='center', ha='left')
    
    # Label downregulated genes (left side)
    for i, (_, row) in enumerate(top_down.iterrows()):
        x_pos = max(row['log2fc'], -x_limit + 0.3)
        y_pos = -np.log10(row['pvals_adj'])
        y_text = y_positions_down[i]
        ax.annotate(row['names'],
                   xy=(x_pos, y_pos),
                   xytext=(-x_limit - 0.5, y_text),
                   fontsize=10, color='darkblue', fontweight='bold',
                   arrowprops=dict(arrowstyle='->', color='darkblue', alpha=0.6, lw=1),
                   va='center', ha='right')
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add legend
    n_up = len(sig_genes[sig_genes['log2fc'] > 0.5])
    n_down = len(sig_genes[sig_genes['log2fc'] < -0.5])
    legend_elements = [
        Patch(facecolor='red', label=f'Upregulated ({n_up})'),
        Patch(facecolor='blue', label=f'Downregulated ({n_down})'),
        Patch(facecolor='lightgray', label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.2, right=0.8)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {filename}')

if __name__ == '__main__':
    # Read DE results
    de_100pe = pd.read_csv('results/tables/DE_100PE_vs_control.csv')
    de_500pe = pd.read_csv('results/tables/DE_500PE_vs_control.csv')
    
    # Create volcano plots
    create_volcano(de_100pe, '100PE vs Control', 'results/figures/10_volcano_100PE.png')
    create_volcano(de_500pe, '500PE vs Control', 'results/figures/11_volcano_500PE.png')
    
    print('Done!')

