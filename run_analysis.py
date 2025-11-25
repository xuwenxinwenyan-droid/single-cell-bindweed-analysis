#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Single-Cell RNA-seq Analysis Pipeline
运行完整的单细胞RNA-seq分析流程
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("Single-Cell RNA-seq Analysis Pipeline")
print("="*70)

# Step 1: Set Data Path
print("\n[Step 1] Setting data path...")
possible_paths = ['./matrix/', 'matrix/', '../matrix/']

DATA_PATH = None
for path in possible_paths:
    if os.path.exists(path):
        samples_found = [d for d in os.listdir(path) 
                         if os.path.isdir(os.path.join(path, d)) 
                         and d in ['control', '100PE', '500PE']]
        if len(samples_found) >= 2:
            DATA_PATH = path
            print(f"[OK] Data found at: {DATA_PATH}")
            print(f"  Samples found: {sorted(samples_found)}")
            break

if DATA_PATH is None:
    print("[ERROR] Could not automatically find data path.")
    print(f"Current working directory: {os.getcwd()}")
    sys.exit(1)

# Step 2: Import packages
print("\n[Step 2] Importing packages...")
try:
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    print(f"[OK] Scanpy version: {sc.__version__}")
    print("[OK] All packages imported!")
except ImportError as e:
    print(f"[ERROR] Missing package: {e}")
    print("Please install: pip install scanpy python-igraph leidenalg numpy pandas matplotlib seaborn")
    sys.exit(1)

# Step 3: Configure Settings
print("\n[Step 3] Configuring settings...")
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

# Create output directory
os.makedirs('./results/figures', exist_ok=True)
os.makedirs('./results/tables', exist_ok=True)

sc.settings.figdir = './results/figures/'
print("[OK] Settings configured")
print("[OK] Output directories created")

# Step 4: Load Data
print("\n[Step 4] Loading data...")
samples = ['control', '100PE', '500PE']
adatas = {}

for sample in samples:
    print(f"\nLoading {sample}...")
    path = os.path.join(DATA_PATH, sample)
    
    if not os.path.exists(path):
        print(f"  [WARNING] Warning: {path} not found, skipping...")
        continue
    
    try:
        import gzip
        from scipy.io import mmread
        
        # Read matrix
        matrix_file = os.path.join(path, 'matrix.mtx.gz')
        features_file = os.path.join(path, 'features.tsv.gz')
        barcodes_file = os.path.join(path, 'barcodes.tsv.gz')
        
        # Read matrix
        with gzip.open(matrix_file, 'rb') as f:
            matrix = mmread(f).T.tocsr()  # Transpose: genes x cells -> cells x genes
        
        # Read features (gene names)
        with gzip.open(features_file, 'rt') as f:
            features = [line.strip().split('\t')[0] for line in f]
        
        # Read barcodes
        with gzip.open(barcodes_file, 'rt') as f:
            barcodes = [line.strip() for line in f]
        
        # Create AnnData object
        import anndata
        adata = anndata.AnnData(X=matrix)
        adata.obs_names = barcodes
        adata.var_names = features
        
        # Make var_names unique
        adata.var_names_make_unique()
        
        # Add metadata
        adata.obs['sample'] = sample
        adata.obs['treatment'] = sample if sample != 'control' else 'control'
        
        adatas[sample] = adata
        print(f"  [OK] {sample}: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    except Exception as e:
        import traceback
        print(f"  [ERROR] Error loading {sample}: {e}")
        traceback.print_exc()
        continue

if len(adatas) < 2:
    print("[ERROR] Need at least 2 samples to continue!")
    sys.exit(1)

# Concatenate
print("\nCombining samples...")
adata_list = [adatas[s] for s in samples if s in adatas]
adata = adata_list[0].concatenate(
    adata_list[1:],
    batch_key='sample',
    batch_categories=[s for s in samples if s in adatas]
)

print(f"\n[OK] Combined: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Save raw counts
adata.layers['counts'] = adata.X.copy()

# Step 5: Quality Control
print("\n[Step 5] Calculating QC metrics...")

# Identify gene types
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')
adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')

# Calculate metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo', 'hb'],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("\n[INFO] QC Summary:")
print(f"  Mean genes/cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"  Mean UMIs/cell: {adata.obs['total_counts'].mean():.0f}")
print(f"  Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")

# QC plots
print("  Generating QC plots...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Quality Control Metrics', fontsize=16, fontweight='bold')

sc.pl.violin(adata, 'total_counts', groupby='sample', ax=axes[0,0], show=False)
axes[0,0].set_title('Total Counts')

sc.pl.violin(adata, 'n_genes_by_counts', groupby='sample', ax=axes[0,1], show=False)
axes[0,1].set_title('N Genes')

sc.pl.violin(adata, 'pct_counts_mt', groupby='sample', ax=axes[0,2], show=False)
axes[0,2].set_title('MT %')

axes[1,0].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], alpha=0.3, s=1)
axes[1,0].set_xlabel('Total counts')
axes[1,0].set_ylabel('N genes')
axes[1,0].set_xscale('log')
axes[1,0].set_yscale('log')

axes[1,1].scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'], alpha=0.3, s=1)
axes[1,1].set_xlabel('Total counts')
axes[1,1].set_ylabel('MT %')
axes[1,1].set_xscale('log')

sample_counts = adata.obs['sample'].value_counts()
axes[1,2].bar(range(len(sample_counts)), sample_counts.values)
axes[1,2].set_xticks(range(len(sample_counts)))
axes[1,2].set_xticklabels(sample_counts.index)
axes[1,2].set_ylabel('Cell count')

plt.tight_layout()
plt.savefig('./results/figures/01_qc_metrics.png', dpi=300, bbox_inches='tight')
plt.close()
print("  [OK] QC plots saved")

# Step 6: Filtering
print("\n[Step 6] Filtering cells and genes...")
print(f"Before filtering: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Filter cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=500)

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)

# Remove high MT% cells
adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()

# Remove outliers
upper_count = np.percentile(adata.obs['total_counts'], 99)
adata = adata[adata.obs['total_counts'] < upper_count, :].copy()

print(f"After filtering: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"\nCells per sample:")
print(adata.obs['sample'].value_counts())

# Step 7: Normalization & Feature Selection
print("\n[Step 7] Normalizing and selecting features...")
print("Normalizing...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

print("Finding highly variable genes...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    batch_key='sample',
    flavor='seurat_v3'
)

print(f"  [OK] {sum(adata.var['highly_variable'])} HVGs selected")

# Plot
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig('./results/figures/02_hvgs.png', dpi=300, bbox_inches='tight')
plt.close()

# Step 8: Scaling & PCA
print("\n[Step 8] Scaling and running PCA...")
print("Scaling and regressing out confounders...")
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

print("Running PCA...")
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Variance ratio plot
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
plt.savefig('./results/figures/03_pca_variance.png', dpi=300, bbox_inches='tight')
plt.close()

# PCA by sample
sc.pl.pca(adata, color='sample', show=False)
plt.savefig('./results/figures/04_pca_sample.png', dpi=300, bbox_inches='tight')
plt.close()

print("[OK] PCA completed")

# Step 9: Neighbors & UMAP
print("\n[Step 9] Computing neighborhood graph and UMAP...")
print("Computing neighborhood graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

print("Computing UMAP...")
sc.tl.umap(adata)

# UMAP plots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sc.pl.umap(adata, color='sample', ax=axes[0], show=False, 
           title='UMAP by Sample', legend_loc='right margin')
sc.pl.umap(adata, color='treatment', ax=axes[1], show=False,
           title='UMAP by Treatment', legend_loc='right margin')

plt.tight_layout()
plt.savefig('./results/figures/05_umap_overview.png', dpi=300, bbox_inches='tight')
plt.close()

print("[OK] UMAP completed")

# Step 10: Clustering
print("\n[Step 10] Clustering...")
print("Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5)

n_clusters = len(adata.obs['leiden'].unique())
print(f"  [OK] Found {n_clusters} clusters")

# UMAP with clusters
sc.pl.umap(adata, color='leiden', legend_loc='on data', 
           title=f'UMAP by Cluster ({n_clusters} clusters)', show=False)
plt.savefig('./results/figures/06_umap_clusters.png', dpi=300, bbox_inches='tight')
plt.close()

# Step 11: Find Marker Genes
print("\n[Step 11] Finding marker genes...")
print("Finding marker genes for each cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Extract and save
marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
marker_genes.to_csv('./results/tables/cluster_markers_top20.csv', index=False)
print("  [OK] Saved: cluster_markers_top20.csv")

# Visualize
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
plt.savefig('./results/figures/07_markers.png', dpi=300, bbox_inches='tight')
plt.close()

# Dotplot
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
plt.savefig('./results/figures/08_markers_dotplot.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nTop 3 markers per cluster:")
print(marker_genes.head(3))

# Step 12: Cell Composition Analysis
print("\n[Step 12] Analyzing cell composition...")
composition = pd.crosstab(adata.obs['leiden'], adata.obs['sample'])
composition_pct = composition.div(composition.sum(axis=0), axis=1) * 100

composition.to_csv('./results/tables/composition_counts.csv')
composition_pct.to_csv('./results/tables/composition_percent.csv')

print("\nComposition (%):\n")
print(composition_pct.round(1))

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

composition.T.plot(kind='bar', stacked=True, ax=axes[0], colormap='tab20')
axes[0].set_xlabel('Sample')
axes[0].set_ylabel('Cell Count')
axes[0].set_title('Cell Composition (Counts)')
axes[0].legend(title='Cluster', bbox_to_anchor=(1.05, 1))

composition_pct.T.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab20')
axes[1].set_xlabel('Sample')
axes[1].set_ylabel('Percentage (%)')
axes[1].set_title('Cell Composition (Percentage)')
axes[1].legend(title='Cluster', bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.savefig('./results/figures/09_composition.png', dpi=300, bbox_inches='tight')
plt.close()

# Step 13: Differential Expression - 100PE vs Control
print("\n[Step 13] Differential Expression: 100PE vs Control")
print("="*70)

if '100PE' in adata.obs['sample'].values and 'control' in adata.obs['sample'].values:
    adata_100pe = adata[adata.obs['sample'].isin(['100PE', 'control'])].copy()
    
    sc.tl.rank_genes_groups(
        adata_100pe,
        groupby='sample',
        groups=['100PE'],
        reference='control',
        method='wilcoxon'
    )
    
    de_100pe = sc.get.rank_genes_groups_df(adata_100pe, group='100PE')
    de_100pe['log2fc'] = np.log2(np.exp(de_100pe['logfoldchanges']))
    de_100pe['significant'] = (de_100pe['pvals_adj'] < 0.05) & (abs(de_100pe['log2fc']) > 0.5)
    
    print(f"\nTotal DEGs: {de_100pe['significant'].sum()}")
    print(f"  Upregulated: {((de_100pe['log2fc'] > 0.5) & (de_100pe['pvals_adj'] < 0.05)).sum()}")
    print(f"  Downregulated: {((de_100pe['log2fc'] < -0.5) & (de_100pe['pvals_adj'] < 0.05)).sum()}")
    
    de_100pe.to_csv('./results/tables/DE_100PE_vs_control.csv', index=False)
    print("\n[OK] Saved: DE_100PE_vs_control.csv")
    
    print("\nTop 10 upregulated:")
    print(de_100pe[de_100pe['log2fc'] > 0].nsmallest(10, 'pvals_adj')[['names', 'log2fc', 'pvals_adj']])
else:
    print("[WARNING] Skipping: 100PE or control sample not found")
    de_100pe = None

# Step 14: Differential Expression - 500PE vs Control
print("\n[Step 14] Differential Expression: 500PE vs Control")
print("="*70)

if '500PE' in adata.obs['sample'].values and 'control' in adata.obs['sample'].values:
    adata_500pe = adata[adata.obs['sample'].isin(['500PE', 'control'])].copy()
    
    sc.tl.rank_genes_groups(
        adata_500pe,
        groupby='sample',
        groups=['500PE'],
        reference='control',
        method='wilcoxon'
    )
    
    de_500pe = sc.get.rank_genes_groups_df(adata_500pe, group='500PE')
    de_500pe['log2fc'] = np.log2(np.exp(de_500pe['logfoldchanges']))
    de_500pe['significant'] = (de_500pe['pvals_adj'] < 0.05) & (abs(de_500pe['log2fc']) > 0.5)
    
    print(f"\nTotal DEGs: {de_500pe['significant'].sum()}")
    print(f"  Upregulated: {((de_500pe['log2fc'] > 0.5) & (de_500pe['pvals_adj'] < 0.05)).sum()}")
    print(f"  Downregulated: {((de_500pe['log2fc'] < -0.5) & (de_500pe['pvals_adj'] < 0.05)).sum()}")
    
    de_500pe.to_csv('./results/tables/DE_500PE_vs_control.csv', index=False)
    print("\n[OK] Saved: DE_500PE_vs_control.csv")
    
    print("\nTop 10 upregulated:")
    print(de_500pe[de_500pe['log2fc'] > 0].nsmallest(10, 'pvals_adj')[['names', 'log2fc', 'pvals_adj']])
else:
    print("[WARNING] Skipping: 500PE or control sample not found")
    de_500pe = None

# Step 15: Volcano Plots
print("\n[Step 15] Creating volcano plots...")

def create_volcano(de_result, title, filename):
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Filter data for better visualization (limit x-axis range)
    x_limit = 5  # Set x-axis limit to [-5, 5]
    
    colors = np.where(
        (de_result['pvals_adj'] < 0.05) & (de_result['log2fc'] > 0.5), 'red',
        np.where((de_result['pvals_adj'] < 0.05) & (de_result['log2fc'] < -0.5), 'blue', 'lightgray')
    )
    
    ax.scatter(de_result['log2fc'], -np.log10(de_result['pvals_adj']),
              c=colors, alpha=0.5, s=15)
    
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(-0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # Set x-axis limits
    ax.set_xlim(-x_limit, x_limit)
    
    # Get significant genes
    sig_genes = de_result[(de_result['pvals_adj'] < 0.05) & (abs(de_result['log2fc']) > 0.5)].copy()
    
    # Top 10 upregulated genes (most significant, log2fc > 0.5)
    top_up = sig_genes[sig_genes['log2fc'] > 0.5].nsmallest(10, 'pvals_adj')
    
    # Top 10 downregulated genes (most significant, log2fc < -0.5)
    top_down = sig_genes[sig_genes['log2fc'] < -0.5].nsmallest(10, 'pvals_adj')
    
    # Label upregulated genes (right side)
    for i, (_, row) in enumerate(top_up.iterrows()):
        x_pos = min(row['log2fc'], x_limit - 0.3)
        y_pos = -np.log10(row['pvals_adj'])
        ax.annotate(row['names'],
                   xy=(x_pos, y_pos),
                   xytext=(x_limit + 0.3, y_pos),
                   fontsize=9, color='red', fontweight='bold',
                   arrowprops=dict(arrowstyle='-', color='red', alpha=0.5, lw=0.5),
                   va='center', ha='left')
    
    # Label downregulated genes (left side)
    for i, (_, row) in enumerate(top_down.iterrows()):
        x_pos = max(row['log2fc'], -x_limit + 0.3)
        y_pos = -np.log10(row['pvals_adj'])
        ax.annotate(row['names'],
                   xy=(x_pos, y_pos),
                   xytext=(-x_limit - 0.3, y_pos),
                   fontsize=9, color='blue', fontweight='bold',
                   arrowprops=dict(arrowstyle='-', color='blue', alpha=0.5, lw=0.5),
                   va='center', ha='right')
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label=f'Upregulated ({len(sig_genes[sig_genes["log2fc"] > 0.5])})'),
        Patch(facecolor='blue', label=f'Downregulated ({len(sig_genes[sig_genes["log2fc"] < -0.5])})'),
        Patch(facecolor='lightgray', label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.85)  # Make room for labels
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

if de_100pe is not None:
    create_volcano(de_100pe, '100PE vs Control', 
                   './results/figures/10_volcano_100PE.png')
if de_500pe is not None:
    create_volcano(de_500pe, '500PE vs Control',
                   './results/figures/11_volcano_500PE.png')

print("[OK] Volcano plots saved")

# Step 15.5: ER Stress Related DEGs Analysis
print("\n[Step 15.5] Analyzing ER stress related differential genes...")
print("="*70)

# Define ER stress related genes
# Classic ER stress key genes
er_stress_genes_classic = [
    # PERK pathway
    'EIF2AK3', 'PERK',  # PERK
    'EIF2S1', 'EIF2A',  # eIF2α
    'ATF4',  # ATF4
    'DDIT3', 'CHOP',  # CHOP
    # IRE1 pathway
    'ERN1', 'IRE1',  # IRE1
    'XBP1',  # XBP1
    # ATF6 pathway
    'ATF6',  # ATF6
    # Chaperones
    'HSPA5', 'GRP78', 'BIP',  # GRP78/BiP
    'HSP90B1', 'GRP94',  # GRP94
    'DNAJB11', 'ERDJ3',  # ERdj3
    'PDIA3', 'P58IPK',  # PDI family
    'CANX',  # Calnexin
    'CALR',  # Calreticulin
    # ERAD pathway
    'EDEM1', 'EDEM2', 'EDEM3',  # EDEM family
    'DERL1', 'DERL2', 'DERL3',  # Derlin family
    'VCP',  # VCP/p97
    # Other ER stress markers
    'HYOU1',  # Hypoxia up-regulated 1
    'DNAJC3', 'P58IPK',  # P58IPK
    'TRAF2',  # TRAF2
    'ERN2', 'IRE1B',  # IRE1β
]

# Try to get Hallmark UPR genes if available (common genes from MSigDB)
er_stress_genes_hallmark = [
    'HSPA5', 'DNAJB11', 'HYOU1', 'PDIA3', 'PDIA4', 'PDIA6',
    'CANX', 'CALR', 'DERL1', 'EDEM1', 'VCP', 'ATF4', 'DDIT3',
    'XBP1', 'ERN1', 'ATF6', 'EIF2AK3', 'EIF2S1', 'ERN2',
    'DNAJC3', 'TRAF2', 'DERL2', 'DERL3', 'EDEM2', 'EDEM3',
    'HSP90B1', 'SEC61A1', 'SEC61A2', 'SEC61B', 'SEC61G',
    'RPN1', 'RPN2', 'OST48', 'DAD1', 'STT3A', 'STT3B',
    'MOGS', 'UGGT1', 'UGGT2', 'MAN1B1', 'MAN1C1'
]

# Combine and deduplicate
er_stress_genes = list(set(er_stress_genes_classic + er_stress_genes_hallmark))
print(f"Defined {len(er_stress_genes)} ER stress related genes")

def filter_er_stress_degs(de_result, er_genes, comparison_name):
    """Filter ER stress related genes from DEG results"""
    if de_result is None or len(de_result) == 0:
        return None
    
    # Convert gene names to uppercase for matching
    de_genes_upper = de_result['names'].str.upper()
    er_genes_upper = [g.upper() for g in er_genes]
    
    # Find matches
    matches = de_genes_upper.isin(er_genes_upper)
    er_stress_degs = de_result[matches].copy()
    
    if len(er_stress_degs) > 0:
        print(f"\n{comparison_name}:")
        print(f"  Found {len(er_stress_degs)} ER stress related DEGs")
        print(f"  Upregulated: {((er_stress_degs['log2fc'] > 0.5) & (er_stress_degs['pvals_adj'] < 0.05)).sum()}")
        print(f"  Downregulated: {((er_stress_degs['log2fc'] < -0.5) & (er_stress_degs['pvals_adj'] < 0.05)).sum()}")
        
        # Show top ER stress genes
        sig_er = er_stress_degs[er_stress_degs['significant']]
        if len(sig_er) > 0:
            print(f"\n  Top significant ER stress DEGs:")
            for idx, row in sig_er.nsmallest(10, 'pvals_adj').iterrows():
                direction = "UP" if row['log2fc'] > 0 else "DOWN"
                print(f"    {direction} {row['names']}: log2FC={row['log2fc']:.2f}, p_adj={row['pvals_adj']:.2e}")
    else:
        print(f"\n{comparison_name}:")
        print(f"  No ER stress related DEGs found")
    
    return er_stress_degs

# Filter ER stress DEGs
er_stress_100pe = None
er_stress_500pe = None

if de_100pe is not None:
    er_stress_100pe = filter_er_stress_degs(de_100pe, er_stress_genes, "100PE vs Control")
    if er_stress_100pe is not None and len(er_stress_100pe) > 0:
        er_stress_100pe.to_csv('./results/tables/ER_stress_DEGs_100PE_vs_control.csv', index=False)
        print("  [OK] Saved: ER_stress_DEGs_100PE_vs_control.csv")

if de_500pe is not None:
    er_stress_500pe = filter_er_stress_degs(de_500pe, er_stress_genes, "500PE vs Control")
    if er_stress_500pe is not None and len(er_stress_500pe) > 0:
        er_stress_500pe.to_csv('./results/tables/ER_stress_DEGs_500PE_vs_control.csv', index=False)
        print("  [OK] Saved: ER_stress_DEGs_500PE_vs_control.csv")

# Create summary table
if (er_stress_100pe is not None and len(er_stress_100pe) > 0) or \
   (er_stress_500pe is not None and len(er_stress_500pe) > 0):
    
    summary_data = []
    
    if er_stress_100pe is not None and len(er_stress_100pe) > 0:
        for _, row in er_stress_100pe.iterrows():
            summary_data.append({
                'gene': row['names'],
                'comparison': '100PE vs Control',
                'log2fc': row['log2fc'],
                'pvals_adj': row['pvals_adj'],
                'significant': row['significant'],
                'direction': 'Up' if row['log2fc'] > 0 else 'Down'
            })
    
    if er_stress_500pe is not None and len(er_stress_500pe) > 0:
        for _, row in er_stress_500pe.iterrows():
            summary_data.append({
                'gene': row['names'],
                'comparison': '500PE vs Control',
                'log2fc': row['log2fc'],
                'pvals_adj': row['pvals_adj'],
                'significant': row['significant'],
                'direction': 'Up' if row['log2fc'] > 0 else 'Down'
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv('./results/tables/ER_stress_DEGs_summary.csv', index=False)
        print("\n[OK] Saved: ER_stress_DEGs_summary.csv")

# Create visualizations
if (er_stress_100pe is not None and len(er_stress_100pe) > 0) or \
   (er_stress_500pe is not None and len(er_stress_500pe) > 0):
    
    # Bar plot: Count of ER stress DEGs
    fig, ax = plt.subplots(figsize=(8, 6))
    
    comparisons = []
    counts_up = []
    counts_down = []
    
    if er_stress_100pe is not None and len(er_stress_100pe) > 0:
        sig_100pe = er_stress_100pe[er_stress_100pe['significant']]
        comparisons.append('100PE\nvs Control')
        counts_up.append(((sig_100pe['log2fc'] > 0.5) & (sig_100pe['pvals_adj'] < 0.05)).sum())
        counts_down.append(((sig_100pe['log2fc'] < -0.5) & (sig_100pe['pvals_adj'] < 0.05)).sum())
    
    if er_stress_500pe is not None and len(er_stress_500pe) > 0:
        sig_500pe = er_stress_500pe[er_stress_500pe['significant']]
        comparisons.append('500PE\nvs Control')
        counts_up.append(((sig_500pe['log2fc'] > 0.5) & (sig_500pe['pvals_adj'] < 0.05)).sum())
        counts_down.append(((sig_500pe['log2fc'] < -0.5) & (sig_500pe['pvals_adj'] < 0.05)).sum())
    
    if comparisons:
        x = np.arange(len(comparisons))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, counts_up, width, label='Upregulated', color='#d62728', edgecolor='black')
        bars2 = ax.bar(x + width/2, counts_down, width, label='Downregulated', color='#2ca02c', edgecolor='black')
        
        # Add value labels on bars
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{int(height)}', ha='center', va='bottom',
                           fontsize=11, fontweight='bold')
        
        ax.set_xlabel('Comparison', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of ER Stress DEGs', fontsize=12, fontweight='bold')
        ax.set_title('ER Stress Related Differential Genes', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(comparisons)
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('./results/figures/ER_stress_barplot.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("[OK] Saved: ER_stress_barplot.png")
    
    # Heatmap: Expression of ER stress genes across samples
    # Get all ER stress genes found in data
    all_er_genes_found = set()
    if er_stress_100pe is not None and len(er_stress_100pe) > 0:
        all_er_genes_found.update(er_stress_100pe['names'].tolist())
    if er_stress_500pe is not None and len(er_stress_500pe) > 0:
        all_er_genes_found.update(er_stress_500pe['names'].tolist())
    
    if all_er_genes_found:
        # Filter genes that exist in adata
        er_genes_in_data = [g for g in all_er_genes_found if g in adata.var_names]
        
        if er_genes_in_data:
            # Calculate mean expression per sample
            er_adata = adata[:, er_genes_in_data].copy()
            
            # Use raw counts if available, otherwise use normalized
            exp_data = []
            for sample in ['control', '100PE', '500PE']:
                if sample in er_adata.obs['sample'].values:
                    sample_mask = er_adata.obs['sample'] == sample
                    if hasattr(adata, 'raw') and adata.raw is not None:
                        sample_exp = adata.raw.X[sample_mask].mean(axis=0)
                    else:
                        sample_exp = er_adata.X[sample_mask].mean(axis=0)
                    if hasattr(sample_exp, 'A1'):
                        exp_data.append(sample_exp.A1)
                    else:
                        exp_data.append(sample_exp)
            
            if exp_data:
                exp_df = pd.DataFrame(
                    np.array(exp_data).T,
                    index=er_genes_in_data,
                    columns=['control', '100PE', '500PE']
                )
                
                # Log transform if needed
                if exp_df.max().max() > 20:
                    exp_df = np.log1p(exp_df)
                
                # Create heatmap
                fig, ax = plt.subplots(figsize=(10, max(8, len(er_genes_in_data) * 0.3)))
                sns.heatmap(exp_df, annot=False, fmt='.2f', cmap='RdYlBu_r',
                           center=0, vmin=-2, vmax=2, ax=ax,
                           cbar_kws={'label': 'Log Normalized Expression'})
                ax.set_title('ER Stress Genes Expression Across Samples', 
                           fontsize=14, fontweight='bold', pad=20)
                ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
                ax.set_ylabel('Gene', fontsize=12, fontweight='bold')
                
                plt.tight_layout()
                plt.savefig('./results/figures/ER_stress_heatmap.png', dpi=300, bbox_inches='tight')
                plt.close()
                print("[OK] Saved: ER_stress_heatmap.png")
    
    # Enhanced volcano plots with ER stress genes highlighted
    def create_volcano_with_er_stress(de_result, er_stress_degs, title, filename):
        if de_result is None or len(de_result) == 0:
            return
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Get ER stress gene names
        er_genes_set = set()
        if er_stress_degs is not None and len(er_stress_degs) > 0:
            er_genes_set = set(er_stress_degs['names'].str.upper())
        
        # Color coding
        colors = []
        sizes = []
        for idx, row in de_result.iterrows():
            gene_upper = row['names'].upper()
            is_sig = (row['pvals_adj'] < 0.05) & (abs(row['log2fc']) > 0.5)
            is_er = gene_upper in er_genes_set
            
            if is_er and is_sig:
                colors.append('orange')  # ER stress significant
                sizes.append(50)
            elif is_er:
                colors.append('yellow')  # ER stress but not significant
                sizes.append(30)
            elif is_sig and row['log2fc'] > 0.5:
                colors.append('red')
                sizes.append(15)
            elif is_sig and row['log2fc'] < -0.5:
                colors.append('blue')
                sizes.append(15)
            else:
                colors.append('gray')
                sizes.append(15)
        
        ax.scatter(de_result['log2fc'], -np.log10(de_result['pvals_adj']),
                  c=colors, s=sizes, alpha=0.6, edgecolors='black', linewidth=0.5)
        
        ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvline(0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvline(-0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
        
        # Label ER stress genes
        if er_stress_degs is not None and len(er_stress_degs) > 0:
            sig_er = er_stress_degs[er_stress_degs['significant']]
            for _, row in sig_er.iterrows():
                ax.annotate(row['names'],
                           xy=(row['log2fc'], -np.log10(row['pvals_adj'])),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=9, alpha=0.9, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
        
        ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
        ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12, fontweight='bold')
        ax.set_title(f'{title}\n(ER stress genes highlighted)', fontsize=14, fontweight='bold')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='orange', edgecolor='black', label='ER stress (significant)'),
            Patch(facecolor='yellow', edgecolor='black', label='ER stress (not sig)'),
            Patch(facecolor='red', edgecolor='black', label='Upregulated'),
            Patch(facecolor='blue', edgecolor='black', label='Downregulated'),
            Patch(facecolor='gray', edgecolor='black', label='Not significant')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    
    if de_100pe is not None:
        create_volcano_with_er_stress(de_100pe, er_stress_100pe, '100PE vs Control',
                                     './results/figures/10_volcano_100PE_ER_stress.png')
    if de_500pe is not None:
        create_volcano_with_er_stress(de_500pe, er_stress_500pe, '500PE vs Control',
                                     './results/figures/11_volcano_500PE_ER_stress.png')

print("\n[OK] ER stress analysis completed")

# Step 16: DEG Overlap Analysis
print("\n[Step 16] Analyzing DEG overlap...")

if de_100pe is not None and de_500pe is not None:
    sig_100pe = set(de_100pe[de_100pe['significant']]['names'])
    sig_500pe = set(de_500pe[de_500pe['significant']]['names'])
    
    common = sig_100pe & sig_500pe
    unique_100pe = sig_100pe - sig_500pe
    unique_500pe = sig_500pe - sig_100pe
    
    print(f"\nDEG Overlap:")
    print(f"  Common: {len(common)}")
    print(f"  Unique to 100PE: {len(unique_100pe)}")
    print(f"  Unique to 500PE: {len(unique_500pe)}")
    
    # Save
    pd.DataFrame({'gene': list(common)}).to_csv('./results/tables/common_DEGs.csv', index=False)
    pd.DataFrame({'gene': list(unique_100pe)}).to_csv('./results/tables/unique_100PE.csv', index=False)
    pd.DataFrame({'gene': list(unique_500pe)}).to_csv('./results/tables/unique_500PE.csv', index=False)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    categories = ['Common', 'Unique\n100PE', 'Unique\n500PE']
    counts = [len(common), len(unique_100pe), len(unique_500pe)]
    colors = ['purple', 'coral', 'mediumseagreen']
    
    bars = ax.bar(categories, counts, color=colors, edgecolor='black', linewidth=1.5)
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom',
                fontsize=14, fontweight='bold')
    
    ax.set_ylabel('Number of DEGs', fontsize=12, fontweight='bold')
    ax.set_title('DEG Overlap', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('./results/figures/12_overlap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("[OK] Overlap analysis saved")
else:
    print("[WARNING] Skipping: Need both DE results for overlap analysis")
    common = set()
    unique_100pe = set()
    unique_500pe = set()

# Step 17: Dose-Response Analysis
print("\n[Step 17] Dose-response analysis...")

if de_100pe is not None and de_500pe is not None:
    merged = pd.merge(
        de_100pe[['names', 'log2fc', 'pvals_adj']],
        de_500pe[['names', 'log2fc', 'pvals_adj']],
        on='names',
        suffixes=('_100PE', '_500PE')
    )
    
    corr = merged['log2fc_100PE'].corr(merged['log2fc_500PE'])
    print(f"\nDose-response correlation: {corr:.3f}")
    
    merged.to_csv('./results/tables/dose_response.csv', index=False)
    
    # Plot
    fig, ax = plt.subplots(figsize=(9, 9))
    
    ax.scatter(merged['log2fc_500PE'], merged['log2fc_100PE'],
              alpha=0.3, s=10, c='lightgray')
    
    sig_both = merged[
        (merged['pvals_adj_100PE'] < 0.05) & (merged['pvals_adj_500PE'] < 0.05)
    ]
    ax.scatter(sig_both['log2fc_500PE'], sig_both['log2fc_100PE'],
              alpha=0.7, s=30, c='red', edgecolor='darkred',
              label=f'Significant in both (n={len(sig_both)})')
    
    lim = max(abs(merged['log2fc_100PE']).max(), abs(merged['log2fc_500PE']).max())
    ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.5, linewidth=1)
    
    ax.axhline(0, color='black', linewidth=0.8, alpha=0.3)
    ax.axvline(0, color='black', linewidth=0.8, alpha=0.3)
    
    ax.set_xlabel('Log2FC (500PE vs Control)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Log2FC (100PE vs Control)', fontsize=12, fontweight='bold')
    ax.set_title(f'Dose-Response Correlation (r={corr:.3f})',
                fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.2)
    
    plt.tight_layout()
    plt.savefig('./results/figures/13_dose_response.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("[OK] Dose-response analysis saved")
else:
    print("[WARNING] Skipping: Need both DE results for dose-response analysis")
    corr = 0.0

# Step 18: Save Processed Data
print("\n[Step 18] Saving processed data...")

adata.write('./results/processed_data.h5ad', compression='gzip')
adata.obs.to_csv('./results/tables/cell_metadata.csv')
adata.var.to_csv('./results/tables/gene_metadata.csv')

print("[OK] All data saved")

# Step 19: Generate Summary Report
print("\n[Step 19] Generating summary report...")

summary = f"""
{'='*70}
ANALYSIS COMPLETE!
{'='*70}

Dataset:
  Total cells analyzed: {adata.shape[0]:,}
  Total genes: {adata.shape[1]:,}
  
Samples:
  control: {(adata.obs['sample'] == 'control').sum():,} cells
  100PE:   {(adata.obs['sample'] == '100PE').sum():,} cells
  500PE:   {(adata.obs['sample'] == '500PE').sum():,} cells

Clustering:
  Clusters identified: {len(adata.obs['leiden'].unique())}

Differential Expression:
"""

if de_100pe is not None:
    summary += f"""  100PE vs control: {de_100pe['significant'].sum()} DEGs
    - Upregulated: {((de_100pe['log2fc'] > 0.5) & (de_100pe['pvals_adj'] < 0.05)).sum()}
    - Downregulated: {((de_100pe['log2fc'] < -0.5) & (de_100pe['pvals_adj'] < 0.05)).sum()}
  
"""

if de_500pe is not None:
    summary += f"""  500PE vs control: {de_500pe['significant'].sum()} DEGs
    - Upregulated: {((de_500pe['log2fc'] > 0.5) & (de_500pe['pvals_adj'] < 0.05)).sum()}
    - Downregulated: {((de_500pe['log2fc'] < -0.5) & (de_500pe['pvals_adj'] < 0.05)).sum()}

"""

if de_100pe is not None and de_500pe is not None:
    summary += f"""DEG Overlap:
  Common: {len(common)}
  Unique to 100PE: {len(unique_100pe)}
  Unique to 500PE: {len(unique_500pe)}

Dose-Response:
  Correlation: {corr:.3f}

"""

summary += f"""Output Files:
  Figures: 13 plots saved
  Tables: 10+ CSV files
  Processed data: processed_data.h5ad

All results saved in: ./results/
{'='*70}
"""

print(summary)

with open('./results/ANALYSIS_SUMMARY.txt', 'w', encoding='utf-8') as f:
    f.write(summary)

print("\n[OK] Summary saved")
print("\n" + "="*70)
print("[SUCCESS] Analysis Complete!")
print("="*70)
print("\nAll results are saved in: ./results/")
print("  - Figures: ./results/figures/")
print("  - Tables: ./results/tables/")
print("  - Summary: ./results/ANALYSIS_SUMMARY.txt")
print("  - Processed data: ./results/processed_data.h5ad")

