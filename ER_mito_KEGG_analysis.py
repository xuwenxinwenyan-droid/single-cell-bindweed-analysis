#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ER Stress-Induced Mitochondrial Damage Analysis
内质网压力引起线粒体损伤的差异基因分析 + KEGG通路富集 + 基因网络

分析内容:
1. 筛选ER stress和线粒体损伤相关的差异基因
2. KEGG通路富集分析
3. 基因-通路网络可视化
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("ER Stress-Induced Mitochondrial Damage Analysis")
print("内质网压力引起线粒体损伤 - 差异基因 & KEGG富集分析")
print("="*70)

# ============================================================
# Step 1: 检查并安装必要的包
# ============================================================
print("\n[Step 1] 检查依赖包...")

def install_package(package):
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])

required_packages = {
    'scanpy': 'scanpy',
    'gseapy': 'gseapy',
    'networkx': 'networkx',
    'numpy': 'numpy',
    'pandas': 'pandas',
    'matplotlib': 'matplotlib',
    'seaborn': 'seaborn',
}

for import_name, pip_name in required_packages.items():
    try:
        __import__(import_name)
    except ImportError:
        print(f"  安装 {pip_name}...")
        install_package(pip_name)

# 导入包
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import networkx as nx
from matplotlib.patches import Patch
from collections import defaultdict

print("[OK] 所有包已导入")

# ============================================================
# Step 2: 设置数据路径和输出目录
# ============================================================
print("\n[Step 2] 设置路径...")

DATA_PATH = './matrix/'
if not os.path.exists(DATA_PATH):
    DATA_PATH = 'matrix/'

# 创建输出目录
os.makedirs('./results/figures', exist_ok=True)
os.makedirs('./results/tables', exist_ok=True)
os.makedirs('./results/KEGG', exist_ok=True)
os.makedirs('./results/networks', exist_ok=True)

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)
sc.settings.figdir = './results/figures/'

print(f"[OK] 数据路径: {DATA_PATH}")
print("[OK] 输出目录已创建")

# ============================================================
# Step 3: 定义ER stress和线粒体损伤相关基因
# ============================================================
print("\n[Step 3] 定义目标基因集...")

# ER stress / UPR (Unfolded Protein Response) 相关基因
ER_STRESS_GENES = [
    # PERK通路 (EIF2AK3 pathway)
    'EIF2AK3', 'EIF2S1', 'EIF2A', 'ATF4', 'DDIT3', 'CHOP', 'GADD34', 'PPP1R15A',
    'ASNS', 'TRIB3', 'WARS', 'CARS', 'GARS', 'YARS', 'SLC7A11', 'CTH',
    
    # IRE1通路 (ERN1 pathway)
    'ERN1', 'XBP1', 'DNAJB9', 'ERDJ4', 'SEC24D', 'BLOC1S1',
    
    # ATF6通路
    'ATF6', 'ATF6B', 'MBTPS1', 'MBTPS2', 'CREBZF',
    
    # ER伴侣蛋白/分子伴侣 (Chaperones)
    'HSPA5', 'BIP', 'GRP78', 'HSP90B1', 'GRP94', 'HYOU1', 'GRP170',
    'DNAJB11', 'DNAJC3', 'P58IPK', 'DNAJC1', 'ERDJ5', 'PPIB',
    
    # 蛋白二硫键异构酶家族 (PDI family)
    'PDIA3', 'PDIA4', 'PDIA6', 'P4HB', 'ERP29', 'ERP44', 'ERP57',
    
    # 钙网蛋白/钙联蛋白 (Calnexin/Calreticulin)
    'CANX', 'CALR',
    
    # ERAD (ER-associated degradation)
    'EDEM1', 'EDEM2', 'EDEM3', 'DERL1', 'DERL2', 'DERL3',
    'VCP', 'VIMP', 'SEL1L', 'HRD1', 'SYVN1', 'UBXN4',
    'OS9', 'ERLEC1', 'HERPUD1', 'HERP',
    
    # Sec61复合物
    'SEC61A1', 'SEC61A2', 'SEC61B', 'SEC61G',
    
    # N-糖基化相关
    'RPN1', 'RPN2', 'DAD1', 'STT3A', 'STT3B', 'DDOST', 'MOGS', 'GANAB',
    'UGGT1', 'UGGT2', 'MAN1B1',
    
    # ER stress诱导凋亡
    'TRAF2', 'ASK1', 'MAP3K5', 'JNK', 'MAPK8', 'MAPK9', 'MAPK10',
    'BCL2', 'BAX', 'BAK1', 'BIM', 'BCL2L11', 'PUMA', 'BBC3',
    'CASP12', 'CASP4', 'CASP3', 'CASP9',
]

# 线粒体功能/损伤相关基因
MITOCHONDRIA_GENES = [
    # 线粒体呼吸链复合物 (Respiratory chain complexes)
    # Complex I (NADH dehydrogenase)
    'NDUFV1', 'NDUFV2', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5',
    'NDUFS6', 'NDUFS7', 'NDUFS8', 'NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFA4',
    'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFA10', 'NDUFA11',
    'NDUFA12', 'NDUFA13', 'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5',
    'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFB10', 'NDUFB11', 'NDUFC1', 'NDUFC2',
    
    # Complex II (Succinate dehydrogenase)
    'SDHA', 'SDHB', 'SDHC', 'SDHD',
    
    # Complex III (Cytochrome bc1)
    'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRB', 'UQCRQ', 'UQCRH', 'UQCR10', 'UQCR11',
    'CYC1', 'CYTB',
    
    # Complex IV (Cytochrome c oxidase)
    'COX4I1', 'COX4I2', 'COX5A', 'COX5B', 'COX6A1', 'COX6A2', 'COX6B1', 'COX6B2',
    'COX6C', 'COX7A1', 'COX7A2', 'COX7B', 'COX7C', 'COX8A', 'MT-CO1', 'MT-CO2', 'MT-CO3',
    
    # Complex V (ATP synthase)
    'ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5F1E',
    'ATP5PB', 'ATP5MC1', 'ATP5MC2', 'ATP5MC3', 'ATP5PD', 'ATP5ME', 'ATP5MF',
    'ATP5PF', 'ATP5MG', 'ATP5IF1', 'MT-ATP6', 'MT-ATP8',
    
    # 线粒体DNA编码基因
    'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-ND6',
    'MT-CYB',
    
    # 线粒体动力学 (Mitochondrial dynamics)
    'MFN1', 'MFN2', 'OPA1', 'DNM1L', 'DRP1', 'FIS1', 'MFF', 'MIEF1', 'MIEF2',
    
    # 线粒体自噬 (Mitophagy)
    'PINK1', 'PRKN', 'PARKIN', 'BNIP3', 'BNIP3L', 'NIX', 'FUNDC1', 'BCL2L13',
    'OPTN', 'NDP52', 'CALCOCO2', 'SQSTM1', 'P62',
    
    # 线粒体膜蛋白/转运
    'TOMM20', 'TOMM22', 'TOMM40', 'TOMM70', 'TIMM17A', 'TIMM17B', 'TIMM23',
    'TIMM44', 'TIMM50', 'VDAC1', 'VDAC2', 'VDAC3', 'SLC25A4', 'ANT1', 'SLC25A5',
    
    # 线粒体应激/UPRmt
    'CLPP', 'LONP1', 'HSPD1', 'HSP60', 'HSPE1', 'HSP10', 'HSPA9', 'mtHSP70',
    'SIRT3', 'SIRT4', 'SIRT5', 'ATF5', 'CHOP', 'C/EBPβ',
    
    # 线粒体凋亡通路
    'CYCS', 'DIABLO', 'SMAC', 'HTRA2', 'OMI', 'ENDOG', 'AIF', 'AIFM1',
    'APAF1', 'CASP9', 'CASP3', 'CASP7',
    
    # 线粒体ROS相关
    'SOD1', 'SOD2', 'CAT', 'GPX1', 'GPX4', 'PRDX3', 'PRDX5', 'TXN2', 'TXNRD2',
    
    # 线粒体生物发生
    'PPARGC1A', 'PGC1A', 'PPARGC1B', 'NRF1', 'GABPA', 'TFAM', 'TFB1M', 'TFB2M',
    'POLG', 'POLG2', 'TWNK', 'MTIF2', 'MTIF3',
    
    # 线粒体钙信号
    'MCU', 'MICU1', 'MICU2', 'MCUR1', 'LETM1',
]

# ER-线粒体交互相关基因 (MAM - Mitochondria-Associated ER Membranes)
ER_MITO_INTERACTION_GENES = [
    # MAM结构蛋白
    'MFN2', 'VAPB', 'PTPIP51', 'RMDN3', 'GRP75', 'HSPA9', 'IP3R', 'ITPR1', 'ITPR2', 'ITPR3',
    'VDAC1', 'VDAC2', 'VDAC3', 'SIGMAR1', 'PACS2',
    
    # 钙信号传导
    'MCU', 'MICU1', 'MICU2', 'SERCA', 'ATP2A1', 'ATP2A2', 'ATP2A3',
    
    # MAM相关脂质代谢
    'ACAT1', 'FACL4', 'ACSL4', 'PSS1', 'PTDSS1', 'PSS2', 'PTDSS2',
    
    # ER-线粒体Ca2+传递
    'CACNA1A', 'CACNA1B', 'RYR1', 'RYR2', 'RYR3',
]

# 合并所有目标基因并去重
ALL_TARGET_GENES = list(set(ER_STRESS_GENES + MITOCHONDRIA_GENES + ER_MITO_INTERACTION_GENES))
print(f"[OK] 总共定义了 {len(ALL_TARGET_GENES)} 个目标基因")
print(f"    - ER stress基因: {len(ER_STRESS_GENES)}")
print(f"    - 线粒体相关基因: {len(MITOCHONDRIA_GENES)}")
print(f"    - ER-线粒体交互基因: {len(ER_MITO_INTERACTION_GENES)}")

# ============================================================
# Step 4: 加载数据
# ============================================================
print("\n[Step 4] 加载数据...")

import gzip
from scipy.io import mmread
import anndata as ad

def load_custom_10x_mtx(path, sample_name):
    """加载自定义格式的10x mtx数据（features只有一列）"""
    # 读取features (gene names)
    features_file = os.path.join(path, 'features.tsv.gz')
    with gzip.open(features_file, 'rt') as f:
        genes = [line.strip() for line in f.readlines()]
    
    # 读取barcodes (cell names)
    barcodes_file = os.path.join(path, 'barcodes.tsv.gz')
    with gzip.open(barcodes_file, 'rt') as f:
        cells = [line.strip() for line in f.readlines()]
    
    # 读取matrix
    matrix_file = os.path.join(path, 'matrix.mtx.gz')
    with gzip.open(matrix_file, 'rb') as f:
        matrix = mmread(f).T.tocsr()  # 转置并转为CSR格式
    
    # 创建AnnData对象
    adata = ad.AnnData(X=matrix)
    adata.var_names = genes
    adata.obs_names = [f"{sample_name}_{c}" for c in cells]
    
    # 处理重复基因名
    adata.var_names_make_unique()
    
    return adata

samples = ['control', '100PE', '500PE']
adatas = {}

for sample in samples:
    print(f"  加载 {sample}...")
    path = os.path.join(DATA_PATH, sample)
    
    if not os.path.exists(path):
        print(f"    [警告] {path} 不存在，跳过...")
        continue
    
    try:
        adata = load_custom_10x_mtx(path, sample)
        adata.obs['sample'] = sample
        adata.obs['treatment'] = sample if sample != 'control' else 'control'
        adatas[sample] = adata
        print(f"    [OK] {sample}: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    except Exception as e:
        print(f"    [错误] 加载 {sample} 失败: {e}")
        import traceback
        traceback.print_exc()
        continue

if len(adatas) < 2:
    print("[错误] 需要至少2个样本！")
    sys.exit(1)

# 合并数据
print("\n  合并样本...")
adata_list = [adatas[s] for s in samples if s in adatas]
adata = adata_list[0].concatenate(
    adata_list[1:],
    batch_key='sample',
    batch_categories=[s for s in samples if s in adatas]
)

print(f"[OK] 合并后: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# 保存原始counts
adata.layers['counts'] = adata.X.copy()

# ============================================================
# Step 5: 质控和预处理
# ============================================================
print("\n[Step 5] 质控和预处理...")

# 计算QC指标
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

print(f"  过滤前: {adata.shape[0]:,} cells")

# 过滤
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()

# 移除异常值
upper_count = np.percentile(adata.obs['total_counts'], 99)
adata = adata[adata.obs['total_counts'] < upper_count, :].copy()

print(f"  过滤后: {adata.shape[0]:,} cells")
print("\n  各样本细胞数:")
print(adata.obs['sample'].value_counts())

# 标准化
print("\n  标准化处理...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# HVG选择
print("  选择高变基因...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='sample', flavor='seurat_v3')

# 只保留高变基因用于降维（但保留raw数据用于差异分析）
adata_hvg = adata[:, adata.var['highly_variable']].copy()

# 降维
print("  PCA和UMAP (仅使用HVGs)...")
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5)

# 将降维结果复制回原始adata
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.obs['leiden'] = adata_hvg.obs['leiden']
adata.obsp['distances'] = adata_hvg.obsp['distances']
adata.obsp['connectivities'] = adata_hvg.obsp['connectivities']
adata.uns['neighbors'] = adata_hvg.uns['neighbors']
adata.uns['umap'] = adata_hvg.uns['umap']
adata.uns['leiden'] = adata_hvg.uns['leiden']
adata.uns['pca'] = adata_hvg.uns['pca']

del adata_hvg  # 释放内存

print("[OK] 预处理完成")

# ============================================================
# Step 6: 差异基因分析
# ============================================================
print("\n[Step 6] 差异基因分析...")

de_results = {}

# 100PE vs Control
if '100PE' in adata.obs['sample'].values and 'control' in adata.obs['sample'].values:
    print("\n  分析 100PE vs Control...")
    adata_100pe = adata[adata.obs['sample'].isin(['100PE', 'control'])].copy()
    sc.tl.rank_genes_groups(adata_100pe, groupby='sample', groups=['100PE'], 
                            reference='control', method='wilcoxon')
    de_100pe = sc.get.rank_genes_groups_df(adata_100pe, group='100PE')
    de_100pe['log2fc'] = np.log2(np.exp(de_100pe['logfoldchanges']))
    de_100pe['significant'] = (de_100pe['pvals_adj'] < 0.05) & (abs(de_100pe['log2fc']) > 0.5)
    de_results['100PE_vs_control'] = de_100pe
    de_100pe.to_csv('./results/tables/DE_100PE_vs_control.csv', index=False)
    print(f"    总DEGs: {de_100pe['significant'].sum()}")

# 500PE vs Control
if '500PE' in adata.obs['sample'].values and 'control' in adata.obs['sample'].values:
    print("\n  分析 500PE vs Control...")
    adata_500pe = adata[adata.obs['sample'].isin(['500PE', 'control'])].copy()
    sc.tl.rank_genes_groups(adata_500pe, groupby='sample', groups=['500PE'], 
                            reference='control', method='wilcoxon')
    de_500pe = sc.get.rank_genes_groups_df(adata_500pe, group='500PE')
    de_500pe['log2fc'] = np.log2(np.exp(de_500pe['logfoldchanges']))
    de_500pe['significant'] = (de_500pe['pvals_adj'] < 0.05) & (abs(de_500pe['log2fc']) > 0.5)
    de_results['500PE_vs_control'] = de_500pe
    de_500pe.to_csv('./results/tables/DE_500PE_vs_control.csv', index=False)
    print(f"    总DEGs: {de_500pe['significant'].sum()}")

# 500PE vs 100PE
if '500PE' in adata.obs['sample'].values and '100PE' in adata.obs['sample'].values:
    print("\n  分析 500PE vs 100PE...")
    adata_pe = adata[adata.obs['sample'].isin(['500PE', '100PE'])].copy()
    sc.tl.rank_genes_groups(adata_pe, groupby='sample', groups=['500PE'], 
                            reference='100PE', method='wilcoxon')
    de_500vs100 = sc.get.rank_genes_groups_df(adata_pe, group='500PE')
    de_500vs100['log2fc'] = np.log2(np.exp(de_500vs100['logfoldchanges']))
    de_500vs100['significant'] = (de_500vs100['pvals_adj'] < 0.05) & (abs(de_500vs100['log2fc']) > 0.5)
    de_results['500PE_vs_100PE'] = de_500vs100
    de_500vs100.to_csv('./results/tables/DE_500PE_vs_100PE.csv', index=False)
    print(f"    总DEGs: {de_500vs100['significant'].sum()}")

print("[OK] 差异分析完成")

# ============================================================
# Step 7: 筛选ER stress和线粒体损伤相关差异基因
# ============================================================
print("\n[Step 7] 筛选ER stress和线粒体损伤相关差异基因...")

def filter_target_genes(de_df, target_genes, category_name):
    """筛选目标基因"""
    de_genes_upper = de_df['names'].str.upper()
    target_genes_upper = [g.upper() for g in target_genes]
    matches = de_genes_upper.isin(target_genes_upper)
    filtered = de_df[matches].copy()
    filtered['category'] = category_name
    return filtered

# 存储所有筛选结果
filtered_results = {}

for comparison, de_df in de_results.items():
    print(f"\n  {comparison}:")
    
    # 筛选各类别基因
    er_stress = filter_target_genes(de_df, ER_STRESS_GENES, 'ER_stress')
    mito = filter_target_genes(de_df, MITOCHONDRIA_GENES, 'Mitochondria')
    er_mito = filter_target_genes(de_df, ER_MITO_INTERACTION_GENES, 'ER-Mito_interaction')
    
    # 合并
    combined = pd.concat([er_stress, mito, er_mito], ignore_index=True)
    combined = combined.drop_duplicates(subset=['names'])
    
    # 统计
    sig_combined = combined[combined['significant']]
    print(f"    ER stress相关: {len(er_stress)} (显著: {(er_stress['significant']).sum()})")
    print(f"    线粒体相关: {len(mito)} (显著: {(mito['significant']).sum()})")
    print(f"    ER-线粒体交互: {len(er_mito)} (显著: {(er_mito['significant']).sum()})")
    print(f"    总计: {len(combined)} (显著: {len(sig_combined)})")
    
    filtered_results[comparison] = combined
    
    # 保存
    combined.to_csv(f'./results/tables/ER_Mito_DEGs_{comparison}.csv', index=False)
    sig_combined.to_csv(f'./results/tables/ER_Mito_DEGs_{comparison}_significant.csv', index=False)

print("\n[OK] 目标基因筛选完成")

# ============================================================
# Step 8: KEGG通路富集分析
# ============================================================
print("\n[Step 8] KEGG通路富集分析...")

# KEGG富集分析函数
def run_kegg_enrichment(gene_list, comparison_name, organism='Human'):
    """运行KEGG富集分析"""
    if len(gene_list) < 5:
        print(f"    [跳过] {comparison_name}: 基因数量不足 ({len(gene_list)})")
        return None
    
    try:
        # 使用gseapy进行KEGG富集
        enr = gp.enrichr(
            gene_list=list(gene_list),
            gene_sets=['KEGG_2021_Human'],
            organism='human',
            outdir=f'./results/KEGG/{comparison_name}',
            cutoff=0.05
        )
        
        results = enr.results
        if len(results) > 0:
            results = results[results['Adjusted P-value'] < 0.1]
            results.to_csv(f'./results/KEGG/{comparison_name}_KEGG_enrichment.csv', index=False)
            print(f"    [OK] {comparison_name}: {len(results)} 个显著通路")
            return results
        else:
            print(f"    [无结果] {comparison_name}")
            return None
    except Exception as e:
        print(f"    [错误] {comparison_name}: {e}")
        return None

kegg_results = {}

for comparison, filtered_df in filtered_results.items():
    print(f"\n  {comparison}:")
    
    # 获取显著差异基因
    sig_genes = filtered_df[filtered_df['significant']]['names'].tolist()
    
    # 上调基因
    up_genes = filtered_df[(filtered_df['significant']) & (filtered_df['log2fc'] > 0.5)]['names'].tolist()
    
    # 下调基因
    down_genes = filtered_df[(filtered_df['significant']) & (filtered_df['log2fc'] < -0.5)]['names'].tolist()
    
    # 分别进行富集分析
    kegg_results[f'{comparison}_all'] = run_kegg_enrichment(sig_genes, f'{comparison}_all')
    kegg_results[f'{comparison}_up'] = run_kegg_enrichment(up_genes, f'{comparison}_up')
    kegg_results[f'{comparison}_down'] = run_kegg_enrichment(down_genes, f'{comparison}_down')

# 同时对所有差异基因（不限于目标基因）进行KEGG富集
print("\n  全基因组KEGG富集分析:")
for comparison, de_df in de_results.items():
    sig_genes = de_df[de_df['significant']]['names'].tolist()
    kegg_results[f'{comparison}_genome_wide'] = run_kegg_enrichment(sig_genes[:500], f'{comparison}_genome_wide')

print("\n[OK] KEGG富集分析完成")

# ============================================================
# Step 9: 创建基因-通路网络
# ============================================================
print("\n[Step 9] 创建基因-通路网络...")

def create_kegg_network(kegg_df, filtered_df, comparison_name, top_n_pathways=15):
    """创建基因-KEGG通路网络"""
    if kegg_df is None or len(kegg_df) == 0:
        print(f"    [跳过] {comparison_name}: 无富集结果")
        return
    
    # 取top通路
    top_pathways = kegg_df.nsmallest(top_n_pathways, 'Adjusted P-value')
    
    # 创建网络
    G = nx.Graph()
    
    # 基因节点颜色映射
    gene_colors = {}
    for _, row in filtered_df.iterrows():
        if row['significant']:
            if row['log2fc'] > 0.5:
                gene_colors[row['names']] = '#e74c3c'  # 红色 - 上调
            else:
                gene_colors[row['names']] = '#3498db'  # 蓝色 - 下调
    
    # 添加节点和边
    pathway_nodes = []
    gene_nodes = set()
    
    for _, pathway in top_pathways.iterrows():
        pathway_name = pathway['Term'].split('_')[0] if '_' in pathway['Term'] else pathway['Term']
        pathway_name = pathway_name[:40] + '...' if len(pathway_name) > 40 else pathway_name
        
        G.add_node(pathway_name, node_type='pathway', pvalue=pathway['Adjusted P-value'])
        pathway_nodes.append(pathway_name)
        
        # 解析基因
        genes = pathway['Genes'].split(';')
        for gene in genes:
            gene = gene.strip().upper()
            if gene in [g.upper() for g in filtered_df['names'].values]:
                original_gene = filtered_df[filtered_df['names'].str.upper() == gene]['names'].values[0]
                G.add_node(original_gene, node_type='gene')
                G.add_edge(pathway_name, original_gene)
                gene_nodes.add(original_gene)
    
    if len(G.nodes()) < 3:
        print(f"    [跳过] {comparison_name}: 网络节点不足")
        return
    
    # 绘制网络
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # 布局
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # 绘制通路节点
    pathway_sizes = [3000 - kegg_df[kegg_df['Term'].str.contains(p.replace('...', ''))]['Adjusted P-value'].values[0] * 10000 
                     if len(kegg_df[kegg_df['Term'].str.contains(p.replace('...', ''))]) > 0 else 2000
                     for p in pathway_nodes if p in G.nodes()]
    
    nx.draw_networkx_nodes(G, pos, nodelist=[n for n in pathway_nodes if n in G.nodes()],
                          node_color='#2ecc71', node_size=pathway_sizes,
                          node_shape='s', alpha=0.8, ax=ax)
    
    # 绘制基因节点
    gene_node_colors = [gene_colors.get(g, '#95a5a6') for g in gene_nodes if g in G.nodes()]
    nx.draw_networkx_nodes(G, pos, nodelist=[n for n in gene_nodes if n in G.nodes()],
                          node_color=gene_node_colors, node_size=800,
                          node_shape='o', alpha=0.8, ax=ax)
    
    # 绘制边
    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='gray', ax=ax)
    
    # 标签
    labels = {n: n for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=7, ax=ax)
    
    # 图例
    legend_elements = [
        Patch(facecolor='#2ecc71', label='KEGG Pathway'),
        Patch(facecolor='#e74c3c', label='Upregulated Gene'),
        Patch(facecolor='#3498db', label='Downregulated Gene'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
    
    ax.set_title(f'Gene-KEGG Pathway Network\n{comparison_name}', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(f'./results/networks/{comparison_name}_KEGG_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 保存网络数据
    nx.write_gexf(G, f'./results/networks/{comparison_name}_KEGG_network.gexf')
    
    print(f"    [OK] {comparison_name}: {len(pathway_nodes)} 通路, {len(gene_nodes)} 基因")

# 为每个比较创建网络
for comparison, filtered_df in filtered_results.items():
    kegg_df = kegg_results.get(f'{comparison}_all')
    create_kegg_network(kegg_df, filtered_df, comparison)

print("\n[OK] 基因网络创建完成")

# ============================================================
# Step 10: 创建综合可视化
# ============================================================
print("\n[Step 10] 创建综合可视化...")

# 1. ER stress和线粒体基因表达热图
print("  创建表达热图...")

# 获取在数据中存在的目标基因
available_genes = [g for g in ALL_TARGET_GENES if g in adata.var_names]
print(f"    数据中找到 {len(available_genes)} 个目标基因")

if len(available_genes) > 10:
    # 计算每个样本的平均表达
    exp_data = []
    for sample in ['control', '100PE', '500PE']:
        if sample in adata.obs['sample'].values:
            sample_mask = adata.obs['sample'] == sample
            if adata.raw is not None:
                sample_exp = adata.raw[:, available_genes].X[sample_mask].mean(axis=0)
            else:
                sample_exp = adata[:, available_genes].X[sample_mask].mean(axis=0)
            if hasattr(sample_exp, 'A1'):
                exp_data.append(sample_exp.A1)
            else:
                exp_data.append(np.array(sample_exp).flatten())
    
    if exp_data:
        exp_df = pd.DataFrame(
            np.array(exp_data).T,
            index=available_genes,
            columns=[s for s in ['control', '100PE', '500PE'] if s in adata.obs['sample'].values]
        )
        
        # 标准化用于热图显示
        exp_df_scaled = (exp_df - exp_df.mean()) / (exp_df.std() + 0.001)
        
        # 选择变化最大的基因
        exp_df_scaled['variance'] = exp_df_scaled.var(axis=1)
        top_var_genes = exp_df_scaled.nlargest(50, 'variance').index.tolist()
        
        # 创建热图
        fig, ax = plt.subplots(figsize=(10, 14))
        heatmap_data = exp_df_scaled.loc[top_var_genes].drop(columns='variance')
        
        sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, 
                   vmin=-2, vmax=2, ax=ax,
                   cbar_kws={'label': 'Z-score'})
        ax.set_title('ER Stress & Mitochondrial Genes Expression\n(Top 50 Variable Genes)', 
                    fontsize=14, fontweight='bold')
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Gene', fontsize=12)
        
        plt.tight_layout()
        plt.savefig('./results/figures/ER_Mito_expression_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("    [OK] 表达热图已保存")

# 2. 差异基因火山图（突出ER-Mito基因）
print("  创建火山图...")

for comparison, de_df in de_results.items():
    filtered_df = filtered_results.get(comparison)
    if filtered_df is None:
        continue
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # 获取目标基因集
    target_genes_upper = set([g.upper() for g in ALL_TARGET_GENES])
    
    # 分类颜色
    colors = []
    sizes = []
    for _, row in de_df.iterrows():
        gene_upper = row['names'].upper()
        is_target = gene_upper in target_genes_upper
        is_sig = row['significant']
        
        if is_target and is_sig and row['log2fc'] > 0.5:
            colors.append('#e74c3c')  # 目标基因-上调
            sizes.append(80)
        elif is_target and is_sig and row['log2fc'] < -0.5:
            colors.append('#3498db')  # 目标基因-下调
            sizes.append(80)
        elif is_target:
            colors.append('#f39c12')  # 目标基因-不显著
            sizes.append(40)
        elif is_sig and row['log2fc'] > 0.5:
            colors.append('#ffcccc')  # 其他-上调
            sizes.append(15)
        elif is_sig and row['log2fc'] < -0.5:
            colors.append('#cce5ff')  # 其他-下调
            sizes.append(15)
        else:
            colors.append('#d3d3d3')  # 不显著
            sizes.append(10)
    
    ax.scatter(de_df['log2fc'], -np.log10(de_df['pvals_adj']),
              c=colors, s=sizes, alpha=0.6, edgecolors='none')
    
    # 阈值线
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(-0.5, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # 标注显著的目标基因
    sig_target = filtered_df[filtered_df['significant']].nsmallest(15, 'pvals_adj')
    for _, row in sig_target.iterrows():
        ax.annotate(row['names'],
                   xy=(row['log2fc'], -np.log10(row['pvals_adj'])),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7))
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12, fontweight='bold')
    ax.set_title(f'Volcano Plot: {comparison}\n(ER Stress & Mitochondrial Genes Highlighted)', 
                fontsize=14, fontweight='bold')
    
    # 图例
    legend_elements = [
        Patch(facecolor='#e74c3c', label='ER/Mito - Upregulated'),
        Patch(facecolor='#3498db', label='ER/Mito - Downregulated'),
        Patch(facecolor='#f39c12', label='ER/Mito - Not Sig'),
        Patch(facecolor='#ffcccc', label='Other - Upregulated'),
        Patch(facecolor='#cce5ff', label='Other - Downregulated'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'./results/figures/volcano_{comparison}_ER_Mito.png', dpi=300, bbox_inches='tight')
    plt.close()

print("    [OK] 火山图已保存")

# 3. KEGG富集条形图
print("  创建KEGG富集条形图...")

for key, kegg_df in kegg_results.items():
    if kegg_df is None or len(kegg_df) == 0:
        continue
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 取top 15通路
    top_kegg = kegg_df.nsmallest(15, 'Adjusted P-value').copy()
    top_kegg['Term_short'] = top_kegg['Term'].apply(lambda x: x[:50] + '...' if len(x) > 50 else x)
    top_kegg = top_kegg.sort_values('Adjusted P-value', ascending=True)
    
    # 条形图
    bars = ax.barh(range(len(top_kegg)), -np.log10(top_kegg['Adjusted P-value']),
                   color=plt.cm.Reds(np.linspace(0.3, 0.9, len(top_kegg))))
    
    ax.set_yticks(range(len(top_kegg)))
    ax.set_yticklabels(top_kegg['Term_short'])
    ax.set_xlabel('-Log10(Adjusted P-value)', fontsize=12, fontweight='bold')
    ax.set_title(f'KEGG Pathway Enrichment\n{key}', fontsize=14, fontweight='bold')
    
    # 添加基因数量标注
    for i, (idx, row) in enumerate(top_kegg.iterrows()):
        gene_count = len(row['Genes'].split(';'))
        ax.text(-np.log10(row['Adjusted P-value']) + 0.1, i, f'n={gene_count}', 
               va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'./results/KEGG/{key}_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()

print("    [OK] KEGG富集条形图已保存")

# 4. 基因分类统计图
print("  创建基因分类统计图...")

fig, axes = plt.subplots(1, len(filtered_results), figsize=(5*len(filtered_results), 6))
if len(filtered_results) == 1:
    axes = [axes]

for idx, (comparison, filtered_df) in enumerate(filtered_results.items()):
    ax = axes[idx]
    
    # 统计各类别
    categories = ['ER_stress', 'Mitochondria', 'ER-Mito_interaction']
    up_counts = []
    down_counts = []
    
    for cat in categories:
        cat_df = filtered_df[filtered_df['category'] == cat]
        up_counts.append(((cat_df['significant']) & (cat_df['log2fc'] > 0.5)).sum())
        down_counts.append(((cat_df['significant']) & (cat_df['log2fc'] < -0.5)).sum())
    
    x = np.arange(len(categories))
    width = 0.35
    
    ax.bar(x - width/2, up_counts, width, label='Upregulated', color='#e74c3c')
    ax.bar(x + width/2, down_counts, width, label='Downregulated', color='#3498db')
    
    ax.set_ylabel('Number of DEGs')
    ax.set_title(f'{comparison}')
    ax.set_xticks(x)
    ax.set_xticklabels(['ER Stress', 'Mitochondria', 'ER-Mito\nInteraction'], rotation=45, ha='right')
    ax.legend()
    
    # 添加数值标注
    for i, (up, down) in enumerate(zip(up_counts, down_counts)):
        if up > 0:
            ax.text(i - width/2, up + 0.5, str(up), ha='center', fontsize=9)
        if down > 0:
            ax.text(i + width/2, down + 0.5, str(down), ha='center', fontsize=9)

plt.suptitle('ER Stress & Mitochondrial Damage Related DEGs', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('./results/figures/ER_Mito_DEGs_summary.png', dpi=300, bbox_inches='tight')
plt.close()
print("    [OK] 基因分类统计图已保存")

# ============================================================
# Step 11: 生成分析报告
# ============================================================
print("\n[Step 11] 生成分析报告...")

report = """
================================================================================
ER Stress-Induced Mitochondrial Damage Analysis Report
内质网压力引起线粒体损伤 - 分析报告
================================================================================

数据概览:
---------
"""

report += f"总细胞数: {adata.shape[0]:,}\n"
report += f"总基因数: {adata.shape[1]:,}\n\n"

for sample in ['control', '100PE', '500PE']:
    if sample in adata.obs['sample'].values:
        report += f"  {sample}: {(adata.obs['sample'] == sample).sum():,} cells\n"

report += f"""
目标基因定义:
------------
  ER stress相关基因: {len(ER_STRESS_GENES)}
  线粒体相关基因: {len(MITOCHONDRIA_GENES)}
  ER-线粒体交互基因: {len(ER_MITO_INTERACTION_GENES)}
  总计: {len(ALL_TARGET_GENES)} (去重后)

差异表达分析结果:
----------------
"""

for comparison, de_df in de_results.items():
    sig_count = de_df['significant'].sum()
    up_count = ((de_df['log2fc'] > 0.5) & (de_df['pvals_adj'] < 0.05)).sum()
    down_count = ((de_df['log2fc'] < -0.5) & (de_df['pvals_adj'] < 0.05)).sum()
    report += f"\n{comparison}:\n"
    report += f"  总DEGs: {sig_count}\n"
    report += f"  上调: {up_count}, 下调: {down_count}\n"

report += """
ER Stress & 线粒体损伤相关DEGs:
-----------------------------
"""

for comparison, filtered_df in filtered_results.items():
    sig_filtered = filtered_df[filtered_df['significant']]
    report += f"\n{comparison}:\n"
    report += f"  总目标DEGs: {len(sig_filtered)}\n"
    
    for cat in ['ER_stress', 'Mitochondria', 'ER-Mito_interaction']:
        cat_df = sig_filtered[sig_filtered['category'] == cat]
        up = (cat_df['log2fc'] > 0.5).sum()
        down = (cat_df['log2fc'] < -0.5).sum()
        report += f"    {cat}: {len(cat_df)} (↑{up}, ↓{down})\n"
    
    # 列出top基因
    if len(sig_filtered) > 0:
        report += "\n  Top 10 显著基因:\n"
        for _, row in sig_filtered.nsmallest(10, 'pvals_adj').iterrows():
            direction = "↑" if row['log2fc'] > 0 else "↓"
            report += f"    {direction} {row['names']}: log2FC={row['log2fc']:.2f}, p_adj={row['pvals_adj']:.2e}\n"

report += """
KEGG通路富集结果:
----------------
"""

for key, kegg_df in kegg_results.items():
    if kegg_df is not None and len(kegg_df) > 0:
        report += f"\n{key}:\n"
        report += f"  显著富集通路数: {len(kegg_df)}\n"
        report += "  Top 5 通路:\n"
        for _, row in kegg_df.nsmallest(5, 'Adjusted P-value').iterrows():
            report += f"    - {row['Term'][:60]}\n"
            report += f"      P_adj={row['Adjusted P-value']:.2e}, Genes: {row['Genes'][:50]}...\n"

report += """
================================================================================
输出文件:
--------
./results/tables/
  - DE_*_vs_*.csv: 完整差异基因列表
  - ER_Mito_DEGs_*.csv: ER stress和线粒体相关差异基因

./results/KEGG/
  - *_KEGG_enrichment.csv: KEGG富集结果
  - *_barplot.png: 富集条形图

./results/networks/
  - *_KEGG_network.png: 基因-通路网络图
  - *_KEGG_network.gexf: 网络文件(可用Gephi/Cytoscape打开)

./results/figures/
  - volcano_*_ER_Mito.png: 火山图
  - ER_Mito_expression_heatmap.png: 表达热图
  - ER_Mito_DEGs_summary.png: DEGs统计汇总

================================================================================
"""

print(report)

with open('./results/ER_Mito_Analysis_Report.txt', 'w', encoding='utf-8') as f:
    f.write(report)

# 保存处理后的数据
print("\n保存处理后的数据...")
adata.write('./results/processed_data.h5ad', compression='gzip')

print("\n" + "="*70)
print("[完成] 所有分析已完成！")
print("="*70)
print("\n结果保存在: ./results/")
print("  - figures/: 可视化图表")
print("  - tables/: 差异基因表格")
print("  - KEGG/: KEGG富集分析结果")
print("  - networks/: 基因-通路网络")
print("  - ER_Mito_Analysis_Report.txt: 分析报告")

