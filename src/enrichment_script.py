import os
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from matplotlib.colors import LinearSegmentedColormap

def load_genes(file_path):
    logger = logging.getLogger()
    try:
        df = pd.read_csv(file_path)
        genes = df['Gene'].dropna().unique().tolist()
        logger.info("[INFO]  %d genes loaded", len(genes))
        return genes
    except Exception as e:
        logger.error("[ERROR]  Failed to load genes: %s", e)
        return []

def fetch_latest_libraries():
    logger = logging.getLogger()
    try:
        libs = gp.get_library_name(organism='human')
        filtered = [
            lib for lib in libs
            if any(x in lib for x in ['GO_Biological_Process', 'GO_Molecular_Function', 'GO_Cellular_Component', 'KEGG'])
            and any(str(y) in lib for y in range(2021, 2026))
        ]
        recent = sorted(filtered, key=lambda x: int(x.split('_')[-1]) if x.split('_')[-1].isdigit() else 0, reverse=True)
        selected = []
        seen = set()
        for lib in recent:
            base = '_'.join(lib.split('_')[:-1])
            if base not in seen:
                selected.append(lib)
                seen.add(base)
        logger.debug("[DEBUG] Enrichment libraries: %s", ', '.join(selected[:4]))
        return selected[:4] if selected else ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human']
    except Exception:
        return ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human']

def run_enrichment(genes, libs):
    logger = logging.getLogger()
    try:
        result = gp.enrichr(gene_list=genes, gene_sets=libs, organism='human', outdir=None)
        return result.results
    except Exception as e:
        logger.error("[ERROR]  Enrichment failed: %s", e)
        return None

def prepare_plot_data(df):
    logger = logging.getLogger()
    df = df[df['Adjusted P-value'] < 0.05].copy()
    if df.empty:
        logger.info("[INFO] ⚠ No significant terms found")
        return None
    df['Gene Count'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
    df['Gene Ratio'] = df['Overlap'].apply(lambda x: eval(x))
    df['-log10(p)'] = -np.log10(df['Adjusted P-value'].replace(0, 1e-10))
    return df.sort_values('-log10(p)', ascending=False).head(20)

def plot_bubble(df, out_path):
    logger = logging.getLogger()
    if df is None or df.empty:
        logger.info("[INFO]  No data to plot")
        return
    plt.figure(figsize=(10, 7))
    bubble = plt.scatter(
        df['Gene Ratio'], df['Term'],
        s=df['Gene Count'] * 40,
        c=df['-log10(p)'],
        cmap='coolwarm',
        edgecolor='black', alpha=0.8
    )
    plt.colorbar(bubble, label='-log10(adj p-value)')
    plt.xlabel('Gene Ratio')
    plt.title('Top Enriched Terms')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logger.info("[INFO]  Bubble plot: %s", out_path)

def main(file_path, output_dir='enrichment_results'):
    logger = logging.getLogger()
    os.makedirs(output_dir, exist_ok=True)
    genes = load_genes(file_path)
    if not genes:
        logger.info("[INFO]  No valid genes")
        return

    libs = fetch_latest_libraries()
    results = run_enrichment(genes, libs)
    if results is None or results.empty:
        logger.info("[INFO]  No enrichment results")
        return

    filtered = prepare_plot_data(results)
    if filtered is None:
        logger.info("[INFO]  No significant terms")
        return

    out_csv = os.path.join(output_dir, 'enrichment_results.csv')
    filtered.to_csv(out_csv, index=False)
    logger.info("[INFO]  Saved: %s", out_csv)
    plot_bubble(filtered, os.path.join(output_dir, 'bubble_plot.png'))