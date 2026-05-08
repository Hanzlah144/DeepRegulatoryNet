import os
import logging
import time
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests_cache
from chembl_webresource_client.new_client import new_client

# Setup caching for API calls
requests_cache.install_cache('chembl_cache', expire_after=86400)

# ================= CONFIG =================
DEFAULT_INPUT_FILE = "output/hub_genes.csv"
OUTPUT_DIR = "output"

OUTPUT_CSV = os.path.join(OUTPUT_DIR, "chembl_drug_gene_interactions.csv")
# New consolidated visualization path
OUTPUT_PLOT_HEATMAP = os.path.join(OUTPUT_DIR, "top_drugs_potency_heatmap.png")

IC50_THRESHOLD = 5000  
TOP_N_GENES = 5  # Number of genes to show in heatmap
TOP_N_DRUGS = 10 # Number of drugs per gene to show in heatmap

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ================= LOGGING =================
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s")
logger = logging.getLogger("ChEMBLPipeline")

# ================= LOAD GENES =================
def load_genes(file_path):
    df = pd.read_csv(file_path)
    if "Gene" not in df.columns:
        raise ValueError("CSV must contain 'Gene' column")
    return df["Gene"].dropna().unique().tolist()

# ================= MAP GENE → TARGET =================
def get_targets(genes):
    target_client = new_client.target
    gene_to_target = {}
    for gene in genes:
        res = target_client.search(gene)
        for r in res:
            if r.get("organism") == "Homo sapiens":
                gene_to_target[gene] = r["target_chembl_id"]
                break
    return gene_to_target

# ================= FETCH BIOACTIVITY =================
def fetch_activities(gene_to_target):
    activity_client = new_client.activity
    records = []

    for gene, target_id in gene_to_target.items():
        logger.info("Fetching activities for %s", gene)
        try:
            acts = activity_client.filter(target_chembl_id=target_id)[:300]
        except Exception as e:
            logger.warning("Failed for %s: %s", gene, str(e))
            continue

        for a in acts:
            if a.get("standard_type") not in ["IC50", "Ki"]:
                continue
            try:
                value = float(a.get("standard_value"))
                if value <= 0 or value > IC50_THRESHOLD:
                    continue
                
                records.append({
                    "Gene": gene,
                    "Drug": a.get("molecule_pref_name") or a.get("molecule_chembl_id"),
                    "Value_nM": value,
                    "pIC50": -np.log10(value * 1e-9)
                })
            except:
                continue

    return pd.DataFrame(records)

# ================= UPDATED VISUALIZATION =================

def plot_labeled_potency_heatmap(df):
    """
    Consolidated Visualization: Ranks top drugs per gene and labels them with names 
    and pIC50 values for actionable insights.
    """
    logger.info("Generating Potency Heatmap with Drug Names...")
    
    # Identify top genes by a combination of potency and data density
    target_metrics = df.groupby('Gene').agg(
        Avg_pIC50=('pIC50', 'mean'),
        Drug_Count=('Drug', 'count')
    ).reset_index()
    target_metrics['Score'] = target_metrics['Avg_pIC50'] * target_metrics['Drug_Count']
    top_genes = target_metrics.sort_values(by='Score', ascending=False).head(TOP_N_GENES)['Gene'].tolist()

    # Filter and rank top 10 drugs per gene
    subset = df[df['Gene'].isin(top_genes)].copy()
    subset = subset.sort_values(['Gene', 'pIC50'], ascending=[True, False])
    subset['Drug_Rank'] = subset.groupby('Gene').cumcount() + 1
    subset = subset[subset['Drug_Rank'] <= TOP_N_DRUGS]

    # Pivot data for the heatmap
    pivoted_values = subset.pivot(index='Gene', columns='Drug_Rank', values='pIC50')
    pivoted_names = subset.pivot(index='Gene', columns='Drug_Rank', values='Drug')

    # Create annotation strings (Name + pIC50)
    annot_matrix = pivoted_names + "\n(" + pivoted_values.round(2).astype(str) + ")"

    plt.figure(figsize=(18, 8))
    sns.heatmap(pivoted_values, 
                annot=annot_matrix, 
                fmt="", 
                cmap='YlGnBu', 
                cbar_kws={'label': 'pIC50 Potency'},
                annot_kws={"size": 9})

    plt.title(f'Top {TOP_N_DRUGS} Potent Drugs per Primary Target Gene', fontsize=16)
    plt.xlabel('Drug Rank (by pIC50)', fontsize=12)
    plt.ylabel('Target Gene', fontsize=12)
    plt.tight_layout()
    
    plt.savefig(OUTPUT_PLOT_HEATMAP, dpi=300)
    plt.close()
    logger.info("Optimized heatmap saved to %s", OUTPUT_PLOT_HEATMAP)

# ================= MAIN =================
def main(hub_genes_path=None, max_genes=None):
    if hub_genes_path is None:
        hub_genes_path = DEFAULT_INPUT_FILE

    if not os.path.exists(hub_genes_path):
        logger.error("Input file %s not found.", hub_genes_path)
        return

    genes = load_genes(hub_genes_path)
    if max_genes:
        genes = genes[:max_genes]

    gene_to_target = get_targets(genes)
    df = fetch_activities(gene_to_target)

    if not df.empty:
        df.to_csv(OUTPUT_CSV, index=False)
        
        # Apply the new consolidated visualization
        plot_labeled_potency_heatmap(df)
        
        logger.info("Pipeline completed successfully. Insights generated.")
    else:
        logger.warning("No interactions found with current filters.")

if __name__ == "__main__":
    main()