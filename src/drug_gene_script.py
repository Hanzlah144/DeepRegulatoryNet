import os
import logging
import time
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt
import requests_cache
from chembl_webresource_client.new_client import new_client

# Setup caching for API calls
requests_cache.install_cache('chembl_cache', expire_after=86400)  # Cache for 1 day

# ================= CONFIG =================
DEFAULT_INPUT_FILE = "output/hub_genes.csv"
OUTPUT_DIR = "output"

OUTPUT_CSV = os.path.join(OUTPUT_DIR, "chembl_drug_gene_interactions.csv")
OUTPUT_PLOT = os.path.join(OUTPUT_DIR, "chembl_top_drugs.png")

IC50_THRESHOLD = 1000   # nM (loose filter)
TOP_N_DRUGS = 15

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ================= LOGGING =================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
)
logger = logging.getLogger("ChEMBLPipeline")

# ================= LOAD GENES =================
def load_genes(file_path):
    df = pd.read_csv(file_path)
    if "Gene" not in df.columns:
        raise ValueError("CSV must contain 'Gene' column")
    genes = df["Gene"].dropna().unique().tolist()
    logger.info("Loaded %d genes", len(genes))
    return genes

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

    logger.info("Mapped %d/%d genes to ChEMBL targets",
                len(gene_to_target), len(genes))
    return gene_to_target

# ================= FETCH BIOACTIVITY =================
def fetch_activities(gene_to_target):
    activity_client = new_client.activity
    records = []

    for gene, target_id in gene_to_target.items():
        logger.info("Fetching activities for %s (%s)", gene, target_id)

        try:
            # Limit to first 100 activities per target to speed up
            acts = activity_client.filter(target_chembl_id=target_id)[:100]
            logger.info("Retrieved %d activities for %s", len(acts), gene)
        except Exception as e:
            logger.warning("Failed to fetch activities for %s: %s", gene, str(e))
            continue

        for a in acts:
            if a.get("standard_type") not in ["IC50", "Ki"]:
                continue

            try:
                value = float(a.get("standard_value"))
            except:
                continue

            if value > IC50_THRESHOLD:
                continue

            records.append({
                "Gene": gene,
                "Drug": a.get("molecule_pref_name"),
                "ChEMBL_ID": a.get("molecule_chembl_id"),
                "Type": a.get("standard_type"),
                "Value_nM": value,
                "Assay_Type": a.get("assay_type"),
                "Confidence": a.get("confidence_score"),
            })

    df = pd.DataFrame(records)
    logger.info("Collected %d filtered interactions", len(df))
    return df

# ================= RANK DRUGS =================
def rank_drugs(df):
    if df.empty:
        return df

    grouped = df.groupby("Drug").agg({
        "Gene": "nunique",
        "Value_nM": "mean"
    }).rename(columns={
        "Gene": "Num_Targets",
        "Value_nM": "Avg_IC50"
    })

    # Ranking score: more targets + lower IC50
    grouped["Score"] = grouped["Num_Targets"] / grouped["Avg_IC50"]

    ranked = grouped.sort_values("Score", ascending=False).reset_index()
    return ranked

# ================= PLOT =================
def plot_top_drugs(ranked_df):
    top = ranked_df.head(TOP_N_DRUGS)

    plt.figure(figsize=(11, 7))

    y_positions = range(len(top))
    sizes = (top["Num_Targets"] / top["Num_Targets"].max()) * 400 + 50

    scatter = plt.scatter(
        top["Score"],
        y_positions,
        s=sizes,
        c=top["Score"],
        cmap="viridis",
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )

    plt.yticks(y_positions, top["Drug"])
    plt.xlabel("Drug Score (Targets / Avg IC50)")
    plt.title("Top Candidate Drugs (ChEMBL Analysis) - Bubble Scatter")
    plt.colorbar(scatter, label="Score")
    plt.grid(axis="x", linestyle="--", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    plt.close()

    logger.info("Plot saved → %s", OUTPUT_PLOT)

# ================= MAIN =================
def main(hub_genes_path=None, max_genes=None):
    if hub_genes_path is None:
        hub_genes_path = DEFAULT_INPUT_FILE

    start = datetime.now()
    logger.info("=== ChEMBL Pipeline Started ===")

    genes = load_genes(hub_genes_path)
    total_genes = len(genes)

    if max_genes is not None and max_genes > 0 and len(genes) > max_genes:
        logger.info("[INFO] Fast mode: limiting drug-gene analysis to top %d of %d genes", max_genes, total_genes)
        genes = genes[:max_genes]
    else:
        logger.info("[INFO] Drug-gene analysis on %d genes", total_genes)

    t0 = time.time()
    gene_to_target = get_targets(genes)
    logger.info("[TIMING] get_targets duration: %.2f seconds", time.time() - t0)

    t1 = time.time()
    df = fetch_activities(gene_to_target)
    logger.info("[TIMING] fetch_activities duration: %.2f seconds", time.time() - t1)

    if df.empty:
        logger.warning("No interactions found")
        return

    df.to_csv(OUTPUT_CSV, index=False)
    logger.info("Saved interactions → %s", OUTPUT_CSV)

    ranked = rank_drugs(df)

    plot_top_drugs(ranked)

    logger.info("Pipeline completed in %s", datetime.now() - start)

