import os
import json
import logging
from datetime import datetime

import pandas as pd
import requests
import matplotlib.pyplot as plt

INPUT_FILE = "output/hub_genes.csv"        
OUTPUT_DIR = "output"
OUTPUT_FULL_CSV = os.path.join(OUTPUT_DIR, "drug_gene_interactions.csv")
OUTPUT_BAR = os.path.join(OUTPUT_DIR, "drug_gene_interaction_barplot.png")

TOP_GENES = 15  

os.makedirs(OUTPUT_DIR, exist_ok=True)

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("DGIdbPipeline")

logger = setup_logging()

def load_gene_names(file_path):
    try:
        df = pd.read_csv(file_path)
        if "Gene" not in df.columns:
            raise ValueError("CSV is missing a 'Gene' column.")
        genes = df["Gene"].dropna().unique().tolist()
        logger.info("Loaded %d unique genes.", len(genes))
        return genes
    except Exception as e:
        logger.error("Failed to read input file: %s", e)
        return []

def query_dgidb(genes):
    url = "https://dgidb.org/api/graphql"
    query = f"""
    {{
      genes(names: {json.dumps(genes)}) {{
        nodes {{
          name
          interactions {{
            drug {{ name conceptId }}
            interactionScore
            interactionTypes {{ type directionality }}
            interactionAttributes {{ name value }}
            publications {{ pmid }}
            sources {{ sourceDbName }}
          }}
        }}
      }}
    }}
    """
    try:
        r = requests.post(url, json={"query": query}, timeout=60)
        if r.status_code != 200:
            logger.error("DGIdb API error (%d): %s", r.status_code, r.text)
            return []
        return r.json().get("data", {}).get("genes", {}).get("nodes", [])
    except Exception as e:
        logger.error("DGIdb request failed: %s", e)
        return []

def extract_interactions(data):
    records = []
    for item in data:
        gene = item.get("name", "")
        for inter in item.get("interactions", []):
            drug_info = inter.get("drug", {})
            records.append({
                "Gene": gene,
                "Drug": drug_info.get("name", ""),
                "Concept_ID": drug_info.get("conceptId", ""),
                "Score": inter.get("interactionScore", 0),
                "Interaction_Types": "; ".join([
                    f"{t.get('type','')} ({t.get('directionality','')})"
                    for t in inter.get("interactionTypes", [])
                ]) or None,
                "Interaction_Attributes": "; ".join([
                    f"{a.get('name','')}={a.get('value','')}"
                    for a in inter.get("interactionAttributes", [])
                ]),
                "Publications_PMIDs": "; ".join([
                    str(p.get("pmid", "")) for p in inter.get("publications", [])
                ]),
                "Sources": "; ".join([
                    s.get("sourceDbName", "") for s in inter.get("sources", [])
                ])
            })
    df = pd.DataFrame(records)
    logger.info("Extracted %d drug–gene interactions.", len(df))
    return df




def export_full_csv(df):
    try:
        df.to_csv(OUTPUT_FULL_CSV, index=False)
        logger.info("Full CSV saved → %s", OUTPUT_FULL_CSV)
    except Exception as e:
        logger.error("Failed to save full CSV: %s", e)


def draw_stacked_bar(df):
    
    stacked_df = df.copy()
    stacked_df["Interaction_Types"] = stacked_df["Interaction_Types"].fillna("Unknown")

    
    gene_interaction_counts = (
        stacked_df.groupby(["Gene", "Interaction_Types"])["Drug"]
        .nunique()
        .reset_index()
    )

    
    pivot_df = (
        gene_interaction_counts.pivot(
            index="Gene", columns="Interaction_Types", values="Drug"
        )
        .fillna(0)
    )

    
    top_genes = (
        stacked_df.groupby("Gene")["Drug"]
        .nunique()
        .sort_values(ascending=False)
        .head(TOP_GENES)
        .index
    )
    pivot_df = pivot_df.loc[pivot_df.index.intersection(top_genes)]

    
    ax = pivot_df.plot(
        kind="bar",
        stacked=True,
        figsize=(12, 7),
        colormap="tab20"
    )
    ax.grid(True, which="both", axis="both", linestyle="--", alpha=0.6)
    plt.title("Drugs per Gene by Interaction Type", fontsize=14, weight="bold")
    plt.xlabel("Genes", fontsize=12)
    plt.ylabel("Number of Drugs", fontsize=12)
    plt.xticks(rotation=75, ha="right")
    plt.legend(title="Interaction Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(OUTPUT_BAR, dpi=300)
    plt.close()
    logger.info("Stacked bar plot saved → %s", OUTPUT_BAR)




def run_pipeline():
    start = datetime.now()
    logger.info("=== DGIdb Drug–Gene Interaction Pipeline Started ===")

    genes = load_gene_names(INPUT_FILE)
    if not genes:
        logger.warning("No valid gene names. Exiting.")
        return

    logger.info("Querying DGIdb for %d genes...", len(genes))
    data = query_dgidb(genes)
    if not data:
        logger.warning("No interaction data retrieved.")
        return

    df = extract_interactions(data)
    if df.empty:
        logger.warning("No interactions to process.")
        return

    export_full_csv(df)
    draw_stacked_bar(df)

    logger.info("Pipeline completed in %s.", datetime.now() - start)




def main(hub_genes_path="hub_genes.csv", output_dir="output"):
    """
    Run DGIdb drug gene interaction pipeline.
    Parameters
    ----------
    hub genes path : str
        Path to CSV file containing a 'Gene' column (hub genes).
    output dir : str
        Directory where results (CSV + plot) will be saved.
    """
    global INPUT_FILE, OUTPUT_DIR, OUTPUT_FULL_CSV, OUTPUT_BAR
    INPUT_FILE = hub_genes_path
    OUTPUT_DIR = output_dir
    OUTPUT_FULL_CSV = os.path.join(OUTPUT_DIR, "drug_gene_interactions.csv")
    OUTPUT_BAR = os.path.join(OUTPUT_DIR, "drug_gene_interaction_barplot.png")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    run_pipeline()