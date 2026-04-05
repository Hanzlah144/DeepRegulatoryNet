import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
from datetime import datetime

logger = logging.getLogger(__name__)

# -----------------------
# MAIN ENRICHMENT FUNCTION
# -----------------------
def run_enrichment_pipeline(genes, temp_dir, output_dir, top_n=30):

    start_time = datetime.now()

    if not genes:
        logger.warning("No genes provided for enrichment")
        return

    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------
    # FIXED LIBRARIES (2024–2025)
    # -----------------------
    libraries = [
        "GO_Biological_Process_2023",   # latest stable GO
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023",
        "KEGG_2021_Human",              # KEGG latest available in Enrichr
        "Reactome_Pathways_2024"        # latest Reactome
    ]

    all_results = []

    for lib in libraries:
        try:
            logger.info(f"[Enrichment] Processing {lib}")

            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=lib,
                organism='human',
                outdir=None
            )

            df = enr.results

            if df.empty:
                continue

            # -----------------------
            # STRICT FILTER
            # -----------------------
            df = df[df['Adjusted P-value'] < 0.01].copy()

            df['Gene Count'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
            df['Term Size'] = df['Overlap'].apply(lambda x: int(x.split('/')[1]))

            df = df[df['Gene Count'] >= 3]

            df['-log10(p)'] = -np.log10(df['Adjusted P-value'].replace(0, 1e-10))

            # -----------------------
            # ENRICHMENT FACTOR
            # -----------------------
            total_genes = 20000
            df['Expected'] = (df['Term Size'] / total_genes) * len(genes)
            df['Enrichment Factor'] = df['Gene Count'] / df['Expected']

            df = df[df['Enrichment Factor'] >= 1.5]

            if df.empty:
                continue

            # -----------------------
            # TOP N TERMS
            # -----------------------
            df = df.sort_values('-log10(p)', ascending=False).head(top_n)

            # Save intermediate
            df.to_csv(
                os.path.join(temp_dir, f"enrichment_raw_{lib}.csv"),
                index=False
            )

            # -----------------------
            # BUBBLE PLOT
            # -----------------------
            plt.figure(figsize=(8, 6))

            plot_df = df.sort_values('-log10(p)', ascending=False)

            scatter = plt.scatter(
                x=plot_df['Enrichment Factor'],
                y=plot_df['Term'],
                s=plot_df['Gene Count'] * 20,   # bubble size
                c=plot_df['-log10(p)'],         # color = significance
                cmap='viridis',
                alpha=0.8
            )

            plt.colorbar(scatter, label='-log10(Adjusted P-value)')

            plt.xlabel("Enrichment Factor")
            plt.ylabel("Terms")
            plt.title(f"Enrichment Bubble Plot - {lib}")

            plt.grid(True)

            plt.tight_layout()

            plt.savefig(
                os.path.join(output_dir, f"enrichment_bubble_{lib}.png"),
                dpi=300
            )

            plt.close()

            all_results.append(df)

        except Exception as e:
            logger.warning(f"{lib} failed: {e}")

    # -----------------------
    # FINAL SUMMARY
    # -----------------------
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        final_df.to_csv(
            os.path.join(output_dir, "enrichment_summary.csv"),
            index=False
        )

    end_time = datetime.now()
    logger.info(
        f"[Enrichment Phase Completed] "
        f"{start_time.strftime('%H:%M')}–{end_time.strftime('%H:%M')}"
    )