import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
from file_loader import FileLoader
from data_grabber import DataGrabber
from data_prepper import DataPrepper
from predictor import Predictor
from mrna_overlap import overlap_mrnas
from network_constructor import construct_circrna_mirna_mrna_network

class AnalysisPipeline:
    def __init__(self, circ_file, mirna_file, deg_file,
                 temp_dir="temp", output_dir="output",
                 model_file="trained models/calibrated_catboost_site_type_model.pkl",
                 encoder_file="trained models/label_encoder.pkl",
                 scaler_file="trained models/robust_scaler.pkl",
                 data_dir=None):
        self.loader = FileLoader(circ_file, mirna_file, deg_file)
        
        self.data_dir = data_dir
        self.grabber = DataGrabber(data_dir or temp_dir)
        if data_dir:
            logging.getLogger().info("[INFO] DataGrabber will use local data directory: %s", data_dir)
        self.prepper = DataPrepper()
        self.predictor = Predictor(model_file, encoder_file, scaler_file)
        self.temp_dir = temp_dir
        self.output_dir = output_dir
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

    def process_single_circ(self, circ):
        logger = logging.getLogger()
        data = self.grabber.fetch(circ)
        if data is None:
            logger.debug(f"[DEBUG] No data for {circ}")
            return None

        features, full_data = self.prepper.clean(data, self.predictor.encoder)
        if features is None:
            logger.debug(f"[DEBUG] Failed to clean data for {circ}")
            return None

        preds, probs, codes = self.predictor.predict(features)
        if preds is None:
            logger.debug(f"[DEBUG] Prediction failed for {circ}")
            return None

        full_data["predicted_site_type"] = preds
        full_data["encoded_prediction"] = codes
        for i, cname in enumerate(self.predictor.class_names):
            full_data[f"prob_{cname}"] = probs[:, i]

        mask = full_data["predicted_site_type"].isin(["7mer-m8", "8mer-1a"])
        filtered = full_data[mask].copy()

        if filtered.empty:
            logger.debug(f"[DEBUG] No strong/medium sites for {circ}")
            return None

        cols = [
            "circ_id", "mirna_id", "TargetScan miRNA predictions_Site Type",
            "TargetScan miRNA predictions_CircRNA Start", "TargetScan miRNA predictions_CircRNA End",
            "TargetScan miRNA predictions_3' pairing", "TargetScan miRNA predictions_local AU",
            "TargetScan miRNA predictions_TA", "TargetScan miRNA predictions_SPS",
            "TargetScan miRNA predictions_context+ score", "TargetScan miRNA predictions_context+ score percentile",
            "predicted_site_type", "encoded_prediction",
            "prob_7mer-1a", "prob_7mer-m8", "prob_8mer-1a"
        ]

        filtered[cols].to_csv(os.path.join(self.temp_dir, f"{circ}_strong_medium_results.csv"), index=False)
        logger.debug(f"[DEBUG] {circ}: Cleaned {len(full_data)} rows → {len(filtered)} sites")

        if logging.getLogger().level == logging.DEBUG:
            breakdown = filtered["predicted_site_type"].value_counts(normalize=True)
            sizes = [breakdown.get("7mer-m8", 0)*100, breakdown.get("8mer-1a", 0)*100]
            labels = ["7mer-m8", "8mer-1a"]
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.pie(sizes, labels=labels, autopct="%1.1f%%", colors=["#66c2a5", "#fc8d62"], startangle=90)
            ax.set_title(f"{circ} Site Type Breakdown")
            plt.savefig(os.path.join(self.temp_dir, f"{circ}_interaction_chart.png"))
            plt.close()

        return filtered

    def process_all_circs(self):
        logger = logging.getLogger()
        results = {}
        for circ in self.loader.get_circs():
            data = self.process_single_circ(circ)
            if data is not None:
                results[circ] = data
        logger.info("[INFO]  Processed %d circRNAs", len(results))
        return results

    def find_strong_hits(self, results):
        logger = logging.getLogger()
        logger.info(" STEP 2: Matching Strong/Medium Hits")
        strong_hits = {}
        all_mirnas = set()
        for circ, df in results.items():
            hits = df[df["predicted_site_type"].isin(["7mer-m8", "8mer-1a"])]["mirna_id"].unique().tolist()
            if hits:
                strong_hits[circ] = hits
                all_mirnas.update(hits)
        logger.info("[INFO]  %d circRNAs with %d total miRNAs", len(strong_hits), len(all_mirnas))
        return strong_hits, all_mirnas

    def match_mirnas(self, results, strong_hits, all_mirnas):
        logger = logging.getLogger()
        input_mirnas = self.loader.get_mirnas()
        matched = all_mirnas.intersection(input_mirnas)
        logger.info("[INFO]  %d matched miRNAs", len(matched))
        matched_data = {}
        for circ, mirnas in strong_hits.items():
            common = set(mirnas).intersection(matched)
            logger.debug(f"[DEBUG] {circ}: {len(common)} matched miRNAs")
            if common:
                df = results[circ]
                subset = df[(df["predicted_site_type"].isin(["7mer-m8", "8mer-1a"])) &
                            (df["mirna_id"].isin(common))]
                matched_data[circ] = subset
                subset.to_csv(os.path.join(self.temp_dir, f"{circ}_strong_medium_matches.csv"), index=False)
                logger.debug(f"[DEBUG] Saved: {circ}_strong_medium_matches.csv")
        return matched_data

    def analyze_mrna_overlap(self):
        logger = logging.getLogger()
        logger.info(" STEP 3: DEG–miRNA mRNA Overlap ===")
        try:
            df = overlap_mrnas(self.loader.deg_path, self.temp_dir, self.output_dir)
            if df.empty:
                logger.info("[INFO]  No overlaps found")
            else:
                logger.info("[INFO]  %d pairs, %d unique genes", len(df), df['gene'].nunique())
            return df
        except Exception as e:
            logger.error("[ERROR]  Overlap failed: %s", e)
            return pd.DataFrame()

    def construct_network(self, results, strong_hits, all_mirnas):
        logger = logging.getLogger()
        logger.info(" STEP 4: Regulatory Network Construction")
        G = construct_circrna_mirna_mrna_network(results, strong_hits, all_mirnas, self.temp_dir, self.output_dir)
        return G