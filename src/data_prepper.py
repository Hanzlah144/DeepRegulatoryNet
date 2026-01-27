import pandas as pd

class DataPrepper:
    def __init__(self):
        self.cols_to_keep = [
            "TargetScan miRNA predictions_CircRNA Mirbase ID",
            "TargetScan miRNA predictions_CircRNA (Top) - miRNA (Bottom) pairing",
            "TargetScan miRNA predictions_CircRNA Start",
            "TargetScan miRNA predictions_CircRNA End",
            "TargetScan miRNA predictions_3' pairing",
            "TargetScan miRNA predictions_local AU",
            "TargetScan miRNA predictions_position",
            "TargetScan miRNA predictions_TA",
            "TargetScan miRNA predictions_SPS",
            "TargetScan miRNA predictions_context+ score",
            "TargetScan miRNA predictions_context+ score percentile"
        ]

        self.numeric_cols = [
            "TargetScan miRNA predictions_CircRNA Start",
            "TargetScan miRNA predictions_CircRNA End",
            "TargetScan miRNA predictions_3' pairing",
            "TargetScan miRNA predictions_local AU",
            "TargetScan miRNA predictions_position",
            "TargetScan miRNA predictions_TA",
            "TargetScan miRNA predictions_SPS",
            "TargetScan miRNA predictions_context+ score",
            "TargetScan miRNA predictions_context+ score percentile"
        ]

    def clean(self, data, encoder):
        if data is None or data.empty:
            return None, None

        def split_line(text):
            text = text.replace("\xa0", " ")
            pieces = text.split()
            circ = pieces[0]
            mirna = "none"
            for piece in pieces:
                if "hsa-miR-" in piece:
                    mirna = piece
                    break
            return circ, mirna

        try:
            circ_ids, mirna_ids = zip(*[split_line(row) for row in data["TargetScan miRNA predictions_CircRNA Mirbase ID"]])
            data["circ_id"] = circ_ids
            data["mirna_id"] = mirna_ids
        except KeyError as e:
            return None, None

        # Convert to numeric
        for col in self.numeric_cols:
            data[col] = pd.to_numeric(data[col], errors="coerce")

        # Extract and fill features
        features = data[self.cols_to_keep].copy()
        features["TargetScan miRNA predictions_CircRNA Mirbase ID"] = features["TargetScan miRNA predictions_CircRNA Mirbase ID"].fillna("missing")
        features["TargetScan miRNA predictions_CircRNA (Top) - miRNA (Bottom) pairing"] = features["TargetScan miRNA predictions_CircRNA (Top) - miRNA (Bottom) pairing"].fillna("missing")
        for col in self.numeric_cols:
            features[col] = features[col].fillna(0)

        return features, data