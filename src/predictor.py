import os
import logging
import pickle

class Predictor:
    FEATURE_COLS = [
        "TargetScan miRNA predictions_CircRNA Start",
        "TargetScan miRNA predictions_CircRNA End",
        "TargetScan miRNA predictions_local AU",
        "TargetScan miRNA predictions_TA",
        "TargetScan miRNA predictions_SPS",
        "TargetScan miRNA predictions_position",
    ]

    def __init__(self,
                 model_file=None,
                 encoder_file=None,
                 scaler_file=None):
        self.logger = logging.getLogger()

        # Default to CatBoost artifacts saved in the project-level
        # "trained models" directory.
        base_dir = os.path.dirname(os.path.dirname(__file__))
        models_dir = os.path.join(base_dir, "trained models")
        if model_file is None:
            model_file = os.path.join(models_dir, "catboost_model.pkl")
        if encoder_file is None:
            encoder_file = os.path.join(models_dir, "label_encoder.pkl")
        if scaler_file is None:
            scaler_file = os.path.join(models_dir, "scaler.pkl")

        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file not found: {model_file}")
        try:
            self.model = pickle.load(open(model_file, 'rb'))
            self.logger.info("[INFO] Loaded model: %s", model_file)
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {model_file}: {e}")

        
        if os.path.exists(encoder_file):
            try:
                self.encoder = pickle.load(open(encoder_file, 'rb'))
                self.class_names = getattr(self.encoder, 'classes_', None)
                self.logger.info("[INFO] Loaded label encoder: %s", encoder_file)
            except Exception as e:
                self.encoder = None
                self.class_names = None
                self.logger.warning("[WARN] Failed to load encoder %s: %s", encoder_file, e)
        else:
            self.encoder = None
            self.class_names = None
            self.logger.warning("[WARN] Encoder file not found: %s", encoder_file)

        
        if os.path.exists(scaler_file):
            try:
                self.scaler = pickle.load(open(scaler_file, 'rb'))
                self.logger.info("[INFO] Loaded scaler: %s", scaler_file)
            except Exception as e:
                self.scaler = None
                self.logger.warning("[WARN] Failed to load scaler %s: %s", scaler_file, e)
        else:
            self.scaler = None
            self.logger.warning("[WARN] Scaler file not found: %s", scaler_file)

    def _prepare(self, features):
        import pandas as pd
        import numpy as np

        if features is None:
            return None

        # Always construct a DataFrame with explicit feature names so that
        # scikit-learn transformers (e.g. StandardScaler) fitted with
        # feature names do not emit warnings.
        if isinstance(features, pd.DataFrame):
            X_df = features.reindex(columns=self.FEATURE_COLS).fillna(0)
        else:
            X_df = pd.DataFrame(features, columns=self.FEATURE_COLS)

        if self.scaler is not None:
            try:
                X_scaled = self.scaler.transform(X_df)
            except Exception as e:
                self.logger.warning("[WARN] Scaler transform failed: %s", e)
                X_scaled = X_df.values
        else:
            X_scaled = X_df.values

        return np.asarray(X_scaled)

    def predict(self, features):
        X = self._prepare(features)
        if X is None:
            return None, None, None

        raw_preds = self.model.predict(X)
        probs = None
        if hasattr(self.model, 'predict_proba'):
            try:
                probs = self.model.predict_proba(X)
            except Exception:
                probs = None

        
        if self.encoder is not None:
            try:
                labels = self.encoder.inverse_transform(raw_preds.ravel())
            except Exception:
                labels = raw_preds
        else:
            labels = raw_preds

        return labels, probs, raw_preds