import os
import logging
import joblib

class Predictor:
    FEATURE_COLS = [
        "TargetScan miRNA predictions_CircRNA Start",
        "TargetScan miRNA predictions_CircRNA End",
        "TargetScan miRNA predictions_local AU",
        "TargetScan miRNA predictions_TA",
        "TargetScan miRNA predictions_context+ score percentile",
    ]

    def __init__(self, model_file="trained models/calibrated_catboost_site_type_model.pkl",
                 encoder_file="trained models/label_encoder.pkl",
                 scaler_file="trained models/robust_scaler.pkl"):
        self.logger = logging.getLogger()
        
        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file not found: {model_file}")
        try:
            self.model = joblib.load(model_file)
            self.logger.info("[INFO] Loaded model: %s", model_file)
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {model_file}: {e}")

        
        if os.path.exists(encoder_file):
            try:
                self.encoder = joblib.load(encoder_file)
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
                self.scaler = joblib.load(scaler_file)
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

        if isinstance(features, pd.DataFrame):
            X = features.reindex(columns=self.FEATURE_COLS).fillna(0).values
        else:
            X = features

        if self.scaler is not None:
            try:
                X = self.scaler.transform(X)
            except Exception as e:
                self.logger.warning("[WARN] Scaler transform failed: %s", e)
        return X

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