import pandas as pd

class FileLoader:
    def __init__(self, circ_path, mirna_path, deg_path):
        self.circ_path = circ_path
        self.mirna_path = mirna_path
        self.deg_path = deg_path

        self.circs = self._load_ids(self.circ_path, "circRNA")
        self.mirnas = self._load_ids(self.mirna_path, "miRNA")
        self.degs = self._load_degs(self.deg_path)

    def _load_ids(self, filepath, kind):
        try:
            with open(filepath, "r") as f:
                ids = {line.strip() for line in f if line.strip()}
            if not ids:
                raise ValueError(f"No {kind} IDs in {filepath}")
            return ids
        except Exception as e:
            print(f"[ERROR]  {kind} file error: {e}")
            raise

    def _load_degs(self, filepath):
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                degs = {line.strip().upper() for line in f if line.strip()}
            if not degs:
                raise ValueError(f"No DEGs in {filepath}")
            return degs
        except Exception as e:
            print(f"[ERROR]  DEG file error: {e}")
            raise

    def get_circs(self):
        return self.circs

    def get_mirnas(self):
        return self.mirnas

    def get_degs(self):
        return self.degs