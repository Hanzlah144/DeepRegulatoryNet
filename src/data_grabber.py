import os
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from io import StringIO

class DataGrabber:
    def __init__(self, save_dir="my_output"):
        self.base_url = "https://circinteractome.nia.nih.gov/api/v2/mirnasearch"
        self.save_dir = save_dir
        self.session = requests.Session()

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

    def fetch(self, circ_id):
        save_path = os.path.join(self.save_dir, f"{circ_id}_targets.xlsx")

        
        if os.path.exists(save_path):
            try:
                df = pd.read_excel(save_path)
                return df
            except Exception:
                pass

        params = {"circular_rna_query": circ_id}
        try:
            resp = self.session.get(self.base_url, params=params, verify=False, timeout=10)
            if resp.status_code != 200:
                return None

            soup = BeautifulSoup(resp.text, "html.parser")
            table = soup.find("table", {"border": "1", "bordercolor": "#006699"})
            if table:
                df = pd.read_html(StringIO(str(table)))[0]
                if isinstance(df.columns, pd.MultiIndex):
                    df.columns = [f"{col[0]}_{col[1]}" for col in df.columns]
                df.to_excel(save_path, index=False)
                return df

            return None

        except Exception:
            return None