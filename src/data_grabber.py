import os
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from io import StringIO
from requests.exceptions import Timeout, ConnectionError, RequestException


class CircInteractomeUnavailableError(RuntimeError):
    """Raised when the NIH CircInteractome server is unreachable or down."""
    pass


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
                # If cached file is corrupt, fall back to fetching from server
                pass

        params = {"circular_rna_query": circ_id}
        try:
            resp = self.session.get(
                self.base_url,
                params=params,
                verify=False,
                timeout=10,
            )
        except (Timeout, ConnectionError) as e:
            # Explicitly signal that the external NIH server is unavailable
            raise CircInteractomeUnavailableError(
                "NIH CircInteractome request timed out or the server is unreachable."
            ) from e
        except RequestException:
            # Other request-related problems (e.g., bad URL, SSL issues)
            return None

        # Treat 5xx errors as server downtime
        if 500 <= resp.status_code < 600:
            raise CircInteractomeUnavailableError(
                f"NIH CircInteractome responded with server error {resp.status_code}."
            )

        if resp.status_code != 200:
            # Non-success but not clear downtime (e.g., 400/404) → behave as no data
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