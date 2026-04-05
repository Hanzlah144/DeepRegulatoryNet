import os
import glob
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from matplotlib_venn import venn2_unweighted
import matplotlib.pyplot as plt
import logging
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import json
from pathlib import Path


warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib_venn")


CACHE_DIR = Path("temp/mirna_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

def _get_cache_path(mirna_name):
    
    return CACHE_DIR / f"{mirna_name.replace('/', '_')}.json"

def _load_from_cache(mirna_name):
    
    cache_path = _get_cache_path(mirna_name)
    if cache_path.exists():
        try:
            with open(cache_path, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return None

def _save_to_cache(mirna_name, targets):
    
    cache_path = _get_cache_path(mirna_name)
    try:
        with open(cache_path, 'w') as f:
            json.dump(targets, f)
    except Exception:
        pass


def query_mirdb_optimized(mirna_name, session=None, max_retries=2, retry_delay=1):
    
    
    cached_targets = _load_from_cache(mirna_name)
    if cached_targets is not None:
        return cached_targets
    
    if session is None:
        session = requests.Session()
        session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
    
    url = "https://mirdb.org/cgi-bin/search.cgi"
    payload = {
        "species": "Human",
        "searchBox": mirna_name,
        "searchType": "miRNA",
        "submitButton": "Go"
    }

    for attempt in range(max_retries + 1):
        try:
            response = session.post(url, data=payload, timeout=15)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')

            rows = soup.find_all('tr')
            targets = []
            for row in rows:
                cols = row.find_all('td')
                if len(cols) > 4:
                    gene_link = cols[4].find('a')
                    if gene_link:
                        targets.append(gene_link.text.strip().upper())

            targets = list(set(targets))  # Deduplicated
            _save_to_cache(mirna_name, targets)  # Cache results
            return targets

        except Exception as e:
            if attempt < max_retries:
                time.sleep(retry_delay)
            else:
                logging.getLogger().warning(f"[WARN] Failed to query miRDB for {mirna_name}: {e}")
                return []

def query_mirdb(mirna_name, max_retries=3, retry_delay=5):
    return query_mirdb_optimized(mirna_name, max_retries=max_retries, retry_delay=retry_delay)

def overlap_mrnas(deg_file, matches_dir='my_output', output_dir='.'):
    logger = logging.getLogger()
    os.makedirs(output_dir, exist_ok=True)

    mirna_files = glob.glob(os.path.join(matches_dir, '*_strong_medium_matches.csv'))
    if not mirna_files:
        raise FileNotFoundError("[ERROR] No match files found")

    mirnas = set()
    for file in mirna_files:
        df = pd.read_csv(file)
        if 'mirna_id' in df.columns:
            mirnas.update(df['mirna_id'].dropna().str.strip())
    mirnas = list(mirnas)

    with open(deg_file, 'r', encoding='utf-8') as f:
        degs = pd.Series(f.readlines()).str.strip().str.upper().drop_duplicates()
    degs = degs[degs != '']
    logger.info("[INFO] DEGs loaded: %d", len(degs))

    
    logger.info("[INFO]  Querying miRDB for %d miRNAs in parallel...", len(mirnas))
    mirna_targets = []
    
    def process_mirna(mirna):
        """Process a single miRNA query"""
        targets = query_mirdb_optimized(mirna)
        if targets:
            return [{'mirna': mirna, 'gene': g} for g in targets]
        return []
    
    
    with ThreadPoolExecutor(max_workers=min(8, len(mirnas))) as executor:
        future_to_mirna = {executor.submit(process_mirna, mirna): mirna for mirna in mirnas}
        
        for future in as_completed(future_to_mirna):
            mirna = future_to_mirna[future]
            try:
                result = future.result()
                mirna_targets.extend(result)
                logger.debug(f"[DEBUG] Processed {mirna}: {len(result)} targets")
            except Exception as e:
                logger.error(f"[ERROR] Failed to process {mirna}: {e}")
    
    logger.info("[INFO]  Completed parallel miRDB queries")

    if not mirna_targets:
        logger.error("[ERROR] No miRNA targets found")
        return pd.DataFrame()

    targets_df = pd.DataFrame(mirna_targets)
    overlapping = targets_df[targets_df['gene'].isin(degs)].drop_duplicates()
    logger.info("[INFO] Overlaps found: %d pairs | %d unique genes", len(overlapping), overlapping['gene'].nunique())

    overlap_path = os.path.join(matches_dir, "overlapping_mrnas.csv")
    overlapping.to_csv(overlap_path, index=False)
    logger.info("[INFO] Saved: %s", overlap_path)

    logger.info("--------------------------------------------------")
    logger.info("[INFO] Generating Venn diagram...")
    plt.figure(figsize=(6, 6))
    venn = venn2_unweighted(
        [set(targets_df['gene']), set(degs)],
        set_labels=(f'miRDB ({len(targets_df["gene"].unique())})', f'DEGs ({len(degs)})'),
        set_colors=("#FF4C4C", "#4C8CFF"),
        alpha=0.7
    )

    plt.title("miRDB–DEG Overlap", fontsize=12, fontweight='bold')
    venn_path = os.path.join(output_dir, "venn_overlap.png")
    plt.savefig(venn_path, dpi=400, bbox_inches='tight')
    plt.close()
    logger.info("[INFO] Venn diagram saved: %s", venn_path)
    logger.info("--------------------------------------------------")

    return overlapping
