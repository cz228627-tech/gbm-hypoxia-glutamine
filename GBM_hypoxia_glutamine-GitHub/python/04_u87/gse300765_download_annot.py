#!/usr/bin/env python3
"""
Download GPL10558 annotation and extract probe-to-gene mapping
for the 13 fate genes. Uses streaming to avoid timeout.
"""
import urllib.request
import gzip
import io
import csv
import pandas as pd
from pathlib import Path
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
OUT = Path(os.environ.get("GBM_U87_DIR", DATA_DIR / "gse300765")).expanduser()
OUT.mkdir(parents=True, exist_ok=True)

# Genes we need
FATE_GENES = ['GLUL','GLS','GLS2','SLC1A5','SLC38A1','SLC38A2',
              'SLC7A5','SLC7A11','GLUD1','GOT1','GOT2','GPT','GPT2']

# GPL10558 text format from GEO
url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558&targ=GPL&form=text&view=data"

print("Downloading GPL10558 annotation (streaming)...")
# The file is large (~100MB), so we stream and filter
probe_to_gene = {}
header = None
gene_col_idx = None
id_col_idx = None
symbol_col_idx = None

req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req, timeout=300) as response:
    for line in response:
        line = line.decode('utf-8', errors='ignore').strip()
        if line.startswith('#') and '=' in line:
            # Header line like: #ID = ...
            parts = line.lstrip('#').split(' = ', 1)
            if len(parts) == 2:
                if header is None:
                    header = []
                header.append(parts[0].strip())
            continue
        if line.startswith('^PLATFORM'):
            continue
        if header and line:
            # Data line
            fields = line.split('\t')
            if id_col_idx is None and header:
                try:
                    id_col_idx = header.index('ID')
                except ValueError:
                    id_col_idx = 0
                # Try different column names for gene symbol
                for col_name in ['GeneSymbol', 'Symbol', 'Gene_Symbol', 'ILMN_GENE']:
                    try:
                        symbol_col_idx = header.index(col_name)
                        break
                    except ValueError:
                        continue
                if symbol_col_idx is None:
                    # Print header for debugging
                    print(f"Header columns: {header[:15]}...")
                    # Try to find any column with 'gene' or 'symbol'
                    for i, h in enumerate(header):
                        if 'gene' in h.lower() or 'symbol' in h.lower():
                            print(f"  Found potential gene column: {h} (index {i})")
                            symbol_col_idx = i
                            break
                    if symbol_col_idx is None:
                        print("  No gene symbol column found!")
                        break

            if id_col_idx is not None and symbol_col_idx is not None:
                if len(fields) > max(id_col_idx, symbol_col_idx):
                    probe_id = fields[id_col_idx]
                    symbol = fields[symbol_col_idx]
                    if symbol and symbol != '':
                        # Check if this symbol is one we need
                        if symbol in FATE_GENES:
                            probe_to_gene[probe_id] = symbol

        if len(probe_to_gene) > 0 and line == '' and header:
            # Might be end of data
            pass

print(f"\nFound {len(probe_to_gene)} probe-to-gene mappings for fate genes:")
for probe, gene in sorted(probe_to_gene.items(), key=lambda x: x[1]):
    print(f"  {probe} -> {gene}")

# Save mapping
df_map = pd.DataFrame(list(probe_to_gene.items()), columns=['ProbeID', 'GeneSymbol'])
df_map.to_csv(OUT / "probe_to_gene_mapping.csv", index=False)
print(f"\nSaved to {OUT / 'probe_to_gene_mapping.csv'}")
