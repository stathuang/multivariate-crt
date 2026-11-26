#!/usr/bin/env python3
"""Explicitly fetch the pinned public inputs, then verify every input hash."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
PIN = "35a3647910792d39e27986b05e07b04b0a6f6461"
BASE = f"https://raw.githubusercontent.com/xiaotong0201/GGMtest/{PIN}/GGM-Test-Version0/Application/"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2034/matrix/GSE2034_series_matrix.txt.gz"


def verify(data, item):
    if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
        raise ValueError(f"Snapshot differs: {item['file']}; the existing file was not replaced.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--check", action="store_true", help="Only verify existing files; no network")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--dia-source", type=Path, help="Exact 2026-10-23 DIA.json supplied locally")
    group.add_argument("--dia-url", help="Author's eventual public URL for that exact DIA.json")
    args = parser.parse_args()
    inputs = json.loads((ROOT / "data/required-inputs.json").read_text())
    for item in inputs:
        name = item["file"]
        dest = args.directory / name
        if dest.exists():
            verify(dest.read_bytes(), item)
            print(f"Verified {name}", flush=True)
            continue
        if args.check:
            raise FileNotFoundError(f"Missing {name}")
        if name == "Stock/DIA.json" and args.dia_source:
            data = args.dia_source.read_bytes()
        else:
            url = (args.dia_url if name == "Stock/DIA.json" else
                   GEO if name == "Breast/GSE2034_series_matrix.txt" else BASE + name)
            if not url:
                raise ValueError("The six public inputs are available. Supply the exact DIA.json with "
                                 "--dia-source or --dia-url; no substitute prices will be downloaded.")
            if not url.startswith("https://"):
                raise ValueError("Downloads require an HTTPS URL")
            print(f"Downloading {name}", flush=True)
            request = urllib.request.Request(url, headers={"User-Agent": "JRSSB-replication/2"})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
            if name == "Breast/GSE2034_series_matrix.txt":
                data = gzip.decompress(data)
        verify(data, item)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("xb") as handle:
            handle.write(data)
        print(f"Saved and verified {name}", flush=True)
    print("All seven exact input snapshots verified.")


if __name__ == "__main__":
    main()
