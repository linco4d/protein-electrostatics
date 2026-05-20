#!/usr/bin/env python3
# fetch_rep_structures.py
import argparse
import csv
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

RCSB_PDB_URL  = "https://files.rcsb.org/download/{pdb}.pdb"
RCSB_CIF_URL  = "https://files.rcsb.org/download/{pdb}.cif"

def read_unique_pdb_ids(tsv_path: Path) -> List[str]:
    if not tsv_path.exists():
        sys.exit(f"[ERR] TSV not found: {tsv_path}")
    ids = set()
    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "pdb_id" not in reader.fieldnames:
            sys.exit(f"[ERR] TSV missing 'pdb_id' column: {tsv_path}")
        for row in reader:
            pid = (row.get("pdb_id") or "").strip()
            if pid:
                ids.add(pid.upper())
    if not ids:
        sys.exit("[ERR] No PDB IDs found in TSV.")
    return sorted(ids)

def http_get(url: str, timeout: int = 60) -> bytes:
    req = Request(url, headers={"User-Agent": "foldseek-fetch/1.0"})
    with urlopen(req, timeout=timeout) as resp:
        return resp.read()

def download_one(pdb_id: str, out_dir: Path, fmt: str, retries: int, delay: float) -> Tuple[str, str, bool, str]:
    """
    Returns (pdb_id, out_path, success, msg)
    """
    out_ext = ".pdb" if fmt == "pdb" else ".cif"
    out_path = out_dir / f"{pdb_id}{out_ext}"
    if out_path.exists() and out_path.stat().st_size > 0:
        return (pdb_id, str(out_path), True, "exists")

    url = (RCSB_PDB_URL if fmt == "pdb" else RCSB_CIF_URL).format(pdb=pdb_id)
    last_err = ""
    for attempt in range(1, retries + 1):
        try:
            blob = http_get(url)
            out_dir.mkdir(parents=True, exist_ok=True)
            with open(out_path, "wb") as w:
                w.write(blob)
            # quick sanity check
            if out_path.stat().st_size < 1000:
                last_err = "file too small (corrupt?)"
                out_path.unlink(missing_ok=True)
                raise IOError(last_err)
            return (pdb_id, str(out_path), True, "downloaded")
        except (HTTPError, URLError, IOError) as e:
            last_err = f"{type(e).__name__}: {e}"
            if attempt < retries:
                time.sleep(delay * attempt)
            else:
                return (pdb_id, str(out_path), False, last_err)
    return (pdb_id, str(out_path), False, last_err or "unknown error")

def write_manifest(manifest_path: Path, rows: List[Tuple[str, str, bool, str]], fmt: str):
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb_id", "path", "ok", "note", "format"])
        for pdb_id, path, ok, note in rows:
            w.writerow([pdb_id, path, "1" if ok else "0", note, fmt])

def main():
    ap = argparse.ArgumentParser(description="Download PDB/mmCIF for representatives in one_rep_per_superfamily.tsv")
    ap.add_argument("--tsv", required=True, help="Path to one_rep_per_superfamily.tsv")
    ap.add_argument("--out", default="structures/reps", help="Output folder for structures")
    ap.add_argument("--fmt", choices=["pdb", "cif"], default="cif", help="Download format (default: cif)")
    ap.add_argument("--retries", type=int, default=3, help="Retries per file")
    ap.add_argument("--delay", type=float, default=1.0, help="Backoff base delay (seconds)")
    ap.add_argument("--workers", type=int, default=8, help="Parallel downloads")
    args = ap.parse_args()

    tsv_path = Path(args.tsv)
    out_dir  = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_ids = read_unique_pdb_ids(tsv_path)
    print(f"[info] Found {len(pdb_ids)} unique PDB IDs in {tsv_path.name}")
    print(f"[info] Downloading as {args.fmt.upper()} into: {out_dir}")

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(download_one, pid, out_dir, args.fmt, args.retries, args.delay) for pid in pdb_ids]
        for fut in as_completed(futs):
            pdb_id, path, ok, note = fut.result()
            results.append((pdb_id, path, ok, note))
            status = "OK " if ok else "ERR"
            print(f"[{status}] {pdb_id} -> {path} ({note})")

    manifest_path = out_dir / "manifest.csv"
    write_manifest(manifest_path, results, args.fmt)
    n_ok = sum(1 for _,_,ok,_ in results if ok)
    n_er = len(results) - n_ok
    print(f"[done] Saved {n_ok} / {len(results)} files. Manifest: {manifest_path}")
    if n_er:
        print("[warn] Some downloads failed; check manifest notes and re-run with higher --retries or different --fmt.")

if __name__ == "__main__":
    main()
