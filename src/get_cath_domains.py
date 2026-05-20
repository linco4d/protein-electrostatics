#!/usr/bin/env python3
"""
CATH domain PDB bulk downloader (Option A).
- Input: CATH domain IDs (e.g., 1cukA01)
- Output: one PDB file per domain in ./cath_domains_pdb/
"""

import os
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ---------- CONFIG ----------
API_BASE = "https://www.cathdb.info/version/latest/api/rest/id"
OUT_DIR = Path("cath_domains_pdb")   # where .pdb files will be saved
DEFAULT_MAX_WORKERS = 16             # tune for your machine/network
DEFAULT_TIMEOUT = 20                 # seconds for each GET
# -----------------------------


def iter_domain_ids_from_cath_domain_list(path: str) -> Iterable[str]:
    """
    Stream domain IDs from cath-domain-list.txt.
    Assumes first token on each non-comment line is the domain ID.
    """
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            toks = line.split()
            if not toks:
                continue
            yield toks[0]  # e.g., '1cukA01'


def make_session() -> requests.Session:
    """
    Create a requests Session with connection pooling and robust retries.
    Retries handle transient HTTP errors (429, 5xx) with backoff.
    """
    s = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1.2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET"},
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=100, pool_maxsize=100)
    s.mount("https://", adapter)
    s.headers.update({"User-Agent": "cath-bulk-downloader/1.0"})
    return s


def looks_like_pdb(text: str) -> bool:
    """
    Very light validation that the response is a plausible PDB.
    Filters HTML error pages / tiny bodies.
    """
    if len(text) < 200:
        return False
    head = text[:800]
    if "DOCTYPE html" in head:
        return False
    return any(tag in head for tag in ("HEADER", "TITLE", "ATOM", "HETATM", "MODEL"))


def fetch_one(session: requests.Session, domain_id: str, timeout: int = DEFAULT_TIMEOUT) -> Tuple[str, str]:
    """
    Download one domain PDB, return (domain_id, status).
    Status is one of: 'ok', 'skip', 'bad-format', 'http-<code>', 'error:<msg>'
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{domain_id}.pdb"

    # Skip if already present (simple resume)
    if out_path.exists() and out_path.stat().st_size > 0:
        return domain_id, "skip"

    url = f"{API_BASE}/{domain_id}.pdb"
    try:
        r = session.get(url, timeout=timeout)
        if r.status_code == 200:
            text = r.text
            if looks_like_pdb(text):
                out_path.write_text(text)
                return domain_id, "ok"
            else:
                return domain_id, "bad-format"
        else:
            return domain_id, f"http-{r.status_code}"
    except Exception as e:
        return domain_id, f"error:{e}"


def download_domains(
    domain_ids,
    max_workers: int = DEFAULT_MAX_WORKERS,
    limit: int | None = None,
) -> None:
    """
    Orchestrate parallel downloads with a progress bar (if tqdm is available).
    """
    ids = list(domain_ids)
    if limit:
        ids = ids[:limit]

    session = make_session()
    successes = 0
    skips = 0
    fails: list[tuple[str, str]] = []

    # Optional: pretty progress if tqdm is present
    try:
        from tqdm import tqdm
        def progress(iterable): return tqdm(iterable, total=len(ids), desc="Downloading domains")
    except Exception:
        def progress(iterable): return iterable

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(fetch_one, session, did): did for did in ids}
        for fut in progress(as_completed(futures)):
            did, status = fut.result()
            if status == "ok":
                successes += 1
            elif status == "skip":
                skips += 1
            else:
                fails.append((did, status))

    print(f"\nDone. ok={successes}, skipped={skips}, failed={len(fails)}")
    if fails:
        print("Sample failures (first 20):")
        for did, st in fails[:20]:
            print(f"  {did}: {st}")


def main():
    """
    CLI usage:
      python get_cath_domains.py <optional:cath-domain-list.txt> <optional:max_workers> <optional:limit>

    Examples:
      python get_cath_domains.py cath-domain-list.txt 24
      python get_cath_domains.py cath-domain-list.txt 24 10000
      # If you already have `domain_data` in memory (same Python session), run:
      #   download_domains([row[0] for row in domain_data], max_workers=24)
    """
    # Parse CLI args
    path = sys.argv[1] if len(sys.argv) > 1 else None
    max_workers = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_MAX_WORKERS
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else None

    if path:
        # Standalone mode: read IDs from the .txt file
        domain_ids = iter_domain_ids_from_cath_domain_list(path)
        download_domains(domain_ids, max_workers=max_workers, limit=limit)
    else:
        # In-session mode: expects `domain_data` (from your parse) to exist in globals
        if "domain_data" not in globals():
            raise SystemExit(
                "No filepath provided and `domain_data` not found. "
                "Run with: python get_cath_domains.py cath-domain-list.txt"
            )
        ids = [row[0] for row in domain_data]
        download_domains(ids, max_workers=max_workers, limit=limit)


if __name__ == "__main__":
    main()

