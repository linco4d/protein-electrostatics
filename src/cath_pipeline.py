#!/usr/bin/env python3
import os
import argparse
import random
import csv
import math
from pathlib import Path
from collections import defaultdict

CATH_DOMAIN_LIST_URL = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt"
CATH_SF_LIST_URL     = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-superfamily-list.txt"

# ----------------------------
# Utilities: download & I/O
# ----------------------------
def download_if_missing(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    os.system(f"wget -c -q '{url}' -O '{dest}'")

def write_tsv(rows, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        header = ["domain_id","pdb_id","chain_id","c","a","t","h","cat","cath"]
        with open(out_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
        return
    header = list(rows[0].keys())
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

# ----------------------------
# Parse CATH lists
# ----------------------------
def parse_cath_domain_list(filepath: Path):
    """
    Returns list of dicts with keys:
      domain_id, pdb_id, chain_id, c, a, t, h, cat, cath (c.a.t.h)
    """
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            toks = line.strip().split()
            if len(toks) < 5:
                continue
            domain_id = toks[0]
            c, a, t, h = toks[1:5]
            pdb_id  = domain_id[:4]
            chain_id = domain_id[4] if len(domain_id) >= 5 and domain_id[4].isalpha() else ""
            cat = f"{c}.{a}.{t}"
            cath = f"{c}.{a}.{t}.{h}"
            rows.append({
                "domain_id": domain_id,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "c": c, "a": a, "t": t, "h": h,
                "cat": cat,
                "cath": cath
            })
    return rows

def parse_cath_superfamily_list(filepath: Path):
    """
    Optional: map superfamilies to descriptions.
    Returns dict: { cath_code: {"c":..., "a":..., "t":..., "h":..., "cat":..., "cath":..., "description":...} }
    """
    mapping = {}
    with open(filepath, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split(None, 6)
            if len(toks) < 4:
                continue
            c, a, t, h = toks[0:4]
            desc = " ".join(toks[4:]).strip() if len(toks) >= 5 else ""
            cat = f"{c}.{a}.{t}"
            cath = f"{c}.{a}.{t}.{h}"
            mapping[cath] = {"c": c, "a": a, "t": t, "h": h, "cat": cat, "cath": cath, "description": desc}
    return mapping

# ----------------------------
# Writers: master & per-CAT (full domain dump)
# ----------------------------
def write_master_and_per_cat(rows, out_dir: Path):
    write_tsv(rows, out_dir / "master_domains.tsv")
    by_cat = defaultdict(list)
    for r in rows:
        by_cat[r["cat"]].append(r)
    for cat, cat_rows in by_cat.items():
        safe = cat.replace(".", "_")
        write_tsv(cat_rows, out_dir / "by_cat" / f"{safe}.tsv")

# ----------------------------
# Old behavior: 1 random domain per C.A.T.H
# ----------------------------
def sample_one_per_superfamily(rows, seed: int = 7):
    random.seed(seed)
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["cath"]].append(r)
    reps = []
    for cath_code, items in grouped.items():
        reps.append(random.choice(items))
    return reps

# ----------------------------
# NEW (1/4): group by superfamily + unique PDBs
# ----------------------------
def group_by_superfamily(rows):
    cath_to_domains = defaultdict(list)
    cath_to_pdbs_set = defaultdict(set)
    for r in rows:
        cath = r["cath"]
        pdb = r["pdb_id"].upper()
        cath_to_domains[cath].append(r)
        if pdb:
            cath_to_pdbs_set[cath].add(pdb)
    cath_to_pdbs = {k: sorted(list(v)) for k, v in cath_to_pdbs_set.items()}
    return cath_to_domains, cath_to_pdbs

# ----------------------------
# NEW (2/4): deterministic even spacing
# ----------------------------
def pick_evenly_spaced(items, k):
    n = len(items)
    if k >= n:
        return items[:]
    if k <= 1:
        return [items[0]]
    idxs = [round(i * (n - 1) / (k - 1)) for i in range(k)]
    seen = set()
    picked = []
    for idx in idxs:
        j = idx
        while j < n and j in seen:
            j += 1
        if j < n:
            seen.add(j)
            picked.append(items[j])
            continue
        j = idx - 1
        while j >= 0 and j in seen:
            j -= 1
        if j >= 0:
            seen.add(j)
            picked.append(items[j])
    return picked[:k]

# ----------------------------
# NEW (3/4): reps = max(1, floor(pdb_count/50))
# ----------------------------
def select_reps_by_pdb_fraction(cath_to_pdbs, method="even", seed=7):
    if method == "random":
        random.seed(seed)

    rep_rows = []
    def sort_key(cath_code: str):
        try:
            return tuple(int(x) for x in cath_code.split("."))
        except Exception:
            return cath_code

    for cath_code in sorted(cath_to_pdbs.keys(), key=sort_key):
        pdbs = cath_to_pdbs[cath_code]
        pdb_count = len(pdbs)
        if pdb_count == 0:
            reps_needed = 0
            reps = []
        else:
            reps_needed = max(1, pdb_count // 50)
            reps_needed = min(reps_needed, pdb_count)
            reps = pick_evenly_spaced(pdbs, reps_needed) if method == "even" else sorted(random.sample(pdbs, reps_needed))
        rep_rows.append({
            "superfamily": cath_code,  # C.A.T.H
            "pdb_count": pdb_count,
            "reps_needed": reps_needed,
            "representatives": ",".join(reps),
            "representatives_list": reps,
        })
    return rep_rows

# ----------------------------
# NEW (4/4): write reps wide + flat
# ----------------------------
def write_representatives_outputs(rep_rows, out_dir: Path, basename="reps_by_pdb_fraction"):
    out_dir.mkdir(parents=True, exist_ok=True)
    wide = out_dir / f"{basename}.tsv"
    flat = out_dir / f"{basename}_flat.tsv"

    with open(wide, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["superfamily", "pdb_count", "reps_needed", "representatives"])
        for r in rep_rows:
            w.writerow([r["superfamily"], r["pdb_count"], r["reps_needed"], r["representatives"]])

    with open(flat, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["superfamily", "pdb_id"])
        for r in rep_rows:
            for pid in r["representatives_list"]:
                w.writerow([r["superfamily"], pid])

    print(f"[done] Reps (wide): {wide}")
    print(f"[done] Reps (flat): {flat}")

# ----------------------------
# NEW: write per-CAT files that ONLY contain the chosen representatives
# ----------------------------
def write_by_cat_reps_from_rep_rows(rep_rows, out_dir: Path):
    """
    Create outputs/by_cat/<C_A_T>.tsv that ONLY list the chosen representatives.
    Each row: cat, superfamily, pdb_id
    """
    by_cat_rows = defaultdict(list)
    for r in rep_rows:
        cath = r["superfamily"]         # C.A.T.H
        cat = ".".join(cath.split(".")[:3])  # C.A.T
        for pid in r["representatives_list"]:
            by_cat_rows[cat].append({
                "cat": cat,
                "superfamily": cath,
                "pdb_id": pid
            })

    by_cat_dir = out_dir / "by_cat"
    by_cat_dir.mkdir(parents=True, exist_ok=True)
    for cat, rows in by_cat_rows.items():
        safe = cat.replace(".", "_")
        path = by_cat_dir / f"{safe}.tsv"
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["cat","superfamily","pdb_id"], delimiter="\t")
            w.writeheader()
            for row in rows:
                w.writerow(row)
    print(f"[done] Reps per CAT TSVs: {by_cat_dir}/<C_A_T>.tsv")

# ----------------------------
# Optional structure prep (kept)
# ----------------------------
def extract_pdb_ids_from_tsv(tsv_path: Path, out_path: Path):
    pdbs = []
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "pdb_id" in row and row["pdb_id"]:
                pdbs.append(row["pdb_id"].upper())
            elif "representatives" in row and row["representatives"]:
                for pid in row["representatives"].split(","):
                    pid = pid.strip().upper()
                    if pid:
                        pdbs.append(pid)
    pdbs = sorted(set(pdbs))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for pid in pdbs:
            f.write(pid + "\n")
    return len(pdbs)

# ----------------------------
# CLI
# ----------------------------
def main():
    p = argparse.ArgumentParser(
        description="Build TSVs for every CATH domain, and select PDB representatives per superfamily."
    )
    p.add_argument("--workdir", default="cath_cache", help="Where to store downloaded inputs and outputs")
    p.add_argument("--domain-list", default="", help="Path to cath-domain-list.txt (optional; will download if missing)")
    p.add_argument("--sf-list", default="", help="Path to cath-superfamily-list.txt (optional; will download if missing)")
    p.add_argument("--seed", type=int, default=7, help="Random seed for sampling")
    p.add_argument("--filter-cat", default="", help="Optional C.A.T filter like '3.40.50' to also emit a dedicated TSV")
    p.add_argument("--extract-pdb-from", default="", help="If set, path to a TSV to extract PDB IDs into a .txt list")
    p.add_argument("--reps-by-pdb-fraction", action="store_true",
                   help="Enable reps per superfamily by rule: max(1, floor(#unique PDBs / 50)).")
    p.add_argument("--rep-mode", choices=["even", "random"], default="even",
                   help="Representative selection method from unique PDBs per superfamily.")
    args = p.parse_args()

    work = Path(args.workdir)
    inputs = work / "inputs"
    outputs = work / "outputs"
    inputs.mkdir(parents=True, exist_ok=True)
    outputs.mkdir(parents=True, exist_ok=True)

    domain_list = Path(args.domain_list) if args.domain_list else inputs / "cath-domain-list.txt"
    sf_list     = Path(args.sf_list)     if args.sf_list     else inputs / "cath-superfamily-list.txt"

    download_if_missing(CATH_DOMAIN_LIST_URL, domain_list)
    download_if_missing(CATH_SF_LIST_URL, sf_list)

    rows = parse_cath_domain_list(domain_list)

    # (Optional) full domain dump + per-CAT full files (not just reps)
    write_master_and_per_cat(rows, outputs)
    if args.filter_cat:
        cat_rows = [r for r in rows if r["cat"] == args.filter_cat]
        safe = args.filter_cat.replace(".", "_")
        write_tsv(cat_rows, outputs / "by_cat" / f"{safe}.tsv")

    # Legacy: one domain per C.A.T.H (not used for reps)
    reps_domains = sample_one_per_superfamily(rows, seed=args.seed)
    write_tsv(reps_domains, outputs / "one_rep_per_superfamily.tsv")

    # Reps per superfamily by unique PDB fraction rule
    if args.reps_by_pdb_fraction:
        cath_to_domains, cath_to_pdbs = group_by_superfamily(rows)
        rep_rows = select_reps_by_pdb_fraction(cath_to_pdbs, method=args.rep_mode, seed=args.seed)
        write_representatives_outputs(rep_rows, outputs, basename="reps_by_pdb_fraction")
        # <<< NEW: by_cat files that contain ONLY the chosen representatives >>>
        write_by_cat_reps_from_rep_rows(rep_rows, outputs)

    print(f"[done] Master TSV:           {outputs/'master_domains.tsv'}")
    print(f"[done] Per-CAT TSVs:         {outputs/'by_cat'/'<C_A_T>.tsv'}")
    print(f"[done] One-per-superfamily:  {outputs/'one_rep_per_superfamily.tsv'}")

if __name__ == "__main__":
    main()
