#!/usr/bin/env python3
import argparse, csv, os, re, sys
from typing import Tuple, Optional, Dict, Set

# Optional: only needed for FASTA extraction
try:
    import gemmi  # pip install gemmi
except Exception:
    gemmi = None

ROSSMANN_PREFIX = "3.40.50"
ROSSMANN_SUPERFAM = "3.40.50.720"
CATH_CODE_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

def parse_cath_domain_id(dom_id: str) -> Tuple[str, str, str]:
    """'1abcA01' -> ('1ABC','A','01'). Chain can be multi-char."""
    pdb = dom_id[:4].upper()
    tail = dom_id[4:]
    m = re.match(r"([A-Za-z]+)(\d.*)", tail)
    if not m:
        letters = ''.join([c for c in tail if c.isalpha()])
        digits  = tail[len(letters):] or "00"
        chain = letters or "_"
        return pdb, chain, digits
    chain, digits = m.group(1), m.group(2)
    return pdb, chain, digits

def parse_cath_domain_list_clf2(path: str):
    """
    CLF 2.0 format: domain_id, then numeric columns:
      Class, Architecture, Topology, Homologous-superfamily, ... (more cols)
    Build 'cath_code' as 'C.A.T.H' from tokens[1:5].
    Yields: (domain_id, cath_code, pdb_id, chain_id)
    """
    with open(path, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            toks = line.strip().split()
            if len(toks) < 5:
                continue
            domain_id = toks[0]
            # Build dotted code from the next 4 numeric columns
            try:
                c, a, t, h = toks[1], toks[2], toks[3], toks[4]
                # ensure numeric (robust to odd whitespace)
                int(c); int(a); int(t); int(h)
                cath_code = ".".join([c, a, t, h])
            except Exception:
                continue
            pdb_id, chain_id, _ = parse_cath_domain_id(domain_id)
            yield domain_id, cath_code, pdb_id, chain_id

def parse_superfamily_names(path: Optional[str]) -> Dict[str, str]:
    """
    Parse cath-superfamily-list.txt into {cath_code: name}.
    File header: 'CATH_ID\tS35_REPS\tDOMAINS\tNAME'
    Skip header by requiring code to look like d.d.d.d
    """
    mapping = {}
    if not path or not os.path.exists(path):
        return mapping
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            toks = re.split(r'\s+', s, maxsplit=3)
            if not toks:
                continue
            code = toks[0]
            if not CATH_CODE_RE.match(code):
                continue  # skip header or malformed
            name = toks[3].strip() if len(toks) >= 4 else ""
            mapping[code] = name
    return mapping

def write_tsv(rows, out_tsv: str):
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "pdb_id","chain_id","domain_id","cath_code",
            "is_rossmann","is_high_conf_superfamily","superfamily_name"
        ])
        for r in rows:
            w.writerow(r)

def extract_chain_sequence_from_mmcif(mmcif_path: str, chain_id: str) -> Optional[str]:
    if gemmi is None:
        return None
    try:
        doc = gemmi.cif.read(mmcif_path)
        st = gemmi.make_structure_from_block(doc.sole_block())
        for model in st:
            target = None
            for c in model:
                if c.name == chain_id:
                    target = c
                    break
            if target is None:
                continue
            seq = []
            for res in target.get_polymer():
                aa = gemmi.find_aa(res)
                if aa is not None:
                    seq.append(gemmi.one_letter_code(aa))
            if seq:
                return "".join(seq)
    except Exception:
        return None
    return None

def build_fasta(hit_pairs: Set[Tuple[str,str]], mmcif_dir: str, out_fasta: str):
    os.makedirs(os.path.dirname(out_fasta) or ".", exist_ok=True)
    n_written = 0
    with open(out_fasta, "w") as fout:
        for (pdb_id, chain_id) in sorted(hit_pairs):
            candidates = [
                os.path.join(mmcif_dir, f"{pdb_id.lower()}.mmcif"),
                os.path.join(mmcif_dir, f"{pdb_id.lower()}.cif"),
                os.path.join(mmcif_dir, f"{pdb_id.upper()}.mmcif"),
                os.path.join(mmcif_dir, f"{pdb_id.upper()}.cif"),
            ]
            mmcif_path = next((p for p in candidates if os.path.exists(p)), None)
            if mmcif_path is None:
                continue
            seq = extract_chain_sequence_from_mmcif(mmcif_path, chain_id)
            if not seq:
                continue
            fout.write(f">{pdb_id}_{chain_id}\n")
            for i in range(0, len(seq), 80):
                fout.write(seq[i:i+80] + "\n")
            n_written += 1
    print(f"Wrote {n_written} chain FASTAs to {out_fasta}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cath-domain-list", required=True, type=str,
                    help="Path to cath-domain-list.txt (CLF 2.0)")
    ap.add_argument("--cath-superfamily-list", type=str, default=None,
                    help="Optional path to cath-superfamily-list.txt (for names)")
    ap.add_argument("--out-tsv", required=True, type=str,
                    help="Output TSV path for Rossmann hits")
    ap.add_argument("--mmcif-dir", type=str, default=None,
                    help="Directory of local mmCIFs; if given with --out-fasta, emits FASTA for hit chains")
    ap.add_argument("--out-fasta", type=str, default=None,
                    help="Optional FASTA path for hit chains (requires --mmcif-dir and gemmi)")
    ap.add_argument("--restrict-to-local", action="store_true",
                    help="If set, TSV includes only entries whose mmCIF file exists in --mmcif-dir")
    args = ap.parse_args()

    sf_names = parse_superfamily_names(args.cath_superfamily_list)

    if not os.path.exists(args.cath_domain_list):
        print(f"ERROR: CATH domain list not found: {args.cath_domain_list}")
        sys.exit(1)

    def entry_available(pdb_id: str) -> bool:
        if not args.mmcif_dir:
            return True
        for fn in (f"{pdb_id.lower()}.mmcif", f"{pdb_id.lower()}.cif",
                   f"{pdb_id.upper()}.mmcif", f"{pdb_id.upper()}.cif"):
            if os.path.exists(os.path.join(args.mmcif_dir, fn)):
                return True
        return False

    rows = []
    hit_pairs: Set[Tuple[str,str]] = set()

    # Parse CLF 2.0 & filter
    total = 0
    class3 = 0
    for domain_id, cath_code, pdb_id, chain_id in parse_cath_domain_list_clf2(args.cath_domain_list):
        total += 1
        if cath_code.startswith("3."):
            class3 += 1
        if not cath_code.startswith(ROSSMANN_PREFIX + "."):
            continue
        if args.restrict_to_local and not entry_available(pdb_id):
            continue
        is_high = int(cath_code == ROSSMANN_SUPERFAM)
        name = sf_names.get(cath_code, "")
        rows.append([pdb_id, chain_id, domain_id, cath_code, 1, is_high, name])
        hit_pairs.add((pdb_id, chain_id))

    write_tsv(rows, args.out_tsv)
    print(f"Scanned {total} domains (class 3: {class3}).")
    print(f"Wrote {len(rows)} Rossmann domain rows to {args.out_tsv} "
          f"({len(hit_pairs)} unique chains).")

    if args.out_fasta:
        if args.mmcif_dir is None:
            print("WARNING: --out-fasta provided but --mmcif-dir is missing; skipping FASTA.")
        elif gemmi is None:
            print("WARNING: gemmi not installed; run: pip install gemmi")
        else:
            build_fasta(hit_pairs, args.mmcif_dir, args.out_fasta)

if __name__ == "__main__":
    main()
