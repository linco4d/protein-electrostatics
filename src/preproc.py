#!/usr/bin/env python3

"""
preproc.py

Mac/HPC-safe preprocessing for the protein electrostatics pipeline.

This script:
1. Cleans a full PDB file using PDBFixer/OpenMM when available.
2. Falls back to a lightweight ATOM-only cleaner if PDBFixer fails or is unavailable.
3. Extracts k-mer residue fragments from the cleaned PDB.
4. Writes each k-mer as its own .pdb file into a fix/ directory.
5. Optionally removes the cleaned full-length PDB after k-mer extraction.

Expected output structure:

outdir/
    <pdb_id>_3mers/
        cleaned/
            <pdb_id>_clean.pdb
        fix/
            <pdb_id>_A_1_3_000000.pdb
            <pdb_id>_A_2_4_000001.pdb
            ...
        pqr_files/
        atompot/

Typical use:

    python3 src/preproc.py \
        --pdb rep_subset/5xus.pdb \
        --outdir /data/users_bigdata/Lincoln.Aftergood.26/protein-analysis/3structs \
        --k 3 \
        --remove_raw

From another script:

    from preproc import preprocess_pdb

    result = preprocess_pdb(
        pdb_file="rep_subset/5xus.pdb",
        outdir="3structs",
        k=3,
        remove_raw=True,
    )

    print(result["fragment_files"])
"""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class ResidueID:
    """
    A stable identifier for a residue in a PDB file.

    PDB fixed-width columns:
    - chain ID: column 22, Python index 21
    - residue sequence number: columns 23-26, Python slice 22:26
    - insertion code: column 27, Python index 26
    """

    chain_id: str
    resseq: str
    icode: str
    resname: str

    def safe_label(self) -> str:
        chain = self.chain_id.strip() or "NA"
        resseq = self.resseq.strip() or "NA"
        icode = self.icode.strip()
        resname = self.resname.strip() or "UNK"

        if icode:
            return f"{chain}_{resname}{resseq}{icode}"

        return f"{chain}_{resname}{resseq}"


def atomic_write_text(path: Path, text: str) -> None:
    """
    Write text atomically.

    This prevents zero-byte or half-written files if the process is interrupted.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
    ) as tmp:
        tmp.write(text)
        tmp_path = Path(tmp.name)

    os.replace(tmp_path, path)


def is_pdb_atom_line(line: str, include_hetero: bool = False) -> bool:
    """
    Return True if line should be treated as a structural atom line.
    """
    if line.startswith("ATOM"):
        return True

    if include_hetero and line.startswith("HETATM"):
        return True

    return False


def parse_residue_id(line: str) -> ResidueID:
    """
    Parse residue identity from a PDB ATOM/HETATM line.
    """
    return ResidueID(
        chain_id=line[21:22],
        resseq=line[22:26],
        icode=line[26:27],
        resname=line[17:20],
    )


def lightweight_clean_pdb(
    pdb_file: str | Path,
    output_file: str | Path,
    include_hetero: bool = False,
) -> str:
    """
    Fallback cleaner.

    This does not repair missing atoms or residues. It simply:
    - keeps ATOM records
    - optionally keeps HETATM records
    - drops waters/solvent by default
    - appends END

    This is useful when PDBFixer is unavailable or fails on a malformed PDB.
    """
    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()
    output_path = Path(os.path.expanduser(str(output_file))).resolve()

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file does not exist: {pdb_path}")

    kept_lines: List[str] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if is_pdb_atom_line(line, include_hetero=include_hetero):
                resname = line[17:20].strip()

                # Drop common waters unless the user explicitly keeps hetero records.
                if resname in {"HOH", "WAT", "H2O"}:
                    continue

                kept_lines.append(line.rstrip("\n"))

    if not kept_lines:
        raise ValueError(f"No ATOM records found after lightweight cleaning: {pdb_path}")

    text = "\n".join(kept_lines) + "\nEND\n"
    atomic_write_text(output_path, text)

    return str(output_path)


def clean_pdb_file(
    pdb_file: str | Path,
    outdir: str | Path,
    keep_water: bool = False,
    add_hydrogens: bool = True,
    ph: float = 7.0,
    use_pdbfixer: bool = True,
    fallback_to_lightweight: bool = True,
) -> str:
    """
    Clean a PDB file.

    Preferred path:
        PDBFixer/OpenMM cleanup.

    Fallback path:
        Lightweight ATOM-only cleanup.

    Parameters
    ----------
    pdb_file:
        Input PDB path.

    outdir:
        Directory where cleaned PDB should be written.

    keep_water:
        Whether PDBFixer should keep crystallographic water.

    add_hydrogens:
        Whether PDBFixer should add hydrogens.

    ph:
        pH used when adding hydrogens.

    use_pdbfixer:
        Whether to attempt PDBFixer cleanup.

    fallback_to_lightweight:
        Whether to fall back to ATOM-only cleanup if PDBFixer fails.

    Returns
    -------
    str
        Path to cleaned PDB file.
    """
    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()
    cleaned_dir = Path(os.path.expanduser(str(outdir))).resolve()
    cleaned_dir.mkdir(parents=True, exist_ok=True)

    pdb_id = pdb_path.stem
    cleaned_path = cleaned_dir / f"{pdb_id}_clean.pdb"

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file does not exist: {pdb_path}")

    if use_pdbfixer:
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile

            print(f"[clean] PDBFixer cleaning {pdb_path.name}")

            fixer = PDBFixer(filename=str(pdb_path))

            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(keepWater=keep_water)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            if add_hydrogens:
                fixer.addMissingHydrogens(ph)

            cleaned_path.parent.mkdir(parents=True, exist_ok=True)

            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=str(cleaned_path.parent),
                delete=False,
            ) as tmp:
                PDBFile.writeFile(fixer.topology, fixer.positions, tmp)
                tmp_path = Path(tmp.name)

            os.replace(tmp_path, cleaned_path)

            if cleaned_path.exists() and cleaned_path.stat().st_size > 0:
                return str(cleaned_path)

            raise RuntimeError(f"PDBFixer wrote empty cleaned file: {cleaned_path}")

        except Exception as exc:
            print(f"[WARN] PDBFixer failed for {pdb_path.name}: {exc}")

            if not fallback_to_lightweight:
                raise

            print(f"[clean] Falling back to lightweight cleaner for {pdb_path.name}")

    return lightweight_clean_pdb(
        pdb_file=pdb_path,
        output_file=cleaned_path,
        include_hetero=False,
    )


def group_residue_lines(
    pdb_file: str | Path,
    include_hetero: bool = False,
) -> List[Tuple[ResidueID, List[str]]]:
    """
    Group PDB ATOM lines by residue, preserving residue order.
    """
    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file does not exist: {pdb_path}")

    grouped: List[Tuple[ResidueID, List[str]]] = []
    current_residue: Optional[ResidueID] = None
    current_lines: List[str] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not is_pdb_atom_line(line, include_hetero=include_hetero):
                continue

            residue_id = parse_residue_id(line)

            if current_residue is None:
                current_residue = residue_id
                current_lines = [line]
                continue

            if residue_id == current_residue:
                current_lines.append(line)
            else:
                grouped.append((current_residue, current_lines))
                current_residue = residue_id
                current_lines = [line]

    if current_residue is not None and current_lines:
        grouped.append((current_residue, current_lines))

    return grouped


def split_residues_by_chain(
    residues: Sequence[Tuple[ResidueID, List[str]]],
) -> Dict[str, List[Tuple[ResidueID, List[str]]]]:
    """
    Split ordered residues into chain-specific lists.
    """
    chains: Dict[str, List[Tuple[ResidueID, List[str]]]] = {}

    for residue_id, lines in residues:
        chain = residue_id.chain_id.strip() or "NA"
        chains.setdefault(chain, []).append((residue_id, lines))

    return chains


def extract_kmer(
    pdb_file: str | Path,
    k: int = 3,
    include_hetero: bool = False,
    max_fragments: Optional[int] = None,
) -> Dict[str, List[str]]:
    """
    Extract residue k-mers from a cleaned PDB file.

    Parameters
    ----------
    pdb_file:
        Cleaned PDB path.

    k:
        Number of consecutive residues per fragment.

    include_hetero:
        Whether to allow HETATM records.

    max_fragments:
        Optional cap for debugging.

    Returns
    -------
    Dict[str, List[str]]
        Mapping from fragment key to PDB lines.

        Example key:
            5xus_A_ALA1_GLY3_000000
    """
    if k <= 0:
        raise ValueError(f"k must be positive, got: {k}")

    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()
    pdb_id = pdb_path.stem.replace("_clean", "")

    residues = group_residue_lines(pdb_path, include_hetero=include_hetero)

    if not residues:
        raise ValueError(f"No residues found in cleaned PDB: {pdb_path}")

    chains = split_residues_by_chain(residues)
    fragments: Dict[str, List[str]] = {}

    counter = 0

    for chain, chain_residues in chains.items():
        if len(chain_residues) < k:
            continue

        for start_idx in range(0, len(chain_residues) - k + 1):
            window = chain_residues[start_idx : start_idx + k]

            start_residue = window[0][0]
            end_residue = window[-1][0]

            fragment_key = (
                f"{pdb_id}_"
                f"{chain}_"
                f"{start_residue.safe_label()}_"
                f"{end_residue.safe_label()}_"
                f"{counter:06d}"
            )

            fragment_lines: List[str] = []

            for _, residue_lines in window:
                fragment_lines.extend(residue_lines)

            # A TER helps downstream tools understand the fragment boundary.
            fragment_lines.append("TER")
            fragment_lines.append("END")

            fragments[fragment_key] = fragment_lines

            counter += 1

            if max_fragments is not None and counter >= max_fragments:
                return fragments

    if not fragments:
        raise ValueError(
            f"No k-mers extracted from {pdb_path}. "
            f"Residues found: {len(residues)}, k={k}"
        )

    return fragments


def write_to_file(
    kmer_dict: Dict[str, List[str]],
    dirname: str | Path,
    out_subdir: str = "fix",
    overwrite: bool = True,
) -> List[str]:
    """
    Write k-mer fragments to PDB files.

    Parameters
    ----------
    kmer_dict:
        Mapping from fragment key to PDB lines.

    dirname:
        Parent output directory for one protein.

    out_subdir:
        Subdirectory for fragment PDBs. Default: fix.

    overwrite:
        Whether to overwrite existing files.

    Returns
    -------
    List[str]
        Paths to written fragment PDB files.
    """
    base_dir = Path(os.path.expanduser(str(dirname))).resolve()
    fix_dir = base_dir / out_subdir
    fix_dir.mkdir(parents=True, exist_ok=True)

    written_files: List[str] = []

    for key, lines in kmer_dict.items():
        output_file = fix_dir / f"{key}.pdb"

        if output_file.exists() and not overwrite:
            written_files.append(str(output_file))
            continue

        text = "\n".join(lines).rstrip() + "\n"
        atomic_write_text(output_file, text)

        if output_file.exists() and output_file.stat().st_size > 0:
            written_files.append(str(output_file))
        else:
            raise RuntimeError(f"Failed to write non-empty fragment file: {output_file}")

    return written_files


def preprocess_pdb(
    pdb_file: str | Path,
    outdir: str | Path,
    k: int = 3,
    remove_raw: bool = False,
    keep_water: bool = False,
    add_hydrogens: bool = True,
    ph: float = 7.0,
    use_pdbfixer: bool = True,
    max_fragments: Optional[int] = None,
    overwrite: bool = True,
) -> Dict[str, object]:
    """
    Full preprocessing pipeline for one PDB.

    Creates:

        outdir/<pdb_id>_<k>mers/
            cleaned/
            fix/
            pqr_files/
            atompot/

    Returns metadata useful for comboscript.py.
    """
    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()
    root_outdir = Path(os.path.expanduser(str(outdir))).resolve()

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file does not exist: {pdb_path}")

    pdb_id = pdb_path.stem
    protein_dir = root_outdir / f"{pdb_id}_{k}mers"
    cleaned_dir = protein_dir / "cleaned"
    fix_dir = protein_dir / "fix"
    pqr_dir = protein_dir / "pqr_files"
    atompot_dir = protein_dir / "atompot"

    for directory in [protein_dir, cleaned_dir, fix_dir, pqr_dir, atompot_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    cleaned_pdb = clean_pdb_file(
        pdb_file=pdb_path,
        outdir=cleaned_dir,
        keep_water=keep_water,
        add_hydrogens=add_hydrogens,
        ph=ph,
        use_pdbfixer=use_pdbfixer,
        fallback_to_lightweight=True,
    )

    fragments = extract_kmer(
        pdb_file=cleaned_pdb,
        k=k,
        include_hetero=False,
        max_fragments=max_fragments,
    )

    fragment_files = write_to_file(
        kmer_dict=fragments,
        dirname=protein_dir,
        out_subdir="fix",
        overwrite=overwrite,
    )

    if remove_raw:
        try:
            Path(cleaned_pdb).unlink()
            cleaned_pdb_removed = True
        except OSError:
            cleaned_pdb_removed = False
    else:
        cleaned_pdb_removed = False

    result: Dict[str, object] = {
        "pdb_id": pdb_id,
        "protein_dir": str(protein_dir),
        "cleaned_pdb": str(cleaned_pdb),
        "cleaned_pdb_removed": cleaned_pdb_removed,
        "fix_dir": str(fix_dir),
        "pqr_dir": str(pqr_dir),
        "atompot_dir": str(atompot_dir),
        "k": k,
        "num_fragments": len(fragment_files),
        "fragment_files": fragment_files,
    }

    return result


def safe_preprocess_pdb(
    pdb_file: str | Path,
    outdir: str | Path,
    log_dir: Optional[str | Path] = None,
    **kwargs,
) -> Dict[str, object]:
    """
    Wrapper that catches errors and optionally writes traceback logs.

    Useful for xargs/GNU parallel/SLURM-style runs where one bad PDB
    should not kill the entire batch.
    """
    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()

    try:
        return preprocess_pdb(
            pdb_file=pdb_path,
            outdir=outdir,
            **kwargs,
        )

    except Exception as exc:
        error_info = {
            "pdb_id": pdb_path.stem,
            "pdb_file": str(pdb_path),
            "success": False,
            "error": str(exc),
        }

        if log_dir is not None:
            log_path = Path(os.path.expanduser(str(log_dir))).resolve()
            log_path.mkdir(parents=True, exist_ok=True)

            traceback_file = log_path / f"{pdb_path.stem}.preproc.traceback.txt"
            atomic_write_text(traceback_file, traceback.format_exc())
            error_info["traceback_file"] = str(traceback_file)

        print(f"[ERROR] preprocessing failed for {pdb_path.name}: {exc}")

        return error_info


def copy_original_pdb(
    pdb_file: str | Path,
    destination_dir: str | Path,
) -> str:
    """
    Optional helper for preserving the original PDB inside the output folder.
    """
    pdb_path = Path(os.path.expanduser(str(pdb_file))).resolve()
    destination = Path(os.path.expanduser(str(destination_dir))).resolve()
    destination.mkdir(parents=True, exist_ok=True)

    output_file = destination / pdb_path.name
    shutil.copy2(pdb_path, output_file)

    return str(output_file)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clean a PDB file and extract residue k-mer PDB fragments."
    )

    parser.add_argument(
        "--pdb",
        required=True,
        help="Input PDB file.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Root output directory.",
    )

    parser.add_argument(
        "--k",
        type=int,
        default=3,
        help="Residue k-mer size. Default: 3.",
    )

    parser.add_argument(
        "--remove_raw",
        action="store_true",
        help="Remove cleaned full-length PDB after fragment extraction.",
    )

    parser.add_argument(
        "--no_pdbfixer",
        action="store_true",
        help="Skip PDBFixer and use lightweight ATOM-only cleaning.",
    )

    parser.add_argument(
        "--no_hydrogens",
        action="store_true",
        help="Do not add hydrogens during PDBFixer cleaning.",
    )

    parser.add_argument(
        "--keep_water",
        action="store_true",
        help="Keep water molecules during PDBFixer cleaning.",
    )

    parser.add_argument(
        "--ph",
        type=float,
        default=7.0,
        help="pH for hydrogen addition. Default: 7.0.",
    )

    parser.add_argument(
        "--max_fragments",
        type=int,
        default=None,
        help="Optional cap on number of k-mer fragments for debugging.",
    )

    parser.add_argument(
        "--log_dir",
        default=None,
        help="Optional directory for traceback logs.",
    )

    args = parser.parse_args()

    result = safe_preprocess_pdb(
        pdb_file=args.pdb,
        outdir=args.outdir,
        log_dir=args.log_dir,
        k=args.k,
        remove_raw=args.remove_raw,
        keep_water=args.keep_water,
        add_hydrogens=not args.no_hydrogens,
        ph=args.ph,
        use_pdbfixer=not args.no_pdbfixer,
        max_fragments=args.max_fragments,
    )

    if result.get("success") is False:
        raise SystemExit(1)

    print("[DONE] preprocessing complete")
    print(f"PDB ID: {result['pdb_id']}")
    print(f"Protein dir: {result['protein_dir']}")
    print(f"Fragments written: {result['num_fragments']}")


if __name__ == "__main__":
    main()
