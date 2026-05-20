#!/usr/bin/env python3
"""
identify_motifs.py

Identify and visualize common vs rare structural motifs (CATH architectures)
using outputs from frequency_analysis.py.

Inputs (one of):
  - --arch-file <path to architecture_distribution.tsv>
  - default: latest 'architecture_distribution_*.tsv' in ../analysis_outputs/

Expected columns in architecture_distribution:
  - architecture
  - count
  - fraction

Outputs (images only, no TSVs):
  - plots/common_architectures_<date>.png
  - plots/rare_architectures_<date>.png
  - plots/common_vs_rare_architectures_<date>.png
"""

import os
import glob
import argparse
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt


# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(SCRIPT_DIR, "..", "analysis_outputs")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")

os.makedirs(PLOT_DIR, exist_ok=True)


# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
def latest_architecture_tsv() -> str | None:
    """Return path to latest architecture_distribution_*.tsv, or None."""
    pattern = os.path.join(ANALYSIS_DIR, "architecture_distribution_*.tsv")
    paths = glob.glob(pattern)
    if not paths:
        return None
    return max(paths, key=os.path.getmtime)


def load_architecture_table(path: str) -> pd.DataFrame:
    """Load architecture_distribution TSV produced by frequency_analysis.py."""
    df = pd.read_csv(path, sep="\t")
    required = {"architecture", "count", "fraction"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")
    df["count"] = df["count"].astype(int)
    df["fraction"] = df["fraction"].astype(float)
    return df


def save_plot(fig, basename: str) -> str:
    """Save a matplotlib figure into plots/ with a date-stamped filename."""
    stamp = datetime.now().strftime("%Y-%m-%d")
    fname = f"{basename}_{stamp}.png"
    path = os.path.join(PLOT_DIR, fname)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved plot] {path}")
    return path


# -------------------------------------------------------------------
# Visualization functions
# -------------------------------------------------------------------
def plot_common(arch_df: pd.DataFrame, top_k: int) -> None:
    """Bar plot of most common architectures."""
    data = (
        arch_df.sort_values("count", ascending=False)
        .head(top_k)
        .copy()
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(10, 0.4 * len(data) + 2))
    ax.barh(data["architecture"].astype(str), data["count"])
    ax.invert_yaxis()  # largest at top
    ax.set_xlabel("Domain count")
    ax.set_ylabel("Architecture")
    ax.set_title(f"Most common CATH architectures (top {top_k})")
    fig.tight_layout()

    save_plot(fig, f"common_architectures_top{top_k}")
    plt.close(fig)


def plot_rare(arch_df: pd.DataFrame, bottom_k: int) -> None:
    """Bar plot of rarest architectures."""
    data = (
        arch_df.sort_values("count", ascending=True)
        .head(bottom_k)
        .copy()
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(10, 0.4 * len(data) + 2))
    ax.barh(data["architecture"].astype(str), data["count"])
    ax.set_xlabel("Domain count")
    ax.set_ylabel("Architecture")
    ax.set_title(f"Rarest CATH architectures (bottom {bottom_k})")
    fig.tight_layout()

    save_plot(fig, f"rare_architectures_bottom{bottom_k}")
    plt.close(fig)


def plot_common_vs_rare(arch_df: pd.DataFrame, top_k: int, bottom_k: int) -> None:
    """Combined bar plot for common vs rare architectures."""
    common = (
        arch_df.sort_values("count", ascending=False)
        .head(top_k)
        .copy()
    )
    rare = (
        arch_df.sort_values("count", ascending=True)
        .head(bottom_k)
        .copy()
    )

    common["group"] = "common"
    rare["group"] = "rare"

    combined = pd.concat([common, rare], ignore_index=True)
    combined["label"] = combined["group"] + ":" + combined["architecture"].astype(str)

    # Sort so rare at bottom, common at top (within group)
    combined = combined.sort_values(["group", "count"], ascending=[True, True])

    fig, ax = plt.subplots(figsize=(10, 0.4 * len(combined) + 2))
    ax.barh(combined["label"], combined["count"])
    ax.set_xlabel("Domain count")
    ax.set_ylabel("Group:Architecture")
    ax.set_title(f"Common (top {top_k}) vs rare (bottom {bottom_k}) CATH architectures")
    fig.tight_layout()

    save_plot(fig, f"common_vs_rare_architectures_top{top_k}_bottom{bottom_k}")
    plt.close(fig)


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Visualize common vs rare structural motifs using architecture_distribution from frequency_analysis.py"
    )
    parser.add_argument(
        "--arch-file",
        type=str,
        default=None,
        help="Path to architecture_distribution TSV (default: latest architecture_distribution_*.tsv in ../analysis_outputs/)",
    )
    parser.add_argument("--top-k", type=int, default=20, help="Number of most common motifs to visualize.")
    parser.add_argument("--bottom-k", type=int, default=20, help="Number of rare motifs to visualize.")

    args = parser.parse_args()

    arch_path = args.arch_file or latest_architecture_tsv()
    if not arch_path:
        raise SystemExit(
            "No architecture_distribution_*.tsv found in ../analysis_outputs/ "
            "and no --arch-file provided. Run frequency_analysis.py first."
        )

    print(f"[using architecture file] {arch_path}")
    arch_df = load_architecture_table(arch_path)

    plot_common(arch_df, top_k=args.top_k)
    plot_rare(arch_df, bottom_k=args.bottom_k)
    plot_common_vs_rare(arch_df, top_k=args.top_k, bottom_k=args.bottom_k)


if __name__ == "__main__":
    main()
