#!/usr/bin/env python3
"""
Frequency analysis of protein architectural classes (CATH "a" level).

- Reads a TSV of CATH domains (one row per domain).
- Computes:
    * Global frequency of each architecture (a).
    * Frequency of (class, architecture) pairs: (c, a).
- Outputs:
    * architecture_distribution.tsv
    * class_architecture_distribution.tsv
- Optionally: bar plot of architecture frequencies.

"""

import os
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------- CONFIG ----------------------------------------

FILE_PATH = "cath_cache/outputs/master_domains.tsv"
OUTPUT_DIR = "analysis_outputs"
PLOT_DIR = os.path.join(OUTPUT_DIR, "plots")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# ------------------------------ HELPERS --------------------------------------


def load_cath_tsv(path: str) -> pd.DataFrame:
    """
    Load the CATH domain table from a TSV file.

    We explicitly set sep="\\t" because many CATH exports are tab-separated.
    """
    df = pd.read_csv(path, sep="\t")

    print("Loaded file:", path)
    print("Columns:", df.columns.tolist())
    print("First few rows:")
    print(df.head(), "\n")

    return df


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure the expected CATH hierarchy columns exist and are of type string.

    Expected minimal columns:
        - 'c' : class
        - 'a' : architecture
        - 't' : topology
        - 'h' : homologous superfamily
    """
    expected = {"c", "a"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns for architecture analysis: {missing}\n"
            f"Available columns: {df.columns.tolist()}"
        )

    # Cast to string to avoid numeric issues
    for col in ["c", "a", "t", "h"]:
        if col in df.columns:
            df[col] = df[col].astype(str)

    return df


def architecture_frequency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute frequency of each architecture ('a') across all domains.

    Returns a DataFrame with:
        architecture (a), count, fraction
    """
    total = len(df)
    counts = (
        df.groupby("a")
          .size()
          .reset_index(name="count")
          .sort_values("count", ascending=False)
    )
    counts["fraction"] = counts["count"] / total
    counts = counts.rename(columns={"a": "architecture"})

    return counts


def class_architecture_frequency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute frequency of each (class, architecture) pair.

    Returns a DataFrame with:
        class (c), architecture (a), count, fraction_within_class
    """
    # Count per (c, a)
    grouped = (
        df.groupby(["c", "a"])
          .size()
          .reset_index(name="count")
    )

    # Total per class
    class_totals = (
        df.groupby("c")
          .size()
          .reset_index(name="class_total")
    )

    merged = grouped.merge(class_totals, on="c", how="left")
    merged["fraction_within_class"] = merged["count"] / merged["class_total"]

    merged = merged.sort_values(["c", "count"], ascending=[True, False])
    merged = merged.rename(columns={"c": "class", "a": "architecture"})

    return merged


def save_table(df: pd.DataFrame, name: str):
    """
    Save a DataFrame as TSV in OUTPUT_DIR with a date-stamped filename.
    """
    stamp = datetime.now().strftime("%Y-%m-%d")
    fname = f"{name}_{stamp}.tsv"
    path = os.path.join(OUTPUT_DIR, fname)
    df.to_csv(path, sep="\t", index=False)
    print(f"[saved table] {path}")


def plot_architecture_bar(counts: pd.DataFrame, top_n: int = 20):
    """
    Simple bar plot of architecture frequencies.

    counts: DataFrame from architecture_frequency()
    """
    data = counts.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.bar(data["architecture"].astype(str), data["count"])
    ax.set_xlabel("Architecture (a)")
    ax.set_ylabel("Number of domains")
    ax.set_title(f"Top {top_n} architectures by domain count")
    ax.set_xticklabels(data["architecture"].astype(str), rotation=90)

    fig.tight_layout()

    # Save plot
    stamp = datetime.now().strftime("%Y-%m-%d")
    fname = f"architecture_distribution_top{top_n}_{stamp}.png"
    path = os.path.join(PLOT_DIR, fname)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved plot] {path}")

    plt.show()


# ----------------------------- MAIN LOGIC ------------------------------------


def main():
    # 1. Load TSV
    df = load_cath_tsv(FILE_PATH)

    # 2. Ensure required columns are present and normalized
    df = normalize_columns(df)

    # 3. Global architecture frequencies
    arch_freq = architecture_frequency(df)
    print("=== Architecture frequency (global) ===")
    print(arch_freq.head(20), "\n")

    save_table(arch_freq, "architecture_distribution")

    # 4. Class × architecture frequencies
    class_arch_freq = class_architecture_frequency(df)
    print("=== Class × architecture frequency ===")
    print(class_arch_freq.head(20), "\n")

    save_table(class_arch_freq, "class_architecture_distribution")

    # 5. Optional visualization
    plot_architecture_bar(arch_freq, top_n=20)


if __name__ == "__main__":
    main()
