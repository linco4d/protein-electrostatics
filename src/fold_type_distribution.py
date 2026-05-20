import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import os
from datetime import datetime
# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
INPUT_CSV = "cath_cache/outputs/master_domains.tsv"   # <--- change to your path

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
df = pd.read_csv("cath_cache/outputs/master_domains.tsv", sep="\t")


# Sanity check
expected_cols = {"domain_id", "pdb_id", "chain_id", "c", "a", "t", "h", "cat", "cath"}
missing = expected_cols - set(df.columns)
if missing:
    print("Warning: missing expected columns:", missing)

# Make sure c, a, t, h are strings to avoid numeric quirks
for col in ["c", "a", "t", "h"]:
    df[col] = df[col].astype(str)

# Optional convenience columns
df["class"] = df["c"]
df["architecture"] = df["a"]
df["topology"] = df["t"]          # "fold level"
df["homol_superfamily"] = df["h"]

#SAVE IMAGES

# Folder where all plots will be saved
PLOT_DIR = "fold_type_plots"
os.makedirs(PLOT_DIR, exist_ok=True)

def save_fig_matplotlib(fig, name: str):
    """
    Save a Matplotlib figure using a standardized filename.
    Example: 'plots/topology_distribution_2025-12-03.png'
    """
    timestamp = datetime.now().strftime("%Y-%m-%d")
    filename = f"{name}_{timestamp}.png"
    path = os.path.join(PLOT_DIR, filename)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved matplotlib figure] {path}")

def save_fig_plotly(fig, name: str):
    """
    Save a Plotly figure both as:
      - PNG  (static)
      - HTML (interactive)
    """
    timestamp = datetime.now().strftime("%Y-%m-%d")
    base = f"{name}_{timestamp}"

    png_path  = os.path.join(PLOT_DIR, base + ".png")
    html_path = os.path.join(PLOT_DIR, base + ".html")

    # Save PNG (requires kaleido: pip install kaleido)
    try:
        fig.write_image(png_path, scale=2)
        print(f"[saved plotly PNG] {png_path}")
    except Exception as e:
        print("[warning] Could not save PNG (install kaleido). Saving HTML only.")

    # Always save interactive HTML
    fig.write_html(html_path)
    print(f"[saved plotly HTML] {html_path}")

# ------------------------------------------------------------------
# 2. Helper: count distributions
# ------------------------------------------------------------------
def level_counts(frame: pd.DataFrame, level: str) -> pd.DataFrame:
    """
    Count number of domains per given level (e.g., 'class', 'architecture', etc.).
    Returns sorted dataframe with columns [level, count].
    """
    out = (
        frame
        .groupby(level)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )
    return out

counts_class  = level_counts(df, "class")
counts_arch   = level_counts(df, "architecture")
counts_topo   = level_counts(df, "topology")         # fold types
counts_homol  = level_counts(df, "homol_superfamily")

print("Class distribution:\n", counts_class.head(), "\n")
print("Architecture distribution:\n", counts_arch.head(), "\n")
print("Topology (fold) distribution:\n", counts_topo.head(), "\n")
print("Homologous superfamily distribution:\n", counts_homol.head(), "\n")

# ------------------------------------------------------------------
# 3. Simple bar plots for each level
# ------------------------------------------------------------------
def plot_bar(counts_df: pd.DataFrame, level: str, top_n: int = 30):
    """
    Bar plot for distribution of a CATH level.
    Saves the figure automatically.
    """
    data = counts_df.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.bar(data[level].astype(str), data["count"])
    ax.set_xticklabels(data[level].astype(str), rotation=90)
    ax.set_xlabel(level)
    ax.set_ylabel("Number of domains")
    ax.set_title(f"Distribution of domains across CATH {level} (top {top_n})")
    fig.tight_layout()

    # Save it
    save_fig_matplotlib(fig, f"{level}_distribution")

    plt.show()


# Few classes => show all
plot_bar(counts_class, "class", top_n=len(counts_class))

# Architectures, topologies, etc.
plot_bar(counts_arch, "architecture", top_n=20)
plot_bar(counts_topo, "topology", top_n=30)

