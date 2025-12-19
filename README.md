# Different_ways_to_measure_cell_specific_expression
This page explains 4 different ways I tried to build a tool to measure cell specific or enriched gene expression with bash/python and HTML. Then I have the tools I built to compare their effectiveness. 

with the tests, celltype_enrichment_v1_4 (below) stands out as the best version so far. I have other methods I used explained after this section

If you want more info on some explanations on celltype_enrichment_v1_4 that are not in it's description, having a look at the section 'Simple Enrichment Scores' might help

#*** FOR ANY OF THESE TO WORK PROPERLY, YOU HAVE TO INSTALL REQUIRED PACKAGES OR LOAD A PYTHON ENVIRONMENT WITH REQUIRED PACKAGES FIRST **

It looks like this on computecanada/alliancecan
```
module load gcc arrow
module load python
python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate
```
use something like following to install packages if they haven't been installed already
```
pip install numpy
```
#
#*****************************************************************************************************************************************

#----------------------------------__________celltype_enrichment_v1_4________-------------------------------------------------------------

#Link to celltype_enrichment_v1_4

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#best-method---celltype_enrichment_v_14

#----------------------------------__________celltype_enrichment_v1_4________-------------------------------------------------------------

#*****************************************************************************************************************************************

#
#
#
Following are the other methods I tried and built this on

#Method 01 - Simple Enrichment Scores

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#method-01---simple-enrichment-scores


#Method 02 - AUROC SCORES

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#method-02---auroc-scores

#Method 03 - Enrichment score with more robust statistical backup systems to filter useful data

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#method-03---enrichment-score-with-more-robust-statistical-backup-systems-to-filter-useful-data

#I wrote universal_plot_maker.py to make interactive plots with any tsv file changing needed fields 

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#universal-interactive-plot-maker

#And finally I have a script that helps to compare the effectiveness between different methods

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#comparing-method-effectiveness

# METHOD 01 - Simple Enrichment Scores 
Simplest way

This method gives you a python script that allows you to to calculate enrichment scores from a tsv file with a list of genes/ gene names and expression values like nCPM. I found this file from Human Protein Atlas

Enrichment score is a simple the ratio between -  the average nCPM in a particular cell type for a particular gene/  average nCPM in all the other cell types for that particular gene

following is the head of the file I obtained from Human Protein Atlas
https://www.proteinatlas.org/humanproteome/single+cell/single+cell+type/data

rna_single_cell_clusters.tsv.zip

rna_single_cell_cluster.tsv
```txt
Gene    Gene name       Tissue  Cluster Cell type       Read count      nCPM
ENSG00000000003 TSPAN6  ovary   c-0     ovarian stromal cells   493     92.5
ENSG00000000003 TSPAN6  ovary   c-1     ovarian stromal cells   529     80.5
ENSG00000000003 TSPAN6  ovary   c-2     ovarian stromal cells   143     52.3
ENSG00000000003 TSPAN6  ovary   c-3     ovarian stromal cells   456     91.4
```
Then I calculated Enrichment Scores for different gene-cell type combinations with the python script below.

You can run it like following
```
python cal_enrich_with_custom_clustor_no.py   --input-file rna_single_cell_cluster.tsv   --output-file my_custom_clustor_enrichment.tsv   --min-clusters 3
```
--input-file gives you the option to select input tsv file

--output-file gives you the option to select output tsv file
In addition to that, this also produces a file with top 100 enrichment values

--min-clusters THIS IS IMPORTANT. This input tsv can have cell clusters mistakenly assigned to wrong cell types.
Therefore I use a minimum number of cell clusters a particular gene-cell type combination should have to be included
in this calculation

Following is the python script

cal_enrich_with_custom_clustor_no.py
```
import pandas as pd
import numpy as np
import argparse

# --- CLI: optional minimum clusters filter + input/output filenames ---
parser = argparse.ArgumentParser(
    description="Compute enrichment with optional filtering and configurable I/O."
)
parser.add_argument(
    "--min-clusters",
    type=int,
    default=None,
    help="Keep only rows where clusters_used >= this integer. If omitted, no filtering is applied."
)
parser.add_argument(
    "--input-file",
    type=str,
    default="rna_single_cell_cluster.tsv",
    help="Path to the input TSV containing single-cell clusters (default: rna_single_cell_cluster.tsv)."
)
parser.add_argument(
    "--output-file",
    type=str,
    default="all_gene_cell_enrichment_data.tsv",
    help="Path for the final enrichment TSV output (default: all_gene_cell_enrichment_data.tsv)."
)
args = parser.parse_args()

# Load the TSV file into cell_cluster_data
file_path = args.input_file
cell_cluster_data = pd.read_csv(file_path, sep="\t")

# --- Aggregation to one row per (Gene × Cell type) with mean nCPM ---
agg_df = (
    cell_cluster_data
    .groupby(['Gene', 'Gene name', 'Cell type'], as_index=False)
    .agg(nCPM=('nCPM', 'mean'))
)

# --- Count how many clusters contributed to each (Gene × Cell type) ---
# Prefer counting unique cluster IDs if a 'Cluster' column exists; otherwise use row count.
if 'Cluster' in cell_cluster_data.columns:
    cluster_counts = (
        cell_cluster_data
        .groupby(['Gene', 'Gene name', 'Cell type'])['Cluster']
        .nunique()
        .reset_index(name='clusters_used')
    )
else:
    cluster_counts = (
        cell_cluster_data
        .groupby(['Gene', 'Gene name', 'Cell type'])
        .size()
        .reset_index(name='clusters_used')
    )

# --- Merge counts into agg_df ---
agg_df = agg_df.merge(cluster_counts, on=['Gene', 'Gene name', 'Cell type'], how='left')

# --- Filter by clusters_used using CLI input (happens immediately after merge) ---
if args.min_clusters is not None:
    before_n = len(agg_df)
    agg_df = agg_df[agg_df['clusters_used'] >= args.min_clusters].copy()
    after_n = len(agg_df)
    dropped_n = before_n - after_n
    print(f"\033[33mFiltered by clusters_used >= {args.min_clusters}. "
          f"Dropped {dropped_n} row(s); {after_n} row(s) remain.\033[0m")

# Rename column nCPM to avg_nCPM
agg_df = agg_df.rename(columns={'nCPM': 'avg_nCPM'})

#checking for non numeric values in nCPM to avoid errors
# Convert to numeric, invalid entries become NaN
numeric_series = pd.to_numeric(agg_df['avg_nCPM'], errors='coerce')
# Find unique non-numeric values
non_numeric_values = agg_df.loc[numeric_series.isna(), 'avg_nCPM'].unique()
# Check and print message
if len(non_numeric_values) > 0:
    print(f"\033[31mNon numeric values found in avg_nCPM: {non_numeric_values}\033[0m")
else:
    print("\033[32mGOOD DATA nCPM. NO Non numeric\033[0m")

#****This is a function to drop rows with genes that do not express in any cell types ****
def drop_genes_with_no_expression(agg_df, expr_col=None, treat_nan_as_zero=False):
    """
    Remove rows for genes whose expression is zero across all cell types.
    Parameters
    ----------
    agg_df : pd.DataFrame
      Aggregated DataFrame with one row per (Gene × Cell type).
    expr_col : str or None
      Column with aggregated expression (defaults to 'avg_nCPM' if present, else 'nCPM').
    treat_nan_as_zero : bool
      If True, NaN values are treated as zero when determining 'no expression'.
      If False, genes with NaN-only expressions will NOT be classified as zero; they are kept.
    Returns
    -------
    filtered_df : pd.DataFrame
      agg_df with rows for zero-expression genes dropped.
    dropped_genes : list
      List of gene identifiers that were removed.
    """
    # Pick expression column
    if expr_col is None:
        expr_col = 'avg_nCPM' if 'avg_nCPM' in agg_df.columns else 'nCPM'
    # Work on a copy
    df = agg_df.copy()
    # Ensure numeric and clean infs
    df[expr_col] = pd.to_numeric(df[expr_col], errors='coerce')
    df[expr_col] = df[expr_col].replace([np.inf, -np.inf], np.nan)
    # Optionally treat NaN as zero for the test
    expr_for_test = df[expr_col].fillna(0) if treat_nan_as_zero else df[expr_col]
    # Compute per-gene max expression across all cell types
    gene_max = expr_for_test.groupby(df['Gene']).max()
    # Genes where max == 0 → zero across all rows
    genes_all_zero = gene_max[gene_max == 0].index.tolist()
    # Filter out those genes
    filtered_df = df[~df['Gene'].isin(genes_all_zero)].copy()
    return filtered_df, genes_all_zero

# ---- Example usage ----
# Run this BEFORE enrichment calculation (recommended), or after if needed.
agg_df_filtered, dropped = drop_genes_with_no_expression(agg_df, expr_col=None, treat_nan_as_zero=False)
print(f"\033[31mDropped {len(dropped)} row(s) with zero expression across all cell types.\033[0m")

# Count unique (non-null) values in 'Cell type'
n_unique = agg_df['Cell type'].dropna().nunique()
# Create a sorted list of unique (non-null) cell types
unique_cell_types = sorted(agg_df['Cell type'].dropna().unique().tolist())
# Save to TSV (one column named 'Cell type')
pd.Series(unique_cell_types, name='Cell type').to_csv('unique_cell_types.tsv', sep='\t', index=False)
print(f"\033[32mNumber of unique cell types: {n_unique}\033[0m")

# overwrite agg_df:
agg_df = agg_df_filtered

# Count unique gene names after dropping genes with no expression
n_unique_genes = agg_df['Gene'].dropna().nunique()
# Print in green
print(f"\033[32mNumber of unique genes remaining: {n_unique_genes}\033[0m")

#****
print("\033[33mCalculating Enrichment Scores....\033[0m")

# --- enrichment score calculation ---
def add_enrichment(
    agg_df: pd.DataFrame,
    gene_col: str = "Gene",
    value_col: str = "nCPM",
    out_col: str = "Enrichment score",
    min_background: float = 1e-3, # minimum allowed background (denominator)
    min_expression: float = 0.0,  # minimum required numerator (nCPM) to compute enrichment
    min_count: int = 2,           # require at least this many cell-type entries per gene
    pseudocount: float | None = None, # optional stabilizer added to background (and/or numerator)
    clip_max: float | None = None     # optional: cap extreme enrichment values
):
    """
    Compute enrichment per row as nCPM / average(nCPM of other cell types for the same gene),
    with safeguards to avoid infinities and noise from tiny denominators.
    """
    df = agg_df.copy()
    # Ensure numeric; coerce invalids to NaN
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    # Group-level sums and counts
    gene_sums = df.groupby(gene_col)[value_col].transform("sum")
    gene_counts = df.groupby(gene_col)[value_col].transform("count")
    # Average of "other" cell types = (sum - current) / (count - 1)
    denom_counts = gene_counts - 1
    avg_other = (gene_sums - df[value_col]) / denom_counts
    # Invalid if there is no "other" (i.e., count <= 1)
    avg_other = avg_other.mask(denom_counts <= 0, np.nan)
    # Optional stabilization via pseudocount on the denominator (and numerator if desired)
    if pseudocount is not None:
        # df[value_col] = df[value_col] + pseudocount  # uncomment to add to numerator
        avg_other = avg_other + pseudocount
    # Enforce minimum background: raise small denominators to min_background
    denom = np.maximum(avg_other, min_background)
    # Enforce minimum expression: rows below threshold are set to NaN
    numer = df[value_col].where(df[value_col] >= min_expression, np.nan)
    # Safe divide (NaN when denom <= 0 or numer is NaN)
    df[out_col] = np.divide(
        numer, denom,
        out=np.full(df.shape[0], np.nan, dtype=float),
        where=(denom > 0)
    )
    # If avg_other itself is NaN (e.g., insufficient counts), keep NaN
    df.loc[avg_other.isna(), out_col] = np.nan
    # Optional: cap extreme ratios
    if clip_max is not None:
        df[out_col] = df[out_col].clip(upper=clip_max)
    return df

# --- Compute enrichment with sensible safeguards ---
agg_df = add_enrichment(
    agg_df,
    gene_col="Gene",
    value_col="avg_nCPM",
    out_col="Enrichment score",
    min_background=1e-3, # lift very small denominators
    min_expression=0.0,  # require >= this to compute enrichment
    min_count=2,         # at least 2 entries per gene
    pseudocount=None,    # try setting to e.g. 0.01 for stabilization
    clip_max=None        # e.g., set to 100 to cap extreme ratios
)

#****
print("\033[33mDone calculating\033[0m")

#******
# --- Choose the expression column (use avg_nCPM if you renamed it) ---
expr_col = 'avg_nCPM' if 'avg_nCPM' in agg_df.columns else 'nCPM'

# --- Build/confirm the single-cell-type flag ---
gene_celltype_counts = agg_df.groupby('Gene')['Cell type'].transform('nunique')
agg_df['single_cell_type_gene'] = (gene_celltype_counts == 1)
min_count = min(gene_celltype_counts)
print(f"\033[33mNumber of cell types for the gene expressed in the least amount of genes: {min_count}\033[0m")

# Show rows for single-cell-type genes ---
single_cell_rows = agg_df[agg_df['single_cell_type_gene']].copy()
if not single_cell_rows.empty:
    n_genes = single_cell_rows['Gene'].nunique()
    print(f"\033[32mFound {n_genes} gene(s) expressed in exactly one cell type.\033[0m")
    print(single_cell_rows.to_string(index=False))
    single_cell_rows.to_csv('single_cell_type_gene_rows.tsv', sep='\t', index=False)
else:
    print("\033[33mNo genes were found that are only expressed in one cell type\033[0m")

#******
# ssort agg_df by enrichment
agg_df = agg_df.sort_values(by="Enrichment score", ascending=False)

# Save agg_df to a TSV file
agg_df.to_csv(args.output_file, sep='\t', index=False)
print(f"\033[33mFiles saved: {args.output_file}\033[0m")

# select top x rows
# Sort by Enrichment score in descending order and take top 100
top100 = agg_df.head(100)
# Save to TSV file
top100.to_csv("top100_enrichment.tsv", sep="\t", index=False)
print("Saved top 100 rows to top100_enrichment.tsv")
```
# Universal interactive plot maker
I made a plot maker that can be used unversally, giving it commands to analyze custom values

I wanted to plot this in an interactive way. This allows user to decide how many values to include in the plot (top x amount of data points from the previous output). This resulting tool lets the user select different cell types and genes to see the selected data points. In addition, you can export a tsv file with the selected values.

This is the script you need to run universal_plot_maker.py.

run_universal_plot_maker_with_options.sh
```
#!/usr/bin/env bash
# Usage:
#   ./run_universal_plot_maker_with_options.sh [overrides...]
# Example:
#   ./run_universal_plot_maker_with_options.sh --file simple_enrich_1_clustor.tsv --out simple_enrich_1_clustor_plot.html

set -euo pipefail

# Resolve the directory of this script so we can find the Python file reliably
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Default arguments (can be overridden by CLI options appended below)
args=(
  --file my_custom_clustor_enrichment.tsv          # REQUIRED: Input data file (TSV/CSV/etc.)
  --out simple_enrich_c3.html                      # Output HTML file name
  --top 35000                                       # Top N rows to plot (default: 100)
  --dedupe-policy mean                             # [max|mean|median|first] aggregation
  --log                                            # Use log scale on numeric axis
  # --linear                                       # Use linear scale instead (mutually exclusive with --log)
  # --horizontal                                   # Horizontal bars (better for long labels)
  --self-contained                                 # Embed Plotly.js for offline HTML
  # --log-digits D2                                # Log-axis minor ticks: D1 (all) or D2 (2 & 5)
  # --lang en-CA                                   # HTML lang attribute (default: en)
  --initial-zoom 100                               # Initial number of bars visible on load
  # --sep $'\t'                                    # Field separator (auto-detected if omitted)
  --x-col "Gene name"                              # Column for X axis (numeric if horizontal)
  --y-col "Enrichment score"                       # Column for Y axis (categorical if horizontal)
  --label-col "Gene"                               # Explicit label column (optional)
  # --value-col Score                              # Explicit numeric value column (optional)
  --group-col "Cell type"                          # Column for color grouping (legend)
  --search-col "Gene name"                         # Column used for search box
  --details "Gene" "Gene name" "Cell type" "clusters_used" "Enrichment score" "single_cell_type_gene"  # Extra columns
)

# Append any CLI overrides so the LAST occurrence of options wins
args+=( "$@" )

# Invoke the Python script with the collected arguments
python3 "${script_dir}/universal_plot_maker.py" "${args[@]}"
```
Then make it an executeble
```
chmod +x run_universal_plot_maker_with_options.sh
```
and then run it changing needed flags. The rest of the flags will be used from run_universal_plot_maker_with_options.sh
```
./run_universal_plot_maker_with_options.sh --file simple_enrich_1_clustor.tsv --out simple_enrich_1_clustor_plot.html
```
This changes input file and output file keeping the rest from run_universal_plot_maker_with_options.sh


following is the python script

universal_plot_maker.py
```py

#!/usr/bin/env python3
"""
Interactive bar chart (table → HTML) with filter/search, details-on-click,
safe duplicate handling, per-group colors, and TSV export.

Highlights:
- ONE bar per (Gene, Cell type) using a combined category internally (no stacking).
- X-axis tick labels show ONLY the gene name, while categories remain unique (Gene | Cell type).
- Group dropdown (filter) and per-bar colors by group.
- Legend hidden by default (use --show-legend to enable).
- Canonical payload keys: __VAL__ (numeric value), __CAT__ (combined category), __TICK__ (gene-only tick label).
- Full-row hover: ALL columns from --details shown on hover (server + client).
- --dedupe-policy pre-collapses duplicates safely and keeps the group column.

Example:
  python universal_plot_maker.py \
    --file all_gene_cell_auroc.tsv \
    --out auroc_plot.html \
    --x-col "Gene name" \
    --y-col "AUROC" \
    --label-col "Gene" \
    --group-col "Cell type" \
    --search-col "Gene name" \
    --top 1000 \
    --initial-zoom 100 \
    --dedupe-policy mean \
    --self-contained
"""

import argparse
import json
import re
import sys
from typing import Optional, List, Tuple, Dict

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

# ------------------------------ #
# Data Loading & Column Resolution
# ------------------------------ #
def _infer_sep(path: str, sep_arg: Optional[str]) -> Optional[str]:
    if sep_arg is not None:
        return sep_arg
    lower = path.lower()
    if lower.endswith(".tsv") or lower.endswith(".tab"):
        return "\t"
    if lower.endswith(".csv"):
        return ","
    return None

def load_table(path: str, sep: Optional[str]) -> pd.DataFrame:
    if sep is None:
        df = pd.read_csv(path, sep=None, engine="python")
    else:
        df = pd.read_csv(path, sep=sep)
    df.columns = [str(c).strip() for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return df

def resolve_columns(
    df: pd.DataFrame,
    x_col: Optional[str],
    y_col: Optional[str],
    label_col: Optional[str],
    value_col: Optional[str],
    group_col: Optional[str],
    search_col: Optional[str],
    orientation: str,
) -> Tuple[str, str, Optional[str], Optional[str]]:
    """Resolve label/value/group/search columns; coerce numeric values; drop NA."""
    def is_num(series) -> bool:
        return pd.api.types.is_numeric_dtype(series)

    numeric_cols = [c for c in df.columns if is_num(df[c])]
    text_cols = [c for c in df.columns if not is_num(df[c])]

    # Explicit x/y given
    if x_col and y_col:
        x_is_num = is_num(df[x_col])
        y_is_num = is_num(df[y_col])
        if orientation == "v":
            label_col_final = x_col if not x_is_num else (label_col or (text_cols[0] if text_cols else x_col))
            value_col_final = y_col if y_is_num else (value_col or (numeric_cols[0] if numeric_cols else y_col))
        else:
            label_col_final = y_col if not y_is_num else (label_col or (text_cols[0] if text_cols else y_col))
            value_col_final = x_col if x_is_num else (value_col or (numeric_cols[0] if numeric_cols else x_col))
    else:
        # Only y given → treat as numeric value
        if y_col and (y_col in df.columns):
            value_col_final = y_col
        else:
            value_col_final = (
                value_col
                or ("AUROC" if "AUROC" in df.columns else None)
                or (numeric_cols[0] if numeric_cols else None)
            )
        label_col_final = (
            label_col
            or ("Gene name" if "Gene name" in df.columns else None)
            or ("Gene" if "Gene" in df.columns else None)
            or (text_cols[0] if text_cols else df.columns[0])
        )

    if value_col_final is None:
        raise ValueError("Could not infer a numeric value column. Provide --value-col or both --x-col/--y-col.")

    # Prefer common names if present
    group_col_final = (
        group_col if (group_col and group_col in df.columns)
        else ("Cell type" if "Cell type" in df.columns else None)
        or ("CellType" if "CellType" in df.columns else None)
    )
    search_col_final = (
        search_col if (search_col and search_col in df.columns)
        else (label_col_final if label_col_final in df.columns else None)
    )

    # Sanitize numeric values for value column
    if df[value_col_final].dtype == object:
        df[value_col_final] = (
            df[value_col_final]
            .astype(str).str.strip()
            .str.replace(",", ".", regex=False)
            .replace({"NA": None, "N/A": None, "null": None, "None": None, "": None})
        )
    df[value_col_final] = pd.to_numeric(df[value_col_final], errors="coerce")
    df.dropna(subset=[value_col_final], inplace=True)
    if df.empty or df[value_col_final].notna().sum() == 0:
        raise ValueError(f"No numeric values found in column '{value_col_final}'.")

    # Ensure labels are usable
    df[label_col_final] = df[label_col_final].astype(str).fillna("").str.strip()
    df = df[df[label_col_final] != ""]

    return label_col_final, value_col_final, group_col_final, search_col_final

# ------------------------------ #
# Duplicate handling (pre-plot)
# ------------------------------ #
def dedupe_for_plot(
    df: pd.DataFrame,
    label_col: str,
    value_col: str,
    group_col: Optional[str],
    policy: str = "error",
    keep_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Collapse duplicates BEFORE plotting to avoid implicit stacking.
    Keys = (label_col, group_col) if group_col present, else (label_col).
    policy: "error" | "max" | "mean" | "median" | "first"
    keep_cols: columns to retain (e.g., --details)
    """
    keys = [label_col] + ([group_col] if group_col else [])
    cols_safe = keys + [value_col]
    if keep_cols:
        cols_safe += [c for c in keep_cols if c in df.columns]

    if not keys:
        return df.loc[:, list(dict.fromkeys(cols_safe))]

    dup_mask = df.duplicated(subset=keys, keep=False)
    if not dup_mask.any():
        return df.loc[:, list(dict.fromkeys(cols_safe))]

    if policy == "error":
        dupe_rows = df.loc[dup_mask, keys + [value_col]].head(10)
        raise ValueError(
            f"Duplicate rows for keys {keys} detected. Example duplicates:\n{dupe_rows}\n"
            f"Use --dedupe-policy [max|mean|median|first] to aggregate before plotting."
        )

    agg_map = {"max": "max", "mean": "mean", "median": "median", "first": "first"}
    reducer = agg_map.get(policy, "max")

    # Aggregate values per key; then reattach desired extras from the first occurrence
    agg_df = (df.groupby(keys, as_index=False).agg(**{value_col: (value_col, reducer)}))
    if keep_cols:
        extras = [c for c in keep_cols if c in df.columns and c not in agg_df.columns]
        if extras:
            rep = df.drop_duplicates(subset=keys, keep="first")[keys + extras]
            agg_df = agg_df.merge(rep, on=keys, how="left")

    # Always only safe columns
    agg_df = agg_df.loc[:, list(dict.fromkeys(cols_safe))]

    # Ensure group_col did not vanish
    if group_col and group_col not in agg_df.columns:
        raise ValueError("Group column disappeared during dedupe; aborting to prevent empty traces.")

    return agg_df

# ------------------------------ #
# Figure Building (single trace; per-bar colors; combined category + gene-only tick labels; full hover)
# ------------------------------ #
def build_fig(
    df: pd.DataFrame,
    label_col_final: str,
    value_col_final: str,
    group_col_final: Optional[str],
    top_n: int = 100,
    use_log: bool = False,
    orientation: str = "v",  # 'v' or 'h'
    log_dtick: Optional[str] = None,
    initial_zoom: Optional[int] = None,
    label_title: Optional[str] = None,
    value_title: Optional[str] = None,
    detail_cols: Optional[List[str]] = None,
    show_legend: bool = False,  # default: hide legend
) -> Tuple[go.Figure, Dict]:
    """
    Build a single-trace bar chart:
      - Each bar == one (label, group) combination.
      - Axis uses "Label | Group" internally, but tick labels show ONLY the label (e.g., Gene).
      - Colors per bar by group_col_final.
      - FULL hover shows ALL fields in detail_cols.
    """
    df_work = df
    if use_log:
        df_work = df_work[df_work[value_col_final] > 0]

    if df_work.empty or df_work[value_col_final].notna().sum() == 0:
        raise ValueError(f"No plottable numeric values found in '{value_col_final}' after cleaning/log filter.")

    # Combined category for uniqueness, but keep plain label for tick text
    if group_col_final:
        combined = (df_work[label_col_final].astype(str).str.strip() + " | " +
                    df_work[group_col_final].astype(str).str.strip())
        df_work = df_work.assign(__label_combined__=combined)
        label_axis_col = "__label_combined__"   # categories
        tick_text_col = label_col_final         # ticktext (only gene)
        page_label_title = f"{label_col_final}" # display title uses only label (gene)
    else:
        label_axis_col = label_col_final
        tick_text_col = label_col_final
        page_label_title = label_col_final

    # Sort by value and select Top N
    df_plot = df_work.sort_values(value_col_final, ascending=False).head(top_n).copy()

    # Categories and ticktext
    categories = df_plot[label_axis_col].astype(str).tolist()
    ticktext = df_plot[tick_text_col].astype(str).tolist()

    # Detail cols
    default_detail_cols = [
        label_col_final, (group_col_final or ""), "Gene name", "Tissue", "clusters_used", "median_nCPM",
        "robust_ratio", value_col_final
    ]
    if detail_cols:
        detail_cols_use = [c for c in detail_cols if c in df_plot.columns]
    else:
        detail_cols_use = [c for c in default_detail_cols if c and c in df_plot.columns]
    required_cols = [label_col_final, value_col_final]
    if group_col_final:
        required_cols.append(group_col_final)
    for rc in required_cols:
        if rc in df_plot.columns and rc not in detail_cols_use:
            detail_cols_use.insert(0, rc)
    seen = set()
    detail_cols_use = [c for c in detail_cols_use if not (c in seen or seen.add(c))]

    # Color map per group
    if group_col_final:
        colors = px.colors.qualitative.Safe
        ctypes = sorted(df_plot[group_col_final].astype(str).unique())
        cmap = {ct: colors[i % len(colors)] for i, ct in enumerate(ctypes)}
        bar_colors = df_plot[group_col_final].astype(str).map(cmap).tolist()
    else:
        cmap = {}
        bar_colors = ["#636EFA"] * len(df_plot)

    # Values
    numeric_vals = df_plot[value_col_final].astype(float).tolist()

    # FULL customdata and hovertemplate (list all detail columns)
    customdata = df_plot[detail_cols_use].values
    hover_lines = []
    if orientation == "v":
        hover_lines.append("**%{x}**")
        hover_lines.append("Value: %{y:.4g}")
    else:
        hover_lines.append("**%{y}**")
        hover_lines.append("Value: %{x:.4g}")
    for idx, col_name in enumerate(detail_cols_use):
        hover_lines.append(f"{col_name}: %{{customdata[{idx}]}}")
    hovertemplate = "<br>".join(hover_lines) + "<extra></extra>"

    # Single trace (no stacking)
    if orientation == "v":
        x_vals = df_plot[label_axis_col].astype(str).tolist()  # combined categories
        y_vals = numeric_vals
        trace = go.Bar(
            x=x_vals, y=y_vals,
            marker_color=bar_colors,
            customdata=customdata,          # full row
            hovertemplate=hovertemplate,    # full hover
            selected={"marker": {"opacity": 1.0}},
            unselected={"marker": {"opacity": 0.5}},
        )
    else:
        x_vals = numeric_vals
        y_vals = df_plot[label_axis_col].astype(str).tolist()  # combined categories
        trace = go.Bar(
            x=x_vals, y=y_vals, orientation="h",
            marker_color=bar_colors,
            customdata=customdata,          # full row
            hovertemplate=hovertemplate,    # full hover
            selected={"marker": {"opacity": 1.0}},
            unselected={"marker": {"opacity": 0.5}},
        )

    fig = go.Figure(data=[trace])

    # Titles
    label_title = label_title or page_label_title           # show only label (Gene)
    value_title = value_title or value_col_final
    page_title = f"{value_title} by {label_title}"

    fig.update_layout(
        title=page_title,
        template="plotly_white",
        height=700,
        bargap=0.2,
        margin=dict(l=80, r=40, t=70, b=130),
        hovermode="closest",
        legend_title_text=(group_col_final or ""),
        showlegend=show_legend,  # default False
        dragmode="select",
        xaxis_title=(value_title if orientation == "h" else label_title),
        yaxis_title=(label_title if orientation == "h" else value_title),
    )

    if orientation == "v":
        fig.update_xaxes(
            title_standoff=10,
            tickangle=-45,
            automargin=True,
            categoryorder="array",
            categoryarray=categories,   # categories (combined)
            tickmode="array",
            tickvals=categories,        # positions: combined
            ticktext=ticktext,          # labels: ONLY gene
        )
        if use_log:
            kwargs = dict(type="log", title=value_title + " (log scale)")
            if log_dtick in ("D1", "D2"):
                kwargs["dtick"] = log_dtick
            fig.update_yaxes(**kwargs)
        else:
            fig.update_yaxes(title_standoff=10, automargin=True)
    else:
        fig.update_yaxes(
            title_standoff=10,
            automargin=True,
            categoryorder="array",
            categoryarray=categories,   # categories (combined)
            tickmode="array",
            tickvals=categories,        # positions: combined
            ticktext=ticktext,          # labels: ONLY gene
        )
        if use_log:
            kwargs = dict(type="log", title=value_title + " (log scale)")
            if log_dtick in ("D1", "D2"):
                kwargs["dtick"] = log_dtick
            fig.update_xaxes(**kwargs)
        else:
            fig.update_xaxes(title_standoff=10, automargin=True)

    # Initial zoom (category range)
    if initial_zoom is not None:
        zoom_n = max(1, min(int(initial_zoom), len(categories)))
        if orientation == "v":
            fig.update_xaxes(range=[-0.5, zoom_n - 0.5])
        else:
            fig.update_yaxes(range=[-0.5, zoom_n - 0.5])

    # Payload rows: include details + canonical keys for value & category & ticktext(gene)
    rows = df_plot[detail_cols_use].to_dict(orient="records")
    cats = df_plot[label_axis_col].astype(str).tolist()
    ticks = df_plot[tick_text_col].astype(str).tolist()
    for i, r in enumerate(rows):
        r["__VAL__"] = numeric_vals[i]   # numeric value used for bars
        r["__CAT__"] = cats[i]           # combined category used for axis positions
        r["__TICK__"] = ticks[i]         # ONLY gene used for tick labels

    payload = {
        "rows": rows,
        "colors": cmap,  # for client-side per-group colors
        "orientation": orientation,
        "detail_cols": detail_cols_use,
        "label_col": label_col_final,
        "value_col": value_col_final,
        "value_key": "__VAL__",
        "cat_key": "__CAT__",            # canonical combined label key (positions)
        "tick_key": "__TICK__",          # canonical tick text key (ONLY gene)
        "group_col": group_col_final,    # may be None
        "label_title": label_title,
        "value_title": value_title,
        "page_title": page_title,
    }
    return fig, payload

# ------------------------------ #
# Functional UI (Group filter + Search + Export + Details + FULL hover)
# ------------------------------ #
DETAILS_SNIPPET = r"""
<div id="controls" style="margin: 0 0 12px 0; font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;">
  <label for="cellTypeSelect"><strong>Group:</strong></label>
  <select id="cellTypeSelect" aria-label="Group filter"><option value="__ALL__">All</option></select>
  <label for="searchBox"><strong>Search:</strong></label>
  <input id="searchBox" type="text" placeholder="Type to filter..." aria-label="Search by label">
  <button id="resetBtn" type="button" aria-label="Reset filters">Reset</button>
  <button id="exportBtn" type="button" aria-label="Export TSV">Export TSV</button>
</div>

<div id="rowDetails" style="font-size: 13px; color: #333;">
  Click a bar to see full row details here. Lasso/box-select points to export only selected.
</div>

<script>
(function() {
  try {
    var NL = String.fromCharCode(10);
    var TAB = String.fromCharCode(9);

    var payloadEl = document.getElementById('__payload__');
    var P = payloadEl ? JSON.parse(payloadEl.textContent || '{}') : null;
    if (!P) return;

    var plotEl = document.querySelector('div.js-plotly-plot');
    if (!plotEl) return;

    var orient = P.orientation; // 'v' or 'h'
    var LABEL = P.label_col;
    var VALUE = P.value_col;
    var VALUE_KEY = P.value_key || '__VAL__';
    var CAT_KEY = P.cat_key || '__CAT__';
    var TICK_KEY = P.tick_key || '__TICK__';
    var GROUP = P.group_col || null;
    var labelTitle = P.label_title || LABEL;
    var valueTitle = P.value_title || VALUE;

    var selectEl = document.getElementById('cellTypeSelect');
    var searchEl = document.getElementById('searchBox');
    var resetEl = document.getElementById('resetBtn');
    var exportEl = document.getElementById('exportBtn');
    var detailsEl = document.getElementById('rowDetails');

    // Populate group dropdown from color map keys (unique group values), if present
    var colorKeys = (P.colors && typeof P.colors === 'object') ? Object.keys(P.colors) : [];
    for (var i = 0; i < colorKeys.length; i++) {
      var g = colorKeys[i];
      var opt = document.createElement('option');
      opt.value = g; opt.textContent = g;
      selectEl.appendChild(opt);
    }

    var allRows = Array.isArray(P.rows) ? P.rows.slice() : [];

    // Selection tracking
    var selectedLabels = [];
    plotEl.on('plotly_selected', function(ev) {
      var pts = (ev && ev.points) ? ev.points : [];
      selectedLabels = pts.map(function(p){ return orient === 'v' ? p.x : p.y; });
    });
    plotEl.on('plotly_deselect', function(){ selectedLabels = []; });

    // Show details on click
    plotEl.on('plotly_click', function(ev) {
      try {
        var p = (ev && ev.points && ev.points[0]) ? ev.points[0] : null;
        if (!p) return;
        var cd = p.customdata || [];
        var cols = P.detail_cols || [];
        if (!cols.length) {
          detailsEl.textContent = 'No detail columns configured.';
          return;
        }
        var html = '<table style="border-collapse:collapse;">';
        for (var i = 0; i < cols.length; i++) {
          var v = cd[i];
          var vv = (typeof v === 'number') ? (Number.isFinite(v) ? v.toPrecision(4) : String(v)) : String(v != null ? v : '');
          html += '<tr><th style="text-align:left;padding:4px 8px;">' + cols[i] + '</th>' +
                  '<td style="padding:4px 8px;">' + vv + '</td></tr>';
        }
        html += '</table>';
        detailsEl.innerHTML = html;
      } catch (e) {
        console.error('click -> details error:', e);
      }
    });

    // Build full hover template listing ALL detail columns
    function makeFullHoverTemplate() {
      var dcols = P.detail_cols || [];
      var lines = [];
      if (orient === 'v') {
        lines.push("**%{x}**");
        lines.push("Value: %{y:.4g}");
      } else {
        lines.push("**%{y}**");
        lines.push("Value: %{x:.4g}");
      }
      for (var i = 0; i < dcols.length; i++) {
        lines.push(dcols[i] + ": %{customdata[" + i + "]}");
      }
      return lines.join("<br>") + "<extra></extra>";
    }

    // Render: apply filters (group + search) and rebuild single trace
    function render() {
      var term = (searchEl.value || '').toLowerCase().trim();
      var groupSel = selectEl.value;

      var rows = allRows.filter(function(r){
        var okGroup = (groupSel === '__ALL__') || !GROUP || (String(r[GROUP]) === groupSel);
        var hay = String(r[LABEL] || '').toLowerCase();
        var okSearch = (term === '') || (hay.indexOf(term) !== -1);
        return okGroup && okSearch;
      });

      var x = [], y = [], customdata = [], categories = [], ticktext = [], colors = [];
      var dcols = P.detail_cols || [];
      for (var i = 0; i < rows.length; i++) {
        var r = rows[i];
        // ALWAYS use canonical combined label for axis positions
        var cat = String(r[CAT_KEY] || '');
        var tick = String(r[TICK_KEY] || '');  // ONLY gene for tick labels
        var val = (r[VALUE_KEY] != null ? r[VALUE_KEY] : r[VALUE]);
        if (typeof val === 'string') {
          var s = val.replace(',', '.'); var vv = parseFloat(s);
          if (!Number.isNaN(vv)) val = vv;
        }
        if (orient === 'v') { x.push(cat); y.push(val); } else { x.push(val); y.push(cat); }
        categories.push(cat);
        ticktext.push(tick);

        var barColor = (GROUP && P.colors && P.colors[String(r[GROUP])]) ? P.colors[String(r[GROUP])] : '#636EFA';
        colors.push(barColor);

        // Push ALL detail columns in original order into customdata
        var cd = [];
        for (var j = 0; j < dcols.length; j++) { cd.push(r[dcols[j]]); }
        customdata.push(cd);
      }

      var trace = {
        type: 'bar',
        orientation: (orient === 'h') ? 'h' : undefined,
        x: x, y: y,
        customdata: customdata,
        marker: { color: colors },
        selected: { marker: { opacity: 1.0 } },
        unselected: { marker: { opacity: 0.5 } },
        hovertemplate: makeFullHoverTemplate(),  // FULL hover after filtering/search
      };

      var baseLayout = plotEl.layout || {};
      var baseX = baseLayout.xaxis || {};
      var baseY = baseLayout.yaxis || {};
      var layout = {
        title: baseLayout.title || P.page_title || (valueTitle + " by " + labelTitle),
        template: baseLayout.template || 'plotly_white',
        hovermode: baseLayout.hovermode || 'closest',
        bargap: (typeof baseLayout.bargap !== 'undefined') ? baseLayout.bargap : 0.2,
        showlegend: false,  // legend hidden
        dragmode: baseLayout.dragmode || 'select',
        margin: baseLayout.margin || {l:80, r:40, t:70, b:130},
        legend: { tracegroupgap: 0 },
        xaxis: Object.assign({}, baseX, {
          title: { text: (orient === 'h') ? (valueTitle || VALUE) : (labelTitle || LABEL) },
          automargin: true,
          title_standoff: 10,
          tickangle: (orient === 'v') ? -45 : baseX.tickangle,
          categoryorder: (orient === 'v') ? 'array' : baseX.categoryorder,
          categoryarray: (orient === 'v') ? categories : baseX.categoryarray,
          tickmode: 'array',
          tickvals: categories,
          ticktext: ticktext,
        }),
        yaxis: Object.assign({}, baseY, {
          title: { text: (orient === 'h') ? (labelTitle || LABEL) : (valueTitle || VALUE) },
          automargin: true,
          title_standoff: 10,
          categoryorder: (orient === 'h') ? 'array' : baseY.categoryorder,
          categoryarray: (orient === 'h') ? categories : baseY.categoryarray,
          tickmode: (orient === 'h') ? 'array' : baseY.tickmode,
          tickvals: (orient === 'h') ? categories : baseY.tickvals,
          ticktext: (orient === 'h') ? ticktext : baseY.ticktext,
        }),
      };

      Plotly.react(plotEl, [trace], layout);
    }

    selectEl.addEventListener('change', render);
    searchEl.addEventListener('input', render);
    resetEl.addEventListener('click', function(){
      selectEl.value = '__ALL__';
      searchEl.value = '';
      selectedLabels = [];
      render();
    });

    exportEl.addEventListener('click', function(){
      // Export TSV (filtered OR selected)
      var term = (searchEl.value || '').toLowerCase().trim();
      var groupSel = selectEl.value;
      var rows = allRows.filter(function(r){
        var okGroup = (groupSel === '__ALL__') || !GROUP || (String(r[GROUP]) === groupSel);
        var hay = String(r[LABEL] || '').toLowerCase();
        var okSearch = (term === '') || (hay.indexOf(term) !== -1);
        return okGroup && okSearch;
      });

      if (selectedLabels && selectedLabels.length > 0) {
        var sset = new Set(selectedLabels.map(String));
        rows = rows.filter(function(r){
          var cat = String(r[CAT_KEY] || '');
          return sset.has(cat);
        });
      }

      var cols = P.detail_cols || (rows.length ? Object.keys(rows[0]) : []);
      var header = cols.join(TAB);
      var lines = rows.map(function(r){
        return cols.map(function(c){ return (r[c] != null ? String(r[c]) : ''); }).join(TAB);
      });
      var tsv = [header].concat(lines).join(NL);
      var blob = new Blob([tsv], { type: 'text/tab-separated-values;charset=utf-8' });
      var a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'export.tsv';
      a.click();
      URL.revokeObjectURL(a.href);
    });

    // Initial render
    render();
  } catch (e) {
    console.error('UI init error:', e);
  }
})();
</script>
"""

# ------------------------------ #
# HTML Saving
# ------------------------------ #
def save_html(fig: go.Figure, out_path: str, payload: dict, self_contained: bool = False, lang_code: str = "en"):
    include_js = True if self_contained else "cdn"
    html = fig.to_html(full_html=True, include_plotlyjs=include_js)

    html = html.lstrip(" \ufeff\r\n\t")
    if not html.startswith("<!DOCTYPE html>"):
        html = "<!DOCTYPE html>\n" + html

    if re.search(r"(?is)<html[^>]*>", html):
        if not re.search(r"(?is)<html[^>]*\blang\s*=", html):
            html = re.sub(r"(?is)<html(\s*)>", f'<html lang="{lang_code}"\\1>', html, count=1)
    else:
        html = f'<!DOCTYPE html>\n<html lang="{lang_code}">\n{html}\n</html>'

    head_snippet = (
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f'<title>{payload.get("page_title") or ((payload.get("value_title") or "Value") + " by " + (payload.get("label_title") or "Label"))}</title>\n'
    )

    if re.search(r"(?is)<head\s*>", html):
        html = re.sub(r"(?is)<head\s*>", "<head>\n" + head_snippet + "\n", html, count=1)
    else:
        html = re.sub(r"(?is)(<html[^>]*>)", r"\1\n<head>\n" + head_snippet + "\n</head>\n", html, count=1)

    payload_json = json.dumps(payload, ensure_ascii=False).replace("</script", "<\\/script")
    data_block = f'<script type="application/json" id="__payload__">{payload_json}</script>\n'

    if re.search(r"(?is)</body>", html):
        html = re.sub(r"(?is)</body>", data_block + DETAILS_SNIPPET + "\n</body>", html, count=1)
    else:
        html = html + "\n" + data_block + DETAILS_SNIPPET + "\n"

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)

# ------------------------------ #
# CLI
# ------------------------------ #
def main():
    ap = argparse.ArgumentParser(
        description="Interactive bar chart (table → HTML) with combined categories, gene-only tick labels, per-group colors, full-row hover, details, search, and TSV export."
    )
    ap.add_argument("--file", "-f", required=True, help="Input table (TSV/CSV/etc.)")
    ap.add_argument("--out", "-o", default="interactive_markers.html", help="Output HTML file")
    ap.add_argument("--top", "-n", type=int, default=100, help="Top N rows by value column")

    # Scale options
    ap.add_argument("--log", action="store_true", help="Use log scale on numeric axis")
    ap.add_argument("--linear", action="store_true", help="Use linear scale on numeric axis")
    ap.add_argument("--log-digits", choices=["D1", "D2"], help="Log-axis minor digits: D1 (all) or D2 (2 & 5)")

    # Orientation
    ap.add_argument("--horizontal", action="store_true", help="Use horizontal bars")

    # Packaging
    ap.add_argument("--self-contained", action="store_true", help="Embed plotly.js for fully offline HTML")
    ap.add_argument("--lang", default="en", help="HTML lang attribute (e.g., 'en', 'en-CA')")

    # Initial viewport
    ap.add_argument("--initial-zoom", type=int, default=100, help="Initial number of bars to show on load")

    # CSV/TSV
    ap.add_argument("--sep", help="Field separator (e.g., '\\t' or ','); otherwise auto-detected")

    # Column flexibility
    ap.add_argument("--x-col", help="Column to plot on the X axis")
    ap.add_argument("--y-col", help="Column to plot on the Y axis")
    ap.add_argument("--label-col", help="Categorical label column")
    ap.add_argument("--value-col", help="Numeric value column")
    ap.add_argument("--group-col", help="Grouping column (used for colors and combined categories)")
    ap.add_argument("--search-col", help="Column used by the search box (defaults to label column)")
    ap.add_argument("--details", nargs="*", help="Extra columns for hover/details/export")

    # Duplicate handling
    ap.add_argument("--dedupe-policy", choices=["error", "max", "mean", "median", "first"], default="error",
                    help="Collapse duplicate (label, group) rows before plotting (default: error).")

    # Legend control (hidden by default)
    ap.add_argument("--show-legend", action="store_true", help="Show the legend (hidden by default)")

    args = ap.parse_args()

    # Scale choice: DEFAULT linear unless --log is provided
    if args.log and args.linear:
        raise SystemExit("Choose either --log or --linear (not both).")
    use_log = True if args.log else False
    orientation = "h" if args.horizontal else "v"

    # Load and resolve columns
    df = load_table(args.file, sep=_infer_sep(args.file, args.sep))
    try:
        label_col_final, value_col_final, group_col_final, _ = resolve_columns(
            df,
            x_col=args.x_col,
            y_col=args.y_col,
            label_col=args.label_col,
            value_col=args.value_col,
            group_col=args.group_col,
            search_col=args.search_col,
            orientation=orientation,
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        raise SystemExit(1)

    # Dedupe BEFORE plotting
    try:
        df = dedupe_for_plot(
            df,
            label_col_final,
            value_col_final,
            group_col_final,
            policy=args.dedupe_policy,
            keep_cols=args.details
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        raise SystemExit(1)

    if df.empty:
        raise SystemExit("No rows to plot after dedupe/filter; check --dedupe-policy or input data.")

    # Axis titles
    label_title = label_col_final
    value_title = value_col_final

    # Build figure & payload (legend off by default)
    fig, payload = build_fig(
        df,
        label_col_final=label_col_final,
        value_col_final=value_col_final,
        group_col_final=group_col_final,
        top_n=args.top,
        use_log=use_log,
        orientation=orientation,
        log_dtick=args.log_digits,
        initial_zoom=args.initial_zoom,
        label_title=label_title,
        value_title=value_title,
        detail_cols=args.details,
        show_legend=args.show_legend,
    )

    save_html(fig, args.out, payload, self_contained=args.self_contained, lang_code=args.lang)
    print(f"Saved: {args.out}")

if __name__ == "__main__":
    main()

```
# METHOD 02 - AUROC SCORES

This calculates a little more complicated measure- auroc score instead enrichment scores

This uses the same input file 

Just like enrichment scores, this calculates auroc scores with a python script.
Because it takes more computing power, I decided to submit a job on computecanada.

you can run it like following

run_auroc.sbatch
```
#!/bin/bash
#SBATCH --job-name=auroc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=logs/auroc_%A_%a.out
#SBATCH --error=logs/auroc_%A_%a.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# ---------------------------
# Environment
# ---------------------------
module load gcc arrow
module load python

python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate

python compute_auroc_cli.py \
  --file rna_single_cell_cluster.tsv \
  --out-all all_gene_cell_auroc.tsv \
  --out-top top100_auroc.tsv \
  --auroc-min 0.92 \
  --median-min 2.0 \
  --ratio-min 3.0 \
  --clusters-min 2 \
  --alpha 0.1
```
--file -input file

--out-all -output file name for all files 

--out-top -output file name for top 100 file

--auroc-min -minimum AUROC to keep. 0.90: Keeps gene–cell-type pairs with strong separation from other cell types in the same tissue.

--median-min -minimum median nCPM in target cell type. 1.0: Avoids keeping “specific but trivially low” expression.

--ratio-min -minimum (median_target+α)/(median_others+α).2.0: Requires the median in the target cell type to be at least ~2× higher than the median of others (with a small pseudocount α=0.1 to stabilize zeros).

--clusters-min -minimum number of clusters in target cell type. 2: Ensures there’s at least minimal within-cell-type replication.

--alpha -sets the pseudocount added when calculating the robust ratio to stabilize division by very small or zero medians, which is important because it prevents inflated ratios and ensures more reliable specificity scoring in sparse


Following is the python script

compute_auroc_cli.py
```
#!/usr/bin/env python3
"""
Compute AUROC specificity per (Gene × Tissue × Cell type) using cluster-level nCPM,
apply thresholds, and save filtered outputs. Designed for the file:
    rna_single_cell_cluster.tsv

Expected columns:
- Gene
- Gene name
- Tissue
- Cluster
- Cell type
- Read count
- nCPM

CLI arguments let you set thresholds, input/output, and pseudocount α for robust ratio.
"""

import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score


def compute_auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """
    Compute AUROC given scores and binary labels (1=target cell type, 0=others).
    Returns NaN if AUROC is not computable (e.g., only one class).
    If all scores are constant, returns 0.5 (no discrimination).
    """
    pos = int(labels.sum())
    n = len(labels)
    if pos == 0 or pos == n:
        return np.nan
    if np.allclose(scores, scores[0]):
        return 0.5
    try:
        return float(roc_auc_score(labels, scores))
    except Exception:
        return np.nan


def main():
    ap = argparse.ArgumentParser(
        description="Compute AUROC per (Gene × Tissue × Cell type) with thresholds.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--file", "-f", default="rna_single_cell_cluster.tsv",
                    help="Input TSV file with cluster-level data.")
    ap.add_argument("--out-all", default="all_gene_cell_auroc.tsv",
                    help="Output TSV for all passing rows after filtering.")
    ap.add_argument("--out-top", default="top100_auroc.tsv",
                    help="Output TSV for top-100 rows after filtering.")
    ap.add_argument("--auroc-min", type=float, default=0.90,
                    help="Minimum AUROC to keep a row.")
    ap.add_argument("--median-min", type=float, default=1.0,
                    help="Minimum median nCPM in target cell type.")
    ap.add_argument("--ratio-min", type=float, default=2.0,
                    help="Minimum robust ratio (median_target+α)/(median_others+α).")
    ap.add_argument("--clusters-min", type=int, default=2,
                    help="Minimum clusters in target cell type.")
    ap.add_argument("--alpha", type=float, default=0.1,
                    help="Pseudocount α for robust ratio stabilization.")
    ap.add_argument("--keep-constant-auc05", action="store_true",
                    help="If set, constant-score cases keep AUROC=0.5; otherwise they are NaN.")
    args = ap.parse_args()

    # -------- Load data --------
    in_path = args.file
    df = pd.read_csv(in_path, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    required_cols = ["Gene", "Gene name", "Tissue", "Cluster", "Cell type", "nCPM"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Input missing required columns: {missing}")

    df["nCPM"] = pd.to_numeric(df["nCPM"], errors="coerce")
    df = df.dropna(subset=["nCPM"])

    # -------- Main loop: per (Tissue × Gene × Cell type) --------
    rows = []
    group_cols = ["Tissue", "Gene"]

    for (tissue, gene), subdf in df.groupby(group_cols, sort=False):
        ctypes = subdf["Cell type"].dropna().unique().tolist()
        if len(ctypes) < 2:
            # No contrast possible
            continue

        ct_medians = subdf.groupby("Cell type")["nCPM"].median().to_dict()
        ct_counts = subdf.groupby("Cell type")["nCPM"].size().to_dict()

        for ct in ctypes:
            labels = (subdf["Cell type"] == ct).astype(int).values
            scores = subdf["nCPM"].values

            # AUROC computation with optional handling of constant scores
            pos = int(labels.sum())
            n = len(labels)
            if pos == 0 or pos == n:
                auc = np.nan
            elif np.allclose(scores, scores[0]):
                auc = 0.5 if args.keep_constant_auc05 else np.nan
            else:
                auc = compute_auroc(scores, labels)

            # Robust ratio
            target_med = ct_medians.get(ct, np.nan)
            other_meds = [m for k, m in ct_medians.items() if k != ct and not pd.isna(m)]
            other_med = np.median(other_meds) if len(other_meds) > 0 else np.nan

            robust_ratio = np.nan
            if not pd.isna(target_med) and not pd.isna(other_med):
                robust_ratio = (target_med + args.alpha) / (other_med + args.alpha)

            single_ct_flag = (len(ctypes) == 1)

            rows.append({
                "Gene": gene,
                "Gene name": subdf["Gene name"].iloc[0],
                "Tissue": tissue,
                "Cell type": ct,
                "clusters_used": int(ct_counts.get(ct, 0)),
                "median_nCPM": float(target_med) if not pd.isna(target_med) else np.nan,
                "robust_ratio": float(robust_ratio) if not pd.isna(robust_ratio) else np.nan,
                "AUROC": float(auc) if not pd.isna(auc) else np.nan,
                "single_cell_type_gene": single_ct_flag
            })

    auroc_df = pd.DataFrame(rows)

    # -------- Apply thresholds --------
    def passes_thresholds(row) -> bool:
        if pd.isna(row["AUROC"]) or pd.isna(row["median_nCPM"]) or pd.isna(row["robust_ratio"]):
            return False
        if row["clusters_used"] < args.clusters_min:
            return False
        if row["AUROC"] < args.auroc_min:
            return False
        if row["median_nCPM"] < args.median_min:
            return False
        if row["robust_ratio"] < args.ratio_min:
            return False
        return True

    filtered_df = auroc_df[auroc_df.apply(passes_thresholds, axis=1)].copy()

    # Sort by AUROC (desc), then robust_ratio, then median_nCPM
    filtered_df = filtered_df.sort_values(
        by=["AUROC", "robust_ratio", "median_nCPM"],
        ascending=[False, False, False]
    )
    # -------- Save outputs --------
    filtered_df.to_csv(args.out_all, sep="\t", index=False)
    filtered_df.head(100).to_csv(args.out_top, sep="\t", index=False)

    print("\033[33mSaved files:\033[0m", args.out_all, ",", args.out_top)
    print(
        f"Thresholds: AUROC ≥ {args.auroc_min}, median_nCPM ≥ {args.median_min}, "
        f"robust_ratio ≥ {args.ratio_min}, clusters_used ≥ {args.clusters_min}, α={args.alpha}"
    )


if __name__ == "__main__":
    main()
```
# METHOD 03 - Enrichment score with more robust statistical backup systems to filter useful data

This method is much more complicated and takes much more computing power. Therefore I decided to split the input file into multiple files and process them simultaniously as a job array.

Therefore, first you have to split this file, keeping same groups together

Group-aware splitting (use --group-by "Gene,Cell type,Tissue" or "Cell type,Tissue"): keeps each unit intact, so within-tissue aggregation and the bootstrap performed in recommend_min_k() are not starved of clusters or mixed across chunks. Partial groups in different files can bias effective_clusters, CV, relative Δ, and rank ρ; grouping prevents that. Global recompute after merging filtered clusters: Your enrichment is defined as a ratio to other cell types per gene; any chunk-wise run sees only a subset of cell types/genes and can inflate/deflate denominators. Merging all *_filtered_clusters.tsv and running enrichment once globally restores the correct background across all cell types. Header verification & manifest: ensures consistent schema and lets you audit exactly which groups landed in each part.

You can run this splitting script like following
```
python split_tsv_robust.py split \
  --input rna_single_cell_cluster.tsv \
  --parts 20 \
  --group-by "Gene,Cell type,Tissue" \
  --outdir splitted_input --prefix part_ --gzip-out \
  --manifest split_manifest.json
```

This is the script used for that

split_tsv_robust.py split
```

#!/usr/bin/env python3
"""
Robust TSV splitter for single-cell cluster data (and similar large TSVs).

Focus: safeguards so splitting does **not** affect downstream enrichment calculations.

Key features:
- Header preserved in every part.
- Two splitting modes:
  1) block (raw even rows per part),
  2) group-aware (keeps complete groups together in one part), with greedy balancing.
- Recommended group key for your enrichment workflow: **Gene, Cell type, Tissue**
  (or at minimum **Cell type, Tissue**) so bootstrap and aggregation remain intact within parts.
- Optional gzip input/output; auto-detect by filename suffix.
- Zero-padded numeric suffix in filenames (01..N, or wider as needed).
- Safeguards: caps parts if total rows < parts; verifies required columns; writes a manifest JSON
  with per-part row counts and group coverage; optional dry-run.
- Merge subcommand to concatenate outputs (keeping one header) and schema checks.

Usage examples:

Split by groups to preserve units relevant to enrichment (recommended):
  python split_tsv_robust.py split \
    --input rna_single_cell_cluster.tsv \
    --parts 20 \
    --group-by "Gene,Cell type,Tissue" \
    --outdir $SCRATCH/splitted_input --prefix part_ --gzip-out

Block-wise split (simple):
  python split_tsv_robust.py split --input rna_single_cell_cluster.tsv --parts 20

Dry-run (plan only, no files):
  python split_tsv_robust.py split --input rna_single_cell_cluster.tsv --parts 20 --group-by "Cell type,Tissue" --dry-run

Merge chunk outputs back together (one header):
  python split_tsv_robust.py merge \
    --pattern "$SCRATCH/enrich_parts/adjusted_*_final_enrichment.tsv" \
    --output merged_final_enrichment.tsv

Two-pass safeguard (recommended end-to-end):
  1) Run array jobs on parts to produce per-chunk "*_filtered_clusters.tsv".
  2) Merge all filtered cluster chunks: `split_tsv_robust.py merge --pattern '.../*_filtered_clusters.tsv' --output merged_filtered_clusters.tsv`
  3) Run **global** enrichment once on the merged filtered clusters to avoid partial baselines:
     `python min_clusters_and_enrichment.py --clusters merged_filtered_clusters.tsv --out-prefix global_adjusted ...`

Required columns (for your enrichment script):
  Gene, Gene name, Tissue, Cluster, Cell type, Read count, nCPM
"""
import argparse
import gzip
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# -------------------------
# Utilities
# -------------------------
def open_maybe_gzip(path: Path, mode: str = 'rt', encoding: Optional[str] = 'utf-8'):
    """Open plain text or gzip based on suffix. Accepts text ('rt','wt') or binary modes."""
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, mode)
    else:
        if 't' in mode:
            return path.open(mode, encoding=encoding)
        return path.open(mode)

def parse_group_by(group_by: Optional[str]) -> List[str]:
    if not group_by:
        return []
    cols = [c.strip() for c in group_by.split(',') if c.strip()]
    return cols

REQUIRED_COLS = ["Gene", "Gene name", "Tissue", "Cluster", "Cell type", "Read count", "nCPM"]

def check_required_columns(header: str, required_cols: List[str]) -> None:
    cols = header.rstrip('\n').split('\t')
    missing = [c for c in required_cols if c not in cols]
    if missing:
        raise SystemExit(f"Input header missing required columns: {missing}. Found: {cols}")

def zero_pad_width(parts: int) -> int:
    return max(2, len(str(parts)))

def make_out_name(outdir: Path, prefix: str, idx: int, pad: int, gzip_out: bool) -> Path:
    suffix = '.tsv.gz' if gzip_out else '.tsv'
    return outdir / f"{prefix}{idx:0{pad}d}{suffix}"

# -------------------------
# Block-wise split
# -------------------------
def split_block(input_path: Path, parts: int, outdir: Path, prefix: str, gzip_out: bool,
                rows_per_part_override: Optional[int], dry_run: bool) -> Dict:
    # First pass: read header + count rows
    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        if not header:
            raise SystemExit("Empty file or missing header")
        total_rows = sum(1 for _ in f)
    if total_rows == 0:
        raise SystemExit("No data rows found (only header)")

    # Safeguard: cap parts to total_rows or compute from override
    if rows_per_part_override and rows_per_part_override > 0:
        parts = math.ceil(total_rows / rows_per_part_override)
    parts = min(parts, total_rows)

    rows_per_part = math.ceil(total_rows / parts)
    pad = zero_pad_width(parts)

    plan = {
        'mode': 'block',
        'parts': parts,
        'rows_per_part': rows_per_part,
        'total_rows': total_rows,
        'outdir': str(outdir),
        'prefix': prefix,
        'files': []
    }

    if dry_run:
        counts = [rows_per_part] * (parts - 1) + [total_rows - rows_per_part * (parts - 1)]
        for i, cnt in enumerate(counts, start=1):
            plan['files'].append({'name': make_out_name(outdir, prefix, i, pad, gzip_out).name,
                                  'rows_including_header': cnt + 1})
        return plan

    outdir.mkdir(parents=True, exist_ok=True)

    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        part_idx = 1
        rows_written_in_part = 0
        out_file_path = make_out_name(outdir, prefix, part_idx, pad, gzip_out)
        out_file = open_maybe_gzip(out_file_path, 'wt')
        out_file.write(header)
        lines_in_current_file = 1
        per_file_counts = []

        for line in f:
            if rows_written_in_part >= rows_per_part and part_idx < parts:
                out_file.close()
                per_file_counts.append((out_file_path.name, lines_in_current_file))
                part_idx += 1
                rows_written_in_part = 0
                out_file_path = make_out_name(outdir, prefix, part_idx, pad, gzip_out)
                out_file = open_maybe_gzip(out_file_path, 'wt')
                out_file.write(header)
                lines_in_current_file = 1
            out_file.write(line)
            rows_written_in_part += 1
            lines_in_current_file += 1
        out_file.close()
        per_file_counts.append((out_file_path.name, lines_in_current_file))

    plan['files'] = [{'name': n, 'rows_including_header': c} for n, c in per_file_counts]
    return plan

# -------------------------
# Group-aware split
# -------------------------
def split_groups(input_path: Path, parts: int, outdir: Path, prefix: str, gzip_out: bool,
                 group_cols: List[str], required_cols: List[str], dry_run: bool) -> Dict:
    # First pass: header + group sizes
    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        if not header:
            raise SystemExit("Empty file or missing header")
        # Verify required columns exist (so downstream enrichment won't fail)
        check_required_columns(header, required_cols)
        cols = header.rstrip('\n').split('\t')
        col_idx = {c: i for i, c in enumerate(cols)}
        for gc in group_cols:
            if gc not in col_idx:
                raise SystemExit(f"--group-by column not found in header: {gc}. Available: {cols}")
        total_rows = 0
        group_counts: Dict[Tuple, int] = {}
        for line in f:
            total_rows += 1
            fields = line.rstrip('\n').split('\t')
            key = tuple(fields[col_idx[c]] for c in group_cols)
            group_counts[key] = group_counts.get(key, 0) + 1
    if total_rows == 0:
        raise SystemExit("No data rows found (only header)")

    # Safeguard: cap parts to min(total_rows, number_of_groups)
    num_groups = len(group_counts)
    parts = min(parts, max(1, num_groups))

    # Greedy bin packing: largest groups first to balance rows across parts
    sorted_groups = sorted(group_counts.items(), key=lambda kv: kv[1], reverse=True)
    bins = [{'rows': 0, 'groups': []} for _ in range(parts)]
    for key, gsize in sorted_groups:
        target = min(range(parts), key=lambda i: bins[i]['rows'])
        bins[target]['rows'] += gsize
        bins[target]['groups'].append(key)

    pad = zero_pad_width(parts)

    plan = {
        'mode': 'group',
        'parts': parts,
        'total_rows': total_rows,
        'num_groups': num_groups,
        'group_cols': group_cols,
        'outdir': str(outdir),
        'prefix': prefix,
        'files': [{'name': make_out_name(outdir, prefix, i+1, pad, gzip_out).name,
                   'assigned_groups': len(bins[i]['groups']),
                   'planned_rows_excluding_header': bins[i]['rows']}
                  for i in range(parts)]
    }

    if dry_run:
        return plan

    outdir.mkdir(parents=True, exist_ok=True)

    # Second pass: route each line to its assigned bin
    writers = []
    for i in range(parts):
        p = make_out_name(outdir, prefix, i+1, pad, gzip_out)
        fh = open_maybe_gzip(p, 'wt')
        writers.append({'path': p, 'fh': fh, 'count': 0})

    group_to_bin: Dict[Tuple, int] = {}
    for i in range(parts):
        for g in bins[i]['groups']:
            group_to_bin[g] = i

    with open_maybe_gzip(input_path, 'rt') as f:
        header = f.readline()
        for w in writers:
            w['fh'].write(header)
            w['count'] = 1
        cols = header.rstrip('\n').split('\t')
        col_idx = {c: i for i, c in enumerate(cols)}
        for line in f:
            fields = line.rstrip('\n').split('\t')
            key = tuple(fields[col_idx[c]] for c in group_cols)
            bin_idx = group_to_bin.get(key)
            if bin_idx is None:
                raise SystemExit(f"Internal error: group {key} not assigned to any bin")
            w = writers[bin_idx]
            w['fh'].write(line)
            w['count'] += 1

    # Close writers and collect counts (include assigned_groups)
    per_file_infos = []
    for i, w in enumerate(writers):
        w['fh'].close()
        per_file_infos.append({
            'name': w['path'].name,
            'rows_including_header': w['count'],
            'assigned_groups': len(bins[i]['groups'])
        })

    plan['files'] = per_file_infos
    return plan

# -------------------------
# Merge helper (concatenate with single header)
# -------------------------
def merge_files(pattern: str, output: Path, gzip_out: bool = False) -> Dict:
    import glob
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No files matched pattern: {pattern}")

    def read_header(fp: str) -> str:
        p = Path(fp)
        with open_maybe_gzip(p, 'rt') as f:
            hdr = f.readline()
            if not hdr:
                raise SystemExit(f"Empty file or missing header: {fp}")
            return hdr

    main_header = read_header(files[0])

    # Verify all headers match exactly
    for fp in files[1:]:
        hdr = read_header(fp)
        if hdr.rstrip('\n') != main_header.rstrip('\n'):
            raise SystemExit(f"Header mismatch between {files[0]} and {fp}\n{main_header}\n!=\n{hdr}")

    out_path = output
    if gzip_out and not str(out_path).endswith('.gz'):
        out_path = Path(str(out_path) + '.gz')

    with open_maybe_gzip(out_path, 'wt') as out:
        out.write(main_header)
        total_lines = 1
        for fp in files:
            p = Path(fp)
            with open_maybe_gzip(p, 'rt') as f:
                _ = f.readline()  # skip header
                for line in f:
                    out.write(line)
                    total_lines += 1

    return {'output': str(out_path), 'lines_including_header': total_lines, 'merged_files': files}

# -------------------------
# CLI
# -------------------------
def main():
    ap = argparse.ArgumentParser(description='Robust TSV splitter/merger with safeguards for enrichment workflows')
    sub = ap.add_subparsers(dest='cmd', required=True)

    sp = sub.add_parser('split', help='Split a large TSV into parts')
    sp.add_argument('--input', '-i', required=True, help='Input TSV file path (.tsv or .tsv.gz)')
    sp.add_argument('--parts', '-n', type=int, default=20, help='Number of parts to split into (default: 20)')
    sp.add_argument('--rows-per-part', type=int, default=None, help='Override rows per part (block mode only)')
    sp.add_argument('--outdir', '-o', default='splitted_input', help='Output directory (default: splitted_input)')
    sp.add_argument('--prefix', default='part_', help='Output filename prefix (default: part_)')
    sp.add_argument('--gzip-out', action='store_true', help='Write gzip-compressed parts (.tsv.gz)')
    sp.add_argument('--group-by', default='', help='Comma-separated column names to keep rows grouped (e.g., "Gene,Cell type,Tissue")')
    sp.add_argument('--require-cols', default=','.join(REQUIRED_COLS), help='Comma-separated required columns to verify in header')
    sp.add_argument('--manifest', default='split_manifest.json', help='Where to write JSON manifest of the split plan/results')
    sp.add_argument('--dry-run', action='store_true', help='Plan only; do not write files')

    mp = sub.add_parser('merge', help='Concatenate TSVs (one header)')
    mp.add_argument('--pattern', required=True, help='Glob pattern for input files to merge')
    mp.add_argument('--output', required=True, help='Output merged file path (.tsv or .tsv.gz)')
    mp.add_argument('--gzip-out', action='store_true', help='Write gzip-compressed output (adds .gz if missing)')

    args = ap.parse_args()

    if args.cmd == 'split':
        input_path = Path(args.input)
        if not input_path.exists():
            raise SystemExit(f"Input file not found: {input_path}")
        outdir = Path(args.outdir)
        required_cols = [c.strip() for c in args.require_cols.split(',') if c.strip()]
        group_cols = parse_group_by(args.group_by)

        # Choose mode
        if group_cols:
            plan = split_groups(input_path, parts=args.parts, outdir=outdir, prefix=args.prefix,
                                gzip_out=args.gzip_out, group_cols=group_cols,
                                required_cols=required_cols, dry_run=args.dry_run)
        else:
            plan = split_block(input_path, parts=args.parts, outdir=outdir, prefix=args.prefix,
                               gzip_out=args.gzip_out, rows_per_part_override=args.rows_per_part, dry_run=args.dry_run)

        # Write manifest
        manifest_path = Path(args.manifest)
        data = {
            'input': str(input_path),
            'mode': plan['mode'],
            'parts': plan['parts'],
            'outdir': plan['outdir'],
            'prefix': plan['prefix'],
            'total_rows': plan.get('total_rows'),
            'num_groups': plan.get('num_groups'),
            'group_cols': plan.get('group_cols'),
            'files': plan['files']
        }
        with manifest_path.open('w', encoding='utf-8') as mf:
            json.dump(data, mf, indent=2)
        print(f"Split {'planned' if args.dry_run else 'completed'}: {plan['parts']} parts to '{outdir}'. Manifest: {manifest_path}")
        for info in plan['files']:
            name = info['name']
            rows = info.get('rows_including_header', info.get('planned_rows_excluding_header', 0) + 1)
            if plan['mode'] == 'group' and 'assigned_groups' in info:
                print(f"{name}: {rows} lines (groups: {info['assigned_groups']})")
            else:
                print(f"{name}: {rows} lines")

    elif args.cmd == 'merge':
        output = Path(args.output)
        res = merge_files(args.pattern, output, gzip_out=args.gzip_out)
        print(f"Merged {len(res['merged_files'])} files into: {res['output']} ({res['lines_including_header']} lines)")

if __name__ == '__main__':
    main()
```
After splitting, you can run the calculating scripts


Just like previous efforts, this calculations also has a python script that does the calculations for you

You can run it like following
```
sbatch enrich_array.sbatch
```

in addition to running it like sbatch with defaults, you also have the option to customize flags at submit time like

```
export B=200
export CV=0.15
export RELDELTA=0.05
export RHO=0.85
export MIN_K=3
export TISSUE_WEIGHTING=weighted
export PSEUDOCOUNT=0.001

sbatch --array=1-${PARTS} enrich_array.sbatch
```


enrich_array.sbatch
```
#!/bin/bash
#SBATCH --job-name=enrich_array
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=20:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/enrich_%A_%a.out
#SBATCH --error=logs/enrich_%A_%a.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# ---------------------------
# Configurable parameters (defaults target the current directory)
# Override at submit-time if needed, e.g.:
#   PARTS=50 INDIR=./chunks PREFIX=chunk_ OUTDIR=./parts_out sbatch enrich_array.sbatch
# ---------------------------
PARTS=${PARTS:-20}                # total number of chunks (must match your split)
INDIR=${INDIR:-./splitted_input}  # directory holding part_XX.tsv(.gz) relative to CWD
PREFIX=${PREFIX:-part_}           # chunk filename prefix (e.g., part_01.tsv.gz)
OUTDIR=${OUTDIR:-./enrich_parts}  # where per-chunk outputs will be written (relative to CWD)

# ---------------------------
# Environment
# ---------------------------
module load gcc arrow
module load python

python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate

mkdir -p "$OUTDIR" logs

# ---------------------------
# Map array index -> chunk file
# ---------------------------
PADDED=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})

# Prefer .tsv.gz, fall back to .tsv if .gz doesn’t exist
INPUT_GZ="${INDIR}/${PREFIX}${PADDED}.tsv.gz"
INPUT_TSV="${INDIR}/${PREFIX}${PADDED}.tsv"
if [ -s "$INPUT_GZ" ]; then
  INPUT="$INPUT_GZ"
elif [ -s "$INPUT_TSV" ]; then
  INPUT="$INPUT_TSV"
else
  echo "[array] Missing input for index ${SLURM_ARRAY_TASK_ID}:"
  echo "  tried: $INPUT_GZ and $INPUT_TSV" >&2
  exit 2
fi

PREFIX_OUT="${OUTDIR}/adjusted_${PADDED}"

echo "[array] task ${SLURM_ARRAY_TASK_ID}/${PARTS}"
echo "[array] input:      $INPUT"
echo "[array] out-prefix: $PREFIX_OUT"

# ---------------------------
# Run enrichment script per chunk
# ---------------------------
python min_clusters_and_enrichment.py \
  --clusters "$INPUT" \
  --out-prefix "$PREFIX_OUT" \
  --B 100 \
  --cv 0.20 \
  --reldelta 0.10 \
  --rho 0.90 \
  --min-k 2 \
  --tissue-weighting weighted \
  --pseudocount 0.01

```
This array only works properly if you submit it as sbatch enrich_array.sbatch. If you use bash enrich_array.sbatch to test this first, you should use something like following
```

# Emulate array task 1 → should target part_01.tsv(.gz)
export SLURM_ARRAY_TASK_ID=1
export PARTS=20
export INDIR=./splitted_input
export PREFIX=part_
export OUTDIR=./enrich_parts

# Optional: show what the script will try
bash -x enrich_array.sbatch | sed -n '1,120p'
```
--B -sets the number of bootstrap/resampling iterations, which is important for testing the stability and robustness of the enrichment results.

--reldelta -Defines the minimum relative change threshold used for convergence or filtering, helping prevent overfitting and ensuring only meaningful improvements are considered.

--rho -Specifies the similarity/correlation cutoff for grouping or merging signals/clusters, which is crucial to avoid combining weakly related patterns and to preserve biological interpretability.

--min-k -Sets the minimum number of clusters (k) to evaluate, ensuring the method doesn’t collapse to trivial one‑cluster solutions and captures real heterogeneity.

--tissue-weighting -Chooses how to weight tissues (e.g., weighted vs. unweighted) so larger or more sampled tissues don’t dominate the analysis, improving fairness across conditions.

--pseudocount -Adds a small constant to counts before ratios/logs, preventing divide‑by‑zero and stabilizing estimates for sparse data.

--clusters — Specifies the number of clusters to generate during analysis, which is important for controlling granularity and ensuring patterns are grouped meaningfully without over- or under-segmentation.

--out-prefix — Sets the prefix for all output files, which is important for organizing results and avoiding overwriting when running multiple analyses.

Following is the python script to run

min_clusters_and_enrichment.py
```py
#!/usr/bin/env python3
"""
min_clusters_and_enrichment.py

Purpose:
  - Use ONLY rna_single_cell_cluster.tsv (columns: Gene, Gene name, Tissue, Cluster, Cell type, Read count, nCPM).
  - Recommend the minimum number of clusters per (Cell type, Tissue) for stable enrichment (weighted by Read count).
  - Apply these replication adjustments (exclude unstable pairs) and recompute proper enrichment values.

Outputs:
  - min_clusters_per_celltype.tsv : recommendations per (Cell type, Tissue) with diagnostics.
  - filtered_clusters.tsv         : cluster-level rows retained for enrichment after applying recommendations.
  - final_enrichment.tsv          : gene × cell type table with adjusted enrichment.

CLI example:
  python min_clusters_and_enrichment.py \
    --clusters rna_single_cell_cluster.tsv \
    --out-prefix adjusted \
    --B 100 --cv 0.20 --reldelta 0.10 --rho 0.90 --min-k 2 \
    --tissue-weighting weighted --pseudocount 0.01
"""

# --- Method notes (short) ---
# Weighted within‑tissue aggregation:
#   Combine cluster nCPM using weights (Read count) so large clusters contribute more than tiny ones.
#   Formula: weighted_mean = Σ(nCPM_i * w_i) / Σ(w_i). Prevents small/noisy clusters skewing the cell type profile.
#
# Weighted bootstrap (example):
#   Sampling clusters with probability ∝ Read count. If clusters have Read counts [1000, 200, 50],
#   then p = [0.8, 0.16, 0.04]. For k=2, most samples include the large cluster + one smaller —
#   resamples reflect the true population.
#
# CV of enrichment (target ≤ 0.20):
#   CV = SD(enrichment) / Mean(enrichment). Values ≤ 0.20 indicate low relative variability → stable.
#   Example: values [10,12,8,10] → mean=10, SD≈1.63 → CV≈0.163 (stable).
#
# Median relative change vs baseline (target ≤ 0.10):
#   median(|E_boot - E_base| / |E_base|). ≤ 10% means resampled enrichment stays close to baseline.
#   Example: base=[10,12,8], boot=[9,11,8] → changes=[0.10,0.083,0] → median=0.083 (stable).
#
# Spearman rank correlation vs baseline (target ≥ 0.90):
#   Checks if gene rank order is preserved. ρ ≥ 0.90 ⇒ very similar ranking.
#   Example: ranks differ slightly → ρ≈0.95 (stable).
# ---

import argparse
import numpy as np
import pandas as pd

# ---------------------------
# Utilities
# ---------------------------

def spearman_corr(x, y):
    """Spearman correlation without SciPy: rank then Pearson."""
    xr = pd.Series(x).rank(method="average").values
    yr = pd.Series(y).rank(method="average").values
    if len(xr) < 2:
        return np.nan
    c = np.corrcoef(xr, yr)
    return c[0, 1]


def effective_clusters(weights):
    """Effective number of clusters given weights (e.g., Read count)."""
    w = np.asarray(weights, dtype=float)
    ss = (w ** 2).sum()
    s = w.sum()
    return (s ** 2 / ss) if ss > 0 else 0.0


# ---------------------------
# Aggregation & Enrichment (HPA-consistent)
# ---------------------------

# Weighted within‑tissue aggregation using Read count as weight

def aggregate_within_tissue(df, expr_col="nCPM", weight_col="Read count"):
    """
    Weighted mean per (Gene, Gene name, Cell type, Tissue) using Read count as weight.
    Also computes clusters_used_tissue and weight_sum_tissue for diagnostics.
    """
    df = df.copy()
    df["_w"] = pd.to_numeric(df[weight_col], errors="coerce").fillna(1.0)
    df["_val"] = pd.to_numeric(df[expr_col], errors="coerce")

    gcols = ["Gene", "Gene name", "Cell type", "Tissue"]

    # Vectorized weighted sum per group (no groupby.apply)
    df["__prod"] = df["_w"] * df["_val"]
    num = df.groupby(gcols)["__prod"].sum()
    den = df.groupby(gcols)["_w"].sum()

    out = (num / den).reset_index(name="avg_nCPM")

    # diagnostics
    out["clusters_used_tissue"] = df.groupby(gcols)["Cluster"].nunique().values
    out["weight_sum_tissue"] = den.values

    # cleanup temp column
    df.drop(columns="__prod", inplace=True)

    return out


def integrate_across_tissues(df_ct, how="weighted"):
    """
    Final cell type profile by averaging tissue-level profiles.
    - 'weighted': weighted mean by weight_sum_tissue
    - 'unweighted': simple mean
    Also returns diagnostics (#tissues, total clusters across tissues, total weight).
    """
    df = df_ct.copy()
    gcols = ["Gene", "Gene name", "Cell type"]
    if how == "weighted" and ("weight_sum_tissue" in df.columns):
        # Vectorized weighted sum across tissues (no groupby.apply)
        df["__prod_ct"] = df["avg_nCPM"] * df["weight_sum_tissue"]
        num = df.groupby(gcols)["__prod_ct"].sum()
        den = df.groupby(gcols)["weight_sum_tissue"].sum()
        out = (num / den).reset_index(name="avg_nCPM")
        df.drop(columns="__prod_ct", inplace=True)
    else:
        out = (df.groupby(gcols, as_index=False)
                 .agg(avg_nCPM=("avg_nCPM", "mean")))

    diag = (df.groupby(gcols, as_index=False)
              .agg(datasets_used=("Tissue", "nunique"),
                   clusters_used=("clusters_used_tissue", "sum"),
                   total_weight=("weight_sum_tissue", "sum")))
    out = out.merge(diag, on=gcols, how="left")
    return out


def add_enrichment(agg_df, gene_col="Gene", value_col="avg_nCPM",
                   out_col="Enrichment score",
                   min_background=1e-3, min_expression=0.0,
                   pseudocount=None):
    """
    Enrichment = value / mean(value in other cell types of the same gene), with safeguards.
    """
    df = agg_df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    sums = df.groupby(gene_col)[value_col].transform("sum")
    counts = df.groupby(gene_col)[value_col].transform("count")
    denom_counts = counts - 1
    avg_other = (sums - df[value_col]) / denom_counts
    avg_other = avg_other.mask(denom_counts <= 0, np.nan)

    if pseudocount is not None:
        avg_other = avg_other + pseudocount
        # Optionally stabilize numerator: df[value_col] = df[value_col] + pseudocount

    denom = np.maximum(avg_other, min_background)
    numer = df[value_col].where(df[value_col] >= min_expression, np.nan)

    df[out_col] = np.divide(numer, denom,
                            out=np.full(df.shape[0], np.nan),
                            where=(denom > 0))
    df.loc[avg_other.isna(), out_col] = np.nan
    return df


# ---------------------------
# Baseline build
# ---------------------------

def build_baseline(cluster_df, expr_col="nCPM", tissue_weighting="weighted", pseudocount=None):
    """
    Build baseline aggregated profiles and enrichment from full dataset.
    Returns:
      within_tissue (DataFrame), across_tissue (DataFrame), baseline_enrich (DataFrame), base_map (dict)
    """
    within_tissue = aggregate_within_tissue(cluster_df, expr_col=expr_col, weight_col="Read count")
    across_tissue = integrate_across_tissues(within_tissue, how=tissue_weighting)
    baseline_enrich = add_enrichment(across_tissue, pseudocount=pseudocount)

    baseline_enrich["key"] = baseline_enrich["Gene"].astype(str) + "||" + baseline_enrich["Cell type"].astype(str)
    base_map = dict(zip(baseline_enrich["key"], baseline_enrich["Enrichment score"]))
    return within_tissue, across_tissue, baseline_enrich, base_map


# ---------------------------
# Leave-One-Tissue-Out (dataset guard)
# ---------------------------

def lodo_stability(cluster_df, base_map, expr_col="nCPM", tissue_weighting="weighted", pseudocount=None):
    """
    Compute per-cell-type max relative change in enrichment when leaving one Tissue out.
    Returns dict: Cell type -> max_relative_change (lower is better).
    """
    df = cluster_df.copy()
    tissues = df["Tissue"].dropna().unique().tolist()
    ct_max_rel = {}

    if not tissues:
        return ct_max_rel

    for t in tissues:
        sub = df[df["Tissue"] != t]
        if sub.empty:
            continue
        within = aggregate_within_tissue(sub, expr_col=expr_col, weight_col="Read count")
        across = integrate_across_tissues(within, how=tissue_weighting)
        enr = add_enrichment(across, pseudocount=pseudocount)
        enr["key"] = enr["Gene"].astype(str) + "||" + enr["Cell type"].astype(str)
        enr["base"] = enr["key"].map(base_map)
        enr = enr.dropna(subset=["Enrichment score", "base"])
        if enr.empty:
            continue

        rel = np.abs(enr["Enrichment score"] - enr["base"]) / (np.abs(enr["base"]) + 1e-9)

        for ct, grp in enr.groupby("Cell type"):
            mx = rel.loc[grp.index].max()
            prev = ct_max_rel.get(ct, 0.0)
            ct_max_rel[ct] = max(prev, mx)

    return ct_max_rel


# ---------------------------
# Weighted bootstrap for min k
# ---------------------------

# Weighted bootstrap stability to recommend minimal clusters k

def recommend_min_k(cluster_df,
                    B=100,
                    cv_thresh=0.20,
                    rel_thresh=0.10,
                    rank_r_thresh=0.90,
                    min_k_rule=2,
                    expr_col="nCPM",
                    tissue_weighting="weighted",
                    pseudocount=None):
    """
    For each (Cell type, Tissue), recommend minimal number of clusters k that yields stable enrichment.
    Sampling is within the (Cell type, Tissue) cluster set. The rest of the dataset stays intact.
    """
    df = cluster_df.copy()
    within_base, across_base, enr_base, base_map = build_baseline(df, expr_col=expr_col,
                                                                  tissue_weighting=tissue_weighting,
                                                                  pseudocount=pseudocount)

    # weights: Read count per (Cell type, Tissue, Cluster)
    wdf = (df.groupby(["Cell type", "Tissue", "Cluster"])['Read count']
             .sum().reset_index())
    weight_map = {(r['Cell type'], r['Tissue'], r['Cluster']): r['Read count'] for _, r in wdf.iterrows()}

    # dataset/tissue counts per cell type
    ct_tissues_used = df.groupby("Cell type")["Tissue"].nunique().to_dict()

    # LODO guard
    lodo = lodo_stability(df, base_map, expr_col=expr_col,
                          tissue_weighting=tissue_weighting,
                          pseudocount=pseudocount)

    results = []

    # Iterate per (Cell type, Tissue)
    for (ct, tissue), clusters in df.groupby(["Cell type", "Tissue"])["Cluster"].unique().items():
        clusters = list(clusters)
        weights = [weight_map.get((ct, tissue, c), 1.0) for c in clusters]
        neff = effective_clusters(weights)
        k_max = len(clusters)
        ds_used = ct_tissues_used.get(ct, np.nan)
        lodo_max = lodo.get(ct, np.nan)

        # Pre-filters
        if k_max < min_k_rule or neff < 2:
            results.append({
                "Cell type": ct,
                "Tissue": tissue,
                "available_clusters": k_max,
                "effective_clusters": round(neff, 3),
                "datasets_used_for_cell_type": ds_used,
                "lodo_max_rel_change": (None if np.isnan(lodo_max) else round(float(lodo_max), 3)),
                "recommended_min_k": np.nan,
                "reason": f"Insufficient replication (k<{min_k_rule} or neff<2)",
                "median_cv": np.nan,
                "median_rel_delta": np.nan,
                "median_rank_corr": np.nan
            })
            continue

        # Dataset guard: prefer ≥2 tissues or LODO max change ≤ 0.15
        dataset_guard = (not np.isnan(ds_used)) and (ds_used >= 2)
        lodo_guard = (not np.isnan(lodo_max)) and (lodo_max <= 0.15)

        # Candidate k
        k_choice = None
        med_cv = med_rel = med_r = np.nan

        # Weighted sampling without replacement within this (ct, tissue)
        w = np.asarray(weights, dtype=float)
        p = (w / w.sum()) if w.sum() > 0 else None

        for k in range(min_k_rule, k_max + 1):
            cvs, rels, rhos = [], [], []

            for b in range(B):
                idxs = np.random.choice(np.arange(k_max), size=k, replace=False, p=p)
                chosen = [clusters[i] for i in idxs]

                # Build a bootstrapped dataset:
                sub_ct_tissue = df[(df["Cell type"] == ct) & (df["Tissue"] == tissue) & (df["Cluster"].isin(chosen))]
                other = df[~((df["Cell type"] == ct) & (df["Tissue"] == tissue))]
                boot_df = pd.concat([other, sub_ct_tissue], ignore_index=True)

                within_b = aggregate_within_tissue(boot_df, expr_col=expr_col, weight_col="Read count")
                across_b = integrate_across_tissues(within_b, how=tissue_weighting)
                enr_b = add_enrichment(across_b, pseudocount=pseudocount)

                # Align with baseline for this cell type only
                e_ct = enr_b[enr_b["Cell type"] == ct].copy()
                e_ct["key"] = e_ct["Gene"].astype(str) + "||" + e_ct["Cell type"].astype(str)
                e_ct["base"] = e_ct["key"].map(base_map)
                e_ct = e_ct.dropna(subset=["Enrichment score", "base"])
                if e_ct.empty:
                    continue

                vals = e_ct["Enrichment score"].values
                mu = np.mean(vals)
                sd = np.std(vals, ddof=1)
                cvs.append(sd / mu if mu > 0 else np.inf)

                rel = np.median(np.abs(vals - e_ct["base"].values) / (np.abs(e_ct["base"].values) + 1e-9))
                rels.append(rel)

                rhos.append(spearman_corr(vals, e_ct["base"].values))

            med_cv = np.median(cvs) if len(cvs) else np.inf
            med_rel = np.median(rels) if len(rels) else np.inf
            med_r = np.median(rhos) if len(rhos) else 0.0

            if (med_cv <= cv_thresh) and (med_rel <= rel_thresh) and (med_r >= rank_r_thresh):
                # If dataset guards fail, you may choose to require one higher k for conservatism
                k_choice = k if (dataset_guard or lodo_guard) else min(k_max, k + 1)
                break

        if k_choice is None:
            results.append({
                "Cell type": ct,
                "Tissue": tissue,
                "available_clusters": k_max,
                "effective_clusters": round(neff, 3),
                "datasets_used_for_cell_type": ds_used,
                "lodo_max_rel_change": (None if np.isnan(lodo_max) else round(float(lodo_max), 3)),
                "recommended_min_k": k_max,
                "reason": "No k met thresholds; using max_k (flag for caution)",
                "median_cv": (None if np.isinf(med_cv) else round(float(med_cv), 3)),
                "median_rel_delta": (None if np.isinf(med_rel) else round(float(med_rel), 3)),
                "median_rank_corr": (None if np.isnan(med_r) else round(float(med_r), 3))
            })
        else:
            results.append({
                "Cell type": ct,
                "Tissue": tissue,
                "available_clusters": k_max,
                "effective_clusters": round(neff, 3),
                "datasets_used_for_cell_type": ds_used,
                "lodo_max_rel_change": (None if np.isnan(lodo_max) else round(float(lodo_max), 3)),
                "recommended_min_k": k_choice,
                "reason": "Meets stability thresholds" + ("" if (dataset_guard or lodo_guard) else " (dataset guard raised k by +1)"),
                "median_cv": round(float(med_cv), 3),
                "median_rel_delta": round(float(med_rel), 3),
                "median_rank_corr": round(float(med_r), 3)
            })

    return pd.DataFrame(results).sort_values(["Cell type", "Tissue"])


# ---------------------------
# Apply adjustments and recompute enrichment
# ---------------------------

# Apply replication filters and recompute final enrichment

def apply_adjustments_and_enrich(cluster_df, recs_df,
                                 expr_col="nCPM",
                                 tissue_weighting="weighted",
                                 pseudocount=None):
    """
    Filter cluster_df to include only (Cell type, Tissue) pairs that meet replication:
      - available_clusters >= recommended_min_k
      - effective_clusters >= 2
      - recommended_min_k is not NaN
    Then recompute within-tissue aggregation, across-tissue integration, and enrichment.
    Returns: filtered_clusters, final_enrichment
    """
    # Merge recommendations on (Cell type, Tissue)
    key_cols = ["Cell type", "Tissue"]
    keep = recs_df.dropna(subset=["recommended_min_k"]).copy()
    keep = keep[(keep["available_clusters"] >= keep["recommended_min_k"]) & (keep["effective_clusters"] >= 2)]

    if keep.empty:
        # No pairs meet replication; return empty results
        return cluster_df.iloc[0:0].copy(), pd.DataFrame(columns=["Gene","Gene name","Cell type","avg_nCPM","datasets_used","clusters_used","total_weight","Enrichment score"]) 

    allowed = keep[key_cols].drop_duplicates()
    filt = cluster_df.merge(allowed, on=key_cols, how="inner")

    # Recompute enrichment on filtered clusters
    within = aggregate_within_tissue(filt, expr_col=expr_col, weight_col="Read count")
    across = integrate_across_tissues(within, how=tissue_weighting)
    final_enr = add_enrichment(across, pseudocount=pseudocount)

    return filt, final_enr


# ---------------------------
# CLI
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Recommend minimum clusters and recompute enrichment using only rna_single_cell_cluster.tsv")
    ap.add_argument("--clusters", required=True, help="Input TSV: rna_single_cell_cluster.tsv")
    ap.add_argument("--out-prefix", default="adjusted", help="Output file prefix")
    ap.add_argument("--B", type=int, default=100, help="Bootstrap replicates")
    ap.add_argument("--cv", type=float, default=0.20, help="CV threshold")
    ap.add_argument("--reldelta", type=float, default=0.10, help="Relative Δ threshold")
    ap.add_argument("--rho", type=float, default=0.90, help="Spearman rank correlation threshold")
    ap.add_argument("--min-k", type=int, default=2, help="Minimum k to consider")
    ap.add_argument("--tissue-weighting", choices=["weighted","unweighted"], default="weighted", help="Across-tissue integration weighting")
    ap.add_argument("--pseudocount", type=float, default=0.01, help="Optional pseudocount for enrichment denominator")

    args = ap.parse_args()

    # Load clusters file
    df = pd.read_csv(args.clusters, sep='\t')
    required = ["Gene", "Gene name", "Tissue", "Cluster", "Cell type", "Read count", "nCPM"]
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise SystemExit(f"Input is missing required columns: {miss}")

    # Recommend min k
    recs = recommend_min_k(df,
                           B=args.B,
                           cv_thresh=args.cv,
                           rel_thresh=args.reldelta,
                           rank_r_thresh=args.rho,
                           min_k_rule=args.min_k,
                           expr_col="nCPM",
                           tissue_weighting=args.tissue_weighting,
                           pseudocount=args.pseudocount)

    # Apply adjustments & recompute enrichment
    filt_clusters, final_enrichment = apply_adjustments_and_enrich(df, recs,
                                                                  expr_col="nCPM",
                                                                  tissue_weighting=args.tissue_weighting,
                                                                  pseudocount=args.pseudocount)

    # Save outputs
    recs_out = f"{args.out_prefix}_min_clusters_per_celltype.tsv"
    filt_out = f"{args.out_prefix}_filtered_clusters.tsv"
    enr_out  = f"{args.out_prefix}_final_enrichment.tsv"

    recs.to_csv(recs_out, sep='\t', index=False)
    filt_clusters.to_csv(filt_out, sep='\t', index=False)
    final_enrichment.to_csv(enr_out, sep='\t', index=False)

    # Console summary
    n_pairs_total = df.groupby(["Cell type","Tissue"]).ngroups
    n_pairs_kept  = filt_clusters.groupby(["Cell type","Tissue"]).ngroups if not filt_clusters.empty else 0
    print(f"Saved: {recs_out}\nSaved: {filt_out}\nSaved: {enr_out}")
    print(f"Pairs total: {n_pairs_total}, kept after adjustments: {n_pairs_kept}")


if __name__ == "__main__":
    main()
```
After the calculation is done for different parts, you should merge this to a single output
```
python split_tsv_robust.py merge   --pattern "enrich_parts/adjusted_*_filtered_clusters.tsv"   --output enrich_parts/merged_filtered_clusters.tsv
```
And then - global enrichment recompute on the merged filtered clusters (to avoid any chunk-wise denominator bias):
```
#!/bin/bash
#SBATCH --job-name=global_enrich
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/global_enrich_%A_%a.out
#SBATCH --error=logs/global_enrich_%A_%a.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load gcc arrow
module load python
python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate

python min_clusters_and_enrichment.py \
  --clusters enrich_parts/merged_filtered_clusters.tsv \
  --out-prefix enrich_parts/global_adjusted \
  --B 100 --cv 0.20 --reldelta 0.10 --rho 0.90 --min-k 2 \
  --tissue-weighting weighted --pseudocount 0.01
```
# Comparing method effectiveness
Finally I use the following R script to compare the effectiveness of different methods with different variables

For this I have 

1) A excel(xlsx) file with information on known genes that are specific in expression to certain cell types that look like following

gene_markers_extended.xlsx
```xlsx
Gene	Cell Type (with extra info)	Cell type (for comparison)	Reference	Specificity Score (human)
MYH7	Cardiomyocytes	cardiomyocytes	Canonical marker, PanglaoDB	
TNNI3	Cardiomyocytes	cardiomyocytes	PanglaoDB, 100% specificity	
TNNT2	Cardiomyocytes	cardiomyocytes	PanglaoDB specificity	
```
2) Following R script saved in the same directory as methods_comparison.txt

```R

###### Match/No Match -> Excel with Conditional Formatting (All Columns Colored)
# + Blank row + Percentage row with color scale for newly added export columns

# ---- Optional: Set working directory to script location if in RStudio ----
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  if (rstudioapi::isAvailable()) {
    try({
      setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    }, silent = TRUE)
  }
}

# ---- Packages ----
# install.packages("openxlsx") # Run once if you don't have it
# install.packages("readxl")   # Run once if you don't have it
library(openxlsx)
library(readxl)

# ---- CONFIG: Change these paths if needed ----
methods_path <- "gene_markers_extended.xlsx" # Excel methods file (.xlsx/.xls)
export_paths <- paste("outputs_to_compare/",
                      list.files("./outputs_to_compare", recursive = TRUE),
                      sep = "")
out_path     <- "methods_comparison.updated.xlsx" # Output Excel file

# ---- General options ----
options(stringsAsFactors = FALSE)

# ---- Helpers ----
normalize_text <- function(x) {
  x2 <- ifelse(is.na(x), "", x)
  x2 <- trimws(x2)
  x2 <- gsub("\\s+", " ", x2) # collapse multiple spaces
  tolower(x2)
}

mk_key <- function(gene, cell) {
  paste(normalize_text(gene), normalize_text(cell), sep = "\n")
}

# Try regex patterns first, then normalized fuzzy matching
find_col <- function(df, patterns) {
  nms <- names(df)
  
  # Regex search
  for (p in patterns) {
    hit <- which(grepl(p, nms, ignore.case = TRUE, perl = TRUE))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  
  # Fuzzy fallback: lowercase, remove punctuation, collapse spaces
  nms_norm <- tolower(gsub("[^a-z0-9]+", " ", nms))
  nms_norm <- trimws(gsub("\\s+", " ", nms_norm))
  lookup <- c("gene", "gene name", "genename",
              "cell type", "celltype",
              "cell type (comparison)", "cell type for comparison")
  for (cand in lookup) {
    idx <- which(nms_norm == cand)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  return(NULL)
}

# Sanitize label (column name) derived from filename
sanitize_label <- function(path) {
  base <- basename(path)
  base_no_ext <- sub("\\.[^.]+$", "", base)
  lab <- gsub("[^A-Za-z0-9]+", "_", base_no_ext) # non-alphanum -> underscore
  lab <- gsub("^_+|_+$", "", lab) # trim boundary underscores
  if (lab == "") lab <- "export"
  lab
}

# Clean carriage returns and whitespace in character columns
clean_text_columns <- function(df) {
  for (i in seq_along(df)) {
    if (is.character(df[[i]])) {
      df[[i]] <- trimws(gsub("\\r", "", df[[i]]))
    }
  }
  df
}

# ---- Read methods file (Excel instead of TSV) ----
methods <- tryCatch(
  readxl::read_excel(methods_path, col_names = TRUE),
  error = function(e) stop("Failed to read methods file: ", e$message)
)

# Remember base (original) columns before appending export result columns
base_cols <- names(methods)

# Identify gene & cell type columns in methods
col_gene_methods <- find_col(methods, c("^gene$", "^gene\\b", "^gene\\s*name$"))
col_cell_methods <- find_col(methods, c("^cell\\s*type\\s*\\(.*comparison.*\\)$",
                                        "^cell\\s*type\\s*\\(for\\s*comparison\\)$",
                                        "^cell\\s*type.*comparison$",
                                        "^cell\\s*type$"))
if (is.null(col_gene_methods) || is.null(col_cell_methods)) {
  stop("Required columns not found in methods file. ",
       "Looked for gene and cell type (comparison). Found columns: ",
       paste(names(methods), collapse = ", "))
}

# ========================= NEW: Drop exact duplicates =========================
# Drop rows where BOTH 'Gene' and 'Cell type (for comparison)' are identical.
dup_mask <- duplicated(mk_key(methods[[col_gene_methods]], methods[[col_cell_methods]]))
n_dup <- sum(dup_mask, na.rm = TRUE)

if (n_dup > 0) {
  # Build a small preview of dropped pairs for the message (up to 5)
  dropped_pairs <- unique(data.frame(
    Gene = methods[[col_gene_methods]][dup_mask],
    `Cell type (for comparison)` = methods[[col_cell_methods]][dup_mask],
    check.names = FALSE
  ))
  preview_n <- min(nrow(dropped_pairs), 5L)
  preview_txt <- paste(
    apply(dropped_pairs[seq_len(preview_n), , drop = FALSE], 1, function(r) {
      sprintf("Gene='%s' | Cell type (for comparison)='%s'", r[1], r[2])
    }),
    collapse = "; "
  )
  
  message(sprintf(
    "Dropped %d duplicate row(s) because 'Gene' and 'Cell type (for comparison)' were duplicated. %s%s",
    n_dup,
    if (preview_n > 0) "Examples: " else "",
    if (preview_n > 0) preview_txt else ""
  ))
  
  # Keep first occurrence; drop the rest
  methods <- methods[!dup_mask, , drop = FALSE]
}
# ==============================================================================

# Build normalized keys for methods rows (on deduplicated data)
methods_keys <- mk_key(methods[[col_gene_methods]], methods[[col_cell_methods]])

# ---- Process each export file ----
for (exp_path in export_paths) {
  exp_df <- tryCatch(
    utils::read.delim(exp_path, check.names = FALSE),
    error = function(e) stop(sprintf("Failed to read export file '%s': %s", exp_path, e$message))
  )
  
  # Identify gene & cell type columns in export
  col_gene_export <- find_col(exp_df, c("^gene\\s*name$", "^gene\\s*name\\b", "^gene$"))
  col_cell_export <- find_col(exp_df, c("^cell\\s*type$", "^cell\\s*type\\b", "^cell$"))
  if (is.null(col_gene_export) || is.null(col_cell_export)) {
    stop(sprintf("Required columns not found in export file: %s. Found columns: %s",
                 exp_path, paste(names(exp_df), collapse = ", ")))
  }
  
  # Build set of keys present in the export
  export_keys <- unique(mk_key(exp_df[[col_gene_export]], exp_df[[col_cell_export]]))
  
  # Compute Match/No Match vector for methods rows
  result_vec <- ifelse(methods_keys %in% export_keys, "Match", "No Match")
  
  # Add/Update a column named after the export file
  label <- sanitize_label(exp_path)
  methods[[label]] <- result_vec
  message(sprintf("Processed %s -> column '%s' added/updated.", exp_path, label))
}

# ---- Clean text columns so conditional rules match reliably ----
methods <- clean_text_columns(methods)

# ---- Write output to Excel with conditional formatting on ALL columns ----
wb <- createWorkbook()
addWorksheet(wb, "Results")
writeData(wb, "Results", methods)

# Styles (Excel standard green/red)
greenStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
redStyle   <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")

n_rows <- nrow(methods)
n_cols <- ncol(methods)

# Apply formatting to rows 2..(n_rows+1) to skip header, across ALL columns
for (col_idx in seq_len(n_cols)) {
  # Color cells that CONTAIN "Match" (case-sensitive)
  conditionalFormatting(
    wb, "Results", cols = col_idx, rows = 2:(n_rows + 1),
    type = "contains", rule = "Match", style = greenStyle
  )
  # Color cells that CONTAIN "No Match"
  conditionalFormatting(
    wb, "Results", cols = col_idx, rows = 2:(n_rows + 1),
    type = "contains", rule = "No Match", style = redStyle
  )
}

# ---- Append empty spacer row and percentage row for newly appended columns ----
# Determine which columns were newly appended (export result columns)
new_cols    <- setdiff(names(methods), base_cols)
new_col_idx <- match(new_cols, names(methods))

if (length(new_col_idx) > 0) {
  # Compute percentage of 'Match' from total of 'Match' + 'No Match' per new column
  matches_counts  <- sapply(new_col_idx, function(ci) sum(methods[[ci]] == "Match",   na.rm = TRUE))
  nomatch_counts  <- sapply(new_col_idx, function(ci) sum(methods[[ci]] == "No Match", na.rm = TRUE))
  totals          <- matches_counts + nomatch_counts
  perc_vals       <- ifelse(totals == 0, NA_real_, matches_counts / totals) # proportions 0..1
  
  # Row indices:
  # - Data rows occupy 2..(n_rows+1) (row 1 is header).
  # - Blank spacer row is (n_rows + 2).
  # - Percentage row is (n_rows + 3).
  blank_row_index <- n_rows + 2
  perc_row_index  <- n_rows + 3
  
  # Write an empty spacer row across all columns using a 1xN matrix (no dimnames)
  spacer_mat <- matrix("", nrow = 1, ncol = n_cols)
  writeData(wb, sheet = "Results", x = spacer_mat,
            startRow = blank_row_index, startCol = 1, colNames = FALSE)
  
  # Build percentage row: values only in new columns; others left as NA (blank in Excel)
  row_vec <- rep(NA_real_, n_cols)
  row_vec[new_col_idx] <- perc_vals
  
  # Write the percentage row as a 1xN matrix (avoids zero-length variable names)
  perc_mat <- matrix(row_vec, nrow = 1, ncol = n_cols)
  writeData(wb, sheet = "Results", x = perc_mat,
            startRow = perc_row_index, startCol = 1, colNames = FALSE)
  
  # Format percentage cells (e.g., 0.00%) only on new columns
  pctStyle <- createStyle(numFmt = "0.00%")
  addStyle(wb, sheet = "Results", style = pctStyle,
           rows = perc_row_index, cols = new_col_idx, gridExpand = TRUE)
  
  # Apply a 3-color scale to the percentage row (lowest red, mid yellow, highest green)
  conditionalFormatting(
    wb, sheet = "Results",
    cols = new_col_idx,
    rows = perc_row_index,
    type = "colorScale",
    style = c("#F8696B", "#FFEB84", "#63BE7B")
  )
}

# Optional: Auto column widths
setColWidths(wb, "Results", cols = 1:n_cols, widths = "auto")

# Save the workbook
saveWorkbook(wb, out_path, overwrite = TRUE)
message("Processing complete!")


```
3) And the output files from previous analysis that looks like following INSIDE A SUBDIRECTORY IN CURRENT WORKING DIRECTORY NAMED "outputs_to_compare"

Example output file from previous analysis. This comparison only checks for the presense of Gene name and Cell type in both data sheets

```
Gene	Gene name	Cell type	avg_nCPM	clusters_used	Enrichment score	single_cell_type_gene
ENSG00000164871	SPAG11B	epididymal principal cells	6432.6	5	87803.51917	FALSE
ENSG00000158874	APOA2	hepatocytes	44653.42857	7	55570.36378	FALSE
ENSG00000228083	IFNA14	pdcs	50.94285714	7	50942.85714	FALSE
ENSG00000213030	CGB8	syncytiotrophoblasts	47.33333333	3	38222.43617	FALSE
ENSG00000286135	ENSG00000286135	epididymal principal cells	245.88	5	35610.49017	FALSE
ENSG00000203970	DEFB110	epididymal principal cells	531.86	5	31273.12245	FALSE
```
# Best method - celltype_enrichment_v_1.4

This is an enhanced version of simple enrichment score method. I have introduced things like batch normalization, Computing Yanai's τ, Applying τ penalization etc.

Before that, I realized some extreme, false expression values in some cell clusters messes up the algorithm. Therefore I had to get rid of (Gene*Cell type) combinations with too much variation between their weighted nCPM values

Following are different ways you can use this filtering script
```

# --- Primary example: row-level robust outlier removal with in-between weighting ---

python filter_weighted_ncpm.py \
  --input rna_single_cell_cluster.tsv \
  --output rna_single_cell_cluster_filtered.tsv \
  --filter-scope row \                 # Drop only offending rows (not whole groups)
  --outlier-method median-mad \        # Robust detection: Median Absolute Deviation (MAD)
  --mad-k 3.0 \                        # MAD threshold (typical: 3.0; lower -> stricter)
  --mad-log \                          # Apply log1p before MAD (stabilizes heavy tails; optional)
  --pair-base alpha \                  # Use in-between weighting: ReadCount^alpha normalized
  --alpha 0.5 \                        # Alpha in [0..1]; 0=raw nCPM, 1=full weighting; 0.5 is a good compromise
  --group-cols Gene "Cell type" \      # Group definition (default: Gene + Cell type)
  --encoding utf-8 \                   # File encoding (default: utf-8)
  --keep-na \                          # Keep rows where an outlier score can't be computed (optional)
  --summary-output rna_single_cell_cluster_summary.tsv \
  --summary-source filtered \          # Write an extra summary TSV from filtered output
  --summary-cols Gene "Cell type"      # Choose


# --- Alternatives & advanced options (enable ONE per line as needed) ---

# 1) Range-based filtering (not robust): use absolute or percent thresholds
#    NOTE: Requires --threshold and --mode; ignores --mad-* flags
# python filter_weighted_ncpm.py \
#   --input rna_single_cell_cluster.tsv \
#   --output rna_single_cell_cluster_filtered_range.tsv \
#   --filter-scope row \               # Or 'group' to drop entire groups
#   --outlier-method range \
#   --threshold 0.10 \                 # If --mode pct: fraction of Row_base group mean (e.g., 0.10 = 10%)
#   --mode pct \                       # 'pct' or 'abs' (abs uses Row_base units)
#   --pair-base wrow \                 # Weighted per-row nCPM (read fraction weighting)
#   --group-cols Gene "Cell type" \
#   --summary-output rna_single_cell_cluster_summary.tsv

# 2) Pair-base choices (choose ONE; affects how per-row values are compared)
# --pair-base row                      # Row_base = nCPM (no read weighting)
# --pair-base wrow                     # Row_base = nCPM * (ReadCount / sum(ReadCount)) (full weighting)
# --pair-base alpha --alpha 0.3        # In-between weighting; tune alpha (e.g., 0.3–0.7)
# --pair-base wterm                    # Row_base = nCPM * ReadCount (product; useful with MAD, often with --mad-log)

# 3) Group-level filtering (drop whole groups that exceed variation threshold)
#    NOTE: Works only with --outlier-method range (+ --threshold/--mode)
# python filter_weighted_ncpm.py \
#   --input rna_single_cell_cluster.tsv \
#   --output rna_single_cell_cluster_filtered_groups.tsv \
#   --filter-scope group \
#   --outlier-method range \
#   --threshold 0.10 \
#   --mode pct \
#   --pair-base wrow \
#   --group-cols Gene "Cell type" \
#   --summary-output rna_single_cell_cluster_summary.tsv \
#   --summary-source full              # Summary based on full dataset (default)

# 4) Encoding alternatives
# --encoding utf-16                    # Use a different file encoding if needed

# 5) Summary file behavior
# --summary-output PATH                # Writes an extra TSV
# --summary-source full                # Summary from full dataset (default)
# --summary-source filtered            # Summary from filtered output
# --summary-cols Gene "Cell type"      # Choose columns that appear in the summary (custom subset)

# 6) Optional summary group selection (include/exclude specific composite keys)
# --summary-group-keys "GENE|CELLTYPE" "GENE2|CELLTYPE2"   # Keys must match the order in --group-cols
# --summary-group-file summary_groups.txt                   # One key per line (supports comments with '#')
# --summary-group-mode include                              # 'include' keeps only those groups; use 'exclude' to remove them

```
Quick reference for flags

--filter-scope group|row: Drop entire groups or only rows that are outliers.

--outlier-method range|median-mad:

range: simple max-difference; requires --threshold and --mode.

median-mad: robust; uses --mad-k (and optional --mad-log). Row only.


--mode abs|pct (range method only): Use absolute units or fraction of group mean.

--pair-base row|wrow|alpha|wterm: How per-row values are formed for comparison:

row: nCPM
wrow: nCPM × (Read / sum(Read))
alpha: nCPM × (Read^α / sum(Read^α)) → recommended compromise
wterm: nCPM × Read (use with MAD; consider --mad-log)


--alpha: 0..1 (only for --pair-base alpha), set to 0.5 as a good starting point.

--mad-k: MAD threshold; typical robust choice 3.0 (lower for stricter).

--mad-log: Apply log1p transform before MAD (helps when values vary widely).

--summary-output: Writes a compact summary TSV.

--summary-source full|filtered: Choose source for summary (default full).


I am using following 
```
python filter_weighted_ncpm.py \
  --input rna_single_cell_cluster.tsv \
  --output rna_single_cell_cluster_filtered_rows_alpha_mad.tsv \
  --filter-scope row \
  --pair-base alpha \
  --alpha 0.5 \
  --outlier-method median-mad \
  --mad-k 3.0 \
  --group-cols Gene "Cell type" \
  --summary-output rna_single_cell_cluster_summary.tsv \
  --summary-source filtered \
  --summary-cols "Gene" "Gene name" "Tissue" "Cell type" "Read count" "nCPM" "Row_mad_score" \
```
This is the script used

filter_weighted_ncpm.py
```py

#!/usr/bin/env python3
"""
Detect abnormal variation between rows within each group (default: Gene × Cell type),
using a configurable per-row base value, then filter groups or rows based on either
range or robust MAD outlier rules. Includes a summary file with optional column and group selection.

Per-row base options (choose with --pair-base):
  - 'row'   : Row_base = nCPM                                (no read weighting)
  - 'wrow'  : Row_base = nCPM * (ReadCount / sum(ReadCount)) (full read weighting)
  - 'alpha' : Row_base = nCPM * (ReadCount^alpha / sum(ReadCount^alpha))  [in-between]
  - 'wterm' : Row_base = nCPM * ReadCount                    (product for median checks)

Range-based variation (for --outlier-method range):
  - Group_variation_abs = max(Row_base) - min(Row_base)
  - Group_variation_pct = Group_variation_abs / mean(Row_base)   (unit-consistent when --mode pct)
  - Row_max_diff_abs = max(|Row_base - min(Row_base)|, |max(Row_base) - Row_base|)  [row scope]
  - Row_max_diff_pct = Row_max_diff_abs / mean(Row_base)         [row scope, --mode pct]

Robust MAD-based outliers (for --outlier-method median-mad, row scope only):
  - Optionally apply log1p scaling to Row_base before MAD with --mad-log
  - Row_mad_score = |Row_base - median(Row_base)| / MAD(Row_base)
  - Drop rows with Row_mad_score > --mad-k

Summary:
  - --summary-cols lets you pick the **columns** in the summary TSV (e.g., Gene "Cell type")
  - --summary-group-keys / --summary-group-file can include/exclude specific groups (optional)

Progress messages are printed at each major step.
"""

import argparse
import sys
import numpy as np
import pandas as pd
from typing import List, Optional, Set

# ---------------------------- Progress logger ----------------------------

def log(msg: str):
    """Print a progress message and flush immediately."""
    print(msg, flush=True)

# ---------------------------- Validation ----------------------------

REQUIRED_COLUMNS = [
    "Gene", "Gene name", "Tissue", "Cluster", "Cell type", "Read count", "nCPM",
]

def validate_columns(df: pd.DataFrame):
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

# ---------------------------- Summary group selection ----------------------------

def load_summary_keys(
    keys_cli: Optional[List[str]],
    keys_file: Optional[str],
) -> Optional[Set[str]]:
    """
    Build a set of composite group keys from CLI list and/or a file.
    Each key is expected to be a pipe-delimited string: 'col1|col2|...'
    Returns None if no keys were provided.
    """
    keys: Set[str] = set()
    if keys_cli:
        for k in keys_cli:
            if k and k.strip():
                keys.add(k.strip())
    if keys_file:
        try:
            with open(keys_file, "r", encoding="utf-8") as f:
                for line in f:
                    s = line.strip()
                    if not s or s.startswith("#"):
                        continue
                    keys.add(s)
        except Exception as e:
            print(f"ERROR: failed to read summary-group-file '{keys_file}': {e}", file=sys.stderr)
            sys.exit(1)
    return keys if keys else None

def build_group_key(df: pd.DataFrame, group_cols: List[str]) -> pd.Series:
    """Create a pipe-delimited composite key series with the columns in order."""
    return df[group_cols].astype(str).agg("|".join, axis=1)

# ---------------------------- Metrics computation ----------------------------

def compute_metrics(
    df: pd.DataFrame,
    group_cols: List[str],
    weight_col: str,
    value_col: str,
    pair_base: str,       # 'row' | 'wrow' | 'alpha' | 'wterm'
    alpha: float,         # used only when pair_base == 'alpha'
    mad_log: bool,        # apply log1p before MAD?
    # Conditional computation flags:
    compute_range: bool,
    compute_range_pct: bool,
    compute_row_range: bool,  # row-level range diffs (for --filter-scope row & range)
    compute_mad: bool,
):
    """
    Compute only the metrics needed for this run based on flags.

    Always computes:
      - Row_base

    Conditionally computes:
      - Group_variation_abs / Group_variation_pct (range)
      - Row_max_diff_abs / Row_max_diff_pct (row-level range)
      - Row_mad_score (median-MAD)
    """
    log("➡️ Step: Ensuring numeric types for Read count and nCPM...")
    df[weight_col] = pd.to_numeric(df[weight_col], errors="coerce")
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    log(f"➡️ Step: Grouping rows by: {group_cols}")
    grp = df.groupby(group_cols, dropna=False)

    # -------- Row_base (safe normalization) --------
    log(f"➡️ Step: Building Row_base using pair-base = '{pair_base}'" + (f" (alpha={alpha})" if pair_base == "alpha" else ""))
    if pair_base == "row":
        df["Row_base"] = df[value_col]

    elif pair_base == "wrow":
        sum_weights = grp[weight_col].transform("sum")
        w_norm = pd.Series(0.0, index=df.index)
        nonzero_sum = sum_weights != 0
        w_norm.loc[nonzero_sum] = (
            df.loc[nonzero_sum, weight_col] / sum_weights.loc[nonzero_sum]
        )
        df["Row_base"] = df[value_col] * w_norm

    elif pair_base == "alpha":
        df["_w_alpha"] = df[weight_col] ** alpha
        sum_w_alpha = grp["_w_alpha"].transform("sum")
        w_norm_alpha = pd.Series(0.0, index=df.index)
        nonzero_alpha = sum_w_alpha != 0
        w_norm_alpha.loc[nonzero_alpha] = (
            df.loc[nonzero_alpha, "_w_alpha"] / sum_w_alpha.loc[nonzero_alpha]
        )
        df["Row_base"] = df[value_col] * w_norm_alpha

    elif pair_base == "wterm":
        df["Row_base"] = df[value_col] * df[weight_col]

    else:
        raise ValueError(f"Unsupported pair_base: {pair_base}")

    # Compute group statistics only when needed
    row_base_mean = None
    min_base = None
    max_base = None

    if compute_range or compute_row_range:
        log("➡️ Step: Computing min/max of Row_base within groups (range)...")
        min_base = grp["Row_base"].transform("min")
        max_base = grp["Row_base"].transform("max")

    if compute_range_pct or compute_row_range:  # need group mean for pct in either case
        log("➡️ Step: Computing group mean of Row_base (for % metrics)...")
        row_base_mean = grp["Row_base"].transform("mean")
        df["Row_base_group_mean"] = row_base_mean

    # -------- Range-based metrics (group scope) --------
    if compute_range:
        log("➡️ Step: Computing group-level range metrics...")
        df["Group_variation_abs"] = max_base - min_base

        if compute_range_pct:
            df["Group_variation_pct"] = pd.NA
            valid_mean = row_base_mean > 0
            df.loc[valid_mean, "Group_variation_pct"] = (
                df.loc[valid_mean, "Group_variation_abs"] / df.loc[valid_mean, "Row_base_group_mean"]
            )

    # -------- Range-based metrics (row scope) --------
    if compute_row_range:
        log("➡️ Step: Computing row-level max-diff metrics (range)...")
        diff_to_min = (df["Row_base"] - min_base).abs()
        diff_to_max = (max_base - df["Row_base"]).abs()
        df["Row_max_diff_abs"] = pd.DataFrame({"a": diff_to_min, "b": diff_to_max}).max(axis=1)

        if compute_range_pct:
            df["Row_max_diff_pct"] = pd.NA
            valid_mean = row_base_mean > 0
            df.loc[valid_mean, "Row_max_diff_pct"] = (
                df.loc[valid_mean, "Row_max_diff_abs"] / df.loc[valid_mean, "Row_base_group_mean"]
            )

    # -------- Robust MAD-based outlier scores --------
    if compute_mad:
        log("➡️ Step: Computing robust MAD outlier scores..." + (" (log1p scaling applied)" if mad_log else ""))
        base_for_mad = df["Row_base"].copy()
        if mad_log:
            base_for_mad = np.log1p(base_for_mad)

        # Build group key
        group_key = df[group_cols].astype(str).agg("|".join, axis=1)

        # Group-wise median
        median_base = base_for_mad.groupby(group_key).transform("median")

        # Group-wise MAD that returns NaN for empty or single-valued groups
        def _mad_series(s: pd.Series) -> float:
            if len(s) == 0:
                return np.nan
            med = s.median()
            abs_dev = (s - med).abs()
            mad = abs_dev.median()
            return np.nan if mad == 0 else mad

        mad_series = base_for_mad.groupby(group_key).transform(_mad_series)

        # Avoid division by zero MAD: score computed only where MAD is non-NaN
        df["Row_mad_score"] = pd.NA
        valid_mad = mad_series.notna()
        df.loc[valid_mad, "Row_mad_score"] = (
            (base_for_mad.loc[valid_mad] - median_base.loc[valid_mad]).abs() / mad_series.loc[valid_mad]
        )

    # Cleanup helper cols
    log("➡️ Step: Cleaning helper columns...")
    drop_cols = [c for c in ["_w_alpha"] if c in df.columns]
    if drop_cols:
        df.drop(columns=drop_cols, inplace=True)

    log("✅ Metrics computation complete.")
    return df

# ---------------------------- CLI parsing ----------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Detect abnormal variation between rows within groups and filter output."
    )
    p.add_argument("--input", "-i", required=True, help="Path to input TSV.")
    p.add_argument("--output", "-o", required=True, help="Path to output TSV.")
    p.add_argument(
        "--threshold", "-t", required=False, type=float,
        help="Variation threshold for 'range' method. If --mode abs, in base units; if --mode pct, as fraction (e.g., 0.10)."
    )
    p.add_argument(
        "--mode", "-m", choices=["abs", "pct"], default="abs",
        help="Comparison mode for 'range' method: 'abs' (absolute) or 'pct' (fraction of Row_base group mean)."
    )
    p.add_argument(
        "--filter-scope", choices=["group", "row"], default="group",
        help="Drop entire groups ('group') or only offending rows ('row')."
    )
    p.add_argument(
        "--pair-base", choices=["row", "wrow", "alpha", "wterm"], default="alpha",
        help="Per-row base: 'row' (nCPM), 'wrow' (weighted by read fraction), "
             "'alpha' (read^alpha normalized), or 'wterm' (nCPM*ReadCount)."
    )
    p.add_argument(
        "--alpha", type=float, default=0.5,
        help="Exponent for --pair-base alpha (0..1). 0=no weighting; 1=full weighting. Default 0.5."
    )
    p.add_argument(
        "--outlier-method", choices=["range", "median-mad"], default="range",
        help="Outlier detection method. 'median-mad' is robust and applies to row-level filtering."
    )
    p.add_argument(
        "--mad-k", type=float, default=3.0,
        help="Threshold k for MAD-based outlier detection (default: 3.0). Rows with Row_mad_score > k are dropped."
    )
    p.add_argument(
        "--mad-log", action="store_true",
        help="Apply log1p to Row_base before MAD to stabilize heavy tails."
    )
    p.add_argument(
        "--group-cols", nargs="+", default=["Gene", "Cell type"],
        help="Columns to define groups (default: Gene and Cell type)."
    )
    p.add_argument("--encoding", default="utf-8", help="File encoding (default: utf-8).")
    p.add_argument(
        "--keep-na", dest="keep_na", action="store_true",
        help="Keep rows/groups where variation cannot be evaluated (e.g., NA in pct or MAD modes)."
    )
    p.add_argument(
        "--summary-output", "-s", default=None,
        help="(Optional) Path to write a summary TSV."
    )
    p.add_argument(
        "--summary-source", choices=["full", "filtered"], default="full",
        help="Source for summary TSV: 'full' dataset (default) or 'filtered' output."
    )
    # Column selection for summary
    p.add_argument(
        "--summary-cols", nargs="+", default=None,
        help="Columns to include in the summary TSV (e.g., Gene \"Cell type\" Tissue). "
             "If omitted, a default set is used."
    )
    # Optional group selection for summary (still available)
    p.add_argument(
        "--summary-group-keys", nargs="+", default=None,
        help="Composite group keys to include/exclude in summary (each key = col1|col2|... per --group-cols order)."
    )
    p.add_argument(
        "--summary-group-file", default=None,
        help="Path to a text file with one composite group key per line (col1|col2|...)."
    )
    p.add_argument(
        "--summary-group-mode", choices=["include", "exclude"], default="include",
        help="If keys are provided, 'include' keeps only those groups, 'exclude' removes those groups."
    )
    return p.parse_args()

# ---------------------------- Main ----------------------------

def main():
    log("🔧 Parsing command-line arguments...")
    args = parse_args()

    # Method-specific validations
    log("🔎 Validating method-specific options...")
    if args.outlier_method == "median-mad" and args.filter_scope != "row":
        print("ERROR: --outlier-method median-mad requires --filter-scope row.", file=sys.stderr)
        sys.exit(1)
    if args.outlier_method == "range" and args.threshold is None:
        print("ERROR: --threshold is required when --outlier-method range.", file=sys.stderr)
        sys.exit(1)
    if args.pair_base == "alpha" and not (0.0 <= args.alpha <= 1.0):
        print("ERROR: --alpha must be between 0 and 1.", file=sys.stderr)
        sys.exit(1)

    # Load
    log(f"📥 Loading input TSV: {args.input}")
    try:
        df = pd.read_csv(args.input, sep="\t", encoding=args.encoding)
    except Exception as e:
        print(f"ERROR: failed to read '{args.input}': {e}", file=sys.stderr)
        sys.exit(1)
    log(f"📊 Loaded {len(df)} rows, {len(df.columns)} columns.")

    # Validate
    log("✅ Validating required columns...")
    try:
        validate_columns(df)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    log("✅ Required columns present.")

    # Decide computations needed
    compute_range = (args.outlier_method == "range")
    compute_range_pct = compute_range and (args.mode == "pct")
    compute_row_range = compute_range and (args.filter_scope == "row")
    compute_mad = (args.outlier_method == "median-mad")

    # Compute metrics (only what's needed)
    log("🧮 Computing metrics (only those needed for the selected method)...")
    df = compute_metrics(
        df=df,
        group_cols=list(args.group_cols),
        weight_col="Read count",
        value_col="nCPM",
        pair_base=args.pair_base,
        alpha=args.alpha,
        mad_log=args.mad_log,
        compute_range=compute_range,
        compute_range_pct=compute_range_pct,
        compute_row_range=compute_row_range,
        compute_mad=compute_mad,
    )

    # Build keep mask
    log(f"⚙️ Building keep mask | Scope={args.filter_scope} | Method={args.outlier_method}...")
    if args.filter_scope == "group":
        # Group-level filtering uses range-based variation only
        comp = df["Group_variation_abs"] if args.mode == "abs" else df["Group_variation_pct"]
        evaluable = comp.notna()
        keep_mask = evaluable & (comp <= args.threshold)
        if args.keep_na:
            keep_mask = keep_mask | (~evaluable)
        evaluable_count = int(evaluable.sum())
        log(f"🔢 Evaluable groups (broadcast to rows): {evaluable_count} rows evaluable.")

    else:  # row-level filtering
        if args.outlier_method == "median-mad":
            comp = df["Row_mad_score"]
            evaluable = comp.notna()
            keep_mask = evaluable & (comp <= args.mad_k)
            if args.keep_na:
                keep_mask = keep_mask | (~evaluable)
            log(f"🔢 Evaluable rows (MAD): {int(evaluable.sum())} rows evaluable.")
        else:  # range
            comp = df["Row_max_diff_abs"] if args.mode == "abs" else df["Row_max_diff_pct"]
            evaluable = comp.notna()
            keep_mask = evaluable & (comp <= args.threshold)
            if args.keep_na:
                keep_mask = keep_mask | (~evaluable)
            log(f"🔢 Evaluable rows (range): {int(evaluable.sum())} rows evaluable.")

    out = df.loc[keep_mask].copy()
    log(f"📉 Filtering applied: kept {len(out)} / {len(df)} rows.")

    # Save filtered output
    log(f"💾 Writing filtered output to: {args.output}")
    try:
        out.to_csv(args.output, sep="\t", index=False, encoding=args.encoding)
    except Exception as e:
        print(f"ERROR: failed to write '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)
    log("✅ Filtered output written.")

    # Optional separate summary TSV (selected columns + optional group filtering)
    if args.summary_output:
        log(f"🧾 Preparing summary TSV: {args.summary_output} (source={args.summary_source})")
        source_df = df if args.summary_source == "full" else out

        # Determine summary columns
        if args.summary_cols:
            summary_cols = list(args.summary_cols)  # user-specified
            log(f"🧩 Using user-selected summary columns: {summary_cols}")
        else:
            summary_cols = [
                "Gene", "Gene name", "Tissue", "Cluster", "Cell type",
                "Read count", "nCPM", "Group_variation_pct",
            ]
            log(f"🧩 Using default summary columns: {summary_cols}")

        # Ensure requested columns exist; if Group_variation_pct requested but not present, create NA
        if "Group_variation_pct" in summary_cols and "Group_variation_pct" not in source_df.columns:
            log("ℹ️ 'Group_variation_pct' requested but not computed; adding as NA to summary source.")
            source_df = source_df.copy()
            source_df["Group_variation_pct"] = pd.NA

        missing = [c for c in summary_cols if c not in source_df.columns]
        if missing:
            print(f"ERROR: summary source missing columns: {missing}", file=sys.stderr)
            sys.exit(1)

        # Optional group selection for summary
        summary_keys = load_summary_keys(args.summary_group_keys, args.summary_group_file)
        if summary_keys is not None:
            log(f"🔎 Applying summary group {args.summary_group_mode}: {len(summary_keys)} key(s)")
            group_key_series = build_group_key(source_df, list(args.group_cols))
            if args.summary_group_mode == "include":
                mask_summary = group_key_series.isin(summary_keys)
            else:  # exclude
                mask_summary = ~group_key_series.isin(summary_keys)
            source_df = source_df.loc[mask_summary]

        # Write summary
        try:
            source_df.loc[:, summary_cols].to_csv(
                args.summary_output, sep="\t", index=False, encoding=args.encoding
            )
            log(f"✅ Summary file written: {args.summary_output}")
        except Exception as e:
            print(f"ERROR: failed to write summary file '{args.summary_output}': {e}", file=sys.stderr)
            sys.exit(1)

    # Console summary
    total = len(df)
    kept = len(out)
    dropped = total - kept
    log("📣 Run complete.")
    log(f"• Grouping columns: {args.group_cols}")
    log(f"• Pair base: {args.pair_base} (alpha={args.alpha if args.pair_base=='alpha' else 'n/a'})")
    log(f"• Scope: {args.filter_scope}")
    if args.outlier_method == "range":
        log(f"• Outlier method: range | Mode: {args.mode} | Threshold: {args.threshold}")
    else:
        log(f"• Outlier method: median-mad | MAD-k: {args.mad_k} | MAD-log: {args.mad_log}")
    log(f"• Rows kept: {kept}/{total} (dropped: {dropped})")
    log(f"• Output: {args.output}")
    if args.summary_output:
        log(f"• Summary: {args.summary_output} (source: {args.summary_source})")

if __name__ == "__main__":
    main()
```
Following is the enrichment script

celltype_enrichment_v1_4.py
```py

#!/usr/bin/env python3
"""
Cell-type–aware aggregation and enrichment with optional weighting
and Yanai's τ (specificity) — v1.4

Includes:
- τ report column (specificity_tau)
- Modes: --specificity-mode off|filter|penalize
- Safe log2 with epsilon + zero masking
- Series-based τ computation to avoid FutureWarning
"""

import argparse
import numpy as np
import pandas as pd
from typing import Tuple, List, Optional

# ---------- Utilities ----------

def coerce_numeric(series: pd.Series) -> Tuple[pd.Series, List]:
    s = pd.to_numeric(series, errors="coerce")
    non_numeric = series.loc[s.isna()].unique().tolist()
    return s, non_numeric


def drop_genes_with_no_expression(
    agg_df: pd.DataFrame,
    expr_col: Optional[str] = None,
    treat_nan_as_zero: bool = False,
) -> Tuple[pd.DataFrame, List[str], int]:
    if expr_col is None:
        expr_col = "avg_nCPM" if "avg_nCPM" in agg_df.columns else "nCPM"
    df = agg_df.copy()
    df[expr_col] = pd.to_numeric(df[expr_col], errors="coerce")
    df[expr_col] = df[expr_col].replace([np.inf, -np.inf], np.nan)
    expr_for_test = df[expr_col].fillna(0) if treat_nan_as_zero else df[expr_col]
    gene_max = expr_for_test.groupby(df["Gene"]).max()
    genes_all_zero = gene_max[gene_max == 0].index.tolist()
    before = len(df)
    filtered_df = df[~df["Gene"].isin(genes_all_zero)].copy()
    after = len(filtered_df)
    rows_removed = before - after
    return filtered_df, genes_all_zero, rows_removed


def add_enrichment(
    agg_df: pd.DataFrame,
    gene_col: str = "Gene",
    value_col: str = "avg_nCPM",
    out_col: str = "Enrichment score",
    min_background: float = 1e-3,
    min_expression: float = 0.0,
    min_count: int = 2,
    pseudocount: Optional[float] = None,
    pseudocount_to_numerator: bool = False,
    clip_max: Optional[float] = None,
) -> pd.DataFrame:
    """Compute enrichment per row = current / mean(other cell types of the same gene)."""
    df = agg_df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    gene_sums = df.groupby(gene_col)[value_col].transform("sum")
    gene_counts = df.groupby(gene_col)[value_col].transform("count")

    denom_counts = gene_counts - 1
    avg_other = (gene_sums - df[value_col]) / denom_counts
    avg_other = avg_other.mask(denom_counts <= 0, np.nan)

    if pseudocount is not None:
        avg_other = avg_other + pseudocount

    numer = df[value_col]
    if pseudocount is not None and pseudocount_to_numerator:
        numer = numer + pseudocount

    denom = np.maximum(avg_other, min_background)
    numer = numer.where(numer >= min_expression, np.nan)

    df[out_col] = np.divide(
        numer, denom,
        out=np.full(df.shape[0], np.nan, dtype=float),
        where=(denom > 0)
    )
    df.loc[avg_other.isna(), out_col] = np.nan

    df[out_col] = df[out_col].where(gene_counts >= min_count, np.nan)

    if clip_max is not None:
        df[out_col] = df[out_col].clip(upper=clip_max)

    # Safe log2: add epsilon + mask zeros to NaN
    EPS = 1e-12
    log2_vals = np.log2(df[out_col].clip(lower=EPS))
    log2_vals = pd.Series(log2_vals, index=df.index).mask(df[out_col] <= 0)
    df["log2_enrichment"] = log2_vals
    return df


def batch_normalize_if_needed(df: pd.DataFrame, value_col: str, batch_col: Optional[str], batch_normalize: str = "none") -> pd.DataFrame:
    """Median-scale the value column per batch to align batch medians with the global median."""
    if not batch_col or batch_normalize == "none" or batch_col not in df.columns:
        return df

    out = df.copy()
    # Global median across all rows
    out[value_col] = pd.to_numeric(out[value_col], errors="coerce")
    global_median = out[value_col].median()
    if pd.isna(global_median) or global_median == 0:
        return out

    # Per-batch medians
    batch_medians = (
        out.groupby(batch_col, dropna=False)[value_col]
           .median()
           .rename("_batch_median")
    )

    # Map per-row scale = global_median / batch_median
    out = out.merge(batch_medians, on=batch_col, how="left")
    scale = np.where((out["_batch_median"].notna()) & (out["_batch_median"] != 0),
                     global_median / out["_batch_median"], 1.0)
    out[value_col] = out[value_col] * scale
    out.drop(columns=["_batch_median"], inplace=True)
    return out
def aggregate_with_celltype(
    df: pd.DataFrame,
    gene_col: str,
    gene_name_col: str,
    cell_type_col: str,
    value_col: str,
    cluster_col: Optional[str],
    weighted: bool,
    weight_col: Optional[str] = "Read count",
    cluster_aggregate: str = "mean",  # used when weighted is off: mean|median
) -> pd.DataFrame:
    """
    Aggregate to one row per (Gene × Gene name × Cell type) across clusters.
    If weighted=True and weight_col exists: weighted mean Σ(nCPM*w)/Σ(w).
    Else: unweighted mean or median across clusters.
    """
    df = df.copy()
    group_cols = [gene_col, gene_name_col, cell_type_col]

    if weighted and weight_col and (weight_col in df.columns):
        w = pd.to_numeric(df[weight_col], errors="coerce").fillna(1.0)
        vals = pd.to_numeric(df[value_col], errors="coerce")
        df["__prod__"] = vals * w
        agg = (
            df.groupby(group_cols, as_index=False)
              .agg(
                  avg_nCPM=("__prod__", "sum"),
                  weight_sum=(weight_col, "sum"),
                  clusters_used=(cluster_col, "nunique") if (cluster_col and cluster_col in df.columns) else (value_col, "count")
              )
        )
        agg["avg_nCPM"] = agg["avg_nCPM"] / agg["weight_sum"].replace(0, np.nan)
        df.drop(columns=["__prod__"], inplace=True)
    else:
        # Unweighted aggregation across clusters
        func = "mean" if cluster_aggregate == "mean" else "median"
        agg = (
            df.groupby(group_cols, as_index=False)
              .agg(
                  avg_nCPM=(value_col, func),
                  clusters_used=(cluster_col, "nunique") if (cluster_col and cluster_col in df.columns) else (value_col, "count")
              )
        )
    return agg

# ---------- Yanai's τ (specificity) ----------

def gene_specificity_tau(sub: pd.DataFrame, expr_col: str = "avg_nCPM") -> float:
    """Compute Yanai's τ for one gene from its avg_nCPM across cell types (0..1)."""
    vals = pd.to_numeric(sub[expr_col], errors="coerce").fillna(0.0).to_numpy()
    if len(vals) == 0:
        return np.nan
    m = vals.max()
    if m <= 0:
        return 0.0
    y = vals / m  # normalize by max
    K = len(vals)
    tau = (np.sum(1.0 - y)) / (K - 1) if K > 1 else 1.0
    return float(tau)


def gene_specificity_tau_series(s: pd.Series) -> float:
    """Yanai's τ from a Series of avg_nCPM values across cell types (0..1)."""
    vals = pd.to_numeric(s, errors="coerce").fillna(0.0).to_numpy()
    if len(vals) == 0:
        return np.nan
    m = vals.max()
    if m <= 0:
        return 0.0
    y = vals / m
    K = len(vals)
    tau = (np.sum(1.0 - y)) / (K - 1) if K > 1 else 1.0
    return float(tau)


def compute_tau(agg_df: pd.DataFrame, gene_col: str = "Gene", expr_col: str = "avg_nCPM") -> pd.DataFrame:
    """Compute τ per gene and merge as a report column specificity_tau."""
    df = agg_df.copy()
    tau_series = df.groupby(gene_col, group_keys=False)[expr_col].apply(gene_specificity_tau_series)
    tau_df = tau_series.rename("specificity_tau").reset_index()
    df = df.merge(tau_df, on=gene_col, how="left")
    return df


def apply_tau_filter(
    agg_df: pd.DataFrame,
    gene_col: str = "Gene",
    min_specificity: Optional[float] = None,
) -> Tuple[pd.DataFrame, int, int]:
    """
    Drop genes with τ < min_specificity. Assumes specificity_tau is present.
    Returns (filtered_df, n_total_genes, n_dropped_genes).
    """
    if min_specificity is None:
        total = agg_df[gene_col].nunique()
        return agg_df, total, 0
    df = agg_df.copy()
    tau_per_gene = df[[gene_col, "specificity_tau"]].drop_duplicates(subset=[gene_col])
    n_total = tau_per_gene[gene_col].nunique()
    keep_genes = tau_per_gene.loc[tau_per_gene["specificity_tau"] >= min_specificity, gene_col]
    filtered_df = df[df[gene_col].isin(keep_genes)].copy()
    n_dropped = n_total - keep_genes.nunique()
    return filtered_df, n_total, n_dropped

# ---------- Top N helpers ----------

def top_n_overall(df: pd.DataFrame, sort_by: str, n: int) -> pd.DataFrame:
    return df.sort_values(by=sort_by, ascending=False).head(n).copy()


def top_n_per_cell_type(df: pd.DataFrame, cell_type_col: str, sort_by: str, n: int) -> pd.DataFrame:
    out_frames = []
    for ct, sub in df.groupby(cell_type_col):
        sub_sorted = sub.sort_values(by=sort_by, ascending=False).head(n).copy()
        sub_sorted["rank_in_cell_type"] = range(1, len(sub_sorted) + 1)
        out_frames.append(sub_sorted)
    if out_frames:
        return pd.concat(out_frames, axis=0, ignore_index=True)
    else:
        return pd.DataFrame(columns=list(df.columns) + ["rank_in_cell_type"])

# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description="Cell-type–aware enrichment with optional weighting and τ report/filter/penalize — v1.4."
    )

    # I/O
    parser.add_argument("--input-file", type=str, default="rna_single_cell_cluster.tsv",
                        help="Path to input TSV containing single-cell cluster data.")
    parser.add_argument("--output-file", type=str, default="celltype_enrichment.tsv",
                        help="Path for the full enrichment output TSV.")
    parser.add_argument("--top-n", type=int, default=100,
                        help="Top N rows to save overall (default: 100).")
    parser.add_argument("--per-cell-type-top-n", type=int, default=20,
                        help="Top N per cell type to export (0 disables).")

    # Column names (defaults match your example)
    parser.add_argument("--gene-col", type=str, default="Gene")
    parser.add_argument("--gene-name-col", type=str, default="Gene name")
    parser.add_argument("--cell-type-col", type=str, default="Cell type")
    parser.add_argument("--cluster-col", type=str, default="Cluster")
    parser.add_argument("--batch-col", type=str, default="Cell type",
                        help="Batch column used for median scaling (default: Cell type).")
    parser.add_argument("--value-col", type=str, default="nCPM")
    parser.add_argument("--weight-col", type=str, default="Read count",
                        help="Weight column used for weighted aggregation (default: Read count).")

    # Weighting options
    parser.add_argument("--weighted", type=str, choices=["on", "off"], default="on",
                        help="Use weighted aggregation across clusters (on/off).")
    parser.add_argument("--cluster-aggregate", type=str, choices=["mean", "median"], default="mean",
                        help="When --weighted off, choose mean or median across clusters.")

    # Filters and safeguards
    parser.add_argument("--min-clusters", type=int, default=None,
                        help="Keep only rows where clusters_used >= this integer.")
    parser.add_argument("--treat-nan-as-zero", action="store_true",
                        help="Treat NaN as zero when deciding 'no expression' genes to drop.")
    parser.add_argument("--min-expr-threshold", type=float, default=0.0,
                        help="Filter rows with avg_nCPM < threshold before enrichment (default: 0).")

    # Enrichment parameters
    parser.add_argument("--min-background", type=float, default=1e-3,
                        help="Minimum denominator for enrichment (default: 1e-3).")
    parser.add_argument("--min-expression", type=float, default=0.0,
                        help="Minimum numerator expression to compute enrichment (default: 0).")
    parser.add_argument("--min-count", type=int, default=2,
                        help="Require >= this many entries per gene for enrichment.")
    parser.add_argument("--pseudocount", type=float, default=None,
                        help="Optional pseudocount added to denominator; set a small value like 0.01.")
    parser.add_argument("--pseudocount-to-numerator", action="store_true",
                        help="Also add pseudocount to numerator.")
    parser.add_argument("--clip-max", type=float, default=None,
                        help="Optional cap on enrichment score (e.g., 100).")
    parser.add_argument("--sort-by", type=str, choices=[
        "Enrichment score",
        "log2_enrichment",
        "Enrichment score (tau penalized)",
        "log2_enrichment_penalized"
    ], default="log2_enrichment",
        help="Column used for sorting outputs.")

    # Batch normalization
    parser.add_argument("--batch-normalize", type=str, choices=["none", "median_scale"],
                        default="median_scale", help="Per-batch normalization method.")

    # τ report/filter/penalize
    parser.add_argument("--specificity-mode", type=str, choices=["off", "filter", "penalize"], default="off",
                        help="How to use Yanai's τ: off, filter (drop genes below threshold), or penalize enrichment.")
    parser.add_argument("--min-specificity", type=float, default=None,
                        help="Threshold for τ (0..1). Used by --specificity-mode filter|penalize. Example: 0.8.")

    args = parser.parse_args()

    # -------- Load input --------
    print("\033[33mLoading input...\033[0m")
    df = pd.read_csv(args.input_file, sep="\t")

    # Sanity check required columns
    required = [args.gene_col, args.gene_name_col, args.cell_type_col, args.value_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Input is missing required columns: {missing}")

    # -------- Coerce numeric on value_col and report non-numeric --------
    df[args.value_col], non_numeric = coerce_numeric(df[args.value_col])
    if non_numeric:
        print(f"\033[31mNon-numeric values in {args.value_col}: {non_numeric}\033[0m")
    else:
        print("\033[32mGOOD DATA nCPM. NO Non numeric\033[0m")

    # -------- Optional batch normalization (per cell type by default) --------
    if args.batch_col and args.batch_col in df.columns:
        print(f"\033[33mApplying batch normalization per '{args.batch_col}'...\033[0m")
        df = batch_normalize_if_needed(
            df, value_col=args.value_col,
            batch_col=args.batch_col,
            batch_normalize=args.batch_normalize
        )
    else:
        if args.batch_col:
            print(f"\033[33mBatch column '{args.batch_col}' not found; skipping normalization.\033[0m")

    # -------- Aggregation across clusters within each (Gene × Cell type) --------
    print("\033[33mAggregating across clusters within each (Gene × Cell type)...\033[0m")
    agg_df = aggregate_with_celltype(
        df,
        gene_col=args.gene_col,
        gene_name_col=args.gene_name_col,
        cell_type_col=args.cell_type_col,
        value_col=args.value_col,
        cluster_col=args.cluster_col if args.cluster_col in df.columns else None,
        weighted=(args.weighted == "on"),
        weight_col=args.weight_col if args.weight_col in df.columns else None,
        cluster_aggregate=args.cluster_aggregate
    )

    # -------- Optional filter by min clusters --------
    if args.min_clusters is not None:
        before_n = len(agg_df)
        agg_df = agg_df[agg_df["clusters_used"] >= args.min_clusters].copy()
        after_n = len(agg_df)
        print(
            f"\033[33mFiltered by clusters_used >= {args.min_clusters}. "
            f"Dropped {before_n - after_n} row(s); {after_n} row(s) remain.\033[0m"
        )

    # -------- Pre-enrichment expression threshold --------
    if args.min_expr_threshold > 0:
        before_n = len(agg_df)
        agg_df = agg_df[agg_df["avg_nCPM"] >= args.min_expr_threshold].copy()
        after_n = len(agg_df)
        print(
            f"\033[33mFiltered rows with avg_nCPM < {args.min_expr_threshold}. "
            f"Dropped {before_n - after_n} row(s); {after_n} row(s) remain.\033[0m"
        )

    # -------- Drop genes with no expression across all cell types --------
    print("\033[33mDropping genes with no expression across all cell types...\033[0m")
    agg_df, dropped_genes, rows_removed = drop_genes_with_no_expression(
        agg_df, expr_col="avg_nCPM", treat_nan_as_zero=args.treat_nan_as_zero
    )
    print(
        f"\033[31mDropped {len(dropped_genes)} gene(s), removing {rows_removed} row(s).\033[0m"
    )

    # -------- Compute τ and merge as a report column --------
    print("\033[33mComputing Yanai's τ per gene...\033[0m")
    agg_df = compute_tau(agg_df, gene_col=args.gene_col, expr_col="avg_nCPM")

    # -------- Specificity mode: filter or penalize --------
    if args.specificity_mode == "filter" and args.min_specificity is not None:
        print(f"\033[33mApplying τ filtering (threshold={args.min_specificity})...\033[0m")
        agg_df, n_total, n_dropped = apply_tau_filter(
            agg_df, gene_col=args.gene_col, min_specificity=args.min_specificity
        )
        print(
            f"\033[33mτ filter: {n_dropped} gene(s) dropped out of {n_total}. "
            f"{agg_df[args.gene_col].nunique()} gene(s) remain.\033[0m"
        )
    elif args.specificity_mode == "filter" and args.min_specificity is None:
        print("\033[33mWARNING: --specificity-mode filter set but --min-specificity not provided; skipping filter.\033[0m")

    # -------- Unique cell types list --------
    n_cell_types = agg_df[args.cell_type_col].dropna().nunique()
    unique_cell_types = sorted(agg_df[args.cell_type_col].dropna().unique().tolist())
    pd.Series(unique_cell_types, name=args.cell_type_col).to_csv("unique_cell_types.tsv", sep="\t", index=False)
    print(f"\033[32mNumber of unique cell types: {n_cell_types}\033[0m")

    # -------- Count unique genes --------
    n_unique_genes = agg_df[args.gene_col].dropna().nunique()
    print(f"\033[32mNumber of unique genes remaining: {n_unique_genes}\033[0m")

    # -------- Enrichment --------
    print("\033[33mCalculating Enrichment Scores...\033[0m")
    agg_df = add_enrichment(
        agg_df=agg_df,
        gene_col=args.gene_col,
        value_col="avg_nCPM",
        out_col="Enrichment score",
        min_background=args.min_background,
        min_expression=args.min_expression,
        min_count=args.min_count,
        pseudocount=args.pseudocount,
        pseudocount_to_numerator=args.pseudocount_to_numerator,
        clip_max=args.clip_max
    )
    print("\033[33mDone calculating.\033[0m")

    # -------- Penalize enrichment by τ (optional) --------
    if args.specificity_mode == "penalize":
        print("\033[33mApplying τ penalization...\033[0m")
        penalty = agg_df["specificity_tau"].clip(lower=0.0, upper=1.0)
        agg_df["Enrichment score (tau penalized)"] = agg_df["Enrichment score"] * penalty
        EPS = 1e-12
        log2p = np.log2(agg_df["Enrichment score (tau penalized)"].clip(lower=EPS))
        log2p = pd.Series(log2p, index=agg_df.index).mask(agg_df["Enrichment score (tau penalized)"] <= 0)
        agg_df["log2_enrichment_penalized"] = log2p

    # -------- Single cell-type genes --------
    gene_celltype_counts = agg_df.groupby(args.gene_col)[args.cell_type_col].transform("nunique")
    agg_df["single_cell_type_gene"] = (gene_celltype_counts == 1)
    min_ct_per_gene = int(gene_celltype_counts.min()) if len(gene_celltype_counts) else 0
    print(
        f"\033[33mMinimum number of cell types per gene (across genes): {min_ct_per_gene}\033[0m"
    )

    single_cell_rows = agg_df[agg_df["single_cell_type_gene"].fillna(False)].copy()
    if not single_cell_rows.empty:
        n_genes_single = single_cell_rows[args.gene_col].nunique()
        print(f"\033[32mFound {n_genes_single} gene(s) expressed in exactly one cell type- according to raw data.\033[0m")
        single_cell_rows.to_csv("single_cell_type_gene_rows.tsv", sep="\t", index=False)
    else:
        print("\033[33mNo genes were found that are only expressed in one cell type.\033[0m")

    # -------- Sort and save full --------
    sort_col = args.sort_by
    agg_df_sorted = agg_df.sort_values(by=sort_col, ascending=False)
    agg_df_sorted.to_csv(args.output_file, sep="\t", index=False)
    print(f"\033[33mSaved full table: {args.output_file}\033[0m")

    # -------- Top-N overall --------
    print("\033[33mSaving top-N overall...\033[0m")
    top_overall = top_n_overall(agg_df_sorted, sort_by=sort_col, n=args.top_n)
    suffix = 'log2' if 'log2' in sort_col else ('penalized' if 'penalized' in sort_col else 'enrichment')
    top_overall_file = f"top{args.top_n}_{suffix}.tsv"
    top_overall.to_csv(top_overall_file, sep="\t", index=False)
    print(f"\033[33mSaved: {top_overall_file}\033[0m")

    # -------- Top-N per cell type (skip if N=0) --------
    if args.per_cell_type_top_n and args.per_cell_type_top_n > 0:
        print("\033[33mSaving top-N per cell type...\033[0m")
        top_pct = top_n_per_cell_type(agg_df, cell_type_col=args.cell_type_col, sort_by=sort_col, n=args.per_cell_type_top_n)
        top_pct_suffix = 'log2' if 'log2' in sort_col else ('penalized' if 'penalized' in sort_col else 'enrichment')
        top_pct_file = f"top_per_cell_type_{args.per_cell_type_top_n}_{top_pct_suffix}.tsv"
        top_pct.to_csv(top_pct_file, sep="\t", index=False)
        print(f"\033[33mSaved: {top_pct_file}\033[0m")
    else:
        print("\033[33mPer-cell-type top-N export disabled (N=0).\033[0m")


if __name__ == "__main__":
    main()
```

Following is the runner for that script (you can change options. I found following options gave me a nice result)

run_celltype_enrichment_v1_4.sh
```sh

#!/usr/bin/env bash
# run_celltype_enrichment_with_options.sh (v1.4 with τ report & penalization)
# Usage:
#   ./run_celltype_enrichment_with_options.sh [overrides...]
# Examples:
#   ./run_celltype_enrichment_with_options.sh --specificity-mode filter --min-specificity 0.8
#   ./run_celltype_enrichment_with_options.sh --specificity-mode penalize --min-specificity 0.8 --sort-by "log2_enrichment_penalized"

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

args=(
  # Path to your input TSV containing: Gene, Gene name, Cell type, Cluster, Read count, nCPM
  --input-file rna_single_cell_cluster.tsv
  # Path for the full enrichment output table (TSV)
  --output-file celltype_enrichment.tsv
  # Number of top rows to export overall (after sorting)
  --top-n 100
  # Number of top rows to export per cell type (set 0 to disable)
  --per-cell-type-top-n 20

  # Column name for gene ID (e.g., ENSG...)
  --gene-col "Gene"
  # Column name for gene symbol/display name
  --gene-name-col "Gene name"
  # Column name for cell type labels
  --cell-type-col "Cell type"
  # Column name for cluster IDs (replicates within a cell type)
  --cluster-col "Cluster"
  # Batch column used for median scaling (default: per cell type)
  --batch-col "Cell type"
  # Column name for expression values to aggregate
  --value-col "nCPM"
  # Column name for weights used in weighted aggregation
  --weight-col "Read count"

  # Toggle weighted aggregation across clusters: on=weighted mean, off=unweighted
  --weighted on
  # If weighted is off, choose how to aggregate clusters: mean or median
  --cluster-aggregate mean

  # Minimum number of clusters required to keep a (Gene × Cell type) row
  --min-clusters 2
  # Drop rows with avg_nCPM below this threshold before enrichment. This helps you to deal with weird values generated by 0 expression 
  --min-expr-threshold 0.00
  # Treat NaN as zero when deciding to drop genes with no expression (comment to disable)
  # --treat-nan-as-zero

  # Floor for the denominator in enrichment to avoid tiny values
  --min-background 1e-3
  # Minimum numerator expression required to compute enrichment
  --min-expression 0.001
  # Minimum number of cell-type entries per gene to compute enrichment
  --min-count 2
  # Add a small constant to stabilize denominators (uncomment to enable)
  --pseudocount 0.001
  # Also add pseudocount to numerator (use with --pseudocount)
  # --pseudocount-to-numerator
  # Cap extremely large enrichment ratios (uncomment to enable)
  # --clip-max 100

  # Sort outputs by: raw/log2 or penalized variants
  --sort-by "log2_enrichment_penalized"

  # Apply median scaling per batch (here: per cell type)
  --batch-normalize median_scale

  # Yanai's τ usage:
  #   off      → only report `specificity_tau` (no filtering/penalization)
  #   filter   → drop genes whose τ < threshold (strict specificity)
  #   penalize → multiply enrichment by τ (keeps genes but down-ranks broad ones)
  --specificity-mode penalize
  # τ threshold (0..1) used by filter/penalize modes; e.g., 0.8 for high specificity
  --min-specificity 0.8
)

# Allow overrides on the CLI (last wins)
args+=("$@")

# Call the Python script
python3 "${script_dir}/celltype_enrichment_v1_4.py" "${args[@]}"
``
```
save this along with celltype_enrichment_v1_4.py and run it like following overriding whatever arguments needed. The arguments you do not use here will default to the arguments provided in run_celltype_enrichment_v1_4.sh

Run celltype_enrichment_v1_4 with the runner, overriding needed arguments
```bash
./run_celltype_enrichment_v1_4.sh \
  --input-file rna_single_cell_cluster_filtered_rows_alpha_mad.tsv \
  --output-file enrichV1_4_3clusters.tsv \
  --min-clusters 3 \
  --specificity-mode penalize \
```
Then you can plot the results with the universal plotter

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#universal-interactive-plot-maker

Save that script and the runner below, make the runner an executable and run it as explained overrriding whatever commands needed

Example usage
```
./run_universal_plot_maker_with_options.sh --file enrichV1_4_3clusters.tsv --out enrichV1_4_3clusters.html
```

run_universal_plot_maker_with_options.sh
```
#!/usr/bin/env bash
# Usage:
#   ./run_universal_plot_maker_with_options.sh [overrides...]
# Example:
#   ./run_universal_plot_maker_with_options.sh --file simple_enrich_1_clustor.tsv --out simple_enrich_1_clustor_plot.html

set -euo pipefail

# Resolve the directory of this script so we can find the Python file reliably
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Default arguments (can be overridden by CLI options appended below)
args=(
  --file enrichV1_4_3clusters.tsv          # REQUIRED: Input data file (TSV/CSV/etc.)
  --out enrichV1_4_3clusters.html                      # Output HTML file name
  --top 35000                                       # Top N rows to plot (default: 100)
  --dedupe-policy mean                             # [max|mean|median|first] aggregation
  # --log                                            # Use log scale on numeric axis
  --linear                                       # Use linear scale instead (mutually exclusive with --log)
  # --horizontal                                   # Horizontal bars (better for long labels)
  --self-contained                                 # Embed Plotly.js for offline HTML
  # --log-digits D2                                # Log-axis minor ticks: D1 (all) or D2 (2 & 5)
  # --lang en-CA                                   # HTML lang attribute (default: en)
  --initial-zoom 100                               # Initial number of bars visible on load
  # --sep $'\t'                                    # Field separator (auto-detected if omitted)
  --x-col "Gene name"                              # Column for X axis (numeric if horizontal)
  --y-col "log2_enrichment_penalized"                       # Column for Y axis (categorical if horizontal)
  --label-col "Gene name"                               # Explicit label column (optional)
  # --value-col Score                              # Explicit numeric value column (optional)
  --group-col "Cell type"                          # Column for color grouping (legend)
  --search-col "Gene name"                         # Column used for search box
  --details "Gene" "Gene name" "Cell type" "clusters_used" "avg_nCPM" "weight_sum" "specificity_tau" "log2_enrichment" "log2_enrichment_penalized" "Enrichment score" # Extra columns
)

# Append any CLI overrides so the LAST occurrence of options wins
args+=( "$@" )

# Invoke the Python script with the collected arguments
python3 "${script_dir}/universal_plot_maker.py" "${args[@]}"
``
```

