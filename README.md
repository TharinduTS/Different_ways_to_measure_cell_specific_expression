# Different_ways_to_measure_cell_specific_expression
This page explains 3 different ways I used to build a tool to measure cell specific or enriched gene expression with bash/python and HTML

#Method 01 - Simple Enrichment Scores

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#method-01---simple-enrichment-scores


#Method 02 - AUROC SCORES

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#method-02---auroc-scores

#Method 03 - Enrichment score with more robust statistical backup systems to filter useful data

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#method-03---enrichment-score-with-more-robust-statistical-backup-systems-to-filter-useful-data

#And finally I have a script that helps to compare the effectiveness between different methods

https://github.com/TharinduTS/Different_ways_to_measure_cell_specific_expression/blob/main/README.md#comparing-method-effectiveness

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
Then I wanted to plot this in an interactive way. This allows user to decide how many values to include in the plot (top x amount
of data points from the previous output). This resulting tool lets the user select different cell types and genes to see the selected data points. In addition, you can export a tsv file with the selected values.

You can run this as
```
python plot_interactive.py \
  -f my_custom_clustor_enrichment.tsv \
  -o simple_enrichment_with_custom_clustor_number.html \
  -n 1000 \
  --log \
  --self-contained \
  --initial-zoom 100 \
  --lang en-CA
```
-f -input_file

-o -output file

-n -number of datapoints to show in the plot

--log -log transforms the enrichment scores to make it easier to view

--self-contained -lets the user use this tool offline

--initial-zoom -decides the number of data points(bars) to show when you open the tool

--lang en-CA -Language

--horizontal -plots horizontal instead default plot

following is the python script7

plot_interactive.py
```py
#!/usr/bin/env python3
"""
Interactive gene enrichment plot with filter/search and TSV export.

- Loads TSV containing gene enrichment results.
- Builds an interactive Plotly bar chart (log or linear, horizontal or vertical).
- Injects an HTML UI with:
  - Cell type filter
  - Gene search
  - Reset button
  - **Export TSV** button (exports filtered view OR selected points)
- Saves a standards-compliant HTML file (with <!DOCTYPE html>, <html lang="...">, and viewport meta).

Usage example:
    python test.py \
      -f all_gene_cell_enrichment_data.tsv \
      -o interactive_markers_testing.html \
      -n 1000 \
      --log \
      --self-contained \
      --horizontal \
      --initial-zoom 100 \
      --lang en-CA
"""

import argparse
import json
import re
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


# ---------- Data Loading ----------
def load_tsv(path: str) -> pd.DataFrame:
    """
    Load the TSV file with required columns and prepare 'Label' for plotting.
    Expects columns:
      "Gene", "Gene name", "Cell type", "avg_nCPM", "clusters_used",
      "Enrichment score", "single_cell_type_gene"
    """
    df = pd.read_csv(path, sep="\t")
    # Clean column names
    df.columns = [c.strip() for c in df.columns]
    expected = [
        "Gene",
        "Gene name",
        "Cell type",
        "avg_nCPM",
        "clusters_used",
        "Enrichment score",
        "single_cell_type_gene",
    ]
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise ValueError(f"Missing expected columns: {missing}")
    # Coerce numeric columns
    for col in ["avg_nCPM", "clusters_used", "Enrichment score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    # Label: prefer Gene name, fallback to Gene
    df["Label"] = df["Gene name"].fillna(df["Gene"])
    # Drop rows missing critical fields (label/enrichment)
    df = df.dropna(subset=["Label", "Enrichment score"])
    return df


# ---------- Figure Building ----------
def build_fig(
    df: pd.DataFrame,
    top_n: int = 100,
    use_log: bool = True,
    orientation: str = "v",  # 'v' or 'h'
    log_dtick: str | None = None,  # "D1" or "D2"
    initial_zoom: int | None = None,  # initial number of bars to show
):
    """
    Build the Plotly figure and construct a JSON payload for client-side filtering/export.
    """
    # Filter out non-positive values for log scale
    if use_log:
        df = df[df["Enrichment score"] > 0]

    # Sort by enrichment and select top N
    df_plot = df.sort_values("Enrichment score", ascending=False).head(top_n).copy()

    # Prepare customdata for click/selection (detailCols order)
    detail_cols = [
        "Gene",
        "Gene name",
        "Cell type",
        "avg_nCPM",
        "clusters_used",
        "Enrichment score",
        "single_cell_type_gene",
    ]
    customdata = df_plot[detail_cols].values

    # Stable color mapping by cell type
    colors = px.colors.qualitative.Safe
    ctypes = sorted(df_plot["Cell type"].astype(str).unique())
    cmap = {ct: colors[i % len(colors)] for i, ct in enumerate(ctypes)}
    bar_colors = df_plot["Cell type"].astype(str).map(cmap)

    # Build bar trace
    if orientation == "v":
        x_vals = df_plot["Label"]
        y_vals = df_plot["Enrichment score"]
        bar = go.Bar(
            x=x_vals,
            y=y_vals,
            customdata=customdata,
            marker_color=bar_colors,
            hovertemplate=(
                "<b>%{x}</b><br>"
                "Enrichment: %{y:.2f}<br>"
                "Cell type: %{customdata[2]}<br>"
                "avg_nCPM: %{customdata[3]}<br>"
                "clusters_used: %{customdata[4]}<br>"
                "<extra></extra>"
            ),
            # Enable visual feedback for selection
            selected={"marker": {"opacity": 1.0}},
            unselected={"marker": {"opacity": 0.5}},
        )
    else:
        x_vals = df_plot["Enrichment score"]
        y_vals = df_plot["Label"]
        bar = go.Bar(
            x=x_vals,
            y=y_vals,
            orientation="h",
            customdata=customdata,
            marker_color=bar_colors,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Enrichment: %{x:.2f}<br>"
                "Cell type: %{customdata[2]}<br>"
                "avg_nCPM: %{customdata[3]}<br>"
                "clusters_used: %{customdata[4]}<br>"
                "<extra></extra>"
            ),
            selected={"marker": {"opacity": 1.0}},
            unselected={"marker": {"opacity": 0.5}},
        )

    fig = go.Figure(data=[bar])

    # Layout
    fig.update_layout(
        title="Gene Enrichment",
        xaxis_title=("Enrichment score" if orientation == "h" else "Gene name"),
        yaxis_title=("Gene name" if orientation == "h" else "Enrichment score"),
        template="plotly_white",
        height=700,
        bargap=0.2,
        margin=dict(l=60, r=30, t=60, b=120),
        hovermode="closest",
        legend_title_text="Cell type",
        dragmode="select",  # Enable box selection by default
    )

    # Axes tweaks
    if orientation == "v":
        fig.update_xaxes(tickangle=-45, automargin=True, categoryorder="array", categoryarray=list(df_plot["Label"]))
        if use_log:
            kwargs = dict(type="log", title="Enrichment score (log scale)")
            if log_dtick in ("D1", "D2"):
                kwargs["dtick"] = log_dtick
            fig.update_yaxes(**kwargs)
        else:
            fig.update_yaxes(title="Enrichment score")
    else:
        fig.update_yaxes(automargin=True, categoryorder="array", categoryarray=list(df_plot["Label"]))
        if use_log:
            kwargs = dict(type="log", title="Enrichment score (log scale)")
            if log_dtick in ("D1", "D2"):
                kwargs["dtick"] = log_dtick
            fig.update_xaxes(**kwargs)
        else:
            fig.update_xaxes(title="Enrichment score")

    # Initial zoom (index-based range on category axis)
    if initial_zoom is not None:
        zoom_n = max(1, min(int(initial_zoom), len(df_plot)))
        if orientation == "v":
            fig.update_xaxes(range=[-0.5, zoom_n - 0.5])
        else:
            fig.update_yaxes(range=[-0.5, zoom_n - 0.5])

    # Build JSON payload for reliable client-side filtering/export
    rows = df_plot[
        [
            "Gene",
            "Gene name",
            "Cell type",
            "avg_nCPM",
            "clusters_used",
            "Enrichment score",
            "single_cell_type_gene",
            "Label",
        ]
    ].to_dict(orient="records")

    payload = {
        "rows": rows,
        "colors": cmap,  # cell type -> color
        "orientation": orientation,
        "detail_cols": detail_cols,  # authoritative order for export/customdata
        "label_col": "Label",
        "enrich_col": "Enrichment score",
    }
    return fig, payload


# ---------- HEAD meta for mobile ----------
HEAD_SNIPPET = """<meta name="viewport" content="width=device-width, initial-scale=1">"""


# ---------- Enhanced HTML (Details + Filters + Search + Export TSV) ----------
# Raw string to preserve JS regex literals like /\s+/g without needing double backslashes
DETAILS_SNIPPET = r"""
<style>
  #controls {
    font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif;
    margin: 12px 0;
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 8px;
  }
  #controls label {
    font-size: 14px;
    color: #333;
  }
  #controls select, #controls input {
    padding: 6px 8px;
    font-size: 14px;
  }
  #controls button {
    padding: 6px 10px;
    font-size: 14px;
    cursor: pointer;
  }
  #count-info {
    font-size: 13px;
    color: #666;
    margin-left: 8px;
  }
  #details-panel {
    font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif;
    margin-top: 12px;
    padding: 12px;
    border: 1px solid #ddd;
    border-radius: 8px;
    background: #fafafa;
    white-space: pre-wrap;
    font-size: 14px;
  }
  #details-panel table {
    border-collapse: collapse;
    width: 100%;
  }
  #details-panel td {
    padding: 4px 8px;
    border-bottom: 1px solid #eee;
    vertical-align: top;
  }
  #details-panel td:first-child {
    font-weight: 600;
    color: #333;
    width: 220px;
  }
  .plotly-graph-div { width: 100% !important; }
</style>

<div id="controls">
  <label>Cell type:
    <select id="celltype-filter"><option value="__all__">All cell types</option></select>
  </label>
  <label style="margin-left:12px;">Gene search:
    <input type="text" id="gene-search" placeholder="e.g., LALBA">
  </label>
  <button id="reset-filters" title="Clear filters">Reset</button>
  <button id="export-tsv" title="Export current view as TSV">Export TSV</button>
  <span id="count-info"></span>
</div>

<div id="details-panel"><em>Click a bar to see full row details here. Lasso/box-select points to export only selected.</em></div>

<!-- The Python side will inject a <script id="chart-data" type="application/json">...</script> before this script -->
<script>
(function(){
  // Wait until the Plotly figure exists and is initialized
  function initOnceReady() {
    var gd = document.querySelector('div.plotly-graph-div') || document.querySelector('.js-plotly-plot');
    var dataTag = document.getElementById('chart-data');
    if (!gd || !gd.data || !gd.data.length || !dataTag) {
      setTimeout(initOnceReady, 50);
      return;
    }

    // Parse embedded JSON payload from Python
    var payload = {};
    try {
      payload = JSON.parse(dataTag.textContent);
    } catch(e) {
      console.error("JSON parse error:", e);
      return;
    }

    var rows = Array.isArray(payload.rows) ? payload.rows : [];
    var colors = payload.colors || {};
    var orientation = payload.orientation || 'v';
    var detailCols = payload.detail_cols || ["Gene","Gene name","Cell type","avg_nCPM","clusters_used","Enrichment score","single_cell_type_gene"];
    var labelCol = payload.label_col || "Label";
    var enrichCol = payload.enrich_col || "Enrichment score";

    // Controls
    var sel = document.getElementById('celltype-filter');
    var searchBox = document.getElementById('gene-search');
    var resetBtn = document.getElementById('reset-filters');
    var exportBtn = document.getElementById('export-tsv');
    var panel = document.getElementById('details-panel');
    var countInfo = document.getElementById('count-info');

    // Number formatter for details panel
    function fmtNumber(val, maxDigits=2) {
      if (val === null || val === undefined || val === "" || isNaN(val)) return String(val ?? "");
      var num = Number(val);
      if (!isFinite(num)) return String(val);
      return num.toLocaleString('en-US', { maximumFractionDigits: maxDigits });
    }

    // Build base arrays from rows (authoritative source)
    var N = rows.length;

    // Populate unique cell types into dropdown
    var cellTypesSet = {};
    for (var i = 0; i < N; i++) {
      var ct = rows[i]["Cell type"];
      if (ct != null) {
        ct = String(ct).trim();
        if (ct.length) cellTypesSet[ct] = true;
      }
    }
    var cellTypes = Object.keys(cellTypesSet).sort();
    cellTypes.forEach(function(ct){
      var opt = document.createElement('option');
      opt.value = ct;
      opt.textContent = ct;
      sel.appendChild(opt);
    });

    // Details panel render function
    function renderDetails(customdata) {
      if (!customdata || !customdata.length) return;
      var html = "<table>";
      for (var i = 0; i < detailCols.length; i++) {
        var key = detailCols[i];
        var val = customdata[i];
        if (key === "avg_nCPM" || key === "Enrichment score") {
          val = fmtNumber(val, 2);
        }
        html += "<tr><td>" + key + "</td><td>" + (val === null ? "" : val) + "</td></tr>";
      }
      html += "</table>";
      panel.innerHTML = html;
      panel.scrollIntoView({ behavior: "smooth", block: "nearest" });
    }

    // Current filtered rows (used for export when no selection)
    var filteredRows = [];
    // Current selected rows (built from plotly_selected)
    var selectedRows = [];

    // Helper: convert rows to TSV using detailCols order
    function rowsToTSV(rowsArr) {
      var cols = detailCols.slice();
      var lines = [];
      // Header
      lines.push(cols.join('\t'));
      // Body
      for (var i = 0; i < rowsArr.length; i++) {
        var r = rowsArr[i];
        var out = [];
        for (var j = 0; j < cols.length; j++) {
          var key = cols[j];
          var val = r[key];
          // Normalize and escape tabs/newlines
          var s = (val === null || val === undefined) ? '' : String(val);
          s = s.replace(/\r?\n/g, ' ').replace(/\t/g, ' ');
          out.push(s);
        }
        lines.push(out.join('\t'));
      }
      return lines.join('\n');
    }

    // Trigger a browser download of a given string as a file
    function downloadTextFile(filename, text) {
      var blob = new Blob([text], { type: 'text/tab-separated-values;charset=utf-8' });
      var url = URL.createObjectURL(blob);
      var a = document.createElement('a');
      a.href = url;
      a.download = filename;
      document.body.appendChild(a);
      a.click();
      setTimeout(function(){
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
      }, 0);
    }

    // Export button: prefer selected rows if present; otherwise use filtered rows
    exportBtn.addEventListener('click', function(){
      var toExport = (selectedRows && selectedRows.length) ? selectedRows : filteredRows;
      if (!toExport || !toExport.length) {
        alert('No rows to export. Use filters/search or select points first.');
        return;
      }
      var tsv = rowsToTSV(toExport);
      // Build a descriptive filename
      var suffix = [];
      if (selectedRows && selectedRows.length) suffix.push('selection');
      var selectedCT = sel.value;
      var q = (searchBox.value || '').trim();
      if (selectedCT && selectedCT !== '__all__') suffix.push(selectedCT.replace(/\s+/g, '_'));
      if (q) suffix.push('q_' + q.replace(/\s+/g, '_'));
      var fname = 'export' + (suffix.length ? '_' + suffix.join('_') : '') + '.tsv';
      downloadTextFile(fname, tsv);
    });

    // Filtering function: update plot + filteredRows cache
    function applyFilter() {
      var selectedCT = sel.value;
      var q = (searchBox.value || "").trim().toLowerCase();

      var fx = [], fy = [], fcd = [], fcol = [];
      filteredRows = []; // reset cache

      for (var i = 0; i < N; i++) {
        var row = rows[i];
        var ct = (row["Cell type"] == null ? "" : String(row["Cell type"]).trim());
        var geneName = String(row["Gene name"] || "").toLowerCase();
        var labelLower = String(row[labelCol] || "").toLowerCase();
        var ctOk = (selectedCT === "__all__") || (selectedCT === ct);
        var qOk = (!q) || geneName.includes(q) || labelLower.includes(q);

        if (ctOk && qOk) {
          filteredRows.push(row); // cache for export

          if (orientation === 'h') {
            fx.push(row[enrichCol]);   // enrichment on x
            fy.push(row[labelCol]);    // labels on y
          } else {
            fx.push(row[labelCol]);    // labels on x
            fy.push(row[enrichCol]);   // enrichment on y
          }
          fcd.push(detailCols.map(function(k){ return row[k]; }));
          fcol.push(colors[ct] || '#636EFA');
        }
      }

      // Restyle robustly
      var update = { x: [fx], y: [fy], customdata: [fcd], marker: [{ color: fcol }] };
      Plotly.restyle(gd, update, [0]);

      // Keep category order aligned with filtered labels
      var relayout = {};
      if (orientation === 'h') {
        relayout['yaxis.categoryorder'] = 'array';
        relayout['yaxis.categoryarray'] = fy;
      } else {
        relayout['xaxis.categoryorder'] = 'array';
        relayout['xaxis.categoryarray'] = fx;
      }
      Plotly.relayout(gd, relayout);
      countInfo.textContent = (fx.length) + " shown of " + N;

      // Clearing selection state: when filter changes, drop previous selection
      selectedRows = [];
    }

    // Wire filter/search/reset
    sel.addEventListener('change', applyFilter);
    searchBox.addEventListener('input', applyFilter);
    resetBtn.addEventListener('click', function(){
      sel.value = "__all__";
      searchBox.value = "";
      applyFilter();
    });

    // Initial count
    countInfo.textContent = N + " shown of " + N;

    // Click → details panel
    var clickHandler = function(evt) {
      var e = evt && evt.points ? evt : (evt && evt.detail ? evt.detail : null);
      if (!e || !e.points || !e.points.length) return;
      var p = e.points[0];
      renderDetails(p.customdata);
    };
    if (typeof gd.on === 'function') gd.on('plotly_click', clickHandler);
    else gd.addEventListener('plotly_click', clickHandler);

    // Selection → build selectedRows (reconstruct objects using detailCols)
    function collectSelectedRows(evt) {
      selectedRows = [];
      if (!evt || !evt.points) return;
      for (var i = 0; i < evt.points.length; i++) {
        var p = evt.points[i];
        var cd = p.customdata;
        if (Array.isArray(cd)) {
          var obj = {};
          for (var j = 0; j < detailCols.length; j++) {
            obj[detailCols[j]] = cd[j];
          }
          selectedRows.push(obj);
        }
      }
    }
    if (typeof gd.on === 'function') gd.on('plotly_selected', collectSelectedRows);
    else gd.addEventListener('plotly_selected', collectSelectedRows);
  }

  // Start readiness check
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initOnceReady);
  } else {
    initOnceReady();
  }
})();
</script>
"""


# ---------- HTML Saving (standards-compliant, embeds JSON payload) ----------
def save_html(fig: go.Figure, out_path: str, payload: dict, self_contained: bool = False, lang_code: str = "en"):
    """
    Save a robust, standards-compliant HTML:
    - Ensures <!DOCTYPE html> at the very beginning (No Quirks Mode).
    - Adds lang="..." on <html>.
    - Injects viewport meta just after <head>.
    - Embeds a JSON payload with rows & colors for client-side filtering and export.
    - Appends DETAILS_SNIPPET before </body> using a callable replacement to avoid backslash parsing.
    """
    include_js = True if self_contained else "cdn"
    html = fig.to_html(full_html=True, include_plotlyjs=include_js)

    # 1) Ensure <!DOCTYPE html> at the very beginning
    html = html.lstrip(" \ufeff\r\n\t")
    if not re.match(r"(?is)^<!doctype\s+html>", html):
        html = "<!DOCTYPE html>\n" + html

    # 2) Ensure <html lang="...">
    if re.search(r"(?is)<html(?:\s[^>]*)?>", html):
        if not re.search(r"(?is)<html[^>]*\blang\s*=", html):
            html = re.sub(r"(?is)<html(\s*)>", f'<html lang="{lang_code}">', html, count=1)
    else:
        html = f"<!DOCTYPE html>\n<html lang=\"{lang_code}\">\n{html}\n</html>"

    # 3) Inject <meta name="viewport"> immediately after <head>
    if re.search(r"(?is)<head\s*>", html):
        html = re.sub(r"(?is)<head\s*>", "<head>\n" + HEAD_SNIPPET + "\n", html, count=1)
    else:
        html = re.sub(
            r"(?is)(<html[^>]*>)",
            r"\1\n<head>\n" + HEAD_SNIPPET + "\n</head>\n",
            html,
            count=1,
        )

    # 4) Embed JSON payload (safe single <script> with application/json)
    json_str = json.dumps(payload, ensure_ascii=False)
    data_block = f'<script id="chart-data" type="application/json">{json_str}</script>\n'

    # Place the data block just before DETAILS_SNIPPET and </body>
    if re.search(r"(?is)</body\s*>", html):
        # Use a callable replacement to avoid backslash interpretation in the replacement string
        html = re.sub(
            r"(?is)</body\s*>",
            lambda m: data_block + DETAILS_SNIPPET + "\n</body>",
            html,
            count=1,
        )
    else:
        # Ensure a body exists, then append our blocks
        if not re.search(r"(?is)<body[^>]*>", html):
            html = re.sub(r"(?is)</head\s*>", "</head>\n<body>\n", html, count=1)
        html = html + "\n" + data_block + DETAILS_SNIPPET + "\n</body>\n</html>"

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(
        description="Make interactive gene enrichment chart (TSV → HTML) with details, filters, search, and TSV export"
    )
    ap.add_argument("--file", "-f", required=True, help="Input TSV file (e.g., all_gene_cell_enrichment_data.tsv)")
    ap.add_argument("--out", "-o", default="interactive_markers.html", help="Output HTML file")
    ap.add_argument("--top", "-n", type=int, default=100, help="Top N genes by enrichment")
    ap.add_argument("--log", action="store_true", help="Use log scale on enrichment axis")
    ap.add_argument("--linear", action="store_true", help="Use linear scale on enrichment axis")
    ap.add_argument("--horizontal", action="store_true", help="Use horizontal bars (better for long labels)")
    ap.add_argument("--self-contained", action="store_true", help="Embed plotly.js for fully offline HTML")
    ap.add_argument("--log-digits", choices=["D1", "D2"], help="Log-axis minor digits: D1 (all) or D2 (2 & 5)")
    ap.add_argument("--lang", default="en", help="HTML lang attribute (e.g., 'en', 'en-CA')")
    ap.add_argument("--initial-zoom", type=int, default=100, help="Initial number of bars to show on load")

    args = ap.parse_args()

    if args.log and args.linear:
        raise SystemExit("Choose either --log or --linear (not both).")

    use_log = True if args.log or (not args.linear) else False
    orientation = "h" if args.horizontal else "v"

    df = load_tsv(args.file)
    fig, payload = build_fig(
        df,
        top_n=args.top,
        use_log=use_log,
        orientation=orientation,
        log_dtick=args.log_digits,
        initial_zoom=args.initial_zoom,
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
  --keep-constant-auc05
```
--file -input file

--out-all -output file name for all files 

--out-top -output file name for top 100 file

--auroc-min -minimum AUROC to keep. 0.90: Keeps gene–cell-type pairs with strong separation from other cell types in the same tissue.

--median-min -minimum median nCPM in target cell type. 1.0: Avoids keeping “specific but trivially low” expression.

--ratio-min -minimum (median_target+α)/(median_others+α).2.0: Requires the median in the target cell type to be at least ~2× higher than the median of others (with a small pseudocount α=0.1 to stabilize zeros).

--clusters-min -minimum number of clusters in target cell type. 2: Ensures there’s at least minimal within-cell-type replication.

--alpha -sets the pseudocount added when calculating the robust ratio to stabilize division by very small or zero medians, which is important because it prevents inflated ratios and ensures more reliable specificity scoring in sparse

--keep-constant-auc05 -Keep AUROC = 0.5 for constant-score groups (instead of NaN):


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
# Comparing method effectiveness
Finally I use the following R script to compare the effectiveness of different methods with different variables

For this I have 

1) A tsv file with information on known genes that are specific in expression to certain cell types that look like following

methods_comparison.txt
```txt
Gene	Cell Type (with extra info)	Cell type (for comparison)	Reference	Specificity Score (human)
MYH7	Cardiomyocytes	cardiomyocytes	"Canonical marker, PanglaoDB"	n/a
TNNI3	Cardiomyocytes	cardiomyocytes	"PanglaoDB, 100% specificity"	n/a
TNNT2	Cardiomyocytes	cardiomyocytes	PanglaoDB specificity	n/a
ALB	Hepatocytes	hepatocytes	Hepatocyte marker article	1
CYP2E1	Hepatocytes (zone-specific)	hepatocytes	Biocompare hepatocyte zonation	0.575
GFAP	Astrocytes	astrocytes	Canonical astrocyte marker	0.344
```
2) Following R script saved in the same directory as methods_comparison.txt

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

###############################################################################
# Match/No Match -> Excel with Conditional Formatting (All Columns Colored)
###############################################################################

# ---- Optional: Set working directory to script location if in RStudio ----
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  if (rstudioapi::isAvailable()) {
    try({
      setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    }, silent = TRUE)
  }
}

# ---- Packages ----
# install.packages("openxlsx")  # Run once if you don't have it
library(openxlsx)

# ---- CONFIG: Change these paths if needed ----
methods_path <- "methods_comparison.txt"           # Tab-delimited methods file
export_paths <- list.files("./outputs_to_compare")                # One or more tab-delimited export files
out_path     <- "methods_comparison.updated.xlsx"  # Output Excel file

# ---- General options ----
options(stringsAsFactors = FALSE)

# ---- Helpers ----
normalize_text <- function(x) {
  x2 <- ifelse(is.na(x), "", x)
  x2 <- trimws(x2)
  x2 <- gsub("\\s+", " ", x2)    # collapse multiple spaces
  tolower(x2)
}

mk_key <- function(gene, cell) {
  paste(normalize_text(gene), normalize_text(cell), sep = "||")
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
  lab <- gsub("[^A-Za-z0-9]+", "_", base_no_ext)  # non-alphanum -> underscore
  lab <- gsub("^_+|_+$", "", lab)                 # trim boundary underscores
  if (lab == "") lab <- "export"
  lab
}

# Clean carriage returns and whitespace in character columns
clean_text_columns <- function(df) {
  for (i in seq_along(df)) {
    if (is.character(df[[i]])) {
      df[[i]] <- trimws(gsub("\r", "", df[[i]]))
    }
  }
  df
}

# ---- Read methods file ----
methods <- tryCatch(
  utils::read.delim(methods_path, check.names = FALSE),
  error = function(e) stop("Failed to read methods file: ", e$message)
)

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

# Build normalized keys for methods rows
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

# Apply formatting to rows 2..(n+1) to skip header, across ALL columns
for (col_idx in seq_len(n_cols)) {
  # Color cells that CONTAIN "Match" (case-sensitive; adjust if you need case-insensitive)
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

# Optional: Auto column widths
setColWidths(wb, "Results", cols = 1:n_cols, widths = "auto")

# Save the workbook
saveWorkbook(wb, out_path, overwrite = TRUE)
message(sprintf("Wrote Excel file with formatting to: %s", out_path))
```
3) And the output files from previous analyziz that looks like following INSIDE A SUBDIRECTORY IN CURRENT WORKING DIRECTORY NAMED "outputs_to_compare"

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

