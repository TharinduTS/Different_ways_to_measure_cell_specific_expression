# Different_ways_to_measure_cell_specific_expression
This page explains 3 different ways I used to build a tool to measure cell specific or enriched gene expression with bash/python and HTML

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

# METHOD 01- Simple Enrichment Scores 
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
# METHOD03- Enrichment score with more robust statistical backup systems to filter useful data

This method is much more complicated and takes much more computing power. Therefore I decided to split the input file into multiple files and process them simultaniously as a job array.

Just like previous efforts, this also has a python script that does the calculations for you

You can run it like following

enrich_array.sbatch
```
#!/bin/bash
#SBATCH --job-name=enrich_array
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-20%5
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
--B -Number of bootstrap tp be used

--reldelta

--rho

--min-k

--tissue-weighting

--pseudocount
