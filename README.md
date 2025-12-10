# Different_ways_to_measure_cell_specific_expression
This page explains 3 different ways I used to build a tool to measure cell specific or enriched gene expression with bash/python and HTML

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
python cal_enrich_with_custom_clustor_no.py   --input-file rna_single_cell_cluster.tsv   --output-file my_custom_clustor_no_enrichment.tsv   --min-clusters 3
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
