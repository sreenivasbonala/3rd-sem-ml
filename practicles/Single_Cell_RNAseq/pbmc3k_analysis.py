import scanpy as sc

# Load PBMC3k dataset
adata = sc.read_10x_mtx(
    "data/filtered_gene_bc_matrices/hg19",
    var_names="gene_symbols",
    cache=False
)

print("PBMC3k dataset loaded successfully!")
print(adata)

print("\nNumber of cells:", adata.n_obs)
print("Number of genes:", adata.n_vars)

# -----------------------------
# Quality Control
# -----------------------------

# Identify mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)

print("\nQC metrics calculated.")

# Display basic QC statistics
print("\nQC statistics:")
print(adata.obs[
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
].describe())

# Filter low-quality cells
adata = adata[
    (adata.obs["n_genes_by_counts"] >= 200) &
    (adata.obs["pct_counts_mt"] < 5)
].copy()

# Filter genes expressed in at least 3 cells
sc.pp.filter_genes(adata, min_cells=3)

print("\nAfter QC filtering:")
print(adata)

print("\nCells remaining:", adata.n_obs)
print("Genes remaining:", adata.n_vars)
# -----------------------------
# Normalization
# -----------------------------

# Save raw counts before normalization
adata.layers["counts"] = adata.X.copy()

# Normalize each cell to 10,000 total counts
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the normalized counts
sc.pp.log1p(adata)

print("\nNormalization completed.")
print("Data normalized to 10,000 counts per cell and log-transformed.")
# -----------------------------
# Highly Variable Gene Selection
# -----------------------------

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat"
)

print("\nHighly variable genes identified.")
print("Number of highly variable genes:",
      adata.var["highly_variable"].sum())

# Display top 10 highly variable genes
hvg_names = adata.var_names[adata.var["highly_variable"]][:10]

print("\nTop highly variable genes:")
for gene in hvg_names:
    print(gene)
# -----------------------------
# Scaling and PCA
# -----------------------------

# Keep only highly variable genes
adata = adata[:, adata.var["highly_variable"]].copy()

# Scale each gene
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata, n_comps=30, svd_solver="arpack")

print("\nPCA completed.")
print("PCA shape:", adata.obsm["X_pca"].shape)

# Show variance explained by first 10 PCs
print("\nVariance explained by first 10 PCs:")
print(adata.uns["pca"]["variance_ratio"][:10])
# -----------------------------
# Nearest Neighbor Graph
# -----------------------------

sc.pp.neighbors(
    adata,
    n_neighbors=10,
    n_pcs=20
)

print("\nNearest-neighbor graph calculated.")
print("Neighbors:", 10)
print("PCs used:", 20)
# -----------------------------
# UMAP Visualization
# -----------------------------

sc.tl.umap(adata, random_state=42)

print("\nUMAP completed.")
print("UMAP shape:", adata.obsm["X_umap"].shape)

# Save UMAP coordinates
umap_df = adata.obsm["X_umap"]

print("\nFirst 5 UMAP coordinates:")
print(umap_df[:5])
# ============================================================
# STEP 9 — Leiden Clustering
# ============================================================

sc.tl.leiden(
    adata,
    resolution=0.5,
    random_state=42
)

print("\nStep 9: Leiden clustering completed.")

print("\nCells per cluster:")
print(adata.obs["leiden"].value_counts().sort_index())


# ============================================================
# STEP 10 — Marker Gene Discovery
# ============================================================

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
)

print("\nStep 10: Marker gene discovery completed.")

# Display top 10 marker genes for each cluster
print("\nTop marker genes for each cluster:")

for cluster in adata.obs["leiden"].cat.categories:
    genes = adata.uns["rank_genes_groups"]["names"][cluster][:10]

    print(f"\nCluster {cluster}:")
    for gene in genes:
        print(gene)


# Save marker genes to CSV
marker_df = sc.get.rank_genes_groups_df(
    adata,
    group=None
)

marker_df.to_csv(
    "marker_genes.csv",
    index=False
)

print("\nMarker genes saved as: marker_genes.csv")


# ============================================================
# STEP 11 — Cell-Type Annotation
# ============================================================

# Known PBMC marker genes
cell_type_markers = {
    "T cells": ["IL7R", "LTB", "IL32", "LST1"],
    "B cells": ["CD79A", "MS4A1", "CD37", "CD79B"],
    "NK cells": ["NKG7", "GNLY", "GZMB", "PRF1"],
    "Monocytes": ["LYZ", "S100A8", "S100A9", "CTSS"],
    "Dendritic cells": ["FCER1A", "CST3", "CLEC10A"],
    "Megakaryocytes": ["PPBP", "PF4"]
}

print("\nStep 11: Cell-type marker reference created.")

print("\nCell-type marker genes:")

for cell_type, markers in cell_type_markers.items():
    print(f"{cell_type}: {', '.join(markers)}")


# ============================================================
# STEP 12 — UMAP Cluster Visualization
# ============================================================

sc.pl.umap(
    adata,
    color="leiden",
    title="PBMC3k Leiden Clusters",
    save="_clusters.png",
    show=False
)

print("\nStep 12: UMAP cluster plot saved.")


# UMAP with selected marker genes
marker_genes_plot = [
    "IL7R",
    "MS4A1",
    "NKG7",
    "LYZ",
    "FCER1A",
    "PPBP"
]

available_markers = [
    gene for gene in marker_genes_plot
    if gene in adata.var_names
]

sc.pl.umap(
    adata,
    color=available_markers,
    title="PBMC3k Marker Gene Expression",
    save="_markers.png",
    show=False
)

print("Marker gene UMAP plot saved.")
