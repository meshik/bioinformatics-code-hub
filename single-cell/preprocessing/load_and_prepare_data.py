import scanpy as sc
import pandas as pd
import scipy.sparse as sp
from pathlib import Path


def adult_human_heart(
    expression_matrix: str | Path,
    metadata_path: str | Path,
    min_cells: int = 3,
    min_genes: int = 200,
    n_hvg: int = 5_000,
    resolution: float = 0.5,
) -> sc.AnnData:
    """
    Parameters
    ----------
    expression_matrix
        CSV with genes as rows and cells as columns.  
        The first column (“...1”) must contain gene IDs.
    metadata_path
        TSV with per-cell metadata; index column “ID”.
    min_cells, min_genes
        Standard Seurat filtering thresholds.
    n_hvg
        Number of highly variable genes to keep.
    resolution
        Leiden clustering resolution (≈ Seurat’s 0.5).

    Returns
    -------
    AnnData
        Fully processed (clusters in `adata.obs['seurat_clusters']`,
        UMAP in `adata.obsm['X_umap']`).
    """
    # ---------- Load expression ----------
    expr_df = pd.read_csv(expression_matrix)
    expr_df = expr_df.set_index(expr_df.columns[0])       # “…1” → index
    X = sp.csr_matrix(expr_df.values)                     # sparse counts
    var = pd.DataFrame(index=expr_df.index)

    # ---------- Load metadata ----------
    obs = pd.read_csv(metadata_path, sep="\t").set_index("ID")

    # ---------- Build AnnData ----------
    adata = sc.AnnData(X=X, var=var, obs=obs)

    # ---------- QC & filtering ----------
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # % mitochondrial reads
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # ---------- Normalise & HVGs ----------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvg, flavor="seurat_v3", subset=True
    )

    # ---------- Scale & PCA ----------
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

    # ---------- Neighbours, clustering, UMAP ----------
    sc.pp.neighbors(adata, n_pcs=10)
    sc.tl.leiden(adata, resolution=resolution, key_added="seurat_clusters")
    sc.tl.umap(adata)

    return adata


# test
if __name__ == "__main__":  # pragma: no cover
    print(Path.cwd())
    adata = adult_human_heart(
        expression_matrix=Path("single-cell/data/GSE109816_normal_heart_umi_matrix.csv.gz"),
        metadata_path=Path("single-cell/data/GSE109816_normal_heart_cell_cluster_info.txt")
    )
    print(adata)
    


    