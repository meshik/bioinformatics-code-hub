# scripts/adult_human_heart.py
from __future__ import annotations
import scanpy as sc
import pandas as pd
import scipy.sparse as sp
from pathlib import Path


def adult_human_heart(
    expression_matrix: str | Path,
    metadata_path: str | Path,
    *,
    min_cells: int = 3,
    min_genes: int = 200,
    n_hvg: int = 5000,
    resolution: float = 0.5,
) -> sc.AnnData:
    """
    Load and preprocess the adult human heart dataset.
    Matches the R Seurat pipeline.
    """
    # ---------- load counts ----------
    expr_df = pd.read_csv(expression_matrix, compression="infer", index_col=0)

    # ---------- load metadata ----------
    metadata = pd.read_csv(metadata_path, sep="\t", index_col="ID")

    # ---------- keep only cells that exist in both ----------
    shared_cells = metadata.index.intersection(expr_df.columns)
    expr_df = expr_df[shared_cells]
    metadata = metadata.loc[shared_cells]

    # ---------- build AnnData ----------
    X = sp.csr_matrix(expr_df.T.values)  # transpose: cells as rows
    adata = sc.AnnData(X=X, obs=metadata, var=pd.DataFrame(index=expr_df.index))

    # ---------- QC ----------
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # ---------- normalization & features ----------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvg, flavor="seurat_v3", subset=True
    )
    adata.raw = adata  # stores log-normalized, HVG-filtered data
    # ---------- scale & dimensionality reduction ----------
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_pcs=10)
    sc.tl.leiden(adata, resolution=resolution, key_added="seurat_clusters")
    sc.tl.umap(adata)

    return adata

# test
if __name__ == "__main__":
    print(Path.cwd())
    adata = adult_human_heart(
        expression_matrix=Path("single-cell/data/GSE109816_normal_heart_umi_matrix.csv.gz"),
        metadata_path=Path("single-cell/data/GSE109816_normal_heart_cell_cluster_info.txt")
    )
    print(adata)
    


    