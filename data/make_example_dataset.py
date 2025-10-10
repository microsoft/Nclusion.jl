import pandas as pd
import scanpy as sc
import numpy as np
results_file = "data/write/pbmc3k.h5ad"

## Assuming you have done this already
#mkdir -p data
# cd data
# mkdir -p write
# test -f pbmc3k_filtered_gene_bc_matrices.tar.gz || curl https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -o pbmc3k_filtered_gene_bc_matrices.tar.gz
# tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# cd ..

adata = sc.read_10x_mtx(
    "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)  # this does nothing, in this specific case
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
adata = adata[
    (adata.obs.n_genes_by_counts < 2500)
    & (adata.obs.n_genes_by_counts > 200)
    & (adata.obs.pct_counts_mt < 5),
    :,
].copy()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sorted_genes = adata.var_names[np.argsort(-np.asarray(adata.X.todense().std(0)).ravel())]
adata = adata[:, sorted_genes[:500]].copy()
adata.X = adata.X.todense()
adata.X = adata.X.astype(np.float32)
adata.write(results_file)
