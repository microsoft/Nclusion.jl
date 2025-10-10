import pandas as pd
import scanpy as sc
import anndata as ad 
import numpy as np

expression_matrix = pd.read_csv("simulated_data/simulated_data.csv")

adata = ad.AnnData(expression_matrix, dtype=np.float32)

condition = np.ones(shape=adata.n_obs)
adata.obs['condition'] = pd.Categorical(condition)

cell_type_labels = pd.read_csv("simulated_data/simulated_labels.csv")['z1']
cell_type_labels.index = cell_type_labels.index.map(str)
adata.obs["cell_type"] = pd.Categorical(cell_type_labels)

adata.write_h5ad('simulated_data/simulated_data.h5ad')

print('Expression Matrix: ', adata.X)
print('Cell IDs: ', adata.obs_names)
print('Gene Names: ', adata.var_names)
print('Condition Labels: ', adata.obs.condition)
print('Manually-curated Cell Type Labels: ', adata.obs.cell_type)
