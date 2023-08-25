import scanpy as sc
import pandas as pd 
import numpy as np 
import anndata as ad 

# GENERAL COMMAND
# python preprocess_tutorial_data.py --path_to_data /path/to/data/ --path_to_save path/to/save --data_name dataname --n_hvgs int

def filter_hvgs(adata, n_hvgs):
    counts = np.transpose(adata.raw.X.astype(np.float64))
    n_umi = adata.obs.total_counts.to_numpy().astype(np.float64)
    scaling_factor = np.float64(1e4)
    psuedo =  np.float64(1)
    norm_X = scaling_factor * (counts / n_umi)
    lognorm_X = np.log(norm_X + psuedo , dtype=np.float64)
    gene_std_vec = np.std(lognorm_X, axis=1, dtype=np.float64)
    sorted_indx = np.argsort(-gene_std_vec)
    topGGenes = int(n_hvgs)
    top_genes = adata.var_names.values[sorted_indx][:topGGenes]
    top_genes_bool = [el in top_genes for  el in adata.var_names.values]
    hvg_name_vec = np.sort(adata.var_names.values[top_genes_bool])
    d = {'genes': hvg_name_vec}
    hvg_name_df = pd.DataFrame(data=d)
    hvg_name_df.to_csv(str(n_hvgs)+'_highly_variable_genes.csv')
    
    adata = adata[:, list(hvg_name_df['genes'].values)]
    raw_mat = ad.AnnData(adata.raw[:, list(hvg_name_df['genes'].values)].X)
    raw_mat.obs_names = adata.obs_names.values
    raw_mat.var_names = adata.var_names.values
    adata.raw = raw_mat

    return adata

def annotate_cell_lines(x):
    if 'MUTZ3' in x['orig_ident']:
        return 'MUTZ3_'+x['cell_type']
    elif 'OCI.AML3' in x['orig_ident']:
        return 'OCI_AML3_'+x['cell_type']
    else:
        return x['cell_type']

def preprocess(adata, data_name, n_hvgs):
    if data_name == 'pbmc':
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)  
        adata = adata[adata.obs['pct_counts_mt'] < 20, :]
        adata = adata[adata.obs['pct_counts_ribo'] > 5, :]
    
    if data_name == 'galen':
        adata = adata[adata.obs['state'] != 'unclear']
        df = adata.obs
        df['cell_type'] = df.apply(lambda x: annotate_cell_lines(x), axis=1)
        adata.obs['cell_type'] = df['cell_type']
        print('cell_types:', adata.obs.cell_type.unique())
          
    adata.raw = adata
    adata.layers["counts"] = adata.raw.X.copy()
    
    sc.pp.normalize_per_cell(adata, scale_factor)
    sc.pp.log1p(adata)
    
    if n_hvgs != None:
        adata = filter_hvgs(adata, n_hvgs)
    
    return adata

def main(argv):
    
    path_to_data = None
    path_to_save = ''
    data_name = None
    n_hvgs = None
    
    try:                                                                     
        opts, args = getopt.getopt(argv, '', ["path_to_data=", "path_to_save=", "data_name=", "n_hvgs="])
    except getopt.GetoptError: 
        sys.exit()
    
    for opt, arg in opts:                                                     
        if opt == '--path_to_data':                                              
            path_to_data = arg 
        elif opt == '--path_to_save':
            path_to_save = arg
        elif opt == '--data_name':
            data_name = arg
        elif opt == '--n_hvgs':
            n_hvgs = int(arg)
    
    if data_name == 'galen':
        path_to_data = "tutorial_data/galenAML_raw.h5ad"
    elif data_name == 'pbmc':
        path_to_data = "tutorial_data/pbmc_raw.h5ad"
    elif data_name == "tissueimmuneatlas"
        path_to_data = "tutorial_data/tissueimmuneatlas_ptd496_raw.h5ad"
    
    adata = sc.read_h5ad(path_to_data)
    
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    adata = preprocess(adata, data_name, n_hvgs)

    adata.write_h5ad(path_to_save)
    
if __name__ == "__main__":
    main(sys.argv[1:]) 
