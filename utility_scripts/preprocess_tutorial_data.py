import scanpy as sc
import pandas as pd 
import numpy as np 
import anndata as ad 
import sys
import getopt

# GENERAL COMMAND
# python preprocess_tutorial_data.py --path_to_data /path/to/data/ --path_to_save path/to/save --data_name dataname --n_hvgs int
def transform_tissue(path_to_data):
    pt_id="D496"
    adata = sc.read_h5ad(path_to_data)
    adata.X = adata.layers['counts']
    data = adata.X
    data = data.todense()
    adata.X = data.A
    adata = adata[adata.obs.Donor == pt_id]
    return adata
    
def concat_pbmc_data(path_to_data):
    cell_types = ["b_cells", "cd14_monocytes", "cd34", "cd4_t_helper", "cd56_nk", "cytotoxic_t", "memory_t", "naive_cytotoxic", "naive_t", "regulatory_t"]

    for i, c in enumerate(cell_types):
        
        cells_i = path_to_data+str(c)+"_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/"
        
        adata = sc.read_10x_mtx(cells_i, 
        var_names='gene_symbols',               
        cache=True) 
        data = adata.X
        data = data.todense()
        adata.X = data.A
        adata.var_names_make_unique()
        adata.obs["cell_type"] = c

        if i == 0:
            concat_adata = adata
        else:
            concat_adata = ad.concat([concat_adata, adata])

    concat_adata.obs_names_make_unique()
    return concat_adata

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
    
    adata = adata[:, list(hvg_name_df['genes'].values)]
    raw_mat = ad.AnnData(adata.raw[:, list(hvg_name_df['genes'].values)].X)
    raw_mat.obs_names = adata.obs_names.values
    raw_mat.var_names = adata.var_names.values
    adata.raw = raw_mat

    return adata

def annotate_raghavan(x):
    if x['cell.type.ID'] == 'Tumor':
        return x['Donor_ID']+'_'+x['cell.type.ID']
    else:
        return x['cell.type.ID']
    
def annotate_vangalen(x):
    if 'MUTZ3' in x['donor']:
        return 'MUTZ3_'+x['CellType']
    elif 'OCI.AML3' in x['donor']:
        return 'OCI_AML3_'+x['CellType']
    else:
        return x['CellType']

def preprocess(adata, data_name, n_hvgs, min_genes, min_cells, pct_counts_mito, pct_counts_ribo, scale_factor):
    if min_genes != None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if min_cells != None:
        sc.pp.filter_genes(adata, min_cells=min_cells)  
    if pct_counts_mito != None:
        adata = adata[adata.obs['pct_counts_mt'] < pct_counts_mito, :]
    if pct_counts_ribo != None:
        adata = adata[adata.obs['pct_counts_ribo'] > pct_counts_ribo, :]
          
    adata.raw = adata
    adata.layers["counts"] = adata.raw.X.copy()
    
    sc.pp.normalize_per_cell(adata, scale_factor)
    sc.pp.log1p(adata)
    
    if n_hvgs != None:
        adata = filter_hvgs(adata, n_hvgs)
    
    return adata

def main(argv):
    species_gene_symbols={'human': {'mt':('MT-', 'MT', 'Mt', 'mt'), 'ribo': ("RPS","RPL", "Rps", "Rpl"), 'hb':("^HB[^(P)]"), 'malat1' : ('MALAT1', 'Malat1', 'malat1')}, 'mouse': {'mt':('Mt-', 'mt-', 'mt', 'Mt'), 'ribo': ("Rps","Rpl", "rps", "rpl"), 'hb': "^hb[^(p)]|^Hb[^(p)]", 'malat1': ('MALAT1', 'Malat1', 'malat1')}}
    path_to_data = None
    path_to_save = ''
    data_name = None
    n_hvgs = None
    pct_counts_mito = None 
    pct_counts_ribo = None
    min_genes = None
    min_cells = None
    scale_factor = 1e4
    species = 'human'
    
    try:                                                                     
        opts, args = getopt.getopt(argv, '', ["path_to_data=", "path_to_save=", "data_name=", "n_hvgs=", "min_genes=", "min_cells=", "pct_counts_mito=", "pct_countd_ribo=", "scale_factor=", "species="])
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
        elif opt == '--min_genes':
            min_genes = int(arg)
        elif opt == '--min_cells':
            min_cells = int(arg)
        elif opt == '--pct_counts_mito':
            pct_counts_mito = int(arg)
        elif opt == '--pct_counts_ribo':
            pct_counts_ribo = int(arg)
        elif opt == '--scale_factor':
            scale_factor = int(arg)
        elif opt == '--species':
            species = arg
    
    if data_name == 'vanGalen':
        counts = pd.read_csv(path_to_data+'galenAML_countsdf.csv', index_col=0)
        metadata = pd.read_csv(path_to_data+'galenAML_metadf.csv', index_col=0)
        adata = ad.AnnData(counts)
        adata.obs['state'] = metadata.loc[:, 'PredictionRefined']
        adata.obs['donor'] = metadata.loc[:, 'orig.ident']
        adata = adata[adata.obs['state'] != 'unclear']
        adata.obs['cell_type'] = metadata.apply(lambda x: annotate_vangalen(x), axis=1)
        
    elif data_name == 'zheng':
       adata = concat_pbmc_data(path_to_data)
       min_genes = 200
       min_cells = 3
       pct_counts_mito = 20
       pct_counts_ribo = 5
       
    elif data_name == "dominguezconde":
        adata = transform_tissue(path_to_data)
        
    elif data_name == "raghavan":
        counts = pd.read_csv(path_to_data+'Biopsy_RawDGE_23042cells.csv', index_col=0)
        metadata = pd.read_cvs(path_to_data+"Biopsy_Metadata_23042cells.csv", index_col=0)
        adata = ad.AnnData(counts.T)
        adata.obs['cell_type'] = metadata.apply(lambda x: annotate_raghavan(x), axis=1)
        adata.obs['donor'] = metadata.loc[:, "Donor_ID"]
        
    else:
        adata = sc.read_h5ad(path_to_data)
    
    condition = np.ones(adata.n_obs)
    adata.obs['condition'] = pd.Categorical(condition)
    
    qc_vars = []
    if species == 'human' or species == 'mouse':
        mt_symbol=species_gene_symbols[species]['mt']
        ribo_symbol=species_gene_symbols[species]['ribo']

        adata.var['mt'] = adata.var_names.str.startswith(mt_symbol) 
        qc_vars.append('mt')

        adata.var['ribo'] = adata.var_names.str.startswith(ribo_symbol)
        qc_vars.append('ribo')
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)
    
    adata = preprocess(adata, data_name, n_hvgs, min_genes, min_cells, pct_counts_mito, pct_counts_ribo, scale_factor)

    adata.write_h5ad(path_to_save)
    
if __name__ == "__main__":
    main(sys.argv[1:]) 
