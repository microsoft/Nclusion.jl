import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import scanpy as sc
import pandas as pd
import getopt
import sys

def make_tsnes(adata, res, figures_prefix):
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    
    adata.obs['inferred_label'] = pd.Categorical(res.loc[:, 'inferred_label_mapped'])
    adata.obs['called_label'] = pd.Categorical(res.loc[:, 'called_label'])
    adata.obs['cell_type'] = pd.Categorical(res.loc[:, "cell_type"])

    sc.pl.tsne(adata, color="cell_type", save=figures_prefix+'_called_cell_type.png')
    sc.pl.tsne(adata, color="inferred_label", save='_'+figures_prefix+'_inferred_label.png')
        


    clusters = res.loc[:, 'inferred_label'].unique()
    color_dict = {True: "#8567ad", False: 'lightgrey'}

    for c in clusters:
        obs_label = 'cluster'+str(c)
        adata.obs[obs_label] = pd.Categorical((adata.obs['inferred_label'] == c))
        sc.pl.tsne(adata, color=obs_label, save=figures_prefix+'_cluster_'+str(c)+'.png', palette=color_dict, title='Cluster '+str(c))



def main(argv):
    path_to_data = '' 
    path_to_labels = ''
    figures_prefix = ''
    translation=''
    
    try:                                                                     
        opts, args = getopt.getopt(argv, '', ["path_to_data=", "path_to_labels=", "figures_prefix="])
        
        for opt, arg in opts:
            if opt == '--path_to_data':
                path_to_data = arg
            elif opt == '--path_to_labels':
                path_to_labels = arg
            elif opt == '--figures_prefix':
                figures_prefix = arg

    except getopt.GetoptError: 
        sys.exit()
        
    adata = sc.read_h5ad(path_to_data)
    
    res = pd.read_csv(path_to_labels, index_col=1)
    
    if figures_prefix == 'vanGalen2019':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/vanGalen_mapping.csv', index_col=0).to_dict()
        cluster_dict = clusters['0']
    elif figures_prefix == 'dominguezconde2022':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/dominguezconde_mapping.csv', index_col=0).to_dict()
        cluster_dict = clusters['0']
    elif figures_prefix == 'raghavan2021':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/raghavan_mapping.csv')
        cluster_dict = clusters['0']
    elif figures_prefix == 'zheng2017':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/zheng_mapping.csv')
        cluster_dict = clusters['0']
    
    res['inferred_label_mapped'] = res.apply(lambda x: cluster_dict[int(x['inferred_label'])], axis=1)
    print(res)
    
    make_tsnes(adata, res, figures_prefix)

if __name__ == '__main__':
    main(sys.argv[1:])