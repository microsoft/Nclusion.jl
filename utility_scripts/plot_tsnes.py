import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import scanpy as sc
import pandas as pd

def make_tsnes(adata, res, figures_prefix):

    sc.tl.tsne(adata)
    
    adata.obs['inferred_label'] = pd.Categorical(res.loc[:, 'inferred_label'])
    adata.obs['called_label'] = pd.Categorical(res.loc[:, 'called_label'])
    adata.obs['cell_type'] = pd.Categorical(res.loc[:, "cell_type"])

    sc.pl.tsne(adata, color="cell_type", save=figures_prefix+'_celltype.png')
    sc.pl.tsne(adata, color="inferred_label", save='_'+model+'_'+dataset+'_'+hvgs+'_inferredlabel.png')
        


    clusters = res.loc[:, 'inferred_label'].unique()
    color_dict = {True: "#8567ad", False: 'lightgrey'}
    print(clusters)
    for c in clusters:
        obs_label = 'cluster'+str(c)
        adata.obs[obs_label] = pd.Categorical((adata.obs['inferred_label'] == c))
        sc.pl.tsne(adata, color=obs_label, save=figures_prefix+'_cluster_'+c+'.png', palette=color_dict, title='Cluster '+c)



def main(argv):
    path_to_data = '' 
    path_to_labels = ''
    figures_prefix = ''
    
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
    
    make_tsnes(adata, res, figures_prefix)

if __name__ == '__main__':
    main(sys.argv[1:])