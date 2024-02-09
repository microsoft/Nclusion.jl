import pandas as pd 
import matplotlib.pyplot as plt
import sys
import getopt
import numpy as np
import decimal
from itertools import permutations
from matplotlib.colors import Normalize
import numpy as np
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize

class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result


def get_pcts(data, purity_score_file, title,outfile,cell_type_ordering = None):#
    
    cell_types = data["cell_type"].unique()
    # print(np.sort(cell_types))
    # if cell_type_ordering is not None:
    #     print(np.sort(cell_type_ordering))
    #     print(np.sort(cell_types)==np.sort(cell_type_ordering))
    #     print(len(cell_types))
    #     print(len(cell_type_ordering))
    inferred_cell_types = data["inferred_label"].unique()
    pct_df_dict = {}
    purity_score_dict = {}
    for c in cell_types:
        df_cell_type = data[data["cell_type"] == c]
        total_cells = len(df_cell_type)
        inferred_dict = {}
        max_inferred = 0
        for i in inferred_cell_types:
            df_inferred_label = df_cell_type[df_cell_type["inferred_label"] == i]
            total_cells_inferred = len(df_inferred_label) 
            pct_inferred = total_cells_inferred/total_cells
            
            if total_cells_inferred > max_inferred:
                max_inferred = total_cells_inferred
                max_inferred_label = i
          
            inferred_dict[i] = pct_inferred
        
        pct_df_dict[c] = inferred_dict
        
        inferred_cluster_size = len(data[data["inferred_label"] == max_inferred_label])
        purity_score_dict[c] = max_inferred/inferred_cluster_size
        
      
    pct_df = pd.DataFrame.from_dict(pct_df_dict, orient='index', columns=inferred_cell_types)
    
    if cell_type_ordering is not None:
        pct_df = pct_df.reindex(index=cell_type_ordering)
    pct_df.to_csv(outfile)
    
    # pct_df = pct_df.reindex(index=['b_cells', 'cd56_nk', 'cd14_monocytes',  'cd34', 'regulatory_t', 'memory_t', 'naive_t', 'naive_cytotoxic','cytotoxic_t', 'cd4_t_helper'])
    
    
    # purity_score_df = pd.DataFrame.from_dict(purity_score_dict, orient='index', columns=['purity_score'])
    # purity_score_df = purity_score_df.reindex(index=['b_cells', 'cd56_nk', 'cd14_monocytes',  'cd34', 'regulatory_t', 'memory_t', 'naive_t', 'naive_cytotoxic', 'cytotoxic_t', 'cd4_t_helper'])
    # purity_score_df.to_csv(purity_score_file)

    return pct_df
    
def alt_purity_score(data, alt_purity_score_file):
    cell_types = data["cell_type"].unique()
    inferred_cell_types = data["inferred_label"].unique()
    
    purity_dict = {}
    for i in inferred_cell_types:
        inferred_df = data[data["inferred_label"] == i]
        max_cell_type = 0
        
        for c in cell_types:
            cell_type_df = inferred_df[inferred_df["cell_type"] == c]
            
            if len(cell_type_df) > max_cell_type:
                max_cell_type = len(cell_type_df)
        
        purity_dict[i] = max_cell_type
        
    sum_cluster_purity = sum(list(purity_dict.values()))
    total_cells = len(data)
    
    alt_purity_score = sum_cluster_purity/total_cells
    
    with open(alt_purity_score_file, "w") as f:
        f.write('Alternative Purity Score: \n')
        f.write(str(alt_purity_score))
        f.close()      
    
            
def optimize_diag(df, col_sum_threshold):  

    col_sum_dict = {}
    for col in df.columns.values:
        col_sum_dict[col] = df.loc[:, col].sum()
    
    sorted_items = sorted(col_sum_dict.items(), key=lambda x: x[1], reverse=True)
    sorted_items = sorted_items[:10] if len(sorted_items) >= 10 else sorted_items
    keep_cols = dict(sorted_items)
    # keep_cols = dict(filter(lambda elem: elem[1] > col_sum_threshold, col_sum_dict.items()))
    print(keep_cols)
    x = len(keep_cols)
    df2 = df.loc[:, keep_cols.keys()]
    
    perm = list(permutations(df2.columns.values, x))
        
    print('found all permutations')
    print(f"Columns Kept (# is {x}): {keep_cols}")


    diag_max = 0
    for p in perm:
        diag_sum = 0
        df_p = df2.loc[:, list(p)]
        df_p = df_p.reindex(columns=list(p))

        for i in range(x):
            diag_sum += df_p.iloc[i, i]
                
        if diag_sum > diag_max:
            diag_max = diag_sum
            optimal_cols = list(p)
    
    print('found optimal order of columns')
    df_optimal = df.loc[:, optimal_cols]
    df_optimal = df_optimal.reindex(columns=optimal_cols)

    if len(df_optimal.columns.values) != len(df.columns.values):
        print('inserting leftover columns')
        df_unplaced_cols = df.drop(columns=optimal_cols)

        max_val_dict = {}
        for col in df_unplaced_cols.columns.values:
            max_val_dict[col] = max(df_unplaced_cols.loc[:, col].values)
            sorted_dict = dict(sorted(max_val_dict.items(), key=lambda kv: kv[1], reverse=True))
            sorted_cols = sorted_dict.keys()
        for col in sorted_cols:
            insert_idx = len(df_optimal.columns.values)   
            df_optimal.insert(insert_idx, col, df_unplaced_cols[col].values)

    return df_optimal   
    
def make_heatmap(pct_table, save_path, title, colormap, display_details):

    pct_array = pct_table.to_numpy()
    fig, ax = plt.subplots()
    norm = MidPointNorm(midpoint=0.1)
    # im = ax.imshow(pct_array, cmap=colormap, norm=norm)
    im = ax.imshow(pct_array, cmap=colormap)
    
    inferred_groups = pct_table.columns.values
    experimental_labels = pct_table.index.values
    ax.set_xticks(np.arange(len(inferred_groups)), labels=np.arange(1, len(inferred_groups)+1), fontsize='medium')
    ax.set_xlabel('Inferred Cell Type Label', fontsize='x-large')
    
    if display_details == True:
        # cbar=plt.colorbar(im, orientation='horizontal', spacing='proportional')
        # cbar.set_ticks([0, 0.0000001,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
        # cbar.set_ticklabels(['0','0.0000001','0.1', '0.2','0.3','0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'])
        
        cbar = plt.colorbar(im, spacing = 'uniform')
        cbar.set_ticks([0 ,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        cbar.set_ticklabels(['0','0.1', '0.2','0.3','0.4', '0.5', '0.6', '0.7', '0.8', '0.9'])
        ax.set_yticks(np.arange(len(experimental_labels)), labels=experimental_labels, fontsize='x-large')
        ax.set_ylabel('Called Cell Type Label', fontsize='x-large')
    else:
        ax.set(ylabel=None)
        ax.tick_params(left=False, labelleft=False)

    ax.set_title(title, fontsize='xx-large')
    fig.tight_layout()
    plt.savefig(save_path)
    fig.clf()


def main(argv):
    input_path = ''
    save_path = 'pct_table_heatmap.png'
    title = ''
    colormap = 'Purples'
    cluster_dict_file = 'mapping_cluster_to_new_label.csv'
    ref_cluster_dict_file = 'mapping_refcluster_to_new_label.csv'
    purity_score_file = 'purity_score_table.csv'
    purity_score_barplot = "purity_score_barplot.png"
    pct_table_file = 'pct_table.csv'
    pct_df_notoptimal = 'pct_df_notoptimal.csv'
    alt_purity_score_file = 'alternative_purity_score.csv'
    col_sum_threshold = ""#0.275#0.15#0.3#
    display_details = False
    output_path_file = ''
    output_path_figure = ''
    output_suffix = ''
    
    dataset_ordering= ""
    

    



    try:
        opts, args = getopt.getopt(argv, '', ["input_path=", "save_path=", "title=", "colormap=", "cluster_dict_file=", "purity_score_file=", "pct_table_file=","pct_df_notoptimal=", "alt_purity_score_file=", "purity_score_barplot=", "col_sum_threshold=", "output_path_file=","output_path_figure=","output_suffix=", "display_details=","dataset_ordering_mode="])#
        print(opts)
                                
        for opt, arg in opts:
            if opt == '--input_path':
                input_path = arg
            elif opt == '--save_path':
                save_path = arg
            elif opt == '--title':
                title = arg
            elif opt == '--colormap':
                colormap = arg
            elif opt == '--cluster_dict_file':
                cluster_dict_file = arg
            elif opt == '--purity_score_file':
                purity_score_file = arg
            elif opt == 'pct_table_file':
                pct_table_file = arg
            elif opt == 'pct_df_notoptimal':
                pct_df_notoptimal = arg
            elif opt == '--alt_purity_score_file':
                alt_purity_score_file = arg
            elif opt == '--purity_score_barplot':
                purity_score_barplot = arg
            elif opt == '--col_sum_threshold':
                if not arg:
                    raise ValueError('empty arg: --col_sum_threshold')
                col_sum_threshold = float(arg)
            elif opt == '--output_path_file':
                output_path_file = arg
            elif opt == '--output_path_figure':
                output_path_figure = arg
            elif opt == '--output_suffix':
                output_suffix = arg
            elif opt == '--display_details':
                display_details = True
            elif opt == '--dataset_ordering_mode':
                dataset_ordering = arg

        
    except getopt.GetoptError: 
        sys.exit()
        print(opts) 

    
    print(output_path_file)
    data = pd.read_csv(input_path, index_col=1)
    print(data)
    

    # files_to_add_suffix = [save_path,cluster_dict_file,ref_cluster_dict_file,purity_score_file,purity_score_barplot,pct_table_file,pct_df_notoptimal,alt_purity_score_file]

    save_path = f"{save_path.split('.')[0]}{output_suffix}.{save_path.split('.')[1]}"
    cluster_dict_file = f"{cluster_dict_file.split('.')[0]}{output_suffix}.{cluster_dict_file.split('.')[1]}"
    ref_cluster_dict_file = f"{ref_cluster_dict_file.split('.')[0]}{output_suffix}.{ref_cluster_dict_file.split('.')[1]}"
    purity_score_file = f"{purity_score_file.split('.')[0]}{output_suffix}.{purity_score_file.split('.')[1]}"
    purity_score_barplot = f"{purity_score_barplot.split('.')[0]}{output_suffix}.{purity_score_barplot.split('.')[1]}"
    pct_table_file = f"{pct_table_file.split('.')[0]}{output_suffix}.{pct_table_file.split('.')[1]}"
    pct_df_notoptimal = f"{pct_df_notoptimal.split('.')[0]}{output_suffix}.{pct_df_notoptimal.split('.')[1]}"
    alt_purity_score_file = f"{alt_purity_score_file.split('.')[0]}{output_suffix}.{alt_purity_score_file.split('.')[1]}"

    label_file_name = input_path.split("/")[-1]
    cell_type_ordering = None
    if "pbmc" in dataset_ordering:
        cell_type_ordering = ['b_cells', 'cd56_nk', 'cd14_monocytes',  'cd34', 'regulatory_t', 'memory_t', 'naive_t', 'naive_cytotoxic','cytotoxic_t', 'cd4_t_helper']
    elif "galenAML" in dataset_ordering:
        cell_type_ordering = np.sort(data["cell_type"].unique())#["HSC","Prog","earlyEry","lateEry","GMP","ProMono","Mono","cDC","pDC","ProB","Plasma","T","CTL","B","NK","HSC−like","Prog−like","GMP−like","ProMono−like","Mono−like","cDC−like","MUTZ3_HSC−like","MUTZ3_Prog−like","MUTZ3_GMP−like","MUTZ3_ProMono−like","MUTZ3_Mono−like","MUTZ3_cDC−like","OCI_AML3_HSC−like","OCI_AML3_Prog−like","OCI_AML3_GMP−like","OCI_AML3_ProMono−like","OCI_AML3_Mono−like","OCI_AML3_cDC−like"]
    elif "pdac" in dataset_ordering:
        cell_type_ordering = np.sort(data["cell_type"].unique())#["Tumor_PANFR0383","Macrophage","T_NK","XCR1_DC","B_Cells","T_Regs","pDC_cell","DC","Endothelial","Plasma_cell","Mesenchymal","Hepatocyte","Tumor_PANFR0489","Tumor_PANFR0526","Tumor_PANFR0543","Tumor_PANFR0545","Tumor_PANFR0557","Tumor_PANFR0562","Tumor_PANFR0489R","Tumor_PANFR0575","Tumor_PANFR0576","Tumor_PANFR0580","Tumor_PANFR0583","Tumor_PANFR0588","Tumor_PANFR0593","Tumor_PANFR0598","Tumor_PANFR0552","Tumor_PANFR0605","Tumor_PANFR0592","Tumor_PANFR0578","Tumor_PANFR0631","Tumor_PANFR0635","Tumor_PANFR0637"]
    elif "tissueimmuneatlas" in dataset_ordering:
        cell_type_ordering = np.sort(data["cell_type"].unique())#["Alveolar macrophages","Tem/emra_CD8","Classical monocytes","Plasma cells","Memory B cells","NK_CD16+","Intestinal macrophages","Tgd_CRTAM+","Tnaive/CM_CD8","Tnaive/CM_CD4","Teffector/EM_CD4","Trm/em_CD8","Tfh","Trm_gut_CD8","MAIT","T_CD4/CD8","Tregs","Mast cells","Trm_Th1/Th17","Naive B cells","NK_CD56bright_CD16−","Trm_Tgd","MNP/B doublets","Cycling","Intermediate macrophages","Plasmablasts","ABCs","Cycling T&NK","DC2","Nonclassical monocytes","MNP/T doublets","ILC3","Erythrophagocytic macrophages","Tnaive/CM_CD4_activated","GC_B (I)","Progenitor","T/B doublets","Megakaryocytes","DC1","migDC","GC_B (II)","pDC","Erythroid","Pre−B"]
    else:
        cell_type_ordering = None

    outfile_pct_df_notoptimal=output_path_file+pct_df_notoptimal
    pct_table = get_pcts(data, purity_score_file, title,outfile_pct_df_notoptimal,cell_type_ordering=cell_type_ordering)
    # pct_table = get_pcts(data, purity_score_file, title,output_path_file)
    print(pct_table)

    outfile_alt=output_path_file+alt_purity_score_file
    print(outfile_alt)
    print(np.unique(data.inferred_label).shape[0])
    if np.unique(data.inferred_label).shape[0] < 50:
        alt_purity_score(data, outfile_alt)
        df_optimal = optimize_diag(pct_table, col_sum_threshold)
        # if "pbmc" in dataset_ordering:
        #     df_optimal = optimize_diag(pct_table, col_sum_threshold)
        # else:
        #     df_optimal = pct_table
        print(df_optimal.columns.values)
        pct_table = pct_table.reindex(columns=df_optimal.columns.values)
        
        
    print(np.unique(data.inferred_label).shape[0])
    cluster_dict = {}
    refcluster_dict = {}
    i = 1
    for col in pct_table.columns.values:
        cluster_dict[col] = i
        i += 1

    i = 1
    for col in pct_table.index.values:
        refcluster_dict[col] = i
        i += 1
    

    cluster_mapping = pd.DataFrame.from_dict(cluster_dict, orient='index')
    outfile_map=output_path_file+cluster_dict_file
    print(outfile_map)
    cluster_mapping.to_csv(outfile_map)  


    refcluster_mapping = pd.DataFrame.from_dict(refcluster_dict, orient='index')
    refoutfile_map=output_path_file+ref_cluster_dict_file
    print(refoutfile_map)
    refcluster_mapping.to_csv(refoutfile_map)  
    
    outfile_pct=output_path_file+pct_table_file
    pct_table = pct_table.rename(columns=cluster_dict)
    print(outfile_pct)
    pct_table.to_csv(outfile_pct) 

    tranlsated_labels = np.array([cluster_dict[el] for el in data.inferred_label.values])
    data["tranlsated_inferred_label"] = tranlsated_labels
    
    outfile_trans=output_path_file+dataset_ordering+'_NCLUSION_mapped_clusters.csv'
    print(outfile_pct)
    data.to_csv(outfile_trans)

    
                            
if __name__ == '__main__':
    main(sys.argv[1:])