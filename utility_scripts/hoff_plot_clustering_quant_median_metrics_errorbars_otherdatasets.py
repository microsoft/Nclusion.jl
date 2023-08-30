import argparse
import numpy as np
import pandas as pd
import anndata as ad
import sys
import os
import matplotlib.pyplot as plt
from IPython.core import ultratb
import itertools
import matplotlib as mpl
sys.excepthook = ultratb.FormattedTB(color_scheme='Linux', call_pdb=False)

from plotnine import (
    ggplot,
    aes,
    geom_col,
    geom_text,
    position_dodge,
    lims,
    theme,
    element_text,
    element_blank,
    element_rect,
    element_line,
    after_scale,
    geom_bar,
    theme_classic,
		scale_fill_manual
)

def parse_args():
    parser=argparse.ArgumentParser(description="Make the metrics median plot for other datasets")
    parser.add_argument("full_data_path")
    parser.add_argument("dataset", default="")
    parser.add_argument("save_file", default="")
    args=parser.parse_args()
    return args

# def basal_condition(df, i):
#     return df[i,:max_score] == df[i,:sc_Basal1] && df[i,:max_score] >= 0.25

def main():

    inputs=parse_args()
    metricsfilename = inputs.full_data_path
    save_file = inputs.save_file
    dataset = inputs.dataset
    if save_file == "":
       save_file_path= "/".join(metricsfilename.split('/')[:-1])+"/"
       filename_list = (metricsfilename.split('/')[-1]).split(".")
       filename_list[-1] = "MedianDiffPlots"
       save_file_namebase="_".join(filename_list)
       save_file=save_file_path+save_file_namebase
    if dataset == "":
        dataset=None

    metricsdf = pd.read_csv(metricsfilename).drop_duplicates()
    metricsdf.method = metricsdf.method.str.split("_").str.join(" ")
    # mean_vec = Float64.(metrics_df.mean)
    # metricsdf.groupby('dataset')['mean'].median()
    group_medians = metricsdf.groupby('dataset')['mean'].transform('mean')
    metricsdf['median_diff'] = metricsdf['mean'] - group_medians
    metricsdf['upper_median_diff'] = metricsdf['upper'] - group_medians
    metricsdf['lower_median_diff'] = metricsdf['lower'] - group_medians
    # speedorderdict = {'seurat nonorm':0, 'NCLUSION':1, 'leiden knn':2, 'scLCA':3, 'scCCESS SIMLR':4, 'SOUP':5}
    # speedorderindex = metricsdf['method'].apply(lambda x: speedorderdict[x])
    # metricsdf.groupby('dataset')['mean'].set_index(speedorderindex.index).sort()
    speedordercategories = ["NCLUSION", "seurat nonorm", "leiden knn", "scLCA", "scCCESS SIMLR", "SOUP"]
    # df["month"] = pd.Categorical(metricsdf.groupby('dataset')["method"], categories = speedordercategories)
    # df.sort_values(by = "month")

    # metricsdf["method"] =  metricsdf["method"].astype("category").cat.reorder_categories(speedordercategories)
    metricsdf["method"] = pd.Categorical(metricsdf["method"], categories=speedordercategories,ordered=True)
    metricsdf.sort_values("method", inplace=True)
    colorsdict={'seurat nonorm':"#BE9C47", 'NCLUSION':"#003462", 'leiden knn':"#71989F", 'scLCA':"#779176", 'scCCESS SIMLR':"#9482A9", 'SOUP':"#CB89B2"}

    allcolorscatlist = list(colorsdict.values())
    allcolorscatlist.append("#999999")
    metricsdf["hexcolors"] =pd.Categorical(metricsdf["method"].apply(lambda x: colorsdict[x]), categories=allcolorscatlist,ordered=False)#metricsdf["method"].apply(lambda x: colorsdict[x]) 

    metricsdf.loc[metricsdf.median_diff <= 0.0,"hexcolors"] = '#999999'

    ss = metricsdf.groupby('dataset')[['method','hexcolors']]
    reindexecolors=list(itertools.chain(*[list(el) for el in ss.indices.values()]))
    metricsdf.hexcolors[reindexecolors]
    # metricsdf.plot(x='dataset',
    #     kind='bar',
    #     stacked=False,
    #     title='Grouped Bar Graph with dataframe')

    # (ggplot(metricsdf, aes(x='dataset', y='median_diff',fill='method',width=.9))
    #     + geom_col(aes(fill='factor(method)'),stat='identity', position='dodge', size=1)
    #     + (scale_fill_manual(values=metricsdf.hexcolors))#allcolorscatlist
    #     + theme_classic())
    metricsdf.pivot(index='dataset', columns='method', values='median_diff').plot(kind='bar', color=allcolorscatlist,legend=None)

    ax = metricsdf.pivot(index='dataset', columns='method', values='median_diff').plot(kind='bar', color=allcolorscatlist,legend=None,rot=0,ylabel='Metric Diff',xlabel='dataset')
    rects = [c for c in ax.get_children() if isinstance(c, mpl.patches.Rectangle)]
    for rect, color in zip(rects, metricsdf.hexcolors): #[reindexecolors]
        rect.set_fc(color)

    # ax.grid(which='minor', axis='x', color='#999999', linestyle='--')
    # ax.grid(which='minor', axis='x', linestyle='--')
    # ax.hlines(y=0.2, xmin=4, xmax=20, linewidth=2, color='r')
    # ax.hlines(y=0.0, linewidth=2, color='r')

    # ax.annotate('Galen 2019',xy=(-0.25, 0.275), xytext=(-0.25, 0.275))
    # ax.annotate('PDAC Biopsy 2021',xy=(0.75, 0.275), xytext=(0.75, 0.275))
    # ax.annotate('Tissue Immune \n (pt D496) 2022',xy=(1.55, 0.275), xytext=(1.55, 0.275))
    ax.set_ylim(-0.35, 0.35)
    # Remove the spines (axes borders)
    ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    ax.set_facecolor('white')
    # Remove ticks and labels
    # ax.tick_params()
    ax.tick_params(bottom=False, labelbottom=False, top=False, right=False, labeltop=False, labelright=False)
    # Show only the y-axis line
    ax.spines['left'].set_linewidth(0.5)# Adjust the line width as desired
    ax.spines['left'].set_color('black')  # Adjust the color as desired 

    num_bars = len(metricsdf["dataset"])
    num_groups = len(pd.unique(metricsdf["dataset"]))
    label_map={"galenAML":'Galen 2019',"pdac_biopsy":'PDAC Biopsy 2021',"tissueimmuneatlas":'Tissue Immune \n (pt D496) 2022'}
    for i, k in enumerate(pd.unique(metricsdf["dataset"])):
        label = label_map[k]
        x_pos = (num_bars / 2) + (i * num_bars)
        ax.text(i, ax.get_ylim()[1], label, ha='center')
    

    ax.annotate('Better Performance', xy=(num_groups-0.5, ax.get_ylim()[1]), xycoords='data',
                xytext=(num_groups+0.5, 0), textcoords='data',
                arrowprops=dict(arrowstyle="->", lw=1.5, color='black'),
                ha='center', va='bottom', color='black', rotation='vertical')
    plt.axhline(y = 0.0, color = 'k', linestyle = '-', linewidth=0.5)
    # plt.axvline(x = -0.490, color = 'k', linestyle = '-', linewidth=0.75)
    plt.axvline(x = 0.5, color = 'gray', linestyle = '--')
    plt.axvline(x = 1.5, color = 'gray', linestyle = '--')
    plt.show()
    
        # + (scale_fill_manual(values=allcolorscatlist)) , color='factor(hexcolors)' + =metricsdf['hexcolors'].astype(str)
    # pd.Categorical(metricsdf["method"].apply(lambda x: colorsdict[x]), categories=colorsdict.values(),ordered=False)

    [4, 0, 3, 2, 1, 5]
    print(save_file)
    # print(dataset)
    if dataset == None:
        dataset="no dataset"
        print(dataset)
    else:
        metrics_df = metricsdf[metricsdf.dataset == dataset_name]
        print(dataset)
    # metadf = pd.read_csv(data_path,index_col=0)
    # new_metadf = metadf.loc[adata.obs_names]
    # ids = new_metadf.index
    # tumor_columns=new_metadf[["Hybrid_Specific1", "sc_Classical1", "sc_Basal1"]]
    # new_metadf['max_score'] = tumor_columns.max(axis=1)
    # new_metadf['cell_labels'] = "Organoid" 
    # new_metadf.loc[np.logical_and(new_metadf['max_score'] == new_metadf["sc_Basal1"], new_metadf['max_score'] >= 0.25), 'cell_labels'] = "Basal" 
    # new_metadf.loc[np.logical_and(new_metadf['max_score'] == new_metadf["sc_Classical1"], new_metadf['max_score'] >= 0.0),'cell_labels'] = "Classical" 
    # new_metadf.loc[np.logical_and(new_metadf['max_score'] == new_metadf["Hybrid_Specific1"], new_metadf['max_score'] >= 0.3),'cell_labels'] = "Hybrid" 
    # new_andata = adata[new_metadf.index.values]
    # new_andata.obs['sample_origin'] = new_andata.obs['cell_type']
    # new_andata.obs['cell_labels'] = new_metadf.cell_labels
    # new_andata.obs['max_score'] = new_metadf.max_score
    # new_andata.obs['sc_Basal1'] = new_metadf.sc_Basal1
    # new_andata.obs['sc_Classical1'] = new_metadf.sc_Classical1
    # new_andata.obs['Hybrid_Specific1'] = new_metadf.Hybrid_Specific1
    # new_andata.obs['cell_type'] = new_metadf.cell_labels
    # new_andata.write_h5ad(save_file)

if __name__ == "__main__":
    main()
