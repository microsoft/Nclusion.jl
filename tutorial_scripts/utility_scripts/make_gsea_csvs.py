import pandas as pd 
import gseapy as gp
import matplotlib.pyplot as plt
import sys
import getopt
import os
import numpy as np
import decimal
from itertools import permutations
from matplotlib.colors import Normalize
import numpy as np
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize
import matplotlib.ticker as ticker
import time
from textwrap import wrap
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering
from itertools import compress

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def setup_axes(width=20,height=4):
    fig,ax = plt.subplots(1,1,figsize=(width,height),dpi=300)
    return fig,ax

def effect_size_calc(popmean,notpopmean,popvar,notpopvar):
    numgenes=len(notpopmean)
    effect_size = (popmean-notpopmean) / np.sqrt(0.5 * (popvar + notpopvar))
    effect_size_sign = np.array(["0" for _ in range(numgenes)])
    effect_size_sign[effect_size >0.0] = "+"
    effect_size_sign[effect_size <0.0] = "-"
    return effect_size, effect_size_sign

def main(argv):
    path_to_pips = ''
    path_to_labels = ''
    path_to_translation = ''
    path_to_data = ''
    path_to_nk = ''
    output_path = ''
    
    try:
        opts, args = getopt.getopt(argv, '', ["path_to_pips=", "path_to_labels=", "path_to_translation=",  "path_to_nk=", "output_path=","path_to_data="])
        print(opts)
                                
        for opt, arg in opts:
            if opt == '--path_to_pips':
                path_to_pips = arg
            if opt == '--path_to_labels':
                path_to_labels = arg
            if opt == '--path_to_translation':
                path_to_translation = arg
            elif opt == '--path_to_nk':
                path_to_nk = arg
            elif opt == '--output_path':
                output_path = arg
            if opt == '--path_to_data':
                path_to_data = arg
 
    except getopt.GetoptError: 
        sys.exit()
        print(opts)   
    print(output_path)
    colorscheme=["#D70000", "#8C3CFF", "#028800", "#00ACC7", "#E7A500", "#FF7FD1", "#6C004F", "#583B00", "#005759", "#15E18D", "#0000DD", "#A2766A", "#BCB7FF", "#C004B9", "#645473", "#790000", "#0774D8", "#739B7D", "#FF7852", "#004B00", "#8F7B01", "#F3007B", "#8FBA00", "#A67BB8", "#5A02A3", "#E3AFAF", "#A03A52", "#A2C8C8", "#9E4B00", "#546745", "#BBC389", "#5F7B88", "#60383C", "#8388FF", "#390000", "#E353FF", "#305382", "#7FCAFF", "#C5668F", "#00816A", "#929EB7", "#CC7407", "#7F2B8E", "#00BEA4", "#2DB152", "#4E33FF", "#00E500", "#FF00CE", "#C85848", "#E59CFF", "#1DA1FF", "#6E70AB", "#C89A69", "#78573B", "#04DAE6", "#C1A3C4", "#FF6A8A", "#BB00FE", "#925380", "#9F0274", "#94A150", "#374425", "#AF6DFF", "#596D00", "#FF3147", "#838057", "#006D2E", "#8956AF", "#5A4AA3", "#773516", "#86C39A", "#5F1123", "#D58581", "#A42918", "#0088B1", "#CB0044", "#FFA056", "#EB4E00", "#6C9700", "#538649", "#755A00", "#C8C440", "#92D370", "#4B9894", "#4D230D", "#61345C", "#8400CF", "#8B0031", "#9F6E32", "#AC8499", "#C63189", "#025438", "#086B84", "#87A8EC", "#6466EF", "#C45DBA", "#019F70", "#815159", "#836F8C", "#B3C0DA", "#B99129", "#FF97B2", "#A793E1", "#698DBE", "#4C5001", "#4802CC", "#61006E", "#456A66", "#9D5743", "#7BACB5", "#CD84BD", "#0054C1", "#7B2F4F", "#FB7C00", "#34C000", "#FF9C88", "#E1B769", "#536177", "#5C3A7C", "#EDA5DA", "#F053A3", "#5D7E69", "#C47750", "#D14868", "#6E00EB", "#1F3400", "#C14104", "#6DD5C2", "#46709F", "#A201C4", "#0A8289", "#AFA601", "#A65C6B", "#FE77FF", "#8B85AE", "#C77FE9", "#9AAB85", "#876CD9", "#01BAF7", "#AF5ED2", "#59512B", "#B6005F", "#7CB66A", "#4985FF", "#00C282", "#D295AB", "#A34BA8", "#E306E3", "#16A300", "#392E00", "#843033", "#5E95AA", "#5A1000", "#7B4600", "#6F6F31", "#335826", "#4D60B6", "#A29564", "#624028", "#45D458", "#70AAD0", "#2E6B4E", "#73AF9E", "#FD1500", "#D8B492", "#7A893B", "#7DC6D9", "#DC9137", "#EC615E", "#EC5FD4", "#E57BA7", "#A66C98", "#009744", "#BA5F22", "#BCAD53", "#88D830", "#873573", "#AEA8D2", "#E38C63", "#D1B1EC", "#37429F", "#3ABEC2", "#669D4D", "#9E0399", "#4E4E7A", "#7B4C86", "#C33531", "#8D6677", "#AA002D", "#7F0175", "#01824D", "#734A67", "#727791", "#6E0099", "#A0BA52", "#E16E31", "#C56A71", "#6D5B96", "#A33C74", "#326200", "#880050", "#335869", "#BA8D7C", "#1959FF", "#919202", "#2C8BD5", "#1726FF", "#21D3FF", "#A490AF", "#8B6D4F", "#5E213E", "#DC03B3", "#6F57CA", "#652821", "#AD7700", "#A3BFF7", "#B58446", "#9738DC", "#B25194", "#7242A3", "#878FD1", "#8A70B1", "#6BAF36", "#5A7AC9", "#C79FFF", "#56841A", "#00D6A7", "#824739", "#11431D", "#5AAB75", "#915B01", "#F64570", "#FF9703", "#E14231", "#BA92CF", "#34584D", "#F8807D", "#913400", "#B3CD00", "#2E9FD3", "#798B9F", "#51817D", "#C136D7", "#EC0553", "#B9AC7E", "#487032", "#849565", "#D99D89", "#0064A3", "#4C9078", "#8F6198", "#FF5338", "#A7423B", "#006E70", "#98843E", "#DCB0C8"]
    godotplotspath = output_path + "go_dotplots"
    violinplotspath = output_path + "violin_dotplots"
    gocsvspath = output_path + "go_csvs"

    # Check if the folder already exists
    if not os.path.exists(godotplotspath):
        # Create the new folder
        os.makedirs(godotplotspath)
        print(f"Folder '{godotplotspath}' created successfully.")
    else:
        print(f"Folder '{godotplotspath}' already exists.")

    if not os.path.exists(violinplotspath):
        # Create the new folder
        os.makedirs(violinplotspath)
        print(f"Folder '{violinplotspath}' created successfully.")
    else:
        print(f"Folder '{violinplotspath}' already exists.")

    if not os.path.exists(gocsvspath):
        # Create the new folder
        os.makedirs(gocsvspath)
        print(f"Folder '{gocsvspath}' created successfully.")
    else:
        print(f"Folder '{gocsvspath}' already exists.")

    data = pd.read_csv(path_to_pips)
    labelsdf = pd.read_csv(path_to_labels)
    translationdict = pd.read_csv(path_to_translation)
    adata = sc.read_h5ad(path_to_data)
    data.set_index('cluster_id', inplace=True)
    nk = pd.read_csv(path_to_nk)
    nk.set_index('cluster_ids', inplace=True)
    Kmax_occupied = np.sum(nk.Nk >0)
    forwardtranslationdict={translationdict.loc[id].iloc[0]:translationdict.loc[id].iloc[1] for id in translationdict.index}
    clusterforwardtranslationdict={"Cluster_"+str(translationdict.loc[id].iloc[0]):"Cluster_"+str(translationdict.loc[id].iloc[1]) for id in translationdict.index}
    reversetranslationdict={translationdict.loc[id].iloc[1]:translationdict.loc[id].iloc[0] for id in translationdict.index}
    clsuterreversetranslationdict={"Cluster_"+str(translationdict.loc[id].iloc[1]):"Cluster_"+str(translationdict.loc[id].iloc[0]) for id in translationdict.index}
    colorsdict={f"Cluster_{k+1}":v  for (k,v) in enumerate(colorscheme)}
    Kmax_all = np.sum(nk.Nk >=0)
    thresh=0.5
    adjustment_weights = (1 - np.sum(data > thresh,axis=0)/Kmax_occupied)
    adj_pips = data * adjustment_weights
    genenames = adj_pips.columns
    top_genes_bool = [ j in genenames for j in adata.raw.var_names]
    col = (adj_pips >= 0.5).any()
    onlysigmat=adj_pips.loc[: , col].transpose().loc[:,nk.Nk >0]
    labelsdict={}
    fdrdict={}
    xlimmax=0
    modulegenesdict={}
    moduleposgenesdict={}
    enrdict={}
    modulelendict={}
    moduleposlendict={}
    zscoresigndict={}
    zscoredict={}
    # num_results=10
    for index, row in adj_pips.iterrows():
        if nk.loc[index].values[0] > 0.0:
            print(index)
            k = int(index.split("_")[-1])
            ingroupcells = labelsdf.cell_id[labelsdf.inferred_label == k]
            outgroupcells = labelsdf.cell_id[labelsdf.inferred_label != k]
            popmean = np.mean(adata[ingroupcells].raw.X,axis=0)[top_genes_bool]
            notpopmean = np.mean(adata[outgroupcells].raw.X,axis=0)[top_genes_bool]
            popvar = np.var(adata[ingroupcells].raw.X,axis=0)[top_genes_bool]
            notpopvar = np.var(adata[outgroupcells].raw.X,axis=0)[top_genes_bool]
            effect_size, effect_size_sign = effect_size_calc(popmean,notpopmean,popvar,notpopvar)
            modulegenes = list(genenames[(row >= 0.5)].values)
            modulegenespos = list(genenames[np.where((row >= 0.5) & (effect_size_sign == "+"))])
            if len(modulegenespos) >= 1:
                zscoredict[index]=effect_size
                zscoresigndict[index]=effect_size_sign
                modulegenesdict[index] = modulegenes
                moduleposgenesdict[index] = modulegenespos
                score_name=index+"_score"
                normalized_score_name=index+"_norm_score"
                sc.tl.score_genes(adata,modulegenespos,score_name=score_name)
                normalized_score = (adata.obs[score_name] - np.min(adata.obs[score_name])) / (np.max(adata.obs[score_name])-np.min(adata.obs[score_name]))
                adata.obs[normalized_score_name]=normalized_score
                moduleposlendict[index] = len(modulegenespos)
                modulelendict[index] = len(modulegenes)
                resultsobj= gp.enrichr(modulegenespos,gene_sets=["GO_Biological_Process_2021"],organism="human",outdir=None) #, "Tabula_Sapiens"
                resultsobj.res2d.Term=resultsobj.res2d.Term.str.split(" \(GO").str[0].str.capitalize()
                enrdict[index]=resultsobj
                labels=resultsobj.res2d.Term
                fdrdict[index]=-np.log10(resultsobj.res2d['Adjusted P-value'])
                currxlimmax=np.max(fdrdict[index])
                labelsdict[index]=[l for l in labels]#['\n '.join(wrap(l, 30)) for l in labels]
                if currxlimmax > xlimmax:
                    xlimmax=currxlimmax  
            time.sleep(20)

    num_results_max = 10
    colormappingsdict={k:colorsdict[clusterforwardtranslationdict[k]]  for k in list(fdrdict.keys())}
    xlimmax = np.ceil(xlimmax + 1)

    for key in list(fdrdict.keys()):
        sig_bool  = list((fdrdict[key] >= -np.log10(0.05)).values)
        cluster_gos = {
            "GO_term": [f"UP {x}" for x, keep in zip(labelsdict[key], sig_bool) if keep],
            "-log10FDR": [x for x, keep in zip(fdrdict[key], sig_bool) if keep],
            "cluster_geneset_size": [moduleposlendict[key] for x, keep in zip(fdrdict[key], sig_bool) if keep],
            "xlimmax": [xlimmax for x, keep in zip(fdrdict[key], sig_bool) if keep],
            "color":[colormappingsdict[key] for x, keep in zip(fdrdict[key], sig_bool) if keep]
        }
        cluster_gos_df = pd.DataFrame(cluster_gos)
        outfile_name=f"{gocsvspath}/{clusterforwardtranslationdict[key]}_GOTerms.csv"
        cluster_gos_df.to_csv(outfile_name) 
                         
if __name__ == '__main__':
    main(sys.argv[1:])