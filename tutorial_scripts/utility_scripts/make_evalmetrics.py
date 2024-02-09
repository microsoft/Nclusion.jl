import pandas as pd 
import sys
import getopt
import numpy as np
from sklearn import metrics
from sklearn.metrics import pairwise_distances
import cvi
import scanpy as sc
from scipy.stats import entropy


def scale_labels(labels):
    new_labels = np.copy(labels)
    if np.min(new_labels) != 0:
        new_labels = new_labels -  np.min(new_labels)
    return new_labels

def species_richness(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    return len(set(data))

def shannon_index(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    probs = [data.count(category) / len(data) for category in set(data)]
    return entropy(probs, base=2)

def shannon_eveness(data):
    ''' from: https://cran.r-project.org/web/packages/tabula/vignettes/diversity.html'''
    if isinstance(data,np.ndarray):
        data = list(data)
    return  shannon_index(data)/ np.log(len(set(data)))

def simpsons_index(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    probs = [data.count(category) / len(data) for category in set(data)]
    return sum(p**2 for p in probs)

def gini_simpson_index(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    probs = [data.count(category) / len(data) for category in set(data)]
    return 1 - sum(p**2 for p in probs)

def berger_parker_index(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    max_frequency = max(data.count(category) for category in set(data))
    return max_frequency / len(data)

def diversity_order_1(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    return np.exp(shannon_index(data))

def diversity_order_2(data):
    if isinstance(data,np.ndarray):
        data = list(data)
    return 1 / simpsons_index(data)

def mcintosh_index(data):
    ''' from: https://cran.r-project.org/web/packages/tabula/vignettes/diversity.html''' 
    if isinstance(data,np.ndarray):
        data = list(data)
    return (len(data) - np.sqrt(np.sum([data.count(category)**2 for category in set(data)]))) / (len(data) - np.sqrt(len(data)))

def mcintosh_eveness(data): 
    ''' from: https://cran.r-project.org/web/packages/tabula/vignettes/diversity.html'''
    if isinstance(data,np.ndarray):
        data = list(data)
    return (len(data) - np.sqrt(np.sum([data.count(category)**2 for category in set(data)]))) / (len(data) -  len(data)/np.sqrt(len(set(data))))

def r1_richness(data):
    ''' from: https://cran.r-project.org/web/packages/tabula/vignettes/diversity.html'''
    if isinstance(data,np.ndarray):
        data = list(data)
    return (len(set(data)) - 1) / np.log(len(data))

def r2_richness(data):
    ''' from: https://cran.r-project.org/web/packages/tabula/vignettes/diversity.html'''
    if isinstance(data,np.ndarray):
        data = list(data)
    return (len(set(data))) / np.log10(len(data))

def r3_richness(data):
    ''' from: https://cran.r-project.org/web/packages/tabula/vignettes/diversity.html'''
    if isinstance(data,np.ndarray):
        data = list(data)
    return (len(set(data))) / np.sqrt(len(data))

def main(argv):
    labels_path = ''
    data_path = ''
    eval_metric = ''
    condition_idx = ''
    data_partition = ''
    method = ''
    dataset = ''
    hvgs = ''
    output_path_file = ''
    output_suffix = ''
    
    # idx,method,dataset,hvgs,Metric,

    



    try:
        opts, args = getopt.getopt(argv, '', ["labels_path=", "data_path=","eval_metric=","condition_idx=", "data_partition=","method=", "dataset=", "hvgs=",  "output_path_file=","output_suffix="])#
        print(opts)
                                
        for opt, arg in opts:
            if opt == '--labels_path':
                labels_path = arg
            elif opt == '--data_path':
                data_path = arg
            elif opt == '--eval_metric':
                eval_metric = arg
            elif opt == '--condition_idx':
                condition_idx = int(arg)
            elif opt == '--data_partition':
                if data_partition != "":
                    data_partition = int(arg)    
                else:
                    data_partition = ""
            elif opt == '--method':
                method = arg
            elif opt == '--dataset':
                dataset = arg
            elif opt == '--hvgs':
                hvgs = int(arg)
            elif opt == '--output_path_file':
                output_path_file = arg
            elif opt == '--output_suffix':
                output_suffix = arg


        
    except getopt.GetoptError: 
        sys.exit()
        print(opts) 

    
    print(output_path_file)
    labels = pd.read_csv(labels_path)
    print(labels)
    #evalution metrics available:
    #ri ari mi ami nmi homogeneity completeness jaccardskl vmeasure fmscore sil_skl ch_skl db_skl ch db gd43 gd53 csil ps wb xb speciesrichness shannonindex shannoneveness simpsonsindex ginisimpsonindex bergerparkerindex diversityorder1 diversityorder2 mcintoshindex mcintosheveness r1richness r2richness r3richness ##rcip

    called_labels = scale_labels(labels.called_label.values)
    inferred_labels = scale_labels(labels.tranlsated_inferred_label.values)
    score = 0.0
    if eval_metric == "ri":
        score = metrics.rand_score(called_labels, inferred_labels)
    elif eval_metric == "ari":
        score = metrics.adjusted_rand_score(called_labels, inferred_labels)
    elif eval_metric == "mi":
        score = metrics.mutual_info_score(called_labels, inferred_labels)
    elif eval_metric == "ami":
        score = metrics.adjusted_mutual_info_score(called_labels, inferred_labels) 
    elif eval_metric == "nmi":
        score = metrics.normalized_mutual_info_score(called_labels, inferred_labels)
    elif eval_metric == "homogeneity":
        score = metrics.homogeneity_score(called_labels, inferred_labels)
    elif eval_metric == "completeness":
        score = metrics.completeness_score(called_labels, inferred_labels)
    elif eval_metric == "jaccardskl":
        score = metrics.jaccard_score(called_labels, inferred_labels, average='macro')
    elif eval_metric == "vmeasure":
        beta=1.0
        score = metrics.v_measure_score(called_labels, inferred_labels, beta=beta)
    elif eval_metric == "fmscore":
        score = metrics.fowlkes_mallows_score(called_labels, inferred_labels)
    elif eval_metric == "silskl":
        distance_metric = 'euclidean'
        eval_metric = f"{eval_metric}-{distance_metric}"
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = metrics.silhouette_score(X, inferred_labels, metric=distance_metric)
    elif eval_metric == "chskl":
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = metrics.calinski_harabasz_score(X, inferred_labels)
    elif eval_metric == "dbskl":
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = metrics.davies_bouldin_score(X, inferred_labels)
    elif eval_metric == "ch":
        my_cvi = cvi.CH()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "db":
        my_cvi = cvi.DB()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "gd43":
        my_cvi = cvi.GD43()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "gd53":
        my_cvi = cvi.GD53()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "csil":
        my_cvi = cvi.cSIL()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "ps":
        my_cvi = cvi.PS()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "rcip":
        my_cvi = cvi.rCIP()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "wb":
        my_cvi = cvi.WB()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "xb":
        my_cvi = cvi.XB()
        adata = sc.read_h5ad(data_path)
        X = adata.X
        score = my_cvi.get_cvi(X, inferred_labels)
    elif eval_metric == "speciesrichness":
        score = species_richness(inferred_labels)
    elif eval_metric == "shannonindex":
        score = shannon_index(inferred_labels)
    elif eval_metric == "shannoneveness":
        score = shannon_eveness(inferred_labels)
    elif eval_metric == "simpsonsindex":
        score = simpsons_index(inferred_labels)
    elif eval_metric == "ginisimpsonindex":
        score = gini_simpson_index(inferred_labels)
    elif eval_metric == "bergerparkerindex":
        score = berger_parker_index(inferred_labels)
    elif eval_metric == "diversityorder1":
        score = diversity_order_1(inferred_labels)
    elif eval_metric == "diversityorder2":
        score = diversity_order_2(inferred_labels)
    elif eval_metric == "mcintoshindex":
        score = mcintosh_index(inferred_labels)
    elif eval_metric == "mcintosheveness":
        score = mcintosh_eveness(inferred_labels)
    elif eval_metric == "r1richness":
        score = r1_richness(inferred_labels)
    elif eval_metric == "r2richness":
        score = r2_richness(inferred_labels)
    elif eval_metric == "r3richness":
        score = r3_richness(inferred_labels)



    results_dict={}
    results_dict["condition_idx"] = [condition_idx]
    if data_partition == "":
        data_partition = "ALL"
    results_dict["data_partition"] = [data_partition]
    results_dict["method"] = [method]
    results_dict["dataset"] = [dataset]
    results_dict["hvgs"] = [hvgs]
    results_dict["eval_metric"] = [eval_metric]
    results_dict["score"] = [score]

    results_df = pd.DataFrame.from_dict(results_dict, orient='columns')
    save_filename = f"{eval_metric}_{method}_{dataset}_{hvgs}hvgs{output_suffix}.csv"
    outfile_results=output_path_file+save_filename
    results_df.to_csv(outfile_results, index=False)

                            
if __name__ == '__main__':
    main(sys.argv[1:])