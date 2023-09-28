import pandas as pd 
import sys
import getopt
import numpy as np
from itertools import permutations

def get_pcts(data):
    
    cell_types = data["cell_type"].unique()
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

    return pct_df
    
    
            
def optimize_diag(df, col_sum_threshold):  

    col_sum_dict = {}
    for col in df.columns.values:
        col_sum_dict[col] = df.loc[:, col].sum()
    keep_cols = dict(filter(lambda elem: elem[1] > col_sum_threshold, col_sum_dict.items()))
    print(keep_cols)
    x = len(keep_cols)
    df2 = df.loc[:, keep_cols]
    
    perm = list(permutations(df2.columns.values, x))
        
    print('found all permutations')

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


def main(argv):
    input_path = ''
    cluster_dict_file = 'mapping_cluster_to_new_label.csv'
    pct_table_file = 'cell_type_distribution_table.csv'
    data_name = ''
    cluster_dict = None
    
    try:
        opts, args = getopt.getopt(argv, '', ["input_path=", "save_path=", "pct_table_file=", "cluster_dict_file=", "data_name="])
        print(opts)
                                
        for opt, arg in opts:
            if opt == '--input_path':
                input_path = arg
            elif opt == '--cluster_dict_file':
                cluster_dict_file = arg
            elif opt == '--pct_table_file':
                pct_table_file = arg
            elif opt == '--col_sum_threshold':
                col_sum_threshold = float(arg)
            elif opt == '--col_sum_threshold':
                col_sum_threshold = float(arg)
            elif opt == '--data_name':
                data_name = arg
    
    except getopt.GetoptError: 
        sys.exit()
        print(opts)   
        
    if data_name == 'vanGalen':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/vanGalen_mapping.csv', index_col=0).to_dict()
        cluster_dict = clusters['0']
    elif data_name == 'dominguezconde':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/dominguezconde_mapping.csv', index_col=0).to_dict()
        cluster_dict = clusters['0']
    elif data_name == 'raghavan':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/raghavan_mapping.csv', index_col=0).to_dict()
        cluster_dict = clusters['0']
    elif data_name == 'zheng':
        clusters = pd.read_csv('tutorial_scripts/cluster_mapping/zheng_mapping.csv', index_col=0).to_dict()
        cluster_dict = clusters['0']

    data = pd.read_csv(input_path, index_col=1)
    
    pct_table = get_pcts(data)
    
    pct_table = pct_table.reindex(columns=list(cluster_dict.keys()))
  
    pct_table = pct_table.rename(columns=cluster_dict)

    pct_table.to_csv(pct_table_file) 
    
                            
if __name__ == '__main__':
    main(sys.argv[1:])