#!/bin/bash
RESULTSPATH=''
PCTTABLE=''
CLUSTERDICT=''
HEATMAP=''
RLIB=NULL
DATANAME=''

for i in "$@"
do
case $i in
    -r=*|--nclusion_results=*)
    RESULTSPATH="${i#*=}"
    ;;
    -t=*|--pct_table=*)
    PCTTABLE="${i#*=}"
    ;;
    -d=*|--cluster_dict=*)
    CLUSTERDICT="${i#*=}"
    ;;
    -m=*|--save_heatmap=*)
    HEATMAP="${i#*=}"
    ;;
    -l=*|--rlib=*)
    RLIB="${i#*=}"
    ;;
    -n=*|--data_name=*)
    DATANAME="${i#*=}"
    ;;
    -h|--help)
    echo "Usage: . preprocess.sh [--nclusion_results STR] [--pct_table STR] [--cluster_dict STR] [--save_heatmap INT]"
    echo "-r | --nclusion_results: path to nclusion cluster assignment results file"
    echo "-t | --pct_table: path to save cell type distribution (.csv) file"
    echo "-d | --cluster_dict: path to file that maps orignal nclusion cluster labels to new labels on heatmap"
    echo "-m | --save_heatmap: path to save the cell distribution heatmap (.pdf)"
    echo "-n | --data_name: name of dataset"
    echo "-l | --rlib: path to r library (Optional)"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation'
    return
    ;;
esac
done

nclusion_results=$RESULTSPATH
pct_table=$PCTTABLE
cluster_dict=$CLUSTERDICT
save_heatmap=$HEATMAP
rlib=$RLIB
data_name=$DATANAME


python tutorial_scripts/utility_scripts/make_heatmaps.py --input_path $nclusion_results --pct_table_file $pct_table --cluster_dict_file $cluster_dict --data_name $data_name

Rscript tutorial_scripts/utility_scripts/make_heatmaps.R $pct_table $save_heatmap $rlib