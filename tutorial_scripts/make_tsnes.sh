#!/bin/bash
DATAPATH=''
RESULTSPATH=''
DATANAME=''

for i in "$@"
do
case $i in
    -d=*|--datapath=*)
    DATAPATH="${i#*=}"
    ;;
    -r=*|--nclusion_results=*)
    RESULTSPATH="${i#*=}"
    ;;
    -n=*|--data_name=*)
    DATANAME="${i#*=}"
    ;;
    -h|--help)
    echo "Usage: . make_tsnes.sh [--datapath STR] [--nclusion_results STR] [--data_name STR]"
    echo "-d | --datapath: path to  data (an .rds file) you want to make count and metadata csvs from"
    echo "-r | --nclusion_results: path to NCLUSION cluster assignment results"
    echo "-n | --data_name: name of dataset to append to plot save names"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation'
    return
    ;;
esac
done

data_path=$DATAPATH
nclusion_results=$RESULTSPATH
data_name=$DATANAME

python tutorial_scripts/utility_scripts/plot_tsnes.py --path_to_data $data_path --path_to_labels $nclusion_results --figures_prefix $data_name