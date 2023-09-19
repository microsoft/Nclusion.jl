#!/bin/bash
DATAPATH=''
RLIB=''
PIPS=''
RESULTSPATH=''
OUTPUTPATH=''
DATANAME=''

for i in "$@"
do
case $i in
    -l=*|--rlib=*)
    RLIB="${i#*=}"
    ;;
    -d=*|--datapath=*)
    DATA="${i#*=}"
    ;;
    -p=*|--pips=*)
    PIPS="${i#*=}"
    ;;
    -r=*|--nclusion_results=*)
    RESULTSPATH="${i#*=}"
    ;;
    -o=*|--output_dir=*)
    OUTPUTPATH="${i#*=}"
    ;;
    -n=*|--data_name=*)
    DATANAME="${i#*=}"
    ;;
    -h|--help)
    echo "Usage: . analyze_results.sh [--rlib STR] [--datapath STR] [--pips STR] [--nclusion_results STR] [--output_dir STR] [--data_name STR]"
    echo "-l | --rlib: path to r library, if applicable (Optional)"
    echo "-d | --datapath: path to data NCLUSION was ran on (.h5ad file)"
    echo "-p | --pips: path to pips identified by NCLUSION"
    echo "-r | --nclusion_results: path to NCLUSION cluster assignment results"
    echo "-o | --output_dir: path to directory where all plots will be saved"
    echo "-n | --data_name: name of dataset for naming plot files"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation'
    return
    ;;
esac
done

data_path=$DATAPATH
rlib=$RLIB
pips=$PIPS
nclusion_results=$RESULTSPATH
output_dir=$OUTPUTPATH
data_name=$DATANAME
echo Rscript tutorial_scripts/utility_scripts/result_plots.R -path_to_data=$data_path -path_to_pips=$pips -path_to_labels=$nclusion_results -outfile_base=$output_dir -dataset_name=$data_name -r_library_path=$rlib

Rscript tutorial_scripts/utility_scripts/result_plots.R -path_to_data=$data_path -path_to_pips=$pips -path_to_labels=$nclusion_results -outfile_base=$output_dir -dataset_name=$data_name -r_library_path=$rlib