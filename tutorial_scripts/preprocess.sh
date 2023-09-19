#!/bin/bash
PATHTODATA=''
PATHTOSAVE=''
DATANAME=''
HVGS=5000

for i in "$@"
do
case $i in
    -d=*|--datapath=*)
    PATHTODATA="${i#*=}"
    ;;
    -s=*|--savepath=*)
    PATHTOSAVE="${i#*=}"
    ;;
    -n=*|--data_name=*)
    DATANAME="${i#*=}"
    ;;
    -g=*|--n_hvgs=*)
    HVGS="${i#*=}"
    ;;
    -h|--help)
    echo "Usage: . preprocess.sh [--datapath STR] [--savepath STR] [--data_name STR] [--n_hvgs INT]"
    echo "-d | --datapath: path to data to be preprocessed"
    echo "-s | --savepath: path to save preprocessed data (as an .h5ad file)"
    echo "-n | --data_name: name of the tutorial dataset to be preprocessed"
    echo "-g | --n_hvgs: number of highly variable genes to extract from data"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation'
    return
    ;;
esac
done

path_to_data = $PATHTODATA
path_to_save = $PATHTOSAVE
data_name = $DATANAME
n_hvgs = $HVGS

python tutorial_scripts/utility_scripts/preprocess_tutorial_data.py --path_to_data $path_to_data --path_to_save $path_to_save --data_name $data_name --n_hvgs $n_hvgs