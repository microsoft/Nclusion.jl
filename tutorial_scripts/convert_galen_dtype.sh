#!/bin/bash
RLIB=NULL
DATA="tutorial_data/vanGalen2019/Seurat_AML.rds"

for i in "$@"
do
case $i in
    -r=*|--rlib=*)
    RLIB="${i#*=}"
    ;;
    -d=*|--datapath=*)
    DATA="${i#*=}"
    ;;
    -h|--help)
    echo "Usage: . convert_galen_dtype.sh [--rlib STR] [--datapath STR]"
    echo "-r | --rlib: path to r library, if applicable (Optional)"
    echo "-d | --datapath: path to  data (an .rds file) you want to make count and metadata csvs from"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation'
    return
    ;;
esac
done

rlib=$RLIB
data=$DATA


Rscript tutorial_scripts/utility_scripts/read_rds.R $rlib $data
