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
    -h=*|--help=*)
    echo "-r | --rlib: path to r library, if applicable"
    echo "-d | --datapath: path to  data (an .rds file) you want to make count and metadata csvs from"
    ;;
    *)
    echo 'Input error. Please use the --help option for input documentation' 
    exit 3 # unknown option
    ;;
esac
done

rlib=$RLIB
data=$DATA
echo rlib = ${rlib}
echo data = ${data}

Rscript read_rds.R $rlib $data
