#!/bin/bash
INPUTFILE=$1
JULIAENV=$2
DATASET=$3
OUTPUTDIR=$4
NUM_THREADS=auto

for i in "$@"
do
case $i in
    -d=*|--datapath=*)
    INPUTFILE="${i#*=}"

    ;;
    -j=*|--julia_env=*)
    JULIAENV="${i#*=}"
    ;;
    -n=*|--data_name=*)
    DATASET="${i#*=}"
    ;;
    -o=*|--output_dir=*)
    OUTPUTDIR="${i#*=}"
    ;;
    -h=*|--help=*)
    echo "Usage: . run_nclusion [--datapath STR] [--julia_env STR] [--dataset_name STR] [--output_dir STR]"
    echo "-d | --datapath: path to preprocessed scRNA-seq dataset (an .h5ad file) to run on NCLUSION"
    echo "-j | --julia_env: path to Julia project environment where NCLUSION is installed"
    echo "-n | --dataset_name: dataset name that will be used to label NCLUSION output files"
    echo "-o | --output_dir: path to the directory where NCLUSION results will be stored"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation' 
    return
    ;;
esac
done

input_file=$INPUTFILE
julia_env=$JULIAENV
dataset_name=$DATASET
output_dir=$OUTPUTDIR

if [ $dataset_name == "zheng2017" ];
then
    a=0.0000001
    b=0.0000001
else
    a=1.0
    b=1.0
fi
k=25
elbo_ep="1.0"
julia --project=${julia_env} --thread=${NUM_THREADS} tutorial_scripts/utility_scripts/run_nclusion.jl  "${input_file}" "${k}" "${a}" "${b}" "12345" "${elbo_ep}" "150" "${dataset_name}" "true" "${output_dir}" 
# julia --project=${julia_env} --thread=${NUM_THREADS} tutorial_scripts/utility_scripts/run_script.jl  "${input_file}" "${k}" "${a}" "${b}" "12345" "${elbo_ep}" "150" "${dataset_name}" "true" "${output_dir}" 
