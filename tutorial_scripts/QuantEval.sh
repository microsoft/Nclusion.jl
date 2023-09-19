#!/bin/bash
PATH_TO_RESULTS=$1
PATH_TO_JULIA_ENV=$2
PATH_TO_DATA=$3
allCvi=(ch db gd43 gd53 csil ps wb xb)

for i in "$@"
do
case $i in
    -r=*|--results_path=*)
    PATH_TO_RESULTS="${i#*=}"

    ;;
    -d=*|--datapath=*)
    PATH_TO_DATA="${i#*=}"
    ;;
    -j=*|--julia_env=*)
    PATH_TO_JULIA_ENV="${i#*=}"
    ;;
    -h|--help)
    echo "Usage: . convert_galen_dtype.sh [--results_path STR] [--datapath STR] [--julia_env STR]"
    echo "-r | --results_path: path to NCLUSION cluster assignment results"
    echo "-d | --datapath: path to  data (an .head file) that NCLUSION was ran on"
    ECHO "J | --julia_env: path to Julia project environment where required tutorial R packages are installed"
    return
    ;;
    *)
    echo 'Input error. Please use the --help flag for input documentation'
    return
    ;;
esac
done

results_path=$PATH_TO_RESULTS
julia_env=$PATH_TO_JULIA_ENV
data_path=$PATH_TO_DATA

echo results_path = $results_path
echo julia_env = $julia_env
echo data_path = $data_path

julia --project=${julia_env} tutorial_scripts/utility_scripts/ClusteringQualityMetrics.jl "${results_path}" "" "" "" "" "" 

for j in ${!allCvi[@]};do 
  METHOD="${allCvi[$j]}"
  echo Running ${METHOD} CVI calculations
  julia --project=${julia_env} tutorial_scripts/utility_scripts/QuantitativeEvaluation.jl "${results_path}" "$data_path" "${METHOD}" "5000" "" "" 
done