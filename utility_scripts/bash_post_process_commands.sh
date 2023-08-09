DATA_ROOT="./data/"

# . bash_quantitative_evaluation ~/nclusion ~/nclusion_preprocessed_data/galenAML2019/1000hvgs_galenAML_preprocessed.h5ad ~/nclusion_figures/galenAML_outputs/leiden_knn/1000hvgs/galenAML_1000hvgs_leiden.csv 1000

MAINDIR=$1 
DATAPATH=$2 
RESPATH=$3
HVG=$4

path_to_nclsn=$MAINDIR 
path_to_data=$DATAPATH
path_to_labels=$RESPATH
n_hvg=$HVG

allCvi=(ch db gd43 gd53 csil ps wb xb)

cd $path_to_nclsn

echo "$path_to_labels"
julia --project=$path_to_nclsn  $path_to_nclsn/baselines/ClusteringQualityMetrics.jl "$path_to_labels" "" "" "" "" "" 
for j in ${!allCvi[@]};do 
    METHOD="${allCvi[$j]}"
    echo Running $METHOD CVI calculations
    julia --project=$path_to_nclsn  $path_to_nclsn/quantitative_evaluation.jl "$path_to_labels"  "$path_to_data" "$METHOD" "$n_hvg" "" "" 
done


