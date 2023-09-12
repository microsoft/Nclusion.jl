PATH_TO_RESULTS=$1
PATH_TO_JULIA_ENV=$2
PATH_TO_DATA=$3
allCvi=(ch db gd43 gd53 csil ps wb xb)

results_path=$PATH_TO_RESULTS
julia_env=$PATH_TO_JULIA_ENV
data_path=$PATH_TO_DATA

julia --project=${julia_env} utility_scripts/ClusteringQualityMetrics.jl "${results_path}" "" "" "" "" "" 

for j in ${!allCvi[@]};do 
  METHOD="${allCvi[$j]}"
  echo Running ${METHOD} CVI calculations
  julia --project=${julia_env} utility_scripts/QuantitativeEvaluation.jl "${results_path}" "$data_path" "${METHOD}" "5000" "" "" 
done