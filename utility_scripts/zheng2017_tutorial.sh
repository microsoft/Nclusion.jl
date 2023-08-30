DATA_ROOT="./data/"

mkdir -p tutorial_outputs/

PATH_TO_NCLUSION=$1

#preprocess 
python preprocess_tutorial_data.py --path_to_data tutorial_data/pbmc_raw.h5ad --path_to_save tutorial_data/5000hvgs_pbmc_preprocessed.h5ad --data_name "pbmc" --n_hvgs 5000
python preprocess_tutorial_data.py --path_to_data tutorial_data/pbmc_raw.h5ad --path_to_save tutorial_data/2500hvgs_pbmc_preprocessed.h5ad --data_name "pbmc" --n_hvgs 2500

#run NCLUSION
cd $PATH_TO_NCLUSION
. /home/v-mahughes/SCoOP-sc/tinker/EXPERIMENT_HOFF_VS/mpu/bash_unbounded_vary_sparsity_lognormHVG_noscale_anndata_mpu.sh "utility_scripts/tutorial_data/5000hvgs_pbmc_preprocessed.h5ad" true "$PATH_TO_NCLUSION" "pbmc" 5000

cd $PATH_TO_NCLUSION
. /home/v-mahughes/SCoOP-sc/tinker/EXPERIMENT_HOFF_VS/mpu/bash_unbounded_vary_sparsity_lognormHVG_noscale_anndata_mpu.sh "utility_scripts/tutorial_data/2500hvgs_pbmc_preprocessed.h5ad" true "$PATH_TO_NCLUSION" "pbmc" 2500

#intrinsic + extrinsic metric evaluation
PATH_TO_RESULTS_5000 = 
allCvi=(ch db gd43 gd53 csil ps wb xb)
cd $PATH_TO_NCLUSION
echo "$PATH_TO_RESULTS"
julia --project=$PATH_TO_NCLUSION  $PATH_TO_NCLUSION/baselines/ClusteringQualityMetrics.jl "$PATH_TO_RESULTS" "" "" "" "" "" 
for j in ${!allCvi[@]};do 
    METHOD="${allCvi[$j]}"
    echo Running $METHOD CVI calculations
    julia --project=$PATH_TO_NCLUSION  $PATH_TO_NCLUSION/quantitative_evaluation.jl "$PATH_TO_RESULTS_5000"  "tutorial_data/5000hvgs_pbmc_preprocessed.h5ad" "$METHOD" "5000" "" "" 
done

PATH_TO_RESULTS_2500 = 
allCvi=(ch db gd43 gd53 csil ps wb xb)
cd $PATH_TO_NCLUSION
echo "$PATH_TO_RESULTS"
julia --project=$PATH_TO_NCLUSION  $PATH_TO_NCLUSION/baselines/ClusteringQualityMetrics.jl "$PATH_TO_RESULTS" "" "" "" "" "" 
for j in ${!allCvi[@]};do 
    METHOD="${allCvi[$j]}"
    echo Running $METHOD CVI calculations
    julia --project=$PATH_TO_NCLUSION  $PATH_TO_NCLUSION/quantitative_evaluation.jl "$PATH_TO_RESULTS_2500"  "tutorial_data/2500hvgs_pbmc_preprocessed.h5ad" "$METHOD" "2500" "" "" 
done

# GENERATE CELL DISTRIBUTION HEATMAPS
python make_heatmaps.py --input_path $PATH_TO_RESULTS_5000 --pct_table_file "tutorial_outputs/pbmc_5000hvgs_NCLUSION_cell_distribution_table.csv" --cluster_dict_file "pbmc_5000hvgs_NCLUSION_cluster_mapping_to_new_label.csv" 
Rscript make_heatmaps.R -f "pbmc_5000hvgs_NCLUSION_cell_distribution_table.csv" -s complex_heatmap.pdf

python make_heatmaps.py --input_path $PATH_TO_RESULTS_5000 --pct_table_file "tutotial_outputs/pbmc_2500hvgs_NCLUSION_cell_distribution_table.csv" --cluster_dict_file "pbmc_2500hvgs_NCLUSION_cluster_mapping_to_new_label.csv" 
Rscript make_heatmaps.R -f "pbmc_2500hvgs_NCLUSION_cell_distribution_table.csv" -s complex_heatmap.pdf

# MAKE TSNES 
python plot_tsnes.py --path_to_data tutorial_data/5000hvgs_pbmc_preprocessed.h5ad --path_to_labels $PATH_TO_RESULTS_5000 --figure_prefix tutorial_outputs/pbmc_5000hvgs_NCLUSION
python plot_tsnes.py --path_to_data tutorial_data/2500hvgs_pbmc_preprocessed.h5ad --path_to_labels $PATH_TO_RESULTS_2500 --figure_prefix tutorial_outputs/pbmc_2500hvgs_NCLUSION

#MAKE HEATMAPS + VIOLIN + GO PLOTS
Rscript nclusion/utility_scripts/result_plots.R -path_to_data="/home/v-mahughes/nclusion_preprocessed_data/galenAML2019/5000hvgs_galenAML_preprocessed.h5ad" -path_to_pips="/home/v-mahughes/RESULTS_81523/galenAML/NCLUSION/5000hvgs_alpha_01/5000G-2023-08-17T193822-pips.csv" -path_to_labels="/home/v-mahughes/RESULTS_81523/galenAML/NCLUSION/5000hvgs_alpha_01/galenAML_scanpy-5000G_nclusion_5000genes-2023-08-17T193822.csv" -outfile_base="/home/v-mahughes/test_plot_script/" -dataset_name='' -r_library_path="/home/v-mahughes/r_packages467"

Rscript nclusion/utility_scripts/result_plots.R "/home/v-mahughes/nclusion_preprocessed_data/galenAML2019/5000hvgs_galenAML_preprocessed.h5ad" "/home/v-mahughes/RESULTS_81523/galenAML/NCLUSION/5000hvgs_alpha_01/5000G-2023-08-17T193822-pips.csv" "/home/v-mahughes/RESULTS_81523/galenAML/NCLUSION/5000hvgs_alpha_01/galenAML_scanpy-5000G_nclusion_5000genes-2023-08-17T193822.csv" "/home/v-mahughes/test_plot_script/" '' "/home/v-mahughes/r_packages"


