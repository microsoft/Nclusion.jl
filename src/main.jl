using Nclusion
using Random
using Statistics, LinearAlgebra, StatsBase, Distributions
using Clustering
using Logging,LoggingExtras
using JLD2
logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;
function flushed_logger(msg;logger=nothing)
    if !isnothing(logger)
        with_logger(logger) do
            @info msg
        end
    end
end

#datafilename1 = "/Users/chibuikemnwizu/Documents/Nclusion.jl/data/write/pbmc3k.h5ad";
# run_nclusion(datafilename1;logger=logger,outdir = "/Users/chibuikemnwizu/Documents/Nclusion.jl")
############################################
############################################
############################################
# Some notes to myself
# 1. Because in scanpy we scale the data but clip values at 10, the data is not centered and scaled perfectly on the columms. Does this cause any issues? May need to look at the impact of clipping on the clusters we infer.
# 2. Recall that scanpy defaults to float32 while julia's default is float64. This shouldnt cause any issues but may be worth noting.
function main(ARGS; logger=nothing)
    datafilename1,KMax,significance_prop,min_number_cells,min_percent_of_genes,max_percent_of_genes,seed,elbo_ep,num_iter,dataset_name,outdir,m_initialization_approach,n_hvgs,time_key,individuals_key,slurm_job_id,change_seeds,add_time, add_individuals, layer_index,update_clusterwise, use_alt_representation, check_cluster_interpretability_bool, gene_set_file_path, samplebased_alpha0,train_h = ARGS
    if !isempty(significance_prop)
        significance_prop = parse(Float64, significance_prop)
    else
        significance_prop=0.5
    end
    if !isempty(min_number_cells)
        min_number_cells = parse(Float64, min_number_cells)
    else
        min_number_cells=1.0
    end
    if !isempty(min_percent_of_genes)
        min_percent_of_genes = parse(Float64, min_percent_of_genes)
    else
        min_percent_of_genes=0.05
    end
    if !isempty(max_percent_of_genes)
        max_percent_of_genes = parse(Float64, max_percent_of_genes)
    else
        max_percent_of_genes=0.75
    end
    if !isempty(KMax)
        KMax = parse(Int64, KMax)
    else
        KMax = 25
    end
    if !isempty(seed)
        seed = parse(Int64, seed)
    else
        seed = 12345
    end
    if !isempty(elbo_ep)
        elbo_ep = parse(Float64, elbo_ep)
    else
        elbo_ep = 10^(-6)
    end
    if !isempty(num_iter)
        num_iter = parse(Int64, num_iter)
    else
        num_iter = 500
    end
    if isempty(dataset_name)
        dataset_name = ""
    end
    if isempty(outdir)
        outdir = ""
    end
    if isempty(m_initialization_approach)
        m_initialization_approach = "rand" # element of ["rand","kpp","kpp+rand","cell_label","cell_label+rand"]
    end
    if isempty(time_key)
        time_key = nothing
    end
    if isempty(individuals_key)
        individuals_key = nothing
    end
    if !isempty(n_hvgs)
        if n_hvgs == "None"
            n_hvgs = nothing
        else
            n_hvgs = parse(Int64, n_hvgs)
        end
    else
        n_hvgs = nothing
    end
    if isempty(slurm_job_id)
        slurm_job_id = "No SLURM Job ID Provided"
    end
    if !isempty(change_seeds)
        change_seeds = parse(Bool, change_seeds)
    else
        change_seeds = false
    end
    if !isempty(add_time)
        add_time = parse(Bool, add_time)
    else
        add_time = false
    end
    if !isempty(add_individuals)
        add_individuals = parse(Bool, add_individuals)
    else
        add_individuals = false
    end
    if !isempty(layer_index)
        layer_index = parse(Int64, layer_index)
    else
        layer_index = 0
    end
    if !isempty(update_clusterwise)
        update_clusterwise = parse(Bool, update_clusterwise)
    else
        update_clusterwise = false
    end
    if isempty(use_alt_representation)
        use_alt_representation = "None" # element of ["None","xmat","factor","pca"]
    end
    if !isempty(check_cluster_interpretability_bool)
        check_cluster_interpretability_bool = parse(Bool, check_cluster_interpretability_bool)
    else
        check_cluster_interpretability_bool = false
    end
    if isempty(gene_set_file_path)
        gene_set_file_path = ""
    end
    if !isempty(samplebased_alpha0)
        samplebased_alpha0 = parse(Bool, samplebased_alpha0)
    else
        samplebased_alpha0 = false
    end
    if !isempty(train_h)
        train_h = parse(Bool, train_h)
    else
        train_h = false
    end
    # datafilename1 =  "/users/cnwizu/data/cnwizu/datasets/KleinLabData/in_vivo/adata_Weinreb2020_in_vivo_gene_filtered_nclusion_preprocessed.h5ad" 
    datafilename1 = "/Users/chibuikemnwizu/Documents/Nclusion.jl/data/write/pbmc3k.h5ad"
    significance_prop=0.5
    min_number_cells=1.0
    min_percent_of_genes=0.05
    max_percent_of_genes=0.75
    KMax = 25
    seed = 12345
    elbo_ep = 10^(-6)
    num_iter = 500
    dataset_name = ""
    outdir = ""
    m_initialization_approach = "rand"
    time_key = nothing
    individuals_key = nothing
    n_hvgs = nothing
    slurm_job_id = "No SLURM Job ID Provided"
    change_seeds = false
    add_time = false
    add_individuals = false
    layer_index = 0
    update_clusterwise = false
    use_alt_representation = "None"
    check_cluster_interpretability_bool = false
    gene_set_file_path = ""
    samplebased_alpha0 = false
    elbo_sign_decrease_max_tolerance=10.0
    train_h = false
    alpha0 = 1.0;
    gamma0 = 1.0; 
    phi1 = 1.0;
    phi2 = 1.0; 
    kappa1 = 1.0;
    kappa2 = 1.0;
    xi1 = 1.0;
    xi2 = 1.0;
    varphi1 = 1.0;
    varphi2 = 1.0;
    script_name = ""
    dataset_name = ""
    report_creation_script_path = ""
    make_abridged_bool_as_string = "True"
    parameter_settings_to_reproduce_run = initialize_ordered_dict(;key_type=String,val_type=String)
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,datafilename1,dataset_name,outdir,script_name,slurm_job_id,KMax,significance_prop,min_number_cells,min_percent_of_genes,max_percent_of_genes,seed,elbo_ep,num_iter,m_initialization_approach,n_hvgs,time_key,individuals_key,change_seeds,add_time,add_individuals, layer_index,update_clusterwise, use_alt_representation, check_cluster_interpretability_bool, gene_set_file_path, samplebased_alpha0,elbo_sign_decrease_max_tolerance,report_creation_script_path,make_abridged_bool_as_string)

    K = KMax;
    Random.seed!(seed);
    iseeds_init = nothing;
    m_mu_init_init = nothing;
    m_mu_init_K = K;
    r_init = nothing;
    h1_init = nothing;
    h2_init = nothing;
    w1_init= nothing;
    w2_init= nothing;
    a_init = nothing;
    b_init = nothing;
    u_init = nothing;
    v_init = nothing;
    define_global_sparsity = true;
    define_global_autocorrelation = true;
    sparsity_based_on_used_representation_dim = false;
    # train_h = false;#true; #
    train_w = true;
    train_ab = true;
    train_uv = true;
    multiple_s_sq_updates = true;
    uv_scalar = 100.0;
    use_std_for_hvg = true;
    center_data_cols=true;
    scale_data_cols=true;
    min_percent_cells = 1/min_number_cells;#round(Int,N*min_percent_cells);
    min_genes_detected_for_mask=10
    standardization_of_used_representation=nothing
    rand_init_inputs = true;
    uniform_theta_init = false;
    num_samples_posterior_samples = 1000
    s_value_thresh=0.05
    burnin = 10;



    anndata_dict1= load_data(datafilename1,seed);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,K,train_h,train_w,train_ab,train_uv,uv_scalar,multiple_s_sq_updates,use_std_for_hvg,center_data_cols,scale_data_cols,min_percent_cells,min_genes_detected_for_mask,standardization_of_used_representation,rand_init_inputs,uniform_theta_init,iseeds_init,m_mu_init_init,num_samples_posterior_samples,s_value_thresh,sparsity_based_on_used_representation_dim,burnin);
    ################################################
    ################################################
    ################################################
    ################################################
    ################################################


    # # ## TEST PARAMATERS THAT ARE REPL FREINDLY
    parameter_settings_to_reproduce_run = initialize_ordered_dict(;key_type=String,val_type=String)
    DATASEED=860#3092#
    datafilename1 = "/Users/chibuikemnwizu/Documents/Nclusion.jl/data/write/pbmc3k.h5ad";#"/users/cnwizu/data/cnwizu/nclusion_manuscript_figure_reproducibility/simulations/revision_round1_20250204_scDesign3_sim/datasets/$(DATASEED)_5K_10000cells_50MarkersPerCluster_zhengPMBC100ksubsample.h5ad";#"/users/cnwizu/data/cnwizu/nclusion+/simulations/scDesign3_sim_zhengPBMC100ksubsample/860_5K_10000cells_50MarkersPerCluster_zhengPMBC100ksubsample.h5ad";#"/users/cnwizu/data/cnwizu/nclusion+/simulations/gmm_simulations/_42_KMax5_J10_K5_N10000_balanced_adata.h5ad"#"/users/cnwizu/data/cnwizu/datasets/KleinLabData/in_vivo/adata_Weinreb2020_in_vivo_gene_filtered_nclusion_preprocessed_subsampled-for-testing_8000cells-equal-celltype-occurance.h5ad"#"/users/cnwizu/data/cnwizu/datasets/KleinLabData/in_vivo/adata_Weinreb2020_in_vivo_gene_filtered_nclusion_preprocessed.h5ad" #"/users/cnwizu/data/cnwizu/datasets/winterTMEandPlasticity2021/manuscript_data/adata_files/Primary_Organoids_RawDGE_metadata_210224_32060Cells_PANFR0575_preprocessed.h5ad" 
    deg = nothing
    dataset_name =""
    gene_set_file_path = ""#"/users/cnwizu/data/cnwizu/datasets/contrasts+perturbations/mh.all.v2023.2.Mm.symbols.gmt"
    significance_prop = 	 0.5 
    min_percent_of_genes = 	 0.05 
    max_percent_of_genes = 	 0.5 # 0.1
    KMax = 	 5;#50;#10;#5;#50;#
    # KMin = 10;#5;#
    seed = 	 2020#12345#
    outdir = "/Users/chibuikemnwizu/Documents/Nclusion.jl/"
    num_iter = 100#1000#10#
    elbo_ep = 1.0e-6
    script_name = ""
    time_key = nothing#"time_label"
    check_cluster_interpretability_bool = false
    individuals_key = nothing
    m_initialization_approach = "rand";#"kmeans";#"kmcen";#""kpp"#"cell_label";#"fcmeans";#"kpp+rand"#"kpp"#
    slurm_job_id = "No SLURM Job ID Provided"
    n_hvgs = 2000
    change_seeds = true
    samplebased_alpha0 = false#true#
    add_time = false
    add_individuals = false
    layer_index = 0
    layer_name = "xmat"#"LdvaeLatentEmbeddingsZ10"#"logcounts"#"scaledata"#"xmat"#"factor"#"factor"#"pca"# nothing#
    is_precomputed_latent_representation =  occursin(lowercase("Embeddings"), lowercase(layer_name)) ?  true  : false
    update_clusterwise = false#true#
    m_mu_init = nothing;
    # use_alt_representation = 
    K = KMax;
    Random.seed!(seed);
    iseeds_init = nothing;
    m_mu_init_init = nothing;
    m_mu_init_K = K;
    r_init = nothing;
    h1_init = nothing;
    h2_init = nothing;
    w1_init= nothing;
    w2_init= nothing;
    a_init = nothing;
    b_init = nothing;
    u_init = nothing;
    v_init = nothing;
    define_global_sparsity = true;
    define_global_autocorrelation = true;
    sparsity_based_on_used_representation_dim = false;
    # train_h = false;#true; #
    train_w = true;
    train_ab = true;
    train_uv = true;
    multiple_s_sq_updates = true;
    uv_scalar = 100.0;
    use_std_for_hvg = true;
    center_data_cols=true;
    scale_data_cols=true;
    min_number_cells = 	 1.0;
    min_percent_cells = 1/min_number_cells;#round(Int,N*min_percent_cells);
    min_genes_detected_for_mask=10
    standardization_of_used_representation=nothing
    rand_init_inputs = true;
    uniform_theta_init = false;
    num_samples_posterior_samples = 1000
    s_value_thresh=0.05
    burnin = 10;
    elbo_sign_decrease_max_tolerance = 10.0;
    K = KMax;
    elbo_sign_change_max = elbo_sign_decrease_max_tolerance;
    nu0 = 0.0;
    sigma_sq_nu = 1e-12;#1.0e12;#1.0;#1.0e12;#
    define_global_sparsity = true;
    define_global_autocorrelation = true;
    sparsity_based_on_used_representation_dim = false;
    train_h = true; #false;#
    eta_update_mode = "Global";#"Genewise";#"Global";#"Clusterwise";#
    sigma_update_mode = "Genewise";#"Global";#"Genewise";#"Clusterwise";#
    lambda_update_mode = "Local";#"Global";#"Genewise";#"Clusterwise";#
    train_w = true;
    train_ab = true;#false;#
    train_uv = true;#false;#
    train_m_nu = true;#false;#
    train_s_sq_nu = true;#false;#
    multiple_s_sq_updates = true;
    uv_scalar = 100.0;#1e24;#
    use_std_for_hvg = true;
    center_data_cols=false;#false;#
    scale_data_cols=false;#false;#
    # min_percent_cells = 1/min_number_cells;#round(Int,N*min_percent_cells);
    min_genes_detected_for_mask=10
    standardization_of_used_representation=nothing#"center_scale"
    rand_init_inputs = true;
    uniform_theta_init = false;
    num_samples_posterior_samples = 1000
    s_value_thresh=0.05
    script_name = ""
    dataset_name = ""
    report_creation_script_path = ""
    make_abridged_bool_as_string = "True"
    burnin = 0;
    remove_small_clusters=true;
    min_number_cells = 	 20;
    min_percent_cells = 1/min_number_cells;
    anndata_dict1= load_data(datafilename1,seed);



    regular_nclusion_note = ""
    time_note1 = ""
    time_note2 = ""
    time_note3 = ""
    individuals_note1 = ""
    individuals_note2 = ""
    individuals_note3 = ""
    if !add_time && !add_individuals
        regular_nclusion_note = ""
    end
    if time_key != nothing && add_time 
        time_note1 = "+"
        time_note2 = "The model uses time information."
        time_note3 = ""
    elseif time_key != nothing && !add_time
        anndata_dict1["obs"][time_key]["codes"] .= 1
        time_note1 = ""
        time_note2 = "The model does runs without time information."# 
        time_note3 = "Time information was provided but not used."# 
    end
    if individuals_key != nothing && add_individuals 
        individuals_note1 = "w/Individuals"
        individuals_note2 = "The model uses individuals information."
        individuals_note3 = ""
    elseif individuals_key != nothing && !add_individuals
        anndata_dict1["obs"][individuals_key]["codes"] .= 1
        individuals_note1 = ""
        individuals_note2 = "The model does runs without individuals information."#
        individuals_note3 = "Individuals information was provided but not used."#
    end
    data_structure_note = "This is $regular_nclusion_note NCLUSION$(time_note1)$(individuals_note1) running. $time_note2 $individuals_note2 $time_note3 $individuals_note3"
    # if haskey(anndata_dict1,"layers")
    #     layer_names = collect(keys(anndata_dict1["layers"]))
    #     if layer_index != 0
    #         layer_name = layer_names[layer_index]
    #         anndata_dict1["X"] = anndata_dict1["layers"][layer_name]
    #     else
    #         layer_name = "None"
    #     end
    # end
    #select_data_representation(anndata_dict1;layer_index=layer_index,layer_name=nothing)
    # highly_variable_genes_bool = get_highly_variable_genes_bool(anndata_dict1;n_hvgs=n_hvgs, use_std=use_std_for_hvg);
    # anndata_dict1 = subset_on_highly_variable_genes_bool(anndata_dict1,highly_variable_genes_bool,layer_index);
    # anndata_dict1 = center_and_scale_data_cols(anndata_dict1;center_cols = center_data_cols,scale_cols = scale_data_cols);
    N = size(anndata_dict1["X"])[2];
    J = size(anndata_dict1["X"])[1];
    
    num_var_feat = J;
    unique_time_id,dataset_used_id,experiment_id = make_ids(dataset_name,J,N);
    gene_names, cell_ids, cell_cluster_labels, time_vec,individuals_vec,cell_cluster_dict = preparing_data(anndata_dict1;time_key=time_key,individuals_key=individuals_key);
    sorting_keys = [(individuals_vec[i],time_vec[i],cell_ids[i]) for i in 1:length(cell_ids)]
    new_order_samples = sortperm(sorting_keys)
    new_order_features = sortperm(gene_names);
    gene_names = gene_names[new_order_features];
    time_vec = time_vec[new_order_samples];
    individuals_vec = individuals_vec[new_order_samples];
    cell_ids = cell_ids[new_order_samples];
    sorted_sorting_keys = sorting_keys[new_order_samples];
    # if !isnothing(cell_cluster_labels)
    #     cell_cluster_labels = cell_cluster_labels[new_order_samples];
    #     if !isnothing(deg) && !isnothing(cell_cluster_dict)
    #         rename!(deg,["$(cell_cluster_dict[el])" for el in names(deg)]);
    #     end
    # end
#   Random.seed!(seed);

    data_input,z,alternative_representation,used_representation_feature_name,layer_name = make_nclusion_inputs(anndata_dict1;time_key=time_key,individuals_key=individuals_key, layer_name=layer_name,gene_set_file_path=gene_set_file_path, min_genes_detected=min_genes_detected_for_mask,standardization_of_used_representation=standardization_of_used_representation,is_precomputed_latent_representation=is_precomputed_latent_representation);
    I_ = length(data_input);
    T = [ length(data_input[i]) for i in 1:I_]
    N_t = [[length(data_input[i][t]) for t in 1:T[i]] for i in 1:I_]
    T_all = sum(T);
    used_representation = hcat([data_input[i][t][n] for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]]...);#anndata_dict1["X"];
    # @test used_representation .== anndata_dict1["obsm"][layer_name][:,alternative_representation[4]]
    if layer_name == "scaledata"
        original_mean = zeros(size(used_representation)[1])#vec(median(used_representation,dims=2))
        # original_mean = vec(mean(anndata_dict1["X"][new_order_features,new_order_samples],dims=2))
        sigma_sq_nu =  1e-12;#
    elseif layer_name == "logcounts"
        used_representation = center_and_scale_matrix_cols(used_representation;center_cols = center_data_cols,scale_cols = scale_data_cols);
        # original_mean = vec(mean(used_representation,dims=2));
        original_mean = zeros(size(used_representation)[1])#vec(median(used_representation,dims=2))
    elseif layer_name == "LdvaeLatentEmbeddingsZ10"
        used_representation = center_and_scale_matrix_cols(used_representation;center_cols = center_data_cols,scale_cols = scale_data_cols);
        # original_mean = vec(mean(used_representation,dims=2));
        original_mean = zeros(size(used_representation)[1])#vec(median(used_representation,dims=2))
    else
        original_mean = vec(mean(used_representation,dims=2));
        used_representation = center_and_scale_matrix_cols(used_representation;center_cols = center_data_cols,scale_cols = scale_data_cols);
    end
    alpha0 = 	1.0 # N#
    gamma0 = 	 1#10#
    phi1 = 	 1.0 
    phi2 = 	 1.0 
    kappa1 = 	1.0# 0.01#<-use for Gene expression space#0.00001#5*0.001*(K*size(used_representation)[1])##5.0#1500# 
    kappa2 = 	1.0#<-use for Gene expression space# 0.1#5*(K*size(used_representation)[1])#10.0#0.5*size(used_representation)[1]#1.0#10*100# 
    xi1 = 	 1.0#<-use for Gene expression space##0.5*(size(used_representation)[1])/0.33 + 1#75/50 * 0.5*(size(used_representation)[1])#
    xi2 = 	 1.0#<-use for Gene expression space#0.001*0.5*(K*size(used_representation)[1])#1.0##0.5*(K*size(used_representation)[1])# 0.5*(K)#
    varphi1 = 	 1.0 
    varphi2 = 	 1.0 
    delta_rsum_ep = 10^(-0)
    delta_rsum_lag = 2
    delta_occymean_ep = 10^(0)#10^(-6)
    delta_occymean_lag = 2

    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run, datafilename1,dataset_name,outdir,script_name,slurm_job_id,KMax,significance_prop,min_percent_of_genes,max_percent_of_genes,seed,elbo_ep,num_iter,m_initialization_approach,n_hvgs,time_key,individuals_key,change_seeds,add_time,add_individuals, layer_index,update_clusterwise,  check_cluster_interpretability_bool, gene_set_file_path, samplebased_alpha0,elbo_sign_decrease_max_tolerance,report_creation_script_path,make_abridged_bool_as_string);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,K,train_h,train_w,train_ab,train_uv,train_m_nu,train_s_sq_nu,uv_scalar,multiple_s_sq_updates,use_std_for_hvg,center_data_cols,scale_data_cols,min_percent_cells,min_genes_detected_for_mask,standardization_of_used_representation,rand_init_inputs,uniform_theta_init,iseeds_init,m_mu_init_init,num_samples_posterior_samples,s_value_thresh,sparsity_based_on_used_representation_dim,eta_update_mode,sigma_update_mode,lambda_update_mode,burnin,nu0,sigma_sq_nu,min_number_cells,remove_small_clusters);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,num_var_feat,unique_time_id,dataset_used_id,experiment_id,delta_rsum_ep,delta_rsum_lag,delta_occymean_ep,delta_occymean_lag);

    
    w1_init= nothing;
    w2_init = nothing;
    
    h1_init = nothing;#[ 1 /(1e6 + 1) ];#
    h2_init = nothing;#[ 1e6 /(1e6 + 1) ];#
    a_init = nothing;#ones(size(used_representation)[1],K);
    b_init = nothing;#ones(size(used_representation)[1],K);
    u_init = nothing;#ones(size(used_representation)[1],K);
    v_init = nothing;#ones(size(used_representation)[1],K);
    # y_init = 1/J ./ ones(size(used_representation)[1],K);
    r_init = nothing;
    y_init =  nothing;#ones(size(used_representation)[1],K); #<--- use for gene expression space
    # m_initialization_approach = "rand"#"ssss"#"kpp"#
    # RMin = Clustering.kmeans(used_representation,KMin);
    # RMax = Clustering.kmeans(used_representation,KMax,init=:kmcen);
    # m_mu_init_init = RMax.centers;#nothing#
    # m_mu_init_K = KMax;#25;#
    # RKmcen = kmeans(used_representation,KMin,init=:kmcen);
    # reference_lables = RKmcen.assignments;
    # no_itr_ = 1;
    # r_init_sim = zeros(K+1,N,no_itr_);
    # for itr_ in 1:no_itr_
    #     RCurr = Clustering.kmeans(used_representation,KMin,init=:kmcen);
    #     curr_assignments = RCurr.assignments;
    #     aligned_curr_assignments = relabel_clusters(curr_assignments, RMax.assignments, KMax)
    #     for k in 1:KMax
    #         r_init_sim[k,aligned_curr_assignments .== k,itr_] .= 1.0;# r_init_sim[k,RMax.assignments .== k] .= 1.0;#
    #     end
    # end


    r_init_string="";
    if isnothing(r_init)
        r_init_string = "r_init initialized to nothing";
    else
        r_init_string = "r_init initialized to something, too large to display here, so check the value";
    end
    if isnothing(a_init) && isnothing(b_init)
        ab_init_string = "a_init and b_init initialized to nothing";
    else
        ab_init_string = "a_init and b_init initialized to something, too large to display here, so check the value";
    end
    if isnothing(u_init) && isnothing(v_init)
        uv_init_string = "u_init and v_init initialized to nothing";
    else
        uv_init_string = "u_init and v_init initialized to something, too large to display here, so check the value";
    end
    m_mu_init,iseeds = get_m_init_from_initialization_approach(m_initialization_approach,used_representation,cell_cluster_labels,m_mu_init_K,seed;iseeds=iseeds_init,m_mu_init=m_mu_init_init);
    if !isnothing(m_mu_init)
        if size(m_mu_init)[2] != K
            m_mu_init = hcat(m_mu_init,randn(size(m_mu_init)[1],K-size(m_mu_init)[2]));
        end
    end

    nu0 = original_mean;#vec(mean(used_representation,dims=2));#zeros(length(median(used_representation,dims=2)));#vec(median(used_representation,dims=2));#
    inputs = initialize_model_parameters(data_input,K,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,seed;num_iter=num_iter,rand_init = rand_init_inputs,change_seeds = change_seeds,uniform_theta_init = uniform_theta_init, m_mu_init = m_mu_init,samplebased_alpha0=samplebased_alpha0,update_clusterwise=update_clusterwise,y_init=y_init,r_init = r_init,h1_init=h1_init,h2_init=h2_init,w1_init=w1_init,w2_init=w2_init,a_init=a_init,b_init=b_init,u_init=u_init,v_init=v_init,eta_update_mode=eta_update_mode,sigma_update_mode=sigma_update_mode,lambda_update_mode=lambda_update_mode,train_h=train_h,train_w=train_w,train_ab=train_ab,train_uv=train_uv,m_nu_init=[el for el in nu0]);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,iseeds,r_init_string,h1_init,h2_init,w1_init,w2_init,ab_init_string,uv_init_string,define_global_sparsity,define_global_autocorrelation);
    filepath = mk_outputs_filepath(outdir,experiment_id,dataset_used_id,unique_time_id);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,filepath);
    mk_outputs_pathname(filepath);
    notes_="$(data_structure_note)__*||*__$(layer_name)"
    summary_file = saving_summary_file(filepath;slurm_job_id=slurm_job_id,script_name=script_name,change_seeds = change_seeds,unique_time_id=unique_time_id,datafilename1=datafilename1,KMax=KMax, seed=seed,num_var_feat=num_var_feat,N=N,elbo_ep=elbo_ep,notes_=notes_,alpha0 = alpha0,gamma0 =gamma0 ,phi1 = phi1,phi2 = phi2,kappa1 = kappa1,kappa2 = kappa2,xi1 =xi1,xi2 = xi2,varphi1 = varphi1,varphi2 = varphi2,significance_prop=significance_prop,min_percent_cells=min_percent_cells,min_number_cells=min_number_cells,min_percent_of_genes=min_percent_of_genes,max_percent_of_genes=max_percent_of_genes,num_iter = num_iter,dataset_name = dataset_name,outdir = outdir,time_key = time_key,individuals_key = individuals_key,m_initialization_approach = m_initialization_approach,update_clusterwise=update_clusterwise, use_alt_representation=layer_name, check_cluster_interpretability_bool=check_cluster_interpretability_bool, gene_set_file_path=gene_set_file_path, samplebased_alpha0=samplebased_alpha0,);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,inputs[:dataparams].I,inputs[:dataparams].T,inputs[:dataparams].J,inputs[:dataparams].N_t,inputs[:dataparams].N,inputs[:dataparams].Jlog,inputs[:dataparams].logpi,inputs[:modelparams].K,inputs[:modelparams].alpha0,inputs[:modelparams].gamma0,inputs[:modelparams].phi1,inputs[:modelparams].phi2,inputs[:modelparams].kappa1,inputs[:modelparams].kappa2,inputs[:modelparams].xi1,inputs[:modelparams].xi2,inputs[:modelparams].varphi1,inputs[:modelparams].varphi2,inputs[:modelparams].nu0,inputs[:modelparams].sigma_sq_nu,inputs[:modelparams].significance_prop,inputs[:modelparams].min_number_cells,inputs[:modelparams].min_percent_cells,inputs[:modelparams].min_percent_of_genes,inputs[:modelparams].max_percent_of_genes,inputs[:modelparams].num_iter,inputs[:modelparams].uniform_theta_init,inputs[:modelparams].rand_init,inputs[:modelparams].change_seeds,inputs[:modelparams].init_seed);
    training_features_dict = nothing;
    #testting_convergence_monitoring(inputs;num_iter=100,use_log=true,eta_update_mode=eta_update_mode,lambda_update_mode=lambda_update_mode,sigma_update_mode=sigma_update_mode)
    elapsed_time = @elapsed begin
        st = time();
        outputs_dict,training_features_dict = cavi(inputs;delta_rsum_ep =delta_rsum_ep,delta_rsum_lag = delta_rsum_lag,delta_occymean_ep = delta_occymean_ep,delta_occymean_lag = delta_occymean_lag,logger=logger,update_clusterwise=update_clusterwise,elbo_sign_change_max=elbo_sign_decrease_max_tolerance,check_cluster_interpretability_bool=check_cluster_interpretability_bool,train_h =train_h,train_w=train_w,train_ab=train_ab,train_uv=train_uv,multiple_s_sq_updates=multiple_s_sq_updates,eta_update_mode=eta_update_mode,sigma_update_mode=sigma_update_mode,lambda_update_mode=lambda_update_mode,burnin=burnin,remove_small_clusters=remove_small_clusters,use_log=true,train_m_nu=train_m_nu,train_s_sq_nu=train_s_sq_nu);
    end
    dt = time() - st;
    _flushed_logger("\t \t Finished Training Model. Model took $dt seconds to run...";logger)
    _flushed_logger("\t \t Final ELBO $(outputs_dict[:LB][end])...";logger)
    _flushed_logger("Model Took a total of $elapsed_time seconds to run";logger)
    outputs_dict[:elapsed_time]=elapsed_time;
    inferred_labels = [argmax(outputs_dict[:r_][i][t][n]) for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]]
    countmap(inferred_labels)
    num_clust = length(unique(inferred_labels))#sum(outputs_dict[:Nk_] .>= min_number_cells);
    if !isnothing(cell_cluster_labels)
        println("ARI: $(Clustering.randindex(cell_cluster_labels,inferred_labels)[1])");
    end
    final_elbo = outputs_dict[:LB][end];
    truncation_value = outputs_dict[:truncation_value];
    is_converged = outputs_dict[:is_converged];
    seed_used = outputs_dict[:seed_used];
    split_seed_used_vec = split(seed_used,", ");
    seed_used_vec = parse.(Int64,split_seed_used_vec[broadcast(!,isempty.(split_seed_used_vec))]);
    outputs_dict[:seed_used_vec] = seed_used_vec;
    reason_for_non_convergence = outputs_dict[:reason_for_non_convergence];
    _flushed_logger( "Number of Cluster $(num_clust)";logger)
    _flushed_logger("Appending to Quick Summary";logger)
    vars = [seed_used,num_clust,final_elbo,elapsed_time,truncation_value,is_converged, reason_for_non_convergence];
    varnames = ["Seed Used During Training", "Number of Cluster","Final ELBO","Elapsed Time","Number of Iterations Run","Convergence Reached", "Reason for Non-Convergence"];
    append_summary(summary_file,vars,varnames);
    outputs_dict[:elapsed_time]=elapsed_time;
    outputs_dict[:unique_time_id] = unique_time_id;
    outputs_dict[:filepath] = filepath;
    outputs_dict[:experiment_id] = experiment_id;
    outputs_dict[:dataset_used_id] = dataset_used_id;
    outputs_dict[:datafilename1] = datafilename1;
    outputs_dict[:KMax] = KMax;
    outputs_dict[:seed] = seed;
    outputs_dict[:gene_names] = gene_names;
    outputs_dict[:cell_ids] = cell_ids;
    outputs_dict[:time_vec] = time_vec;
    outputs_dict[:individuals_vec] = individuals_vec;
    if !isnothing(cell_cluster_labels)
        outputs_dict[:cell_cluster_labels] = cell_cluster_labels;
    else
        outputs_dict[:cell_cluster_labels] = "nan";
    end
    outputs_dict[:significance_prop] = significance_prop;
    outputs_dict[:min_number_cells] = min_number_cells;
    outputs_dict[:min_percent_cells] = min_percent_cells;
    outputs_dict[:min_percent_of_genes] = min_percent_of_genes;
    outputs_dict[:max_percent_of_genes] = max_percent_of_genes;
    outputs_dict[:alpha0] = alpha0;
    outputs_dict[:gamma0] = gamma0;
    outputs_dict[:phi1] = phi1;
    outputs_dict[:phi2] = phi2;
    outputs_dict[:kappa1] = kappa1;
    outputs_dict[:kappa2] = kappa2;
    outputs_dict[:xi1] = xi1;
    outputs_dict[:xi2] = xi2;
    outputs_dict[:varphi1] = varphi1;
    outputs_dict[:varphi2] = varphi2;
    outputs_dict[:num_var_feat] = num_var_feat;
    # outputs_dict[:N] = N
    outputs_dict[:elbo_ep] = elbo_ep;
    # outputs_dict[:notes_] = notes_
    outputs_dict[:num_iter] = num_iter;
    outputs_dict[:dataset_name] = dataset_name;
    outputs_dict[:outdir] = outdir;
    outputs_dict[:time_key] = time_key;
    outputs_dict[:m_initialization_approach] = m_initialization_approach;
    outputs_dict[:layer_index] = layer_index;
    # outputs_dict[:use_alt_representation] = use_alt_representation;    
    # Nkplus1 = sum(sum.(r))
    # outputs_dict[:Nkplus1_] = Nkplus1
    outputs_dict[:used_representation_feature_name] = used_representation_feature_name;
    outputs_dict[:individuals_key] = individuals_key;
    r = outputs_dict[:r_];
    # [argmax(r[i][t][n]) for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]];
    z_argmax =     deepcopy(inferred_labels);
    outputs_dict[:z_argmax] = [z_argmax];
    outputs_dict = subset_on_occupied_clusters(outputs_dict);
    outputs_dict = size_reorder_clusters(outputs_dict);
    outputs_dict[:sorting_keys] = sorted_sorting_keys;
    outputs_dict[:new_order_samples] = new_order_samples;
    outputs_dict[:new_order_features] = new_order_features;
    outputs_dict = save_embeddings(anndata_dict1,filepath;outputs_dict = outputs_dict, logger = logger,unique_time_id=unique_time_id,new_order_samples=new_order_samples);
    outputs_dict[:gene_factor_loadings] = nothing;
    if !isnothing(alternative_representation) && !isnothing(alternative_representation[1])
        outputs_dict[:gene_factor_loadings] = alternative_representation[1];
        representation_filename = "$filepath/used_representation.csv";
        used_representation_df = DataFrame(permutedims(used_representation),:auto);
        rename!(used_representation_df,used_representation_feature_name);
        CSV.write(representation_filename,used_representation_df);
    end
    outputs_dict[:parameter_settings_to_reproduce_run] = parameter_settings_to_reproduce_run;
    # if !isnothing(alternative_representation)
    #     outputs_dict = make_new_embeddings(outputs_dict,used_representation,used_representation_feature_name,seed; samples_as_rows=false);
    # end
    # posterior_summaries(outputs_dict, anndata_dict1,seed;num_samples = 1000,s_value_thresh=0.05)
    flushed_logger("Calculating Posterior Summaries...";logger)
    temp_outputs_dict,temp_z_argmax,temp_posterior_gene_summaries_df,temp_cluster_summaries_df,temp_cell_summaries_df = posterior_summaries(outputs_dict, used_representation,used_representation_feature_name,seed;num_samples = num_samples_posterior_samples,s_value_thresh=s_value_thresh,return_data_frames = false,save_data_frames = true)
    if !isnothing(temp_outputs_dict)
        outputs_dict = temp_outputs_dict
    end
    if !isnothing(temp_z_argmax)
        z_argmax = temp_z_argmax
        outputs_dict[:z_argmax] = [z_argmax];
    end
    filename = "$filepath/output.jld2"
    flushed_logger("Saving Outputs...";logger)
    jldsave(filename,true;outputs_dict=outputs_dict)
    trainingfilename = "$filepath/train_history.jld2"
    if !isempty(training_features_dict)
        flushed_logger("Saving Training Features...";logger)
        jldsave(trainingfilename,true;outputs_dict=training_features_dict)
    end
    runParameterSettingsfilename = "$filepath/_RunParameterSettings_$unique_time_id.csv"
    runParameterSettingsdf = DataFrame(parameter_name=parameter_settings_to_reproduce_run.keys,parameter_value=[parameter_settings_to_reproduce_run[key] for key in parameter_settings_to_reproduce_run.keys]);
    CSV.write(runParameterSettingsfilename,runParameterSettingsdf);
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("Attempting to create a report using $report_creation_script_path...";logger)
    # report_creation_status = generate_report(report_creation_script_path, String(chop(filepath)), datafilename1, make_abridged_bool_as_string);
    report_creation_status = submit_slurm_job_to_generate_report(report_creation_script_path, String(chop(filepath)), datafilename1, make_abridged_bool_as_string);
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("Report Creation Status:\n\n\n $report_creation_status \n\n\n";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("**********************************************************************";logger)
    flushed_logger("Finishing Script...";logger)
end
# main(ARGS;logger=logger)
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS;logger=logger)
end