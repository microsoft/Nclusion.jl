
"""
        recursive_flatten(x::AbstractArray)
    This function takes an arbitrarily nested set of vectors and recursively flattens them into one 1-D vector
    ```math

    ```
"""
function recursive_flatten(x::AbstractArray)
    if any(a->typeof(a)<:AbstractArray, x)#eltype(x) <: Vector
        recursive_flatten(vcat(x...))
    else
        return x
    end
end

"""
        outermelt(val,num_repeats)
    This function recursively performs an outer melt of a vector input
"""
function outermelt(val,num_repeats)
    melt = nothing
    if typeof(val) <: Vector && eltype(val) <: Number
        melt = repeat(val, outer = num_repeats)
    elseif typeof(val) <: Number
        val = [val]
        melt = repeat(val, outer = num_repeats)
    end
    return melt
end

"""
       innermelt(val,num_repeats)
    This function recursively performs an intter melt of a vector input
"""
function innermelt(val,num_repeats) 
    melt = nothing
    if typeof(num_repeats) <: Vector
        # println("Condition 1")
        melt = innermelt.(val,num_repeats)
        melt = recursive_flatten(melt)
    else
        if typeof(val) <: Vector && eltype(val) <: Number
            # println("Condition 2")
            melt = repeat(val, inner = num_repeats)
        elseif typeof(val) <: Number
            # println("Condition 3")
            val = [val]
            melt = repeat(val, inner = num_repeats)
        end
    end
    return melt
end

"""
    name(arg...)
This macro turns a string (or list of strings) into a symbol type.
"""    
macro name(arg...)
    x = string(arg)
    quote
        $x
    end
end

"""
    name(arg)
    This macro turns a string (or list of strings) into a symbol type.
"""      
macro name(arg)
    x = string(arg)
    quote
        $x
    end
end

"""
    naming_vec(arg_str_list)
    This function parses the string on commas (,) into a list of strings with a colon appended to the front.
```math

```
"""  
function naming_vec(arg_str_list)
    arg_str_list_trunc = chop(arg_str_list,head=1);
    arg_str_vec = split(arg_str_list_trunc,", ");
    num_var = length(arg_str_vec)
    str_var_vec = Vector{String}(undef, num_var)
    for i in 1:num_var
        el = arg_str_vec[i]
        if el[1] == ':'
            str_var_vec[i] = el[2:end] 
        else
            str_var_vec[i] = el[1:end] 
        end
    
    end
    return str_var_vec
end

"""
    addToDict!(dict,key_array,val_array)
Adds a set of values to a previously initialized dictionary
"""  
function addToDict!(dict,key_array,val_array)
    num_var = length(key_array)
    for i in 1:num_var
        key = key_array[i]
        val = val_array[i]
        dict[key] = val 
    end
    dict
end

"""
    addToOrderedDict!(ordered_dict,key_array,val_array)
    Adds a set of values to a previously initialized ordered dictionary
"""      
function addToOrderedDict!(ordered_dict,key_array,val_array)
    num_var = length(key_array)
    for i in 1:num_var
        key = key_array[i]
        val = val_array[i]
        ordered_dict[key] = val 
    end
    ordered_dict
end

"""
    get_unique_time_id()
This function generates a unique ID based on the current system date and time.

"""  
function get_unique_time_id()
    datetimenow = Dates.now(Dates.UTC)
    now_str = string(datetimenow)
    now_str = string(split(now_str, ".")[1])
    r = ":"
    return replace(now_str,r => "" )
end


"""
    setup_experiment_tag(experiment_filename)
Creates an experiment tag 
```math

```
"""      
function setup_experiment_tag(experiment_filename)
    return "EXPERIMENT_$experiment_filename"
end
function load_data(datafilename1,seed)
    fid1 = h5open(datafilename1,"r")
    anndata_dict1 = read(fid1)
    return anndata_dict1
end

function preparing_data(anndata_dict1)
    gene_names = anndata_dict1["var"]["_index"]
    cell_ids = anndata_dict1["obs"]["_index"]
    cell_cluster_labels = anndata_dict1["obs"]["cell_type"]["codes"]
    cell_cluster_labels = cell_cluster_labels .+ 1

    return gene_names, cell_ids, cell_cluster_labels
end


function make_nclusion_inputs(anndata_dict1;T=1)
    x_mat = anndata_dict1["X"]
    _, _, cell_cluster_labels=preparing_data(anndata_dict1)
    N = size(x_mat)[2]
    G = size(x_mat)[1]
    x_input = Vector{Vector{Vector{Float64}}}(undef,T)
    z = Vector{Vector{Int}}(undef,T)
    C_t = Vector{Float64}(undef,T)
    C_t[1] = N
    x_input[1] = [Float64.(collect(col)) for col in eachcol(x_mat)]
    z[1] = Int.(collect(cell_cluster_labels))
    return x_input,z
end

function select_cells_hvgs(x_mat,num_var_feat,num_cnts,gene_ids,cell_cluster_labels,scale_factor,N;T=1,chosen_cells=nothing)
    cell_intersect_bool = nothing 
    if !isnothing(chosen_cells)
        cell_intersect = Set(chosen_cells)
        cell_intersect_bool = [in(el,cell_intersect) for  el in collect(1:N)]
    else
        cell_intersect_bool = [true for  el in collect(1:N)]
    end
    x_temp = Vector{Vector{Vector{Float64}}}(undef,T)
    x_temp[1] = [Float64.(collect(col)) for col in eachcol(x_mat[:,cell_intersect_bool])]
    numi_temp = Vector{Vector{Float64}}(undef,T)
    numi_temp[1] = Float64.(collect(num_cnts))
    log_norm_x = lognormalization(x_temp;scaling_factor=scale_factor,pseudocount=1.0,numi=numi_temp)
    log_norm_xmat =  hcat(vcat(log_norm_x...)...)
    gene_std_vec = vec(std(log_norm_xmat, dims=2))
    sorted_indx  = sortperm(gene_std_vec,rev=true)
    top_genes = gene_ids[sorted_indx][1:num_var_feat]
    top_genes_bool = [in(el,Set(top_genes)) for  el in gene_ids]


    x = Vector{Vector{Vector{Float64}}}(undef,T)
    z_true = Vector{Vector{Int}}(undef,T)
    C_t = Vector{Float64}(undef,T)
    numi = Vector{Vector{Float64}}(undef,T)

    C_t[1] = sum(cell_intersect_bool)
    x[1] = [Float64.(collect(col)) for col in eachcol(x_mat[top_genes_bool,cell_intersect_bool])]
    z_true[1] = Int.(collect(cell_cluster_labels))[cell_intersect_bool]
    numi[1] = Float64.(collect(num_cnts))[cell_intersect_bool]
    return x,z_true,numi,top_genes,C_t
end


function initialize_model_parameters(x_input,KMax,alpha1,gamma1;sparsity_lvl=nothing,phi1=1.0,num_iter=1000,rand_init = false, uniform_theta_init = true,mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing)
    if typeof(KMax) <: AbstractFloat
        KMax = Int(round(KMax))
    end
    T = length(x_input)
    G = length(x_input[1][1])
    C_t = [length(el) for el in x_input]
    N = sum(C_t)
    if !isnothing(sparsity_lvl)
        ηk = sparsity_lvl
    else
        ηk = 1/G
    end

    α0,γ0,ϕ0 = alpha1,gamma1,phi1
    mk_hat_init = init_mk_hat!(mk_hat_init,x_input,KMax,G;rand_init = rand_init);
    v_sq_k_hat_init = init_v_sq_k_hat_vec!(v_sq_k_hat_init,KMax,G;rand_init = rand_init, lo=0,hi=1);
    λ_sq_init = init_λ_sq_vec!(λ_sq_init,G;rand_init = rand_init, lo=0,hi=1) ;
    σ_sq_k_init = init_σ_sq_k_vec!(σ_sq_k_init,KMax,G;rand_init = rand_init, lo=0,hi=1);
    gk_hat_init,hk_hat_init = init_ghk_hat_vec!(gk_hat_init,hk_hat_init,KMax;rand_init = rand_init, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2);
    st_hat_init = init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    d_hat_init = init_d_hat_vec!(d_hat_init,KMax,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, gk_hat_init = gk_hat_init, hk_hat_init= hk_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,KMax,T,C_t;rand_init = rand_init)
    yjk_init = init_yjk_vec!(yjk_init,G,KMax;rand_init = rand_init)
    float_type=eltype(x_input[1][1])
    Tk = Vector{Float64}(undef,KMax+1);
    
    cellpop = [CellFeatures(t,i,KMax,x_input[t][i]) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:KMax];
    dataparams = DataFeatures(x_input);
    conditionparams = [ConditionFeatures(t,KMax,T;float_type=float_type) for t in 1:T];
    geneparams = [GeneFeatures(j) for j in 1:G];
    modelparams = ModelParameterFeatures(x_input,KMax,ηk,α0,γ0,ϕ0,num_iter,uniform_theta_init,rand_init);
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init);

    input_str_list = @name cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,Tk;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,Tk];
    inputs = OrderedDict()
    addToDict!(inputs,input_key_list,input_var_list);

    
    return inputs
end


function run_cavi(inputs;elbo_ep = 10^(-0),logger=nothing)
    num_iter = inputs[:modelparams].num_iter
    KMax = inputs[:modelparams].K
    elbologger = ElboFeatures(1,KMax,num_iter) 

    _flushed_logger("Maximum KMax initialized at $KMax";logger)

    inputs[:elbolog] = elbologger
    elapsed_time = @elapsed begin
        _flushed_logger("\t Model running now...";logger)
        st = time()
        outputs_dict = cavi(inputs;elbo_ep = elbo_ep);
    end
    dt = time() - st
    elbo_, rtik_, yjk_hat_, mk_hat_, v_sq_k_hat_, σ_sq_k_hat_, var_muk_, Nk_, gk_hat_, hk_hat_, ak_hat_, bk_hat_,x_hat_, x_hat_sq_, d_hat_t_, c_tt_prime_, st_hat_, λ_sq_,per_k_elbo_,ηk_, Tk_, is_converged, truncation_value, ηk_trend_ = (; outputs_dict...);
    _flushed_logger("\t \t Finished Training Model. Model took $dt seconds to run...";logger)
    _flushed_logger("\t \t Final ELBO $(elbo_[end])...";logger)
    _flushed_logger("Model Took a total of $elapsed_time seconds to run";logger)
    outputs_dict[:elapsed_time]=elapsed_time
    return outputs_dict
end

function save_pips(pip,gene_names;unique_time_id="",filepath="")
    KMax =length(pip)
    G = length(pip[1])
    pip_mat = permutedims(hcat(pip...))
    pip_mat = hcat(["Cluster_$el" for el in 1:KMax],pip_mat)
    col_names = vcat("cluster_id",gene_names)
    pip_df  = DataFrame(pip_mat, :auto);
    rename!(pip_df,Symbol.(col_names));
    CSV.write(filepath*"$(G)G-"*unique_time_id*"-pips.csv",  pip_df)

end

function save_Nk(rtik_;unique_time_id="",filepath="")
    KMax = length(rtik_[1][1])
    N = length(rtik_[1])
    T = length(rtik_)
    cluster_name = ["Cluster_$el" for el in 1:KMax]
    suffix = "-Nk.csv"
    Nk = sum(sum.(rtik_))
    filename = filepath*unique_time_id*suffix
    Nk_mat = hcat(cluster_name,Nk)
    Nk_df  = DataFrame(Nk_mat, :auto);
    col_names = vcat("cluster_ids","Nk")
    rename!(Nk_df,Symbol.(col_names));
    CSV.write(filename,  Nk_df)

end

function summarize_parameters(outputs_dict_vec,elapsed_time,final_elbo_vec,elbo_vec,rtik_vec,yjk_vec,perK_elbo_vec,delta_t_vec,nk_perL_vec,ηk_vec)
    L = length(outputs_dict_vec)
    importance_weights = norm_weights(final_elbo_vec)
    weighted_yjk_vec = importance_weights .* yjk_vec;
    weighted_rtik_vec = importance_weights .* rtik_vec;
    weighted_elbo_vec = importance_weights .* elbo_vec;
    weighted_elapsed_time = importance_weights .* delta_t_vec;
    mk_hat_L = [el[:mk_hat_] for el in outputs_dict_vec];
    v_sq_k_hat_L = [el[:v_sq_k_hat_] for el in outputs_dict_vec];
    σ_sq_k_hat_L = [el[:σ_sq_k_hat_] for el in outputs_dict_vec];
    Nk_L = [el[:Nk_] for el in outputs_dict_vec];
    d_hat_t_L = [el[:d_hat_t_] for el in outputs_dict_vec];
    c_tt_prime_L = [el[:c_tt_prime_] for el in outputs_dict_vec];
    st_hat_L = [el[:st_hat_] for el in outputs_dict_vec];
    λ_sq_L = [el[:λ_sq_] for el in outputs_dict_vec];
    weighted_mk_vec = importance_weights .* mk_hat_L;
    weighted_v_sq_k_vec = importance_weights .* v_sq_k_hat_L;
    weighted_σ_sq_k_vec = importance_weights .* σ_sq_k_hat_L;
    weighted_Nk_vec = importance_weights .* Nk_L;
    weighted_d_vec = importance_weights .* d_hat_t_L;
    weighted_c_tt_prime_vec = importance_weights .* c_tt_prime_L;
    weighted_st_vec = importance_weights .* st_hat_L;
    weighted_λ_sq_vec = importance_weights .* λ_sq_L;
    old_weighted_elbo_vec = deepcopy(weighted_elbo_vec)
    maxLen = maximum(length.(old_weighted_elbo_vec))
    for l in 1:L
        curr_len = length(old_weighted_elbo_vec[l])
        diff_ = maxLen - curr_len
        if !iszero(diff_)
            weighted_elbo_vec[l] = vcat(old_weighted_elbo_vec[l],zeros(diff_))
        end
    end
    mean_elbo = sum(weighted_elbo_vec)
    pip = sum(weighted_yjk_vec)
    mean_rtik = sum(weighted_rtik_vec)
    mean_mk= sum(weighted_mk_vec)
    mean_v_sq_k = sum(weighted_v_sq_k_vec)
    mean_σ_sq_k = sum(weighted_σ_sq_k_vec)
    mean_Nk = sum(weighted_Nk_vec)
    mean_d = sum(weighted_d_vec)
    mean_c_tt_prime = sum(weighted_c_tt_prime_vec)
    mean_st = sum(weighted_st_vec)
    mean_λ_sq = sum(weighted_λ_sq_vec) 


    return mean_elbo,pip,mean_rtik,mean_mk,mean_v_sq_k,mean_σ_sq_k,mean_Nk,mean_d,mean_c_tt_prime,mean_st,mean_λ_sq
end
function _flushed_logger(msg;logger=nothing)
    if !isnothing(logger)
        with_logger(logger) do
            @info msg
        end
    end
end

function make_ids(dataset,G,N)
    unique_time_id = get_unique_time_id()
    dataset_used_id = "$(dataset)_$(G)HVGs-$(N)N"
    experiment_filename = "nclusion_$(dataset)"
    experiment_id = setup_experiment_tag(experiment_filename)

    return unique_time_id,dataset_used_id,experiment_id
end

function mk_outputs_filepath(outdir,experiment_id,dataset_used,unique_time_id)
    filepath ="$outdir/outputs/$experiment_id/DATASET_$dataset_used/$unique_time_id/"
    return filepath
end
function mk_outputs_pathname(filepath)
    mkpath(filepath)
end

function saving_summary_file(filepath;unique_time_id="",datafilename1="",alpha1="",gamma1="",KMax="", seed="",num_var_feat="",N="",elbo_ep="",notes_="")
    
    summary_file = filepath*"_QuickSummary_"*unique_time_id*".txt"
    vars = [datafilename1,alpha1,gamma1,KMax, seed, "mpu (Output from maddy;s pipeline)", num_var_feat,N,true,false,false,false,elbo_ep,notes_]
    varnames = ["Filename","alpha","gamma","KMax","Seed","Which Method did I use to select HVGs", "How many Genes were in this anaysis", "How many Cells were in this anaysis","(True/False) I used all cells in this anaysis",  "(True/False) I manually had to standardize the data and did not use the steps in the QC pipeline","(True/False) I started with raw counts","(True/False) I used all genes","Elbo intolerance threshold","Notes"]



    # results = h5open(results_filename, "w")

    # datafilename1 = "/mnt/e/cnwizu/Playground/SCoOP-sc/data/pbmc/labelled_cells/pure_pbmc/pure_pbmc_preprocessed.h5ad"

    # @info "Saving Quick Summary..."
    
    open(summary_file, "w") do f
        for i in eachindex(vars)
            write(f,"$(varnames[i]) = \t $(vars[i]) \n")
        end
    end
    return summary_file
end
function save_embeddings(anndata_dict1,filepath;logger = nothing,unique_time_id="")
    G = size(anndata_dict1["X"])[1]
    N = size(anndata_dict1["X"])[2]
    tsne_data = nothing
    pca_data = nothing
    umap_data = nothing
    # @info "Getting TSNE Transform..."
    _flushed_logger("Getting TSNE Transform...";logger)
    
    if haskey(anndata_dict1["obsm"],"X_tsne")
        tsne_data =  permutedims(anndata_dict1["obsm"]["X_tsne"])
        tsne_data = tsne_data
        tsne_data_df  = DataFrame(tsne_data, :auto);
        ncols = size(tsne_data)[2]
        rename!(tsne_data_df,Symbol.(["TSNE_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-tsne_coordinates.csv",  tsne_data_df)
    elseif N <= 5000
        tsne_data = get_tsne_transform(x_input);
        tsne_data_df  = DataFrame(tsne_data, :auto);
        ncols = size(tsne_data)[2]
        rename!(tsne_data_df,Symbol.(["TSNE_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-tsne_coordinates.csv",  tsne_data_df)
    end

    
    _flushed_logger("Getting PCA Transform...";logger)
    if haskey(anndata_dict1["obsm"],"X_pca")
        pca_data =  permutedims(anndata_dict1["obsm"]["X_pca"])
        pca_data = pca_data
        pca_data_df  = DataFrame(pca_data, :auto);
        ncols = size(pca_data)[2]
        rename!(pca_data_df,Symbol.(["PC_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-pca_coordinates.csv",  pca_data_df)
    end


    _flushed_logger("Getting UMAP Transform..."; logger)
    
    if haskey(anndata_dict1["obsm"],"X_umap")
        umap_data =  permutedims(anndata_dict1["obsm"]["X_umap"])
        umap_data = umap_data
        umap_data_df  = DataFrame(umap_data, :auto);
        ncols = size(umap_data)[2]
        rename!(umap_data_df,Symbol.(["UMAP_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-umap_coordinates.csv",  umap_data_df)
    end

end

function make_labels(x_input,anndata_dict1,z_argmax)
    T = length(x_input)
    timepoint_map = [t * ones(Int,Int(length(x_input[t]))) for t in 1:T]
    timepoint_map = recursive_flatten(timepoint_map)
    # cluster_remap = Dict(v => k for (k,v) in cluster_map)
    cell_ids = anndata_dict1["obs"]["_index"]
    cell_cluster_labels = anndata_dict1["obs"]["cell_type"]["codes"]
    cell_cluster_labels = cell_cluster_labels .+ 1
    z = Vector{Vector{Int}}(undef,T)
    z[1] = Int.(collect(cell_cluster_labels))
    cluster_map = OrderedDict(k => v for (k,v) in enumerate(anndata_dict1["obs"]["cell_type"]["categories"]))
    cell_type_vec = [cluster_map[el] for el in recursive_flatten(z)]

    cluster_results_df = DataFrame(condtion=timepoint_map,cell_id = cell_ids,cell_type = cell_type_vec, called_label = recursive_flatten(z), inferred_label = recursive_flatten(z_argmax))
    
    return cluster_results_df
    
end
function save_labels(cluster_results_df;dataset_used="",G="",unique_time_id="",filepath="")
    CSV.write(filepath*"$(dataset_used)_nclusion-"*unique_time_id*".csv",  cluster_results_df)
end
function run_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = nothing)
    Random.seed!(seed)
    _flushed_logger("Loading data and metadata...";logger)
    anndata_dict1= load_data(datafilename1,seed)
    N = size(anndata_dict1["X"])[2]
    G = size(anndata_dict1["X"])[1]
    num_var_feat = G
    _flushed_logger("Number of cells: $(N) ...";logger)
    _flushed_logger("Number of Genes: $(G) ...";logger)
    unique_time_id,dataset_used_id,experiment_id = make_ids(dataset,G,N)

    _flushed_logger("Preparing Dataset ...";logger)
    gene_names, cell_ids, cell_cluster_labels = preparing_data(anndata_dict1)
    x_input,z = make_nclusion_inputs(anndata_dict1)
    

    filepath = mk_outputs_filepath(outdir,experiment_id,dataset_used_id,unique_time_id)
    _flushed_logger("Preparing saving Directory at $filepath ...";logger)
    mk_outputs_pathname(filepath)


    _flushed_logger("Saving Quick Summary...";logger)
    notes_=""
    summary_file = saving_summary_file(filepath;unique_time_id=unique_time_id,datafilename1=datafilename1,alpha1=alpha1,gamma1=gamma1,KMax=KMax, seed=seed,num_var_feat=num_var_feat,N=N,elbo_ep=elbo_ep,notes_=notes_)
    save_embeddings(anndata_dict1,filepath;logger = logger,unique_time_id=unique_time_id)

    

    _flushed_logger("Initializing Model parameters...";logger)
    inputs = initialize_model_parameters(x_input,KMax,alpha1,gamma1);

    _flushed_logger("Starting Variational Inference";logger)
    outputs_dict = run_cavi(inputs;elbo_ep=elbo_ep,logger=logger)

    
    elbo_, rtik_, pip, _ = (; outputs_dict...);

    z_argmax = [argmax.(el) for el in  rtik_];
    cluster_results_df = make_labels(x_input,anndata_dict1,z_argmax)
    
    outputs_dict[:z_argmax] = z_argmax
    outputs_dict[:cluster_results_df] = cluster_results_df
    outputs_dict[:unique_time_id] = unique_time_id
    outputs_dict[:filepath] = filepath
    

    num_clust = length(unique(recursive_flatten(z_argmax)))
    final_elbo = elbo_[end]
    elapsed_time = outputs_dict[:elapsed_time]

    _flushed_logger( "Number of Cluster $(num_clust)";logger)


    _flushed_logger("Appending to Quick Summary";logger)
    vars = [num_clust,final_elbo,elapsed_time]
    varnames = ["Number of Cluster","Final ELBO","Elapsed Time"]
    append_summary(summary_file,vars,varnames)

    _flushed_logger("Saving PIPs";logger)
    save_pips(pip,gene_names;unique_time_id=unique_time_id,filepath=filepath)
    _flushed_logger("Saving Nk";logger)
    save_Nk(rtik_;unique_time_id=unique_time_id,filepath=filepath)

    _flushed_logger("Saving Cluster Memberships";logger)
    save_labels(cluster_results_df;dataset_used=dataset_used_id,G=G,unique_time_id=unique_time_id,filepath=filepath)


    _flushed_logger("Calculating Extrinsic Metrics";logger)
    clustering_quality_metrics(cluster_results_df,filepath);

    return outputs_dict
end
function append_summary(summary_file,vars,varnames)
    open(summary_file, "a") do f
        for i in eachindex(vars)
            write(f,"$(varnames[i]) = \t $(vars[i]) \n")
        end
    end
end

function create_results_dict(run_,function_name)
    results_dict = OrderedDict{Symbol,Vector{Union{String,Int,Float64}}}()
    results_dict[:name] = [function_name]
    results_dict[:num_alloc] = [run_.allocs]
    results_dict[:memory] = [run_.memory]
    results_dict[:times] = run_.times
    results_dict[:avg_time] = [mean(run_.times)]
    results_dict[:med_time] = [median(run_.times)]
    results_dict[:max_time] = [maximum(run_.times)]
    results_dict[:min_time] = [minimum(run_.times)]
    results_dict[:num_samples] =[ run_.params.samples]
    results_dict[:num_evals] = [run_.params.evals]
    return results_dict
end
function benchmark_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = nothing)
    Random.seed!(seed)
    _flushed_logger("Loading data and metadata...";logger)
    anndata_dict1= load_data(datafilename1,seed)
    N = size(anndata_dict1["X"])[2]
    G = size(anndata_dict1["X"])[1]
    num_var_feat = G
    _flushed_logger("Number of cells: $(N) ...";logger)
    _flushed_logger("Number of Genes: $(G) ...";logger)
    unique_time_id,dataset_used_id,experiment_id = make_ids(dataset,G,N)

    _flushed_logger("Preparing Dataset ...";logger)
    gene_names, cell_ids, cell_cluster_labels = preparing_data(anndata_dict1)
    x_input,z = make_nclusion_inputs(anndata_dict1)
    

    filepath = mk_outputs_filepath(outdir,experiment_id,dataset_used_id,unique_time_id)
    _flushed_logger("Preparing saving Directory at $filepath ...";logger)
    mk_outputs_pathname(filepath)

    _flushed_logger("Saving Quick Summary...";logger)
    notes_=""
    summary_file = saving_summary_file(filepath;unique_time_id=unique_time_id,datafilename1=datafilename1,alpha1=alpha1,gamma1=gamma1,KMax=KMax, seed=seed,num_var_feat=num_var_feat,N=N,elbo_ep=elbo_ep,notes_=notes_)
    
    _flushed_logger("Initializing Model parameters...";logger)
    inputs = initialize_model_parameters(x_input,KMax,alpha1,gamma1);

    
    num_iter = inputs[:modelparams].num_iter
    elbologger = ElboFeatures(1,KMax,num_iter)
    inputs[:elbolog] = elbologger

    _flushed_logger("Benchmarking NCLUSION";logger)
    run_ = @benchmark outputs_dict = cavi($inputs;elbo_ep = $elbo_ep);
    results_dict = create_results_dict(run_,"NCLUSION")
    results_dict[:filepath] = [filepath]
    results_dict[:unique_time_id] = [unique_time_id]
    results_dict[:G] = [G]
    results_dict[:N] =[N]

    return results_dict
end