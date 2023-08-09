ENV["GKSwstype"] = "100"

using Logging,LoggingExtras


logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;

function flushed_logger(logger,msg)
    with_logger(logger) do
        @info msg
    end
end




# @info "Using $(Threads.nthreads()) thread(s)...."
# @info "Loading Packages...."
flushed_logger(logger,"Using $(Threads.nthreads()) thread(s)....")
flushed_logger(logger,"Loading Packages....")
using Random
using Distributions
using Flux
# using Flux, Turing
# using Turing.Variational
# using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, LinearAlgebra, HypothesisTests
using Test
# import Debugger
using CSV,DataFrames
using JSON, JSON3
using Dates
using TSne, MultivariateStats, Clustering
using LaTeXStrings, TypedTables, PrettyTables
# using RCall
using Gnuplot, Colors, ColorSchemes
using SpecialFunctions
using Optim
using BenchmarkTools
using DataFramesMeta
using Profile
using JLD2,FileIO
using Hyperopt
using OrderedCollections
using HDF5
using ClusterValidityIndices
using StaticArrays
using BlackBoxOptim

using Pkg
using Distributed

curr_dir = ENV["PWD"]
src_dir = "/hoff_vs/src/"
env_location = "/home/v-mahughes/miniconda3/envs/nclusion/bin/python"
ENV["PYTHON"] = env_location
Pkg.build("PyCall")
Pkg.build("RCall")
using RCall
using PyCall

# @info "Loading NCLUSION Modules...."
flushed_logger(logger,"Loading NCLUSION Modules....")

include(curr_dir*src_dir*"PlottingUtils.jl")
using .PlottingUtils
include(curr_dir*src_dir*"MathUtils.jl")
using .MathUtils
include(curr_dir*src_dir*"MCMCInferenceUtils.jl")
using .MCMCInferenceUtils
include(curr_dir*src_dir*"DataGenUtils.jl")
using .DataGenerationUtils
include(curr_dir*src_dir*"DataPreprocessingUtils.jl")
using .DataPreprocessingUtils
include(curr_dir*src_dir*"MetricsUtils.jl")
using .MetricsUtils
include(curr_dir*src_dir*"IOUtils.jl")
using .IOUtils
include(curr_dir*src_dir*"VariationalInferenceUtils.jl")
using .VariationalInferenceUtils
include(curr_dir*src_dir*"TidyUtils.jl")
using .TidyUtils


include(curr_dir*src_dir*"test_viFunctions.jl")
include(curr_dir*src_dir*"baselines_RCall.jl")
include(curr_dir*src_dir*"viDevelopment_ModelVersions.jl")

# @everywhere using Distributed
# addprocs(1)

flushed_logger(logger,"Setting Plotting Settings....")
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        # @info "Enabling Headless Plotting...."
        flushed_logger(logger,"\t Enabling Headless Plotting....")
        
        Gnuplot.options.gpviewer = false
    else
        Gnuplot.options.gpviewer = true
    end
else
    Gnuplot.options.gpviewer = true
end

if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        # @info "Setting Plotting enviroment variables...."
        flushed_logger(logger,"\t Setting Plotting enviroment variables....")
        to_display=false
    else
        to_display=true
    end
else
    to_display=true
end
function main(ARGS)
    filename,datafilename1,key_str,num_var_feat,scale_fator,seed = ARGS

    if !isempty(num_var_feat)
        num_var_feat = parse(Int64, num_var_feat)
    else
        num_var_feat = 2500
    end
    if !isempty(scale_fator)
        scale_fator = parse(Int64, scale_fator)
    else
        scale_fator = 10000
    end

    if !isempty(seed)
        seed = parse(Int64, seed)
    else
        seed = 12345
    end
    if isempty(key_str)
        key_str = "ch"
    end

    # Get Save Filepath from results file
    filepath = join(split(filename, "/")[1:end-1],"/")*"/"
    
    # # Get HVG (using seurat-v3 method)
    # new_datafilename_vec = split(datafilename1,".")
    # new_datafilename_vec[1]= new_datafilename_vec[1]*"_$(num_var_feat)-seurat_v3"
    # new_datafilename = join(new_datafilename_vec,".")
    # if isfile(new_datafilename)
    #     datafilename1 = new_datafilename
    # else
    #     get_hv_genes(datafilename1,num_var_feat);
    #     datafilename1 = new_datafilename;
    # end
    
    # open anndata file  
    flushed_logger(logger, "Loading Data...")

    Random.seed!(seed)
    fid1 = h5open(datafilename1,"r")
    anndata_dict1 = read(fid1)
    # Get raw counts
    x_mat = deepcopy(anndata_dict1["raw"]["X"])
    N = size(x_mat)[2]


    # Transform Data
    gene_ids = deepcopy(anndata_dict1["var"]["_index"])
    cell_ids = deepcopy(anndata_dict1["obs"]["_index"])
    cell_cluster_labels = deepcopy(anndata_dict1["obs"]["cell_type"]["codes"])
    unique_clusters_labels =  [join(uppercasefirst.(split(el,"_")), " " ) for el in deepcopy(anndata_dict1["obs"]["cell_type"]["categories"])]
    KCalled = length(unique_clusters_labels)
    cluster_label_dict = OrderedDict(unique_clusters_labels[i] => i for i in 1:KCalled)
    rev_cluster_label_dict = OrderedDict(i =>  unique_clusters_labels[i] for i in 1:KCalled)
    cell_labels = [rev_cluster_label_dict[el+1] for el in cell_cluster_labels]
    if all(isa.(cell_labels,Number))
        if iszero(minimum(cell_labels))
            cell_labels = cell_labels .+ 1
        end
    end
    cell_cluster_labels = cell_cluster_labels .+ 1
    num_cnts = deepcopy(anndata_dict1["obs"]["total_counts"])
    x = nothing
    z = nothing
    C_t = nothing
    numi = nothing
    top_genes_bool = nothing
    top_genes = nothing
    top_genes = gene_ids
    top_genes_bool = [true for  el in gene_ids]
    T=1
    x_temp = Vector{Vector{Vector{Float64}}}(undef,T)
    x_temp[1] = [Float64.(collect(col)) for col in eachcol(x_mat)]
    numi_temp = Vector{Vector{Float64}}(undef,T)
    numi_temp[1] = Float64.(collect(num_cnts))
    flushed_logger(logger, "Log normalizing data...")
    log_norm_x = raghavan2021_lognormalization(x_temp;scaling_factor=scale_fator,pseudocount=1.0,numi=numi_temp)
    log_norm_xmat =  hcat(vcat(log_norm_x...)...)
    gene_std_vec = vec(std(log_norm_xmat, dims=2))
    
    sorted_indx  = sortperm(gene_std_vec,rev=true)
    flushed_logger(logger, "Getting HVGs...")
    topGGenes = num_var_feat#sum(sorted_gene_std_vec.>=1.0)
    top_genes = gene_ids[sorted_indx][1:topGGenes]
    top_genes_bool = [in(el,Set(top_genes)) for  el in gene_ids]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    z = Vector{Vector{Int}}(undef,T)
    C_t = Vector{Float64}(undef,T)
    numi = Vector{Vector{Float64}}(undef,T)
    num_cells = N#5000
    x_mat_to_use = nothing
    x_mat_to_use = x_mat
    chosen_cells = randperm(N)[1:num_cells]
    cell_intersect = Set(chosen_cells)
    cell_intersect_bool = [in(el,cell_intersect) for  el in collect(1:N)]
    C_t[1] = sum(cell_intersect_bool)
    
    x[1] = [Float64.(collect(col)) for col in eachcol(x_mat_to_use[top_genes_bool,cell_intersect_bool])]
    z[1] = Int.(collect(cell_cluster_labels))[cell_intersect_bool]
    numi[1] = Float64.(collect(num_cnts))[cell_intersect_bool]


    G = length(x[1][1])
    x_to_use = raghavan2021_lognormalization(x;scaling_factor=scale_fator,pseudocount=1.0,numi=numi);
    x_input=x_to_use;
    unique_clusters = collect(1:KCalled);
    cell_names = cell_ids[cell_intersect_bool]
    gene_means = vec(mean(hcat(vcat(x_input...)...),dims=2))
    gene_names  = gene_ids;
    # num_var_feat = G
    xmat = hcat(vcat(x_input...)...)

    
    
    
    
 
    key = Symbol(key_str)
    flushed_logger(logger, "Loading Cluster Memberships")
    cluster_results_df =CSV.read(filename,DataFrame,header=true)
    z_post_s = [cluster_results_df[!,end]]
    flushed_logger(logger, "Saving Cluster Validity Indices")
    posteriorSummarizationsDict = calc_time_invariant_CVI(xmat,z_post_s,key);
    unique_time_id = get_unique_time_id()
    cvi_df = DataFrame(posteriorSummarizationsDict[String(key)], :auto);
    rename!(cvi_df,Symbol.(["idx", "mean", "std","upper","lower"]));
    CSV.write(filepath*"$(uppercase(String(key)))_Summary_AllTimepoints_$(num_var_feat)genes-$unique_time_id.csv",  cvi_df)
    
end

function calc_time_invariant_CVI_summarization(cvi; conf_level=0.95)



    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    cvi_values = cvi


    cvi_means = mean(cvi_values)
    cvi_std = std(cvi_values)
    # cvi_conf_level = t_test(cvi_values; conf_level=conf_level)
    # cvi_quantiles_mat = reduce(hcat, [cvi_conf_level...])
    if length(cvi_values) > 1
        cvi_conf_level = t_test(cvi_values; conf_level=conf_level)
        cvi_quantiles_mat = reduce(hcat, [cvi_conf_level...])
    else
        cvi_quantiles_mat = [cvi_means, cvi_means]
    end

    cvi_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    cvi_summary[1, 1] = 1
    cvi_summary[1, 2] = cvi_means
    cvi_summary[1, 3] = cvi_std
    cvi_summary[1, 4] = cvi_quantiles_mat[1]
    cvi_summary[1, 5] = cvi_quantiles_mat[2]


    return cvi_summary
end

function calc_time_invariant_CVI(xmat,z_post_s)
    cvi_func_dict = OrderedDict(:ch => CH(), :db =>DB(), :gd43 => GD43(), :gd53 => GD53()) #:csil => cSIL(), , :ps=> PS(), :rcip=> rCIP(), :wb => WB(),:xb => XB() 
    metrics_list = setup_metrics_list()
    posteriorSummarizationsDict = OrderedDict()
    num_posterior_samples = length(z_post_s)
    for key in keys(cvi_func_dict)
        flushed_logger(logger,"\t Calculating $key metric now...")
        cvi = cvi_func_dict[key]
        cvi_criterion_value_b_vec = Vector{Float64}(undef, num_posterior_samples)
        elsaped_time = @elapsed begin
            for s in 1:num_posterior_samples 
                z_s = z_post_s[s]
                cvi_b = cvi 
                z_infer_vec = vcat(z_s...)
                cvi_criterion_value_b_vec[s] = get_cvi!(cvi_b, xmat, z_infer_vec)
            end
        end
        flushed_logger(logger,"\t Finished Calculating $key metric. Metric calculations took $elsaped_time seconds to run...")
        posteriorSummarizationsDict[String(key)] = calc_time_invariant_CVI_summarization(cvi_criterion_value_b_vec;conf_level=0.95)
    end
    return posteriorSummarizationsDict
end
function calc_time_invariant_CVI(xmat,z_post_s,key)
    cvi_func_dict = OrderedDict(:ch => CH(), :db =>DB(), :gd43 => GD43(), :gd53 => GD53(),:csil => cSIL(),  :ps=> PS(), :rcip=> rCIP(), :wb => WB(),:xb => XB() ) # 
    posteriorSummarizationsDict = OrderedDict()
    num_posterior_samples = length(z_post_s)
    flushed_logger(logger,"\t Calculating $key metric now...")
    cvi = cvi_func_dict[key]
    cvi_criterion_value_b_vec = Vector{Float64}(undef, num_posterior_samples)
    elsaped_time = @elapsed begin
        for s in 1:num_posterior_samples 
            z_s = z_post_s[s]
            cvi_b = cvi 
            z_infer_vec = vcat(z_s...)
            cvi_criterion_value_b_vec[s] = get_cvi!(cvi_b, xmat, z_infer_vec)
        end
    end
    flushed_logger(logger,"\t Finished Calculating $key metric. Metric calculations took $elsaped_time seconds to run...")
    posteriorSummarizationsDict[String(key)] = calc_time_invariant_CVI_summarization(cvi_criterion_value_b_vec;conf_level=0.95)
    return posteriorSummarizationsDict
end

function distributed_calc_time_invariant_CVI(p,cvi_func_dict,xmat,z_post_s)
    # cvi_func_dict = OrderedDict(:ch => CH(), :csil => cSIL(), :db =>DB(), :gd43 => GD43(), :gd53 => GD53()) #, :ps=> PS(), :rcip=> rCIP(), :wb => WB(),:xb => XB() 
    
    num_posterior_samples = length(z_post_s)
    key_vec = collect(keys(cvi_func_dict))
    key = key_vec[p]
    cvi = cvi_func_dict[key]
    cvi_criterion_value_b_vec = Vector{Float64}(undef, num_posterior_samples)
    # @info "Calculating $key metric now..."
    flushed_logger(logger,"\t Calculating $key metric now...")

    elsaped_time = @elapsed begin
        for s in 1:num_posterior_samples 
            z_s = z_post_s[s]
            cvi_b = cvi 
            z_infer_vec = vcat(z_s...)
            cvi_criterion_value_b_vec[s] = get_cvi!(cvi_b, xmat, z_infer_vec)
        end
    end
    # @info "Finished Calculating $key metric. Metric calculations took $elsaped_time seconds to run..."
    flushed_logger(logger,"\t Finished Calculating $key metric. Metric calculations took $elsaped_time seconds to run...")
    return p, String(key), calc_time_invariant_CVI_summarization(cvi_criterion_value_b_vec;conf_level=0.95)
end

function clustering_quality_metrics(results_df,filepath_;beta_ = 1.0, conf_lvl = 0.95)
       #Get Called assignments
       called_assignments_vec = results_df[:,end-1]

       #Get Inferred assignments
       inferred_assignments_vec = results_df[:,end]
   
       #Tranform data into format needed for metric functions
       conditionpoints = sort(unique(results_df[:,1]))
       conditionpoints_dict =countmap(results_df[:,1])
       N_t = [conditionpoints_dict[el] for el in conditionpoints]
       ranges_ = get_timeranges(N_t)
       called_assignments = [called_assignments_vec[st:en] for (st,en) in ranges_]
       inferred_assignments = [inferred_assignments_vec[st:en] for (st,en) in ranges_]
       T = length(conditionpoints);
   
    #    @info "Calculating Metrics..."
        flushed_logger(logger,"\t Calculating Metrics...")

       #Get Metrics on a per time/condition basis
       ari_vov,_,_,_ = getRandIndices(called_assignments,[inferred_assignments])
       nmi_vov = getNMI(called_assignments,[inferred_assignments])
       vmeasure_vov = getVmeasure(called_assignments,[inferred_assignments];beta=beta_)
       varinfo_vov = getVarInfo(called_assignments,[inferred_assignments])
       jaccard_vov = getJaccardSimilarity(called_assignments,[inferred_assignments])
   
       #Get overall Metric
       ari_invar = time_invariant_ari(called_assignments,[inferred_assignments])
       nmi_invar =time_invariant_nmi(called_assignments,[inferred_assignments])
       vmeasure_invar =time_invariant_vmeasure(called_assignments,[inferred_assignments],;beta=beta_)
       varinfo_invar =time_invariant_varinfo(called_assignments,[inferred_assignments])
       jaccard_invar =time_invariant_jaccard(called_assignments,[inferred_assignments])
   
    #    @info "Saving Metrics..."
       flushed_logger(logger,"\t Saving Metrics...")
       #Save ARI as CSV
       ari_invar_summary = calc_time_invariant_ARI_summarization(ari_invar;conf_level=conf_lvl);
       ari_invar_summary_df = DataFrame(ari_invar_summary,:auto);
       rename!(ari_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/ARI_Summary_AllTimepoints.csv", ari_invar_summary_df) 
       ari_var_summary =calc_time_variant_ARI_summarization(ari_vov,T;conf_level=conf_lvl)
       ari_var_summary_df = DataFrame(ari_var_summary,:auto);
       rename!(ari_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/ARI_Summary_PerTimepoints.csv", ari_var_summary_df) 
   
       #Save NMI as CSV
       nmi_invar_summary = calc_time_invariant_NMI_summarization(nmi_invar;conf_level=conf_lvl)
       nmi_invar_summary_df = DataFrame(nmi_invar_summary,:auto);
       rename!(nmi_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/NMI_Summary_AllTimepoints.csv", nmi_invar_summary_df) 
       nmi_var_summary = calc_time_variant_NMI_summarization(nmi_vov,T;conf_level=conf_lvl)
       nmi_var_summary_df = DataFrame(nmi_var_summary,:auto);
       rename!(nmi_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/NMI_Summary_PerTimepoints.csv", nmi_var_summary_df) 
   
       #Save VMeasure as CSV
       vmeasure_invar_summary = calc_time_invariant_Vmeasure_summarization(vmeasure_invar;conf_level=conf_lvl)
       vmeasure_invar_summary_df = DataFrame(vmeasure_invar_summary,:auto);
       rename!(vmeasure_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/VMeasure_Summary_AllTimepoints.csv", vmeasure_invar_summary_df) 
       vmeasure_var_summary =calc_time_variant_Vmeasure_summarization(vmeasure_vov,T;conf_level=conf_lvl)
       vmeasure_var_summary_df = DataFrame(vmeasure_var_summary,:auto);
       rename!(vmeasure_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/VMeasure_Summary_PerTimepoints.csv", vmeasure_var_summary_df) 
   
       #Save VarInfo as CSV
       varinfo_invar_summary = calc_time_invariant_VarInfo_summarization(varinfo_invar;conf_level=conf_lvl)
       varinfo_invar_summary_df = DataFrame(varinfo_invar_summary,:auto);
       rename!(varinfo_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/VarInfo_Summary_AllTimepoints.csv", varinfo_invar_summary_df) 
       varinfo_var_summary = calc_time_variant_VarInfo_summarization(varinfo_vov,T;conf_level=conf_lvl)
       varinfo_var_summary_df = DataFrame(varinfo_var_summary,:auto);
       rename!(varinfo_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/VarInfo_Summary_PerTimepoints.csv", varinfo_var_summary_df) 
   
       #Save Jaccard as CSV
       jaccard_invar_summary = calc_time_invariant_Jaccard_summarization(jaccard_invar;conf_level=conf_lvl)
       jaccard_invar_summary_df = DataFrame(jaccard_invar_summary,:auto);
       rename!(jaccard_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/Jaccard_Summary_AllTimepoints.csv", jaccard_invar_summary_df)
       jaccard_var_summary = calc_time_variant_Jaccard_summarization(jaccard_vov,T;conf_level=conf_lvl)
       jaccard_var_summary_df = DataFrame(jaccard_var_summary,:auto);
       rename!(jaccard_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/Jaccard_Summary_PerTimepoints.csv", jaccard_var_summary_df) 
   
      
end


main(ARGS)