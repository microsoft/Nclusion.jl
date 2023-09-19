module DataPreprocessingUtils
    using Random
    using Distributions
    using Turing
    using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, Combinatorics
    using Test
    import Debugger
    using CSV,DataFrames
    using RCall
    using Gnuplot
    using OrderedCollections
    using HDF5
    using Pkg
    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"
    # env_location = curr_dir*"/nclsn/bin/python"
    # ENV["PYTHON"] = env_location
    # Pkg.build("PyCall")
    # using PyCall

    include(curr_dir*src_dir*"IOUtils.jl")
    using .IOUtils
    include(curr_dir*src_dir*"TidyUtils.jl")
    using .TidyUtils
    include(curr_dir*src_dir*"DataGenUtils.jl")
    using .DataGenerationUtils


           #From synDataPreprocess.jl
    export getCustomersAtTimepoint,
           getTrueClusterMembershipDict, 
           genData_vec2dict,

           #From This file
           raghavan2021_lognormalization,
           raghavan2021_lognormalization_depracated,
           raghavan2021_centerAndScale,
           raghavan2021_filterCells,
           raghavan2021_filterCellsWithMetadata,
           raghavan2021_filterGenes,
           raghavan2021_AddLabels!,
           mode_shift,

           basal_condition,
           classical_condition,
           hybrid_condition,
           get_cell_class_label,
           sum_gene_counts,
           subsample_on_label,
           get_tumor_stats_column_indices,
           append_first_col_indx!,
           get_tumor_stats,
           get_max_score,
           append_maxscore!,
           insertCellLabels!,
           extract_histbin_infomation,
           get_celltype_1feature_metadata,
           get_metadata_at_timepoint,
           get_timepoint_labels,
           get_biopsy_labels,

           R_lognorm_FeatSelect_Scale_noMetaDataGenenames,
           R_lognorm_FeatSelect_Scale,
           R_VariableFeatSelect,
           get_hv_genes,

           generate_processed_raghavan2021_data,
           subsample_tuning_set,
           generate_subsample_data_processed_raghavan2021_data,
           generate_processed_nonmalignant_raghavan2021_data,
           generate_processed_malignant_pt_raghavan2021_data,
           generate_processed_simulated_copulacorrelatedcounts_data,
           generate_subsample_data_processed_simulated_data,
           generate_processed_simulated_copulacorrelatedcounts_data,
           generate_processed_1condition_anndata_data,
           generate_subsample_data_processed_1condition_anndata_data
           
    
    function raghavan2021_lognormalization_depracated(x;scaling_factor=10_000,pseudocount=1.0)
        T = length(x)
        C_t = [length(el) for el in x]
        G = length(x[1][1])
        x_transformed = Vector{Vector{Vector{Float64}}}(undef,T)
        for t in 1:T
            cells = C_t[t]
            normed_val = x[t] ./  sum.(x[t])
            scale_normed = scaling_factor .* normed_val
            transformed_shifted_scale_normed = Vector{Vector{Float64}}(undef,cells)
            for i in 1:cells
                transformed_shifted_scale_normed[i] = log.(scale_normed[i] .+ pseudocount)
            end
            x_transformed[t] =  transformed_shifted_scale_normed
        end
        return x_transformed
    end
    function raghavan2021_lognormalization(x;scaling_factor=10_000,pseudocount=1.0,numi=nothing)
        T = length(x)
        C_t = [length(el) for el in x]
        G = length(x[1][1])
        x_transformed = Vector{Vector{Vector{Float64}}}(undef,T)
        for t in 1:T
            cells = C_t[t]
            transformed_shifted_scale_normed = Vector{Vector{Float64}}(undef,cells)
            for i in 1:cells
                if isnothing(numi)
                    normed_val = x[t][i] ./  sum(x[t][i])
                else
                    normed_val = x[t][i] ./ numi[t][i]
                end
                scale_normed = scaling_factor .* normed_val
            
                transformed_shifted_scale_normed[i] = log.(scale_normed .+ pseudocount)
            end
            x_transformed[t] =  transformed_shifted_scale_normed
        end
        return x_transformed
    end

    function raghavan2021_centerAndScale(x_transformed)
        T = length(x_transformed)
        center_and_scale(cell,mu_x,std_x) = (cell .- mu_x) ./ std_x
        x_transformed_stacked = collect(Base.Iterators.flatten(x_transformed))
        x_mean = mean(x_transformed_stacked,dims = 1)[1]
        x_std = std(x_transformed_stacked)
        x_std[x_std .== 0.0] .= 1.0
        x_transformed = [map(cell -> center_and_scale(cell, x_mean,x_std) ,x_transformed[t]) for t in 1:T]
    end
    function mode_shift(x)
        T = length(x)
        C_t = C_t =[length(el) for el in x]
        xmat = permutedims(reduce(hcat,reduce(vcat,x)));
        modes_x = [mode(col) for col in eachcol(xmat)];
        newx = [[ x[t][i] .- modes_x  for i in 1:C_t[t]] for t in 1:T]
    end
    function raghavan2021_filterCells(data_df,metadata_df;lowerbound_genes=400,detected_UMIs=1000,pct_mito=0.5,upperbound_genes=8000)
        # < 400 genes detected; < 1,000 UMIs; > 50% mitochondrial reads
        # further trimmed to remove cells with > 8,000 genes which represent outliers and likely doublet cells. We also removed genes that were not detected in at least 50 cells.
        original_cell_count = size(data_df)[2]
        conditional_ = (lowerbound_genes .< metadata_df.nGene .< upperbound_genes) .& (metadata_df[!, "percent.mito"] .< pct_mito) .& (metadata_df.nUMI .> detected_UMIs)
        conditional_names = names(data_df)[conditional_];
        subset_data = data_df[!, conditional_names];
        subset_metadata = metadata_df[conditional_,1:end];
        new_cell_count = size(subset_data)[2]
        num_filtered_cells = original_cell_count - new_cell_count
        println("Filtered $(num_filtered_cells) cells")
        return subset_data,subset_metadata
    end
    function raghavan2021_filterCellsWithMetadata(data_df,metadata_df;lowerbound_genes=400,detected_UMIs=1000,pct_mito=0.5,upperbound_genes=8000)
        # < 400 genes detected; < 1,000 UMIs; > 50% mitochondrial reads
        # further trimmed to remove cells with > 8,000 genes which represent outliers and likely doublet cells. We also removed genes that were not detected in at least 50 cells.
        original_cell_count = size(data_df)[2]
        cells_with_metadata = [in(el,Set(metadata_df[:,1])) for el in names(data_df)];
        data_df =  data_df[!, cells_with_metadata];
        conditional_ = (lowerbound_genes .< metadata_df.nGene .< upperbound_genes) .& (metadata_df[!, "percent.mito"] .< pct_mito) .& (metadata_df.nUMI .> detected_UMIs)
        conditional_names = names(data_df)[conditional_];
        subset_data = data_df[!, conditional_names];
        subset_metadata = metadata_df[conditional_,1:end];
        new_cell_count = size(subset_data)[2]
        num_filtered_cells = original_cell_count - new_cell_count
        println("Filtered $(num_filtered_cells) cells")
        return subset_data,subset_metadata
    end
    function raghavan2021_filterGenes(data_df,gene_names_df;detected_genes=50)
        # < 400 genes detected; < 1,000 UMIs; > 50% mitochondrial reads
        # further trimmed to remove cells with > 8,000 genes which represent outliers and likely doublet cells. We also removed genes that were not detected in at least 50 cells.
        original_gene_count = size(data_df)[1]
        pseudobulk_gene_expression =  [sum(r) for r in eachrow(data_df)]
        lowly_expessed_genes_conditional_ = pseudobulk_gene_expression .> detected_genes
        subset_df = data_df[lowly_expessed_genes_conditional_,1:end]
        subset_genes_names = gene_names_df[lowly_expessed_genes_conditional_,1:end]
        new_gene_count = size(subset_df)[1]
        num_filtered_genes = original_gene_count - new_gene_count
        println("Filtered $(num_filtered_genes) genes")
        return subset_df, subset_genes_names
    end
    function raghavan2021_AddLabels!(metadata_df;tumor_stats_columns = ["Hybrid_Specific1", "sc_Classical1", "sc_Basal1"])
        tumor_stats_columns = get_tumor_stats_column_indices(metadata_df;tumor_stats_columns = tumor_stats_columns);
        append_first_col_indx!(tumor_stats_columns);
        tumor_stats = get_tumor_stats(metadata_df, tumor_stats_columns);
        max_score = get_max_score(tumor_stats);
        append_maxscore!(tumor_stats,max_score);
        append_maxscore!(metadata_df,max_score);
        insertCellLabels!(tumor_stats)
        insertCellLabels!(metadata_df,tumor_stats)
        return metadata_df,tumor_stats
    end
    
    function tidy_center_scale(xmat)
        T = length(unique(xmat[:,1]))
        N = length(unique(xmat[:,2]))
        G = length(unique(xmat[:,3]))
        nrows = size(xmat)[1]
        ncols = size(xmat)[2]
        gene_ids = collect(1:G)
        scaled_data_vec = Vector{Float64}(undef,nrows)
        center_scaling = n -> StatsBase.standardize(ZScoreTransform,map(Float64,n))
        for j in 1:G
            gene_query = isequal_condition(j)
            gene_vals = filter_on(xmat,3,4,gene_query)
            scaled_data_vec[gene_query(xmat[:,3])] = center_scaling(gene_vals)
        end
        new_xmat = Matrix{Union{Float64,Int}}(undef,nrows,ncols+1)
        new_xmat[:,1:ncols] =  xmat
        new_xmat[:, end] = scaled_data_vec
        return new_xmat
    end
    
    

    function basal_condition(df, i)
        return df[i,:max_score] == df[i,:sc_Basal1] && df[i,:max_score] >= 0.25
    end
    
    function classical_condition(df, i)
        return df[i,:max_score] == df[i,:sc_Classical1] && df[i,:max_score] >= 0.0
    end
    
    function hybrid_condition(df, i)
        return df[i,:max_score] == df[i,:Hybrid_Specific1] && df[i,:max_score] >= 0.3
    end
    
    function get_cell_class_label(df)
        num_rows = size(df)[1]
        cell_labels = [basal_condition(df, i) ? "Basal" :
         classical_condition(df, i) ? "Classical" :
          hybrid_condition(df, i) ? "Hybrid" :
           "Organoid" for i in 1:num_rows]
        return cell_labels
    end
    # function get_cell_class_label(df)
    #     num_rows = size(df)[1]
    #     cell_labels = []
    #     for i in 1:num_rows
    #         if basal_condition(df, i) 
    #             push!(cell_labels,"Basal")\
    #             continue
    #         elseif classical_condition(df, i)
    #             push!(cell_labels,"Classical")
    #             continue
    #         elseif hybrid_condition(df, i)
    #             push!(cell_labels,"Hybrid")
    #             continue
    #         else
    #             push!(cell_labels,"Organoid")
    #             continue
    #         end
    #     end
    #     return cell_labels
    # end
    function sum_gene_counts(x)
        return [sum(x[i]) for i in 1:length(x)]
    end
    
    function subsample_on_label(meta_df, label)
        return meta_df[(meta_df.cell_labels .== label), :]
    end
    
    function get_tumor_stats_column_indices(metadata_df;tumor_stats_columns = ["Hybrid_Specific1", "sc_Classical1", "sc_Basal1"])
        num_cols = length(names(metadata_df))
        columnnames_indices = collect(1:num_cols)
        included_cols =  [el in tumor_stats_columns for el in names(metadata_df)]
        return columnnames_indices[included_cols]
    end
    
    function append_first_col_indx!(tumor_stats_columns)
        insert!(tumor_stats_columns,1,1)
    end
    
    function get_tumor_stats(metadata_df, tumor_stats_columns)
        tumor_stats = metadata_df[!,tumor_stats_columns]
        # num_rows = size(tumor_stats)[1]
        # max_score = [maximum(tumor_stats[i,2:end]) for i in 1:num_rows]
        # insertcols!(subsampled_metadata, length(names(subsampled_metadata)) + 1, :max_score => max_score)
        return tumor_stats
    end
    
    function get_max_score(tumor_stats)
        num_rows = size(tumor_stats)[1]
        max_score = [maximum(tumor_stats[i,2:end]) for i in 1:num_rows]
        # insertcols!(subsampled_metadata, length(names(subsampled_metadata)) + 1, :max_score => max_score)
        return max_score
    end
    function  append_maxscore!(df,max_score)
        insertcols!(df, length(names(df)) + 1, :max_score => max_score)
    end
    
    function insertCellLabels!(df)
        cell_labels = get_cell_class_label(df)
        insertcols!(df, length(names(df)) + 1, :cell_labels => cell_labels)
    end
    function insertCellLabels!(df1,df2)
        cell_labels = get_cell_class_label(df2)
        insertcols!(df1, length(names(df1)) + 1, :cell_labels => cell_labels)
    end
    
    function get_metadata_at_timepoint(tp,query_metadata, reference_metadata)
        metadata_to_use = query_metadata[reference_metadata[!,:orig_ident] .== sort(unique(reference_metadata[!,:orig_ident]))[tp],:]
        return metadata_to_use
    end
    function get_celltype_1feature_metadata(metadata_to_use, label,feature_column;filter_label_function_dict=nothing)
        if isnothing(filter_label_function_dict)
            isbasal(label::String) = label == "Basal"
            isclassical(label::String) = label == "Classical"
            ishybrid(label::String) = label == "Hybrid"
            isorganoid(label::String) = label == "Organoid"
            function_list = [isbasal,isclassical,ishybrid,isorganoid]
            keys_list = ["Basal","Classical","Hybrid","Organoid"]
            filter_label_function_dict = Dict()
            for i in 1:length(keys_list)
                key = keys_list[i]
                filter_function = function_list[i]
                filter_label_function_dict[key] = filter_function
            end
        end
        celltype_filter_function = filter_label_function_dict[label] 
        celltype_val = Float64.(Matrix(filter(:cell_labels => celltype_filter_function , metadata_to_use[!,["cell_labels",feature_column]]))[:,2])
        return celltype_val
    
    end
    function extract_histbin_infomation(metadata_to_use, unique_clusters,feature_column;filter_label_function_dict=nothing,binsize = 0.1)
        KCalled = length(unique_clusters)
        celltype_histbins = Vector{Vector{Float64}}(undef,KCalled)
        celltype_histcounts = Vector{Vector{Int64}}(undef,KCalled)
        lowest_feature_val = minimum(metadata_to_use[!,feature_column])
        highest_feature_val = maximum(metadata_to_use[!,feature_column])
        for l in 1:KCalled
            label  = unique_clusters[l]
            celltype_val = get_celltype_1feature_metadata(metadata_to_use, label,feature_column;filter_label_function_dict=filter_label_function_dict)
            if length(celltype_val) == 1
                celltype_histbins[l] = collect(lowest_feature_val:0.1:highest_feature_val)
                celltype_histcounts[l] = (celltype_histbins[l] .- 0.5*binsize) .< celltype_val[1] .<= (celltype_histbins[l] .+ 0.5*binsize)
            elseif length(celltype_val) ==0
                celltype_histbins[l] = collect(lowest_feature_val:0.1:highest_feature_val)
                celltype_histcounts[l] = zeros(Float64,length(celltype_histbins[l]))
            elseif length(unique(celltype_val)) ==1
                celltype_histbins[l] = collect(lowest_feature_val:0.1:highest_feature_val)
                celltype_histcounts[l] = (celltype_histbins[l] .- 0.5*binsize) .< unique(celltype_val)[1] .<= (celltype_histbins[l] .+ 0.5*binsize)
            else
                celltype_hist = hist(celltype_val, bs=binsize)
                celltype_histbins[l] = celltype_hist.bins
                celltype_histcounts[l] = celltype_hist.counts
            end        
        end
        max_counts = maximum(vcat(celltype_histcounts...))
        return celltype_histbins,celltype_histcounts,max_counts
    end

    function get_timepoint_labels(z)
        T= length(z)
        C_t = [length(el) for el in z]
        return [[t for i in 1:C_t[t]] for t in 1:T]
    end
    function get_biopsy_labels(z)
        T= length(z)
        C_t = [length(el) for el in z]
        return [[ t==1 ? true : false for i in 1:C_t[t]] for t in 1:T]
    end

    function generate_processed_raghavan2021_data(datafilename,metadatafilename,genenamesfilename,unique_time_id;lowerbound_genes=400,detected_UMIs=1000,pct_mito=0.5,upperbound_genes=8000,detected_genes=50,tumor_stats_columns = ["Hybrid_Specific1", "sc_Classical1", "sc_Basal1"],num_var_feat = 5000,only_var_feat_bool = true,scale_fator = 10000,seed = 1989)
        Random.seed!(seed)
        data_prep = CSV.read(datafilename,DataFrame);
        metadata_prep = CSV.read(metadatafilename,DataFrame);
        gene_names_prep = CSV.read(genenamesfilename,DataFrame);
        data_prep,metadata_prep = raghavan2021_filterCellsWithMetadata(data_prep,metadata_prep;lowerbound_genes=lowerbound_genes,detected_UMIs=detected_UMIs,pct_mito=pct_mito,upperbound_genes=upperbound_genes);
        data_prep,gene_names_prep = raghavan2021_filterGenes(data_prep,gene_names_prep;detected_genes=detected_genes);
        metadata_prep,tumor_stats = raghavan2021_AddLabels!(metadata_prep;tumor_stats_columns = tumor_stats_columns);


        x_int_lognorm_postp,x_int_scaled_postp,metadata,gene_names_lognorm_postp,gene_names_scaled_postp = R_lognorm_FeatSelect_Scale(data_prep,metadata_prep,gene_names_prep;num_var_feat = num_var_feat, scale_fator = scale_fator,only_var_feat=only_var_feat_bool,seed=seed); 

        arg_str_list_prep = @name datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool;
        key_list_prep = naming_vec(arg_str_list_prep);
        var_list_prep = [datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool];

        scaled_data = DataFrame(x_int_scaled_postp, :auto);
        lognorm_data = DataFrame(x_int_lognorm_postp, :auto);
        rename!(scaled_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        rename!(lognorm_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        


        data = scaled_data;
        data_mat = x_int_scaled_postp;
        gene_names = gene_names_scaled_postp;
        dropmissing!(metadata)
        cells_not_missing = [in(el,Set(metadata[:,"Column1"])) for el in names(data)];
        data =  data[!, cells_not_missing];

        G,T,C_t,cell_timepoints_index = generate_modelStatisticsFromDataFrame(data);
        KCalled,unique_clusters_labels,remap_trueMembership,cluster_map = generate_LabelInfoFromMetaDataFrame(metadata);
        x,z,x_scaled,x_lognorm = generate_modelInputsFromDataFrames(data,metadata;scaled_data = scaled_data,lognorm_data = lognorm_data);
        x_used_label_dict = OrderedDict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        N = sum(C_t)
        unique_clusters = collect(1:KCalled)

        x_to_use = x;#x_scaled;
        true_z = z;
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ sum(C_t);
        
        x_used_label = "x_scaled";
        x_used_label_dict[x_used_label] = true;

        scaling_function = "Seurat(R)"
        data_generating_function = "generate_processed_raghavan2021_data"
        file = split(datafilename, "/")[end]
        ptid = split(metadata_prep[1,1],"_")[occursin.("PANFR",split(metadata_prep[1,1],"_"))][1]

        collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]
        
        
        data_verbal_summary = "The data used was taken from the patient with id: $ptid.\n G = $G, T = $T, N_t = $C_t ( $(sum(C_t)) cells total), KTrue = $KCalled. \n Number of Highly Variably genes was $num_var_feat and we $(only_var_feat_bool ? " DID " : " DID NOT" ) exclusively use these genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])). \nData was centered and scaled using the function: scaling_function. The data was preprocessed using the function: $data_generating_function   \n The preprocessing was as follow: cells were filtered with less than $lowerbound_genes genes detected, greater than $upperbound_genes genes dectected, greater than $detected_UMIs UMIs detected, or more the $pct_mito mitochondrial genes detected. Genes were filtered if they did not appear in at least $detected_genes cells. Labels were assigned using the procedure outline in Raghavan et. al. 2021. $scaling_function was used to log normalize (scale factor of $scale_fator and pseudocount of 1) and center and scale the data. The final data set used was the $(collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]) output from the $data_generating_function  \n The seed used is $seed . "

        
        arg_str_list_prep = @name G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name ptid, file, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [ptid, file, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    function generate_processed_nonmalignant_raghavan2021_data(datafilename,metadatafilename,unique_time_id;lowerbound_genes=400,detected_UMIs=1000,pct_mito=0.5,genenamesfilename = "",upperbound_genes=8000,detected_genes=50,num_var_feat = 5000,only_var_feat_bool = true,scale_fator = 10000,seed = 1989)
        # metadatafilename =  "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/NormalCells_metadata_210305.csv" 
        # datafilename = "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/NormalCells_RawDGE_210305.csv"#
        Random.seed!(seed)
        data_df = CSV.read(datafilename,DataFrame);
        metadata_df = CSV.read(metadatafilename,DataFrame);
        gene_names_prep = nothing;
        if !isempty(genenamesfilename)
            gene_names_prep = CSV.read(genenamesfilename,DataFrame);
        else
            gene_names_prep = data_df[!,1];
        end
        data_df = data_df[!,2:end];
        cells_ptID = [el[2] for el in split.(names(data_df),"_")];

        metadata_cells_not_0580 = [el != "PANFR0580" for el in metadata_df.ID];
        data_cells_not_0580 = [el != "PANFR0580" for el in cells_ptID ];
        metadata_prep = deepcopy(metadata_df[metadata_cells_not_0580,:]);
        data_prep_ = deepcopy(data_df[!,data_cells_not_0580]);
        cells_not_missing1 = [in(el,Set(metadata_prep.Column1)) for el in names(data_prep_)];
        data_prep = data_prep_[!, cells_not_missing1];
        data_prep,metadata_prep = raghavan2021_filterCellsWithMetadata(data_prep,metadata_prep;lowerbound_genes=lowerbound_genes,detected_UMIs=detected_UMIs,pct_mito=pct_mito,upperbound_genes=upperbound_genes);
        data_prep,gene_names_prep = raghavan2021_filterGenes(data_prep,gene_names_prep;detected_genes=detected_genes);
        metadata_prep[!,:cell_labels] = metadata_prep[!,"cell.types"];
        dropmissing!(metadata_prep)
        cells_not_missing2 = [in(el,Set(metadata_prep[:,"Column1"])) for el in names(data_prep)];
        data_prep =  data_prep[!, cells_not_missing2];
        metadata_prep = metadata_prep[cells_not_missing2, :];
        cells_not_0580 = [el != "PANFR0580" for el in metadata_prep.ID];
        data_prep =  data_prep[!, cells_not_0580];
        metadata_prep = metadata_prep[cells_not_0580, :];


        x_int_lognorm_postp,x_int_scaled_postp,metadata,gene_names_lognorm_postp,gene_names_scaled_postp = R_lognorm_FeatSelect_Scale(data_prep,metadata_prep,gene_names_prep;num_var_feat = num_var_feat, scale_fator = scale_fator,only_var_feat=only_var_feat_bool,seed=seed); 

        arg_str_list_prep = @name datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool;
        key_list_prep = naming_vec(arg_str_list_prep);
        var_list_prep = [datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool];

        scaled_data = DataFrame(x_int_scaled_postp, :auto);
        lognorm_data = DataFrame(x_int_lognorm_postp, :auto);
        rename!(scaled_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        rename!(lognorm_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        


        data = scaled_data;
        data_mat = x_int_scaled_postp;
        gene_names = gene_names_scaled_postp;
        dropmissing!(metadata)
        cells_not_missing = [in(el,Set(metadata[:,"Column1"])) for el in names(data)];
        data =  data[!, cells_not_missing];

        G,T,C_t,cell_timepoints_index = generate_modelStatisticsFromDataFrame(data);
        T = 1
        C_t = [sum(C_t)]
        cell_timepoints_index = vcat(cell_timepoints_index...)

        KCalled,unique_clusters_labels,remap_trueMembership,cluster_map = generate_LabelInfoFromMetaDataFrame(metadata);
        x,z,x_scaled,x_lognorm = generate_modelInputsFromDataFrames(data,metadata;scaled_data = scaled_data,lognorm_data = lognorm_data);
        x_used_label_dict = OrderedDict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        N = sum(C_t)
        unique_clusters = collect(1:KCalled)
        x = [vcat(x...)]
        x_lognorm = [vcat(x_lognorm...)]
        x_scaled = [vcat(x_scaled...)]
        z = [vcat(z...)]
        x_used_label_dict = Dict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        x_to_use = x;#x_scaled;
        true_z = z;
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ N;
        
        x_used_label = "x_scaled";
        x_used_label_dict[x_used_label] = true;

        scaling_function = "Seurat(R)"
        data_generating_function = "generate_processed_nonmalignant_raghavan2021_data"
        file = split(datafilename, "/")[end]
        metafile = split(metadatafilename, "/")[end]

        collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]
        
        
        data_verbal_summary = "The data used was taken from the patient with id: Non-Malignant Cells.\n G = $G, T = $T, N_t = $C_t ( $(sum(C_t)) cells total), KCalled = $KCalled. \n Number of Highly Variably genes was $num_var_feat and we $(only_var_feat_bool ? " DID " : " DID NOT" ) exclusively use these genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])). \nData was centered and scaled using the function: scaling_function. The data was preprocessed using the function: $data_generating_function   \n The preprocessing was as follow: cells were filtered with less than $lowerbound_genes genes detected, greater than $upperbound_genes genes dectected, greater than $detected_UMIs UMIs detected, or more the $pct_mito mitochondrial genes detected. Genes were filtered if they did not appear in at least $detected_genes cells. Labels were assigned using the procedure outline in Raghavan et. al. 2021. $scaling_function was used to log normalize (scale factor of $scale_fator and pseudocount of 1) and center and scale the data. The final data set used was the $(collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]) output from the $data_generating_function  \n The seed used is $seed . "

        
        arg_str_list_prep = @name G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);
        ptid = "NonMalignant"

        arg_str_list_summary = @name ptid, file,metafile, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [ptid, file,metafile, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    function generate_processed_malignant_pt_raghavan2021_data(datafilename,metadatafilename,unique_time_id,ptid;lowerbound_genes=400,detected_UMIs=1000,pct_mito=0.5,genenamesfilename = "",upperbound_genes=8000,detected_genes=50,num_var_feat = 5000,only_var_feat_bool = true,scale_fator = 10000,seed = 1989)
        # metadatafilename =  "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/Primary_Organoids_metadata_210224_32060Cells.csv" # "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/complete_MetaData_70170cells.csv"
        # datafilename = "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/Primary_Organoids_RawDGE_210224_32060Cells.csv"#
        # ptid = "PANFR0489"#"PANFR0575"#
        Random.seed!(seed)
        data_df = CSV.read(datafilename,DataFrame);
        metadata_df = CSV.read(metadatafilename,DataFrame);
        gene_names_prep = nothing;
        if !isempty(genenamesfilename)
            gene_names_prep = CSV.read(genenamesfilename,DataFrame);
        else
            gene_names_prep = data_df[!,1];
        end
       
        data_df = data_df[!,2:end];
        cells_ptID = [el[2] for el in split.(names(data_df),"_")];
        
        metadata_prep = deepcopy(metadata_df[metadata_df.ID .==ptid,:]);
        data_prep_ = deepcopy(data_df[!,cells_ptID .==ptid]);
        cells_not_missing1 = [in(el,Set(metadata_prep.Column1)) for el in names(data_prep_)];
        data_prep = data_prep_[!, cells_not_missing1];
        data_prep,metadata_prep = raghavan2021_filterCellsWithMetadata(data_prep,metadata_prep;lowerbound_genes=lowerbound_genes,detected_UMIs=detected_UMIs,pct_mito=pct_mito,upperbound_genes=upperbound_genes);
        data_prep,gene_names_prep = raghavan2021_filterGenes(data_prep,gene_names_prep;detected_genes=detected_genes);
        metadata_prep,tumor_stats = raghavan2021_AddLabels!(metadata_prep;tumor_stats_columns = ["Hybrid_Specific1", "sc_Classical1", "sc_Basal1"]);
        # metadata_prep[!,:cell_labels] = metadata_prep[!,"cell.types"];
        dropmissing!(metadata_prep)
        cells_not_missing2 = [in(el,Set(metadata_prep[:,"Column1"])) for el in names(data_prep)];
        data_prep =  data_prep[!, cells_not_missing2];
        metadata_prep = metadata_prep[cells_not_missing2, :];
        cells_not_0580 = [el !== "PANFR0580" for el in metadata_prep.ID];
        data_prep =  data_prep[!, cells_not_0580];
        metadata_prep = metadata_prep[cells_not_0580, :];



        x_int_lognorm_postp,x_int_scaled_postp,metadata,gene_names_lognorm_postp,gene_names_scaled_postp = R_lognorm_FeatSelect_Scale(data_prep,metadata_prep,gene_names_prep;num_var_feat = num_var_feat, scale_fator = scale_fator,only_var_feat=only_var_feat_bool,seed=seed); 

        arg_str_list_prep = @name datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool;
        key_list_prep = naming_vec(arg_str_list_prep);
        var_list_prep = [datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool];

        scaled_data = DataFrame(x_int_scaled_postp, :auto);
        lognorm_data = DataFrame(x_int_lognorm_postp, :auto);
        rename!(scaled_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        rename!(lognorm_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        


        data = scaled_data;
        data_mat = x_int_scaled_postp;
        gene_names = gene_names_scaled_postp;
        dropmissing!(metadata)
        cells_not_missing = [in(el,Set(metadata[:,"Column1"])) for el in names(data)];
        data =  data[!, cells_not_missing];

        G,T,C_t,cell_timepoints_index = generate_modelStatisticsFromDataFrame(data);

        KCalled,unique_clusters_labels,remap_trueMembership,cluster_map = generate_LabelInfoFromMetaDataFrame(metadata);
        x,z,x_scaled,x_lognorm = generate_modelInputsFromDataFrames(data,metadata;scaled_data = scaled_data,lognorm_data = lognorm_data);
        x_used_label_dict = OrderedDict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        N = sum(C_t)
        unique_clusters = collect(1:KCalled)
        x_used_label_dict = Dict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        x_to_use = x;#x_scaled;
        true_z = z;
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ N;
        
        x_used_label = "x_scaled";
        x_used_label_dict[x_used_label] = true;

        scaling_function = "Seurat(R)"
        data_generating_function = "generate_processed_malignant_pt_raghavan2021_data"
        file = split(datafilename, "/")[end]
        metafile = split(metadatafilename, "/")[end]

        collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]
        
        
        data_verbal_summary = "The data used was taken from the patient with id: $ptid.\n G = $G, T = $T, N_t = $C_t ( $(sum(C_t)) cells total), KCalled = $KCalled. \n Number of Highly Variably genes was $num_var_feat and we $(only_var_feat_bool ? " DID " : " DID NOT" ) exclusively use these genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])). \nData was centered and scaled using the function: scaling_function. The data was preprocessed using the function: $data_generating_function   \n The preprocessing was as follow: cells were filtered with less than $lowerbound_genes genes detected, greater than $upperbound_genes genes dectected, greater than $detected_UMIs UMIs detected, or more the $pct_mito mitochondrial genes detected. Genes were filtered if they did not appear in at least $detected_genes cells. Labels were assigned using the procedure outline in Raghavan et. al. 2021. $scaling_function was used to log normalize (scale factor of $scale_fator and pseudocount of 1) and center and scale the data. The final data set used was the $(x_used_label) output from the $data_generating_function  \n The seed used is $seed . "

        
        arg_str_list_prep = @name G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name ptid, file, metafile, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [ptid, file, metafile, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    function generate_subsample_data_processed_raghavan2021_data(datadf,metadatadf,gene_names_prep,unique_time_id,subsetbool_vec,datafilename,metadatafilename,genenamesfilename;is_tuning_set = true,lowerbound_genes=400,detected_UMIs=1000,pct_mito=0.5,upperbound_genes=8000,detected_genes=50,tumor_stats_columns = ["Hybrid_Specific1", "sc_Classical1", "sc_Basal1"],num_var_feat = 5000,only_var_feat_bool = true,scale_fator = 10000,seed = 1989)
        # Random.seed!(seed)
        # data_prep = CSV.read(datafilename,DataFrame);
        # metadata_prep = CSV.read(metadatafilename,DataFrame);
        # gene_names_prep = CSV.read(genenamesfilename,DataFrame);
        # data_prep,metadata_prep = raghavan2021_filterCells(data_prep,metadata_prep;lowerbound_genes=lowerbound_genes,detected_UMIs=detected_UMIs,pct_mito=pct_mito,upperbound_genes=upperbound_genes);
        # data_prep,gene_names_prep = raghavan2021_filterGenes(data_prep,gene_names_prep;detected_genes=detected_genes);
        # metadata_prep,tumor_stats = raghavan2021_AddLabels!(metadata_prep;tumor_stats_columns = tumor_stats_columns);
        data_prep = datadf[!,subsetbool_vec]
        metadata_prep = metadatadf[subsetbool_vec,:]

        x_int_lognorm_postp,x_int_scaled_postp,metadata,gene_names_lognorm_postp,gene_names_scaled_postp = R_lognorm_FeatSelect_Scale(data_prep,metadata_prep,gene_names_prep;num_var_feat = num_var_feat, scale_fator = scale_fator,only_var_feat=only_var_feat_bool,seed=seed); 

        arg_str_list_prep = @name datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool;
        key_list_prep = naming_vec(arg_str_list_prep);
        var_list_prep = [datafilename, seed, unique_time_id, num_var_feat, scale_fator,only_var_feat_bool];

        scaled_data = DataFrame(x_int_scaled_postp, :auto);
        lognorm_data = DataFrame(x_int_lognorm_postp, :auto);
        rename!(scaled_data,Symbol.(metadata[!,:Column1]), makeunique=true);
        rename!(lognorm_data,Symbol.(metadata[!,:Column1]), makeunique=true);



        data = scaled_data;
        data_mat = x_int_scaled_postp;
        gene_names = gene_names_scaled_postp;

        dropmissing!(metadata)
        cells_not_missing = [in(el,Set(metadata[:,"Column1"])) for el in names(data)];
        data =  data[!, cells_not_missing];

        G,T,C_t,cell_timepoints_index = generate_modelStatisticsFromDataFrame(data);
        KCalled,unique_clusters_labels,remap_trueMembership,cluster_map = generate_LabelInfoFromMetaDataFrame(metadata);
        x,z,x_scaled,x_lognorm = generate_modelInputsFromDataFrames(data,metadata;scaled_data = scaled_data,lognorm_data = lognorm_data);
        x_used_label_dict = OrderedDict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        N = sum(C_t)
        unique_clusters = collect(1:KCalled)

        x_to_use = x;#x_scaled;
        true_z = z;
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ sum(C_t);
        
        x_used_label = "x_scaled";
        x_used_label_dict[x_used_label] = true;

        scaling_function = "Seurat(R)"
        data_generating_function = "generate_subsample_data_processed_raghavan2021_data"
        file = split(datafilename, "/")[end]
        metafile = split(datafilename, "/")[end]
        ptid = split(metadata_prep[1,1],"_")[occursin.("PANFR",split(metadata_prep[1,1],"_"))][1]

        collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]
        subset_prefix = ""
        if is_tuning_set
            subset_prefix = "This is the TUNING SET. "
        else
            subset_prefix = "This is the TESTING SET. "
        end
        
        data_verbal_summary = subset_prefix*"The data used was taken from the patient with id: $ptid.\n G = $G, T = $T, N_t = $C_t ( $(sum(C_t)) cells total), KTrue = $KCalled. \n Number of Highly Variably genes was $num_var_feat and we $(only_var_feat_bool ? " DID " : " DID NOT" ) exclusively use these genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])). \nData was centered and scaled using the function: scaling_function. The data was preprocessed using the function: $data_generating_function   \n The preprocessing was as follow: cells were filtered with less than $lowerbound_genes genes detected, greater than $upperbound_genes genes dectected, greater than $detected_UMIs UMIs detected, or more the $pct_mito mitochondrial genes detected. Genes were filtered if they did not appear in at least $detected_genes cells. Labels were assigned using the procedure outline in Raghavan et. al. 2021. $scaling_function was used to log normalize (scale factor of $scale_fator and pseudocount of 1) and center and scale the data. The final data set used was the $(collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]) output from the $data_generating_function  \n The seed used is $seed . "

        
        arg_str_list_prep = @name G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,remap_trueMembership,cluster_map,x_used_label_dict,unique_clusters, datafilename,metadatafilename,genenamesfilename,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_lognorm,gene_names,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name ptid, file,metafile, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [ptid, file,metafile, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, num_var_feat, only_var_feat_bool,data_generating_function,scaling_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end

    
    function generate_processed_simulated_copulacorrelatedcounts_data(KTrue,C_t,G,percent_important,high_de,high_nonde,unique_time_id;corr_feat = true,one_corr_mat = true,mix_prob =nothing,same_prob_t = true,dynamic = false,copula_epsilon=1,random_locations=false,low_de=0,random_deg_values=false,low_nonde=0,random_nondeg_values=false,used_x= "x",seed = 1989)
        # metadatafilename =  "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/Primary_Organoids_metadata_210224_32060Cells.csv" # "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/complete_MetaData_70170cells.csv"
        # datafilename = "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/Primary_Organoids_RawDGE_210224_32060Cells.csv"#
        # ptid = "PANFR0489"#"PANFR0575"#
        Random.seed!(seed)
        unique_clusters = collect(1:KTrue)
        unique_clusters_labels = ["Cluster $k" for k in 1:KTrue]
        gene_names = ["Gene $j" for j in 1:G];
        μk = [zeros(G) for k in 1:KTrue];
        T = length(C_t);
        N = sum(C_t)
        num_important_features = generate_number_of_DEG(G,percent_important);
        deg_locations = generate_DEG_locations(G,KTrue,num_important_features;random_locations=random_locations);
        μk = generate_DEG_values!(μk,deg_locations,high_de;low_de=low_de,random_deg_values=random_deg_values);
        μk = generate_nonDEG_values!(μk,deg_locations,high_nonde;low_nonde=low_nonde,random_nondeg_values=random_nondeg_values);
        Σk= generate_rand_covariance_matrix(G, KTrue;corr_feat = corr_feat,one_corr_mat = one_corr_mat);
        cell_timepoints_index = tidy_get_timeranges(C_t);
        x,z,π_ = fake_mvCountsViaCopulas_corrdata_for_testing(C_t,KTrue,μk ,Σk; mix_prob =mix_prob,same_prob_t = same_prob_t,dynamic = dynamic,epsilon=copula_epsilon);
        x_scaled = raghavan2021_centerAndScale(x);
        x_shifted = mode_shift(x_scaled) ;
        x_used_ = OrderedDict( k => v for (k,v) in zip(["x","x_scaled","x_shifted"],[x,x_scaled,x_shifted]));
        x_used_data_gen_func = OrderedDict( k => v for (k,v) in zip(["x","x_scaled","x_shifted"],["fake_mvCountsViaCopulas_corrdata_for_testing ","fake_mvCountsViaCopulas_corrdata_for_testing  ==> raghavan2021_centerAndScale ","fake_mvCountsViaCopulas_corrdata_for_testing  ==> raghavan2021_centerAndScale ==> mode_shift"]));




        arg_str_list_prep = @name seed, unique_time_id, KTrue,C_t,G,percent_important,high_de,high_nonde,corr_feat,one_corr_mat,mix_prob ,same_prob_t,dynamic,copula_epsilon,random_locations,low_de,random_deg_values,low_nonde,random_nondeg_values;
        key_list_prep = naming_vec(arg_str_list_prep);
        var_list_prep = [seed, unique_time_id, KTrue,C_t,G,percent_important,high_de,high_nonde,corr_feat,one_corr_mat,mix_prob ,same_prob_t,dynamic,copula_epsilon,random_locations,low_de,random_deg_values,low_nonde,random_nondeg_values];

        

        


        x_used_label_dict = OrderedDict( key => false for key in ["x","x_scaled","x_shifted"]);
        x_used_label_dict[used_x] = true;

        scaling_function = "raghavan2021_centerAndScale"
        shifting_function = "mode_shift"
        data_generating_function = x_used_data_gen_func[used_x]

        
        collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in  ["x","x_scaled","x_shifted"]]][1]
        x_to_use=x_used_[used_x]
        
        data_verbal_summary = "The data used was simulated.\n G = $G, T = $T, N_t = $C_t ( $(N) cells total), KCalled = $KTrue. \n Number of driver genes was $num_important_features and we used all genes for analysis. This is $percent_important of the genes that are important. This location of the differentially expressed genes $(random_locations ? " WAS " : " WAS NOT" )  random.The High of the differentially expressed genes are $high_de and the low was set to $low_de. We  $(random_deg_values ? " DID " : " DID NOT" ) have random values for the differentially expressed genes.   The High of the non-differentially expressed genes are $high_nonde and the low was set to $low_nonde. The copula_epsilon was set to $copula_epsilon. We  $(random_nondeg_values ? " DID " : " DID NOT" ) have random values for the non-differentially expressed genes. We  $(corr_feat ? " DID " : " DID NOT" ) have Correlation across the genes. We  $(one_corr_mat ? " DID " : " DID NOT" ) have the same correlation matrix across the clusters. The mixture weights $(random_locations ? " WERE " : " WERE NOT" )  uniform.  The mixture weights $(dynamic ? " DID " : " DID NOT" )  change across time points. That is the final number of genes used was $(length(x_to_use[1][1])).The final data set used was the $(used_x) output from the $data_generating_function.  \n The centering and scaling function: $scaling_function. The shifting_function function: $shifting_function. But the data was preprocessed using the series of functions: $data_generating_function   \n The preprocessing was as follow: descrete counts were transformed using the procedure outline in Raghavan et. al. 2021 set to default values.  \n The seed used is $seed . "

        

        arg_str_list_prep = @name seed, unique_time_id, KTrue,G,T,C_t,N,μk,Σk,deg_locations,cell_timepoints_index,percent_important,high_de,high_nonde,corr_feat,one_corr_mat,mix_prob ,same_prob_t,dynamic,copula_epsilon,random_locations,low_de,random_deg_values,low_nonde,random_nondeg_values,unique_clusters_labels,unique_clusters,used_x,x_used_label_dict,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [seed, unique_time_id, KTrue,G,T,C_t,N,μk,Σk,deg_locations,cell_timepoints_index,percent_important,high_de,high_nonde,corr_feat,one_corr_mat,mix_prob ,same_prob_t,dynamic,copula_epsilon,random_locations,low_de,random_deg_values,low_nonde,random_nondeg_values,unique_clusters_labels,unique_clusters,used_x,x_used_label_dict,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_shifted,gene_names,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_shifted,gene_names,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name seed, G, KTrue, T, C_t,N,unique_clusters_labels,data_generating_function,scaling_function,shifting_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [seed, G, KTrue, T, C_t,N,unique_clusters_labels,data_generating_function,scaling_function,shifting_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    function generate_subsample_data_processed_simulated_data(x_input,z_input,unique_time_id,subsetbool_vec;cell_timepoints_index=nothing,used_x="x_shifted",is_tuning_set = true,seed = 1989)
        Random.seed!(seed)
        # data_prep = CSV.read(datafilename,DataFrame);
        # metadata_prep = CSV.read(metadatafilename,DataFrame);
        # gene_names_prep = CSV.read(genenamesfilename,DataFrame);
        # data_prep,metadata_prep = raghavan2021_filterCells(data_prep,metadata_prep;lowerbound_genes=lowerbound_genes,detected_UMIs=detected_UMIs,pct_mito=pct_mito,upperbound_genes=upperbound_genes);
        # data_prep,gene_names_prep = raghavan2021_filterGenes(data_prep,gene_names_prep;detected_genes=detected_genes);
        # metadata_prep,tumor_stats = raghavan2021_AddLabels!(metadata_prep;tumor_stats_columns = tumor_stats_columns);
    
        T = length(x_input)
        G = length(x_input[1][1])
        C_t = length.(x_input)
        if isnothing(cell_timepoints_index)
            cell_timepoints_index = tidy_get_timeranges(C_t);
        end
        timepoint_index = [[st,en] for (st,en) in  cell_timepoints_index]
        x = Vector{Vector{Vector{eltype(x_input[1][1])}}}(undef,T)
        z = Vector{Vector{eltype(z_input[1][1])}}(undef,T)
        for t in 1:T
            subset_cells_indx = subsetbool_vec[timepoint_index[t][1]:timepoint_index[t][2]]
            num_cells_ = sum(subsetbool_vec[timepoint_index[t][1]:timepoint_index[t][2]])
            # x[t] = Vector{Vector{eltype(x_input[1][1])}}(undef,num_cells_)
            # z[t] = Vector{eltype(z_input[1][1])}(undef,num_cells_)
            # Should check that num_cells_ <= C_t[t]
            cells_indx = collect(1:C_t[t])[subset_cells_indx]
            x[t] = x_input[t][cells_indx]
            z[t] = z_input[t][cells_indx]
            # for i in 1:C_t[t]
            #     if subset_cells_indx[i]
            #         x[t][i] = x_input[t][i]
            #         z[t][i] = z_input[t][i]
            #     end
            # end
        end
        
        # data_prep = datadf[!,subsetbool_vec]
        # metadata_prep = metadatadf[subsetbool_vec,:]


        x_scaled = raghavan2021_centerAndScale(x);
        x_shifted = mode_shift(x_scaled) ;
        x_used_label_dict = OrderedDict( key => false for key in ["x","x_scaled","x_shifted"]);
        x_used_ = OrderedDict( k => v for (k,v) in zip(["x","x_scaled","x_shifted"],[x,x_scaled,x_shifted]));
        x_used_data_gen_func = OrderedDict( k => v for (k,v) in zip(["x","x_scaled","x_shifted"],["From X Input ","From X Input  ==> raghavan2021_centerAndScale ","From X Input  ==> raghavan2021_centerAndScale ==> mode_shift"]));
        N_t = length.(x)
        N = sum(N_t)
        x_to_use=x_used_[used_x];
        x_used_label_dict[used_x] = true;
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ sum(N_t);
        KCalled = unique(vcat(z...))
        unique_clusters = KCalled


    
        scaling_function = "raghavan2021_centerAndScale"
        shifting_function = "mode_shift"
        data_generating_function = x_used_data_gen_func[used_x]
        unique_clusters_labels = ["Cluster $k" for k in KCalled]


        subset_prefix = ""
        if is_tuning_set
            subset_prefix = "This is the TUNING SET. "
        else
            subset_prefix = "This is the TESTING SET. "
        end
        
        data_verbal_summary = subset_prefix*" \n The seed used is $seed ."
        #"The data used was simulated.\n G = $G, T = $T, N_t = $C_t ( $(N) out of $(sum(C_t))  cells total), KCalled = $KCalled. \n We used all genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])).The final data set used was the $(used_x) output from the $data_generating_function.  \n The centering and scaling function: $scaling_function. The shifting_function function: $shifting_function. But the data was preprocessed using the series of functions: $data_generating_function   \n The preprocessing was as follow: descrete counts were transformed using the procedure outline in Raghavan et. al. 2021 set to default values." 

        
        arg_str_list_prep = @name G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,subsetbool_vec,x_used_label_dict,unique_clusters,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,subsetbool_vec,x_used_label_dict,unique_clusters,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_shifted,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_shifted,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name  seed,  KCalled, T, C_t,N,unique_clusters_labels,data_generating_function,scaling_function,shifting_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [seed, KCalled, T, C_t,N,unique_clusters_labels,data_generating_function,scaling_function,shifting_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    function subsample_tuning_set(cell_labels;min_pct_subset= 0.3)

        N = length(cell_labels)
        min_cells = ceil(Int, N * min_pct_subset)
        π__ = [countmap(cell_labels)[k] for k in keys(sort(countmap(cell_labels)))] ./ N;
        amnt_per_label = ceil.(Int,π__ .* min_cells)
        unique_labels = sort(unique(cell_labels))
        num_unique_labels = length(unique_labels)
        state_linindx_vec = [ (state,indx) for (state,indx) in zip(cell_labels,collect(1:N))]
        sample_indx = []
        for i in 1:num_unique_labels
            # println(i)
            label = unique_labels[i]
            lin_label_bool = label .== first.(state_linindx_vec)#[ for el in state_linindx_vec]
            s = sample(state_linindx_vec[lin_label_bool],amnt_per_label[i], replace = false)
            push!(sample_indx,s)
        end
        subsampled_cells_idnx = last.(recursive_flatten(sample_indx))
        subsampled_cells_bool = [in(el, Set(subsampled_cells_idnx)) for el in last.(state_linindx_vec)]
        not_subsampled_cells_bool = broadcast(!,subsampled_cells_bool)
        return Bool.(subsampled_cells_bool),Bool.(not_subsampled_cells_bool)
    end

    function generate_processed_1condition_anndata_data(datafilename,unique_time_id,dataset_id;scale_fator = 10000,seed = 1989)
        # metadatafilename =  "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/NormalCells_metadata_210305.csv" 
        # datafilename = "/mnt/e/cnwizu/Playground/datasets/winterTMEandPlasticity2021/NormalCells_RawDGE_210305.csv"#
        Random.seed!(seed)
        fid = h5open(datafilename,"r")
        anndata_dict = read(fid)
        x_mat = anndata_dict["X"] 
        N = size(x_mat)[2]
        G = size(x_mat)[1]
        gene_ids = anndata_dict["var"]["_index"]
        cell_ids = anndata_dict["obs"]["_index"]
        cell_labels = anndata_dict["obs"]["cell_type"]["codes"]
        if iszero(minimum(cell_labels))
            cell_labels = cell_labels .+ 1
        end
        unique_clusters_labels = anndata_dict["obs"]["cell_type"]["categories"]
        KCalled = length(unique(cell_labels))
        num_cnts = anndata_dict["obs"]["total_counts"]
        x = nothing
        z = nothing
        C_t = nothing
        numi = nothing
        if haskey(anndata_dict["obs"],"conditions")
            T = length(unique(anndata_dict["obs"]["conditions"]))
            ordered_condition = OrderedDict(t=> sort(unique(anndata_dict["obs"]["conditions"]))[t] for t in 1:T)
            # Dict( t=> sort(unique(anndata_dict["obs"]["cell_type"]["codes"]))[t] for t in 1:T)
            condition_indx = OrderedDict(t => [] for t in 1:T)
            for i in 1:N
                curr_condition = anndata_dict["obs"]["conditions"][i]
                condition_order = ordered_condition[curr_condition]
                push!(condition_indx[condition_order],i)
            end
            x = Vector{Vector{Vector{Float64}}}(undef,T)
            z = Vector{Vector{Int}}(undef,T)
            numi = Vector{Vector{Float64}}(undef,T)
            C_t = Vector{Float64}(undef,T)
            for t in 1:T
                C_t[t] = length(condition_indx[t])
                x[t] = [Float64.(collect(col)) for col in eachcol(x_mat[:,condition_indx[t]])]
                z[t] = Int.(collect(cell_labels[condition_indx[t]]))
                numi[t] = Float64.(collect(num_cnts[condition_indx[t]]))
            end
        else
            T = 1
            x = Vector{Vector{Vector{Float64}}}(undef,T)
            z = Vector{Vector{Int}}(undef,T)
            C_t = Vector{Float64}(undef,T)
            numi = Vector{Vector{Float64}}(undef,T)
            C_t[1] = N
            x[1] = [Float64.(collect(col)) for col in eachcol(x_mat)]
            z[1] = Int.(collect(cell_labels))
            numi[1] = Float64.(collect(num_cnts))
        end
        
        x_to_use = mode_shift(raghavan2021_centerAndScale(raghavan2021_lognormalization(x;scaling_factor=scale_fator,pseudocount=1.0,numi=numi))) 
        

        unique_clusters = collect(1:KCalled)

        x_used_label_dict = Dict( key => false for key in ["x_scaled","x_rawcounts","x_lognorm"]);
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ N;
        
        x_used_label = "x_scaled";
        x_used_label_dict[x_used_label] = true;

        scaling_function = "raghavan2021_centerAndScale"
        data_generating_function = "preprocessQC.py"
        file = split(datafilename, "/")[end]


        
        
        only_var_feat_bool=true
        data_verbal_summary = "The data used was taken from the dataset: $dataset_id.\n G = $G, T = $T, N_t = $C_t ( $(sum(C_t)) cells total), KCalled = $KCalled. \n Number of Highly Variably genes was $G and we $(only_var_feat_bool ? " DID " : " DID NOT" ) exclusively use these genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])). \nData was centered and scaled using the function: scaling_function. The data was preprocessed using the function: $data_generating_function   \n The preprocessing was as follows Maddy's pipeline. Labels were assigned prior to this step using the preprocessing pipeline. $scaling_function was used to log normalize (scale factor of $scale_fator and pseudocount of 1) and center and scale the data. We then recentered the data such that the mode oof the ddataset was 0. The final data set used was the $(collect(keys(x_used_label_dict))[[x_used_label_dict[key] for key in ["x_scaled","x_rawcounts","x_lognorm"]]][1]) output from the $data_generating_function  \n The seed used is $seed . "

        
        arg_str_list_prep = @name G,T,C_t,N,KCalled,unique_clusters_labels,x_used_label_dict,unique_clusters, datafilename,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,KCalled,unique_clusters_labels,x_used_label_dict,unique_clusters, datafilename,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,π_,gene_ids,cell_ids,cell_labels,unique_clusters_labels;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,π_,gene_ids,cell_ids,cell_labels,unique_clusters_labels];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name dataset_id, file, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, data_generating_function,scaling_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [dataset_id, file, seed, G, KCalled, T, C_t,N,unique_clusters_labels,scale_fator, data_generating_function,scaling_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    function generate_subsample_data_processed_1condition_anndata_data(x_input,z_input,unique_time_id,subsetbool_vec;cell_timepoints_index=nothing,used_x="x_shifted",is_tuning_set = true,seed = 1989)
        Random.seed!(seed)
        # data_prep = CSV.read(datafilename,DataFrame);
        # metadata_prep = CSV.read(metadatafilename,DataFrame);
        # gene_names_prep = CSV.read(genenamesfilename,DataFrame);
        # data_prep,metadata_prep = raghavan2021_filterCells(data_prep,metadata_prep;lowerbound_genes=lowerbound_genes,detected_UMIs=detected_UMIs,pct_mito=pct_mito,upperbound_genes=upperbound_genes);
        # data_prep,gene_names_prep = raghavan2021_filterGenes(data_prep,gene_names_prep;detected_genes=detected_genes);
        # metadata_prep,tumor_stats = raghavan2021_AddLabels!(metadata_prep;tumor_stats_columns = tumor_stats_columns);
    
        T = length(x_input)
        G = length(x_input[1][1])
        C_t = length.(x_input)
        if isnothing(cell_timepoints_index)
            cell_timepoints_index = tidy_get_timeranges(C_t);
        end
        timepoint_index = [[st,en] for (st,en) in  cell_timepoints_index]
        x = Vector{Vector{Vector{eltype(x_input[1][1])}}}(undef,T)
        z = Vector{Vector{eltype(z_input[1][1])}}(undef,T)
        for t in 1:T
            subset_cells_indx = subsetbool_vec[timepoint_index[t][1]:timepoint_index[t][2]]
            num_cells_ = sum(subsetbool_vec[timepoint_index[t][1]:timepoint_index[t][2]])
            # x[t] = Vector{Vector{eltype(x_input[1][1])}}(undef,num_cells_)
            # z[t] = Vector{eltype(z_input[1][1])}(undef,num_cells_)
            # Should check that num_cells_ <= C_t[t]
            cells_indx = collect(1:C_t[t])[subset_cells_indx]
            x[t] = x_input[t][cells_indx]
            z[t] = z_input[t][cells_indx]
            # for i in 1:C_t[t]
            #     if subset_cells_indx[i]
            #         x[t][i] = x_input[t][i]
            #         z[t][i] = z_input[t][i]
            #     end
            # end
        end
        
        # data_prep = datadf[!,subsetbool_vec]
        # metadata_prep = metadatadf[subsetbool_vec,:]


        x_scaled = raghavan2021_centerAndScale(x);
        x_shifted = mode_shift(x_scaled) ;
        x_used_label_dict = OrderedDict( key => false for key in ["x","x_scaled","x_shifted"]);
        x_used_ = OrderedDict( k => v for (k,v) in zip(["x","x_scaled","x_shifted"],[x,x_scaled,x_shifted]));
        x_used_data_gen_func = OrderedDict( k => v for (k,v) in zip(["x","x_scaled","x_shifted"],["From X Input ","From X Input  ==> raghavan2021_centerAndScale ","From X Input  ==> raghavan2021_centerAndScale ==> mode_shift"]));
        N_t = length.(x)
        N = sum(N_t)
        x_to_use=x_used_[used_x];
        x_used_label_dict[used_x] = true;
        π_ = [countmap(vcat(z...))[k] for k in keys(sort(countmap(vcat(z...))))] ./ sum(N_t);
        KCalled = unique(vcat(z...))
        unique_clusters = KCalled


    
        scaling_function = "raghavan2021_centerAndScale"
        shifting_function = "mode_shift"
        data_generating_function = x_used_data_gen_func[used_x]
        unique_clusters_labels = ["Cluster $k" for k in KCalled]


        subset_prefix = ""
        if is_tuning_set
            subset_prefix = "This is the TUNING SET. "
        else
            subset_prefix = "This is the TESTING SET. "
        end
        
        data_verbal_summary = subset_prefix*" \n The seed used is $seed ."
        #"The data used was simulated.\n G = $G, T = $T, N_t = $C_t ( $(N) out of $(sum(C_t))  cells total), KCalled = $KCalled. \n We used all genes for analysis. That is the final number of genes used was $(length(x_to_use[1][1])).The final data set used was the $(used_x) output from the $data_generating_function.  \n The centering and scaling function: $scaling_function. The shifting_function function: $shifting_function. But the data was preprocessed using the series of functions: $data_generating_function   \n The preprocessing was as follow: descrete counts were transformed using the procedure outline in Raghavan et. al. 2021 set to default values." 

        
        arg_str_list_prep = @name G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,subsetbool_vec,x_used_label_dict,unique_clusters,unique_time_id,data_verbal_summary;
        key_list_prep = Symbol.(naming_vec(arg_str_list_prep));
        var_list_prep = [G,T,C_t,N,cell_timepoints_index,KCalled,unique_clusters_labels,subsetbool_vec,x_used_label_dict,unique_clusters,unique_time_id,data_verbal_summary];

        data_description_dict = OrderedDict()
        addToDict!(data_description_dict,key_list_prep,var_list_prep);

        arg_str_list_data = @name x_to_use,x,z,x_scaled,x_shifted,π_;
        key_list_data= Symbol.(naming_vec(arg_str_list_data));
        var_list_data = [x_to_use,x,z,x_scaled,x_shifted,π_];

        data_dict = OrderedDict()
        addToDict!(data_dict,key_list_data,var_list_data);


        arg_str_list_summary = @name  seed,  KCalled, T, C_t,N,unique_clusters_labels,data_generating_function,scaling_function,shifting_function;
        key_list_summary = Symbol.(naming_vec(arg_str_list_summary));
        var_list_summary = [seed, KCalled, T, C_t,N,unique_clusters_labels,data_generating_function,scaling_function,shifting_function];

        data_summary_dict = OrderedDict()
        addToDict!(data_summary_dict,key_list_summary,var_list_summary);

        return data_description_dict,data_dict,data_summary_dict
    end
    include(curr_dir*src_dir*"synDataPreprocess.jl")
    include(curr_dir*src_dir*"preprocessing_RCall.jl")
    include(curr_dir*src_dir*"preprocessing_PyCall.jl")
end