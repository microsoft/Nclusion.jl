module IOUtils
    using DataFrames
    using JSON, JSON3
    using Dates
    using StatsBase
    using OrderedCollections
    using CSV
    using LaTeXStrings, TypedTables, PrettyTables
    using JLD2,FileIO
    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"

    include(curr_dir*src_dir*"TidyUtils.jl")
    using .TidyUtils


    export @name, 
           naming_vec,
           set_current_value,
           addToDict!,
           addToOrderedDict!,
           generate_modelStatisticsFromDataFrame,
           generate_LabelInfoFromMetaDataFrame,
           generate_modelInputsFromDataFrames,
           getExperimentVariablesString, 
           get_unique_time_id,
           nested_dict_to_json,
           nested_dict_to_json3_with_name,
           append_to_json,generate_json_filename, 
           generate_filenameBase, 
           generate_filenameBase_closure, 
           generate_dirname, 
           json_to_dict,json3_to_dictstruct, 
           makeRelabellingOutputDirectories,
           makeGenericOutputDirectories,
           makeTidyOutputDirectories,
           saveInitializationFile,
           create_summarization_file,
           create_benchmark_summarization_file,
           save_run


        
    macro name(arg...)
        x = string(arg)
        quote
            $x
        end
    end
    macro name(arg)
        x = string(arg)
        quote
            $x
        end
    end
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
    function set_current_value(prev_iteration_values, val_to_update_every_itr)
        
        function extract(d)
            expr = quote end
            for (k, v) in d
                if k in val_to_update_every_itr
                    push!(expr.args, :($(Symbol(k)) = $v))
                end
            end
            eval(expr)
            return
        end
        
        
        extract(prev_iteration_values)
        return β,π_, λ,parameterized_liklihood, γ, θa, θb, conc, ρ
    end
    function addToDict!(dict,key_array,val_array)
        num_var = length(key_array)
        for i in 1:num_var
            key = key_array[i]
            val = val_array[i]
            dict[key] = val 
        end
        dict
    end
    function addToOrderedDict!(ordered_dict,key_array,val_array)
        num_var = length(key_array)
        for i in 1:num_var
            key = key_array[i]
            val = val_array[i]
            ordered_dict[key] = val 
        end
        ordered_dict
    end
    

    function generate_modelStatisticsFromDataFrame(data)
        cell_names = names(data)
        cell_timepoints = [el[1] for el in split.(cell_names,"_")];
        unique_timepoints = sort(unique(cell_timepoints))
        cell_timepoints_counts = countmap(cell_timepoints)
        total_cell_count = size(data)[2]
        G = size(data)[1]
        T = length(unique_timepoints)
        C_t = [cell_timepoints_counts[key] for key in sort(collect(keys(cell_timepoints_counts)))]
        timepoint_map = Dict( st => t for (st,t) in zip(sort(unique_timepoints),collect(1:T)) )
        cell_timepoints_index = [[] for t in 1:T]
        for indx in 1: total_cell_count
            tp = cell_timepoints[indx]
            t_indx = timepoint_map[tp]
            push!(cell_timepoints_index[t_indx],indx)
        end
        return G,T,C_t,cell_timepoints_index
    end
    function generate_LabelInfoFromMetaDataFrame(metadata)
        unique_clusters = sort(unique(vec(Matrix(metadata[!,["cell_labels"]]))))
        total_cell_count = size(metadata)[1]; 
        KCalled = length(unique_clusters)
        KTrue = KCalled
        cluster_map = Dict( k => new_k for (k,new_k) in zip(sort(unique_clusters),collect(1:KTrue)) );
        remap_trueMembership = [cluster_map[vec(Matrix(metadata[!,["cell_labels"]]))[i]] for i in 1: total_cell_count];
        return KCalled,unique_clusters,remap_trueMembership,cluster_map
    end
    function generate_modelInputsFromDataFrames(data,metadata;scaled_data = nothing,lognorm_data = nothing)
        G,T,C_t,cell_timepoints_index = generate_modelStatisticsFromDataFrame(data);
        KCalled,unique_clusters,remap_trueMembership,cluster_map = generate_LabelInfoFromMetaDataFrame(metadata);
        x = Vector{Vector{Vector{Float64}}}(undef,T)
        z = Vector{Vector{Int64}}(undef,T)
        x_scaled = nothing
        x_lognorm = nothing
        if !isnothing(scaled_data)
            x_scaled = Vector{Vector{Vector{Float64}}}(undef,T);
        end
        if !isnothing(lognorm_data)
            x_lognorm = Vector{Vector{Vector{Float64}}}(undef,T);
        end
    
        for t in 1:T
            cells = C_t[t]
            z_t = Vector{Int64}(undef,cells)
            x_t = Vector{Vector{Float64}}(undef,cells)
            if !isnothing(scaled_data)
                x_scaled_t = Vector{Vector{Float64}}(undef,cells);
            end
            if !isnothing(lognorm_data)
                x_lognorm_t = Vector{Vector{Float64}}(undef,cells);
            end
            for i in 1:cells
                lin_indx = cell_timepoints_index[t][i]
                z_t[i] = remap_trueMembership[lin_indx]
                # x_t[i] = subset_data[:,lin_indx]'
                x_t[i] = data[!,lin_indx]
                
                if !isnothing(scaled_data)
                    x_scaled_t[i] = scaled_data[!,lin_indx];
                end
                if !isnothing(lognorm_data)
                    x_lognorm_t[i] = lognorm_data[!,lin_indx];
                end
            end
            x[t] = x_t
            z[t] = z_t
            if !isnothing(scaled_data)
                x_scaled[t] = x_scaled_t;
            end
            if !isnothing(lognorm_data)
                x_lognorm[t] = x_lognorm_t;
            end
        end
        
        return x,z,x_scaled,x_lognorm
    end
    
    function getExperimentVariablesString(experiment_params_to_save,parametersDict)
        param_str = "_"
        for el in experiment_params_to_save
            if haskey(parametersDict,el)
                param_str = param_str * el * string(parametersDict[el]) * "_"
            else
                param_str = param_str * el * string(el) * "_"
            end
        end

        return param_str[1:end-1]
    end


    function get_unique_time_id()
        datetimenow = Dates.now(Dates.UTC)
        now_str = string(datetimenow)
        now_str = string(split(now_str, ".")[1])
        r = ":"
        return replace(now_str,r => "" )
    end

    function nested_dict_to_json(d,job_name_stem,unique_time_id, to_save=true)
        json_stringdata = JSON.json(d)
        if to_save
            # unique_time_id = get_unique_time_id() # from playing_with_datetimes.jl
            filename = job_name_stem*unique_time_id*".json"
            open(filename, "w") do f
                write(f, json_stringdata)
            end
        end
        return json_stringdata
    end
    function nested_dict_to_json3_with_name(d,filename)
        # json_stringdata = JSON.json(d)
        open(filename, "w") do io
            JSON3.write(io, d)
        end
    end
    function append_to_json(json_filename,new_data_dict)
        json_obj = JSON3.read(open(json_filename,"r"))
        json_string = JSON3.write(json_obj)
        mod_json_string = json_string[1:end-1]
        new_data_json_string = JSON3.write(new_data_dict)
        mod_new_data_json_string = new_data_json_string[2:end]
        appended_json_string = mod_json_string *","*mod_new_data_json_string 
        appended_json_obj = JSON3.read(appended_json_string)
        open(json_filename, "w") do io
            JSON3.write(io, appended_json_obj)
        end
    end

    function generate_json_filename(job_name_stem,param_str,unique_time_id)
        # unique_time_id = get_unique_time_id()
        filename = job_name_stem*param_str*"."*unique_time_id*".json"
        
        return filename
    end
    function generate_filenameBase(x,job_name_stem,param_str,unique_time_id)
        # unique_time_id = get_unique_time_id()
        filename = job_name_stem*"_"*x*param_str*"."*unique_time_id
        
        return filename
    end


    generate_filenameBase_closure(job_name_stem,param_str,unique_time_id,filepath,file_ext) = (plotID) -> filepath * generate_filenameBase(plotID ,job_name_stem,param_str,unique_time_id) * file_ext
    
    function saveInitializationFile(initialization_filename,parametersDict)
        vals_to_print = sort(collect(keys(parametersDict)))
        open(initialization_filename, "w") do f
            for i in 1:length(vals_to_print)
                var_ =  vals_to_print[i]
                symb_ = var_
                val_ = parametersDict[symb_]
                val_str_ = string(val_)
                line = var_ * ": \t " * val_str_ * "\n"
                write(f,line)
            end
        end
    end
    
    function generate_dirname(dir_name,job_name_stem,param_str)
        filename = dir_name*job_name_stem*param_str
        
        return filename
    end

    function json_to_dict(filename)
        dict = filename |> open |> JSON.parse
        
    end
    function json3_to_dictstruct(filename)
        dict = filename |> open |> JSON3.read;
        
    end
    function makeRelabellingOutputDirectories(curr_dir, topDirName,job_name_stem,param_str,unique_time_id)

        if !isdir(curr_dir *"/outputs/")
            mkdir(curr_dir *"/outputs/")
        end
        dir_path = "/outputs/"*topDirName
        topDirPath = curr_dir * dir_path 
        if !isdir(topDirPath)
            mkdir(topDirPath)
        end
        dir_name = generate_dirname(topDirPath,job_name_stem,param_str)
        if !isdir(dir_name)
            mkdir(dir_name)
        end
        subdir_name = dir_name*"/"*unique_time_id
        if !isdir(subdir_name)
            mkdir(subdir_name)
        end
        filepath = subdir_name*"/"
        thin_original_filepath = filepath*"thinned-original-chain"*"/"
        if !isdir(thin_original_filepath)
            mkdir(thin_original_filepath)
        end
        thin_relabel_filepath = filepath*"thinned-relabelled-chain"*"/"
        if !isdir(thin_relabel_filepath)
            mkdir(thin_relabel_filepath)
        end
        return dir_path,topDirPath,dir_name,subdir_name,filepath,thin_original_filepath,thin_relabel_filepath
    end

    function makeGenericOutputDirectories(curr_dir, topDirName,job_name_stem,param_str,unique_time_id;inferenceModel_dir = nothing)

        if !isdir(curr_dir *"/outputs/")
            mkdir(curr_dir *"/outputs/")
        end
        dir_path = "/outputs/"*topDirName
        topDirPath = curr_dir * dir_path 
        if !isdir(topDirPath)
            mkdir(topDirPath)
        end
        if !isnothing(inferenceModel_dir)
            topDirPath_wInferenceModel = topDirPath*inferenceModel_dir
            if !isdir(topDirPath_wInferenceModel)
                mkdir(topDirPath_wInferenceModel)
            end
            topDirPath = topDirPath_wInferenceModel
        end
        dir_name = generate_dirname(topDirPath,job_name_stem,param_str)
        if !isdir(dir_name)
            mkdir(dir_name)
        end
        subdir_name = dir_name*"/"*unique_time_id
        if !isdir(subdir_name)
            mkdir(subdir_name)
        end
        filepath = subdir_name*"/"

        return dir_path,topDirPath,dir_name,subdir_name,filepath
    end

    function makeTidyOutputDirectories(curr_dir, topDirName,unique_time_id)

        if !isdir(curr_dir *"/outputs/")
            mkdir(curr_dir *"/outputs/")
        end
        dir_path = "/outputs/"*topDirName
        topDirPath = curr_dir * dir_path 
        if !isdir(topDirPath)
            mkdir(topDirPath)
        end
        dir_name = generate_dirname(topDirPath,"Runs","")
        if !isdir(dir_name)
            mkdir(dir_name)
        end
        subdir_name = dir_name*"/"*unique_time_id
        if !isdir(subdir_name)
            mkdir(subdir_name)
        end
        filepath = subdir_name*"/"

        return dir_path,topDirPath,dir_name,subdir_name,filepath
    end

    function create_summarization_file(filepath,unique_time_id,modeltype;ari_RandIndices_summary_invar=nothing,ari_RandIndices_summary_var=nothing,nmi_summary_invar=nothing,nmi_summary_var=nothing,vmeasure_summary_invar=nothing,vmeasure_summary_var=nothing,verbal_summary=nothing,data_verbal_summary=nothing,data_summary_dict=nothing,elbo_summary_dict=nothing,parameter_summary_dict=nothing,run_stats_summary_dict=nothing,tidy_avg_counts=nothing,conflvl=.95)
        filedir = filepath*"Runs/"*unique_time_id*"/"
        sumfiledir = filepath*"SummarizationFiles/"
        summarization_filename = unique_time_id*".txt"
        conf = set_pt_conf(tf = tf_markdown, alignment = :c);
        files = [filedir*summarization_filename,sumfiledir*summarization_filename]
        if !isdir(filepath)
            mkdir(filepath)
        end
        if !isdir(filedir)
            mkdir(filedir)
        end
        if !isdir(sumfiledir)
            mkdir(sumfiledir)
        end
        for filename in files
            open(filename, "w") do f
                # line = var_ * ": \t " * val_str_ * "\n"
                str_prefix ="SUMMARY : \t "
                write(f,"###!$(unique_time_id) \n")
                write(f,"RUN ID : \t $(unique_time_id) \n")
                write(f,"Model Used : \t $(modeltype) \n")
                if !isnothing(verbal_summary)
                    write(f,"WRITTEN SUMMARY : \t "*verbal_summary*"  \n")
                    write(f,"\n")
                end
                if !isnothing(data_verbal_summary)
                    write(f,"DATA SUMMARY : \t "*data_verbal_summary*"  \n")
                    write(f,"\n")
                end
                if !isnothing(data_summary_dict)
                    write(f,str_prefix*"DATA PARAMETERS \n")
                    pretty_table(f,data_summary_dict; header= (["DATA PARAM.", "DATA PARAM. VALUE."],["",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    write(f,"\n")
                end
                if !isnothing(ari_RandIndices_summary_invar) && !isnothing(ari_RandIndices_summary_var)
                    write(f,str_prefix*"ARI METRIC \n")
                    if !isnothing(ari_RandIndices_summary_invar)
                        pretty_table(f,ari_RandIndices_summary_invar; header= (["", "ARI (mean)", "ARI (std)", "ARI (low)", "ARI (upper)"],["","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                        write(f,"\n")
                    end
                    if !isnothing(ari_RandIndices_summary_var)
                        pretty_table(f,ari_RandIndices_summary_var; header= (["Time Point", "ARI (mean)", "ARI (std)", "ARI (low)", "ARI (upper)"],["T","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                        write(f,"\n")
                        write(f,"\n")
                    end
                end
                if !isnothing(nmi_summary_invar) && !isnothing(nmi_summary_var)
                    write(f,str_prefix*"NMI METRIC \n")
                    if !isnothing(nmi_summary_invar)
                        pretty_table(f,nmi_summary_invar; header= (["", "NMI (mean)", "NMI (std)", "NMI (low)", "NMI (upper)"],["","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                        write(f,"\n")
                    end
                    if !isnothing(nmi_summary_var)
                        pretty_table(f,nmi_summary_var; header= (["Time Point", "NMI (mean)", "NMI (std)", "NMI (low)", "NMI (upper)"],["T","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                        write(f,"\n")
                        write(f,"\n")
                    end
                end
    
                if !isnothing(vmeasure_summary_invar) && !isnothing(vmeasure_summary_var)
                    write(f,str_prefix*"VMEASURE METRIC \n")
                    if !isnothing(vmeasure_summary_invar)
                        pretty_table(f,vmeasure_summary_invar; header= (["", "VMeasure (mean)", "VMeasure (std)", "VMeasure (low)", "VMeasure (upper)"],["","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    end
                    if !isnothing(vmeasure_summary_var)
                        pretty_table(f,vmeasure_summary_var; header= (["Time Point", "VMeasure (mean)", "VMeasure (std)", "VMeasure (low)", "VMeasure (upper)"],["T","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                        write(f,"\n")
                        write(f,"\n")
                    end
                end
                if !isnothing(elbo_summary_dict)
                    write(f,str_prefix*"ELBO \n")
                    pretty_table(f,elbo_summary_dict; header= (["ELBO SUMMARY PARAM.", "ELBO SUMMARY VALUE."],["",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    write(f,"\n")
                end
                if !isnothing(parameter_summary_dict)
                    write(f,str_prefix*"PARAMETERS \n")
                    pretty_table(f,parameter_summary_dict; header= (["PARAMETER", "PARAMETER VALUE"],["",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    write(f,"\n")
                end
                if !isnothing(run_stats_summary_dict)
                    write(f,str_prefix*"INFERENCE RUN \n")
                    pretty_table(f,run_stats_summary_dict; header= (["RUN STATISTIC", "VALUE"],["",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    write(f,"\n")
                end
                if !isnothing(tidy_avg_counts)
                    write(f,str_prefix*"AVERAGE POSTERIOR CLUSTER COUNTS \n")
                    pretty_table(f,tidy_avg_counts; header= (["Timepoint", "Called K", "Inferred K", "Average Counts"],["T","KCalled","KMax",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    write(f,"\n")
                end
            end
        end
        

        summary_runner_file = filepath*"QuickRunSummaries"*".txt"
        if !isfile(summary_runner_file)
            open(summary_runner_file, "w") do f
                write(f,"###! \n")
            end
        end
        new_addition_tmp_file = filepath*"tmp.txt"
        open(new_addition_tmp_file, "w") do f
            # line = var_ * ": \t " * val_str_ * "\n"
            str_prefix ="SUMMARY : \t "
            write(f,"###!$(unique_time_id) \n")
            write(f,"RUN ID : \t $(unique_time_id) \n")
            write(f,"Model Used : \t $(modeltype) \n")
            write(f,"Summarization File Location : \t $(sumfiledir*summarization_filename) \n")
            write(f,"Summarization Directory Location : \t $(filedir) \n")
            if !isnothing(verbal_summary)
                write(f,"WRITTEN SUMMARY : \t "*verbal_summary*"  \n")
                write(f,"\n")
            end
            if !isnothing(elbo_summary_dict)
                write(f,str_prefix*"ELBO \n")
                pretty_table(f,elbo_summary_dict; header= (["ELBO SUMMARY PARAM.", "ELBO SUMMARY VALUE."],["",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                write(f,"\n")
                write(f,"\n")
            end
            if !isnothing(ari_RandIndices_summary_invar)
                write(f,str_prefix*"ARI METRIC \n")
                pretty_table(f,ari_RandIndices_summary_invar; header= (["", "ARI (mean)", "ARI (std)", "ARI (low)", "ARI (upper)"],["","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                write(f,"\n")
            end
            if !isnothing(nmi_summary_invar)
                write(f,str_prefix*"NMI METRIC \n")
                pretty_table(f,nmi_summary_invar; header= (["", "NMI (mean)", "NMI (std)", "NMI (low)", "NMI (upper)"],["","","","$(conflvl*100)% Level","$(conflvl*100)% Level"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                write(f,"\n")
            end
            if !isnothing(run_stats_summary_dict)
                write(f,str_prefix*"INFERENCE RUN \n")
                pretty_table(f,run_stats_summary_dict; header= (["RUN STATISTIC", "VALUE"],["",""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
                write(f,"\n")
                write(f,"\n")
            end
            if !isnothing(data_verbal_summary)
                write(f,"DATA SUMMARY : \t "*data_verbal_summary*"  \n")
                write(f,"\n")
            end
            write(f,"############################################ END OF RUN #######################################################\n")
        end
        str_buffer_array = []
        tmp_array = open(new_addition_tmp_file)
        new_lines = readlines(tmp_array)
        close(tmp_array)
        str_buffer_array = vcat(new_lines,str_buffer_array)
        summary_runner_array = open(summary_runner_file,"r")
        old_lines = readlines(summary_runner_array)
        close(summary_runner_array)
        str_buffer_array = vcat(str_buffer_array,old_lines)
        open(summary_runner_file, "w") do f
            for el in str_buffer_array
                write(f,el*"\n")
            end
        end
        rm(new_addition_tmp_file)
    end
    function save_run(filepath,unique_time_id;outputs=nothing,parametersDict=nothing,initializationDict=nothing,posteriorSummarizationsDict=nothing,dataDescriptionDict=nothing)
        println("Saving Results")
        filename = filepath*unique_time_id*"/"*unique_time_id*".jld2"
        jldsave(filename; outputs, parametersDict, initializationDict,posteriorSummarizationsDict,dataDescriptionDict, compress=true)
    end

    function create_benchmark_summarization_file(filepath,unique_time_id,modeltype,results_dict,verbal_summary;is_summary_dict=false)
        filedir = filepath*"Runs/"*unique_time_id*"/"
        sumfiledir = filepath*"SummarizationFiles/"
        summarization_filename = unique_time_id*".txt"
        conf = set_pt_conf(tf = tf_markdown, alignment = :c);
        files = [filedir*summarization_filename,sumfiledir*summarization_filename]
        if !isdir(filepath)
            mkdir(filepath)
        end
        if !isdir(filedir)
            mkdir(filedir)
        end
        if !isdir(sumfiledir)
            mkdir(sumfiledir)
        end
        function remove_times_from_results_dict(r_dict,key)
            cp_dict = deepcopy(r_dict[key])
            delete!(cp_dict,:times)
            return cp_dict
        end
        function reformat_results_dict(r_dict)
            new_r_dict = OrderedDict( key => r_dict[key][1]  for key in keys(r_dict))
            return new_r_dict
        end
        for filename in files
            open(filename, "w") do f
                # line = var_ * ": \t " * val_str_ * "\n"
                str_prefix ="SUMMARY : \t "
                write(f,"###!$(unique_time_id) \n")
                write(f,"RUN ID : \t $(unique_time_id) \n")
                write(f,"Model Used : \t $(modeltype) \n")
                if !isnothing(verbal_summary)
                    write(f,"WRITTEN SUMMARY : \t "*verbal_summary*"  \n")
                    write(f,"\n")
                end
                result_keys = collect(keys(results_dict))
                write(f,str_prefix*"$(String(result_keys[1])) \n")
                pretty_table(f,results_dict[result_keys[1]];backend = Val(:text),tf = tf_markdown, alignment = :c)
                write(f,"\n")
                write(f,"\n")
                for key in result_keys[2:end]
                    r_dict = remove_times_from_results_dict(results_dict,key)
                    # cp_dict = deepcopy(results_dict[key])
                    # delete!(cp_dict,:times)
                    if is_summary_dict
                        cp_dict = r_dict
                    else
                        cp_dict = reformat_results_dict(r_dict)
                    end
                    write(f,str_prefix*"$(String(key)) \n")
                    pretty_table(f,cp_dict;backend = Val(:text),tf = tf_markdown, alignment = :c)
                    write(f,"\n")
                    write(f,"\n")
                end
            end
        end
        

        summary_runner_file = filepath*"QuickRunSummaries"*".txt"
        if !isfile(summary_runner_file)
            open(summary_runner_file, "w") do f
                write(f,"###! \n")
            end
        end
        new_addition_tmp_file = filepath*"tmp.txt"
        open(new_addition_tmp_file, "w") do f
            # line = var_ * ": \t " * val_str_ * "\n"
            str_prefix ="SUMMARY : \t "
            write(f,"###!$(unique_time_id) \n")
            write(f,"RUN ID : \t $(unique_time_id) \n")
            write(f,"Model Used : \t $(modeltype) \n")
            write(f,"Summarization File Location : \t $(sumfiledir*summarization_filename) \n")
            write(f,"Summarization Directory Location : \t $(filedir) \n")
            if !isnothing(verbal_summary)
                write(f,"WRITTEN SUMMARY : \t "*verbal_summary*"  \n")
                write(f,"\n")
            end
            result_keys = collect(keys(results_dict))
            write(f,str_prefix*"$(String(result_keys[1])) \n")
            pretty_table(f,results_dict[result_keys[1]];backend = Val(:text),tf = tf_markdown, alignment = :c)
            write(f,"\n")
            write(f,"\n")
            max_avg_time_key = result_keys[2:end][sortperm([results_dict[key][:avg_time][1] for key in result_keys[2:end]], rev= true)][1]
            write(f,str_prefix*"MAXIMUM AVG. TIME--> $(results_dict[max_avg_time_key][:name][1]) \n")
            pretty_table(f,results_dict[max_avg_time_key][:avg_time];header= (["AVG. TIME"],["ns"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
            write(f,"\n")
            write(f,"\n")
            max_med_time_key = result_keys[2:end][sortperm([results_dict[key][:med_time][1] for key in result_keys[2:end]], rev= true)][1]
            write(f,str_prefix*"MAXIMUM MEDIAN TIME--> $(results_dict[max_med_time_key][:name][1]) \n")
            pretty_table(f,results_dict[max_med_time_key][:med_time];header= (["MEDIAN TIME"],["ns"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
            write(f,"\n")
            write(f,"\n")
            max_num_alloc_key = result_keys[2:end][sortperm([results_dict[key][:num_alloc][1] for key in result_keys[2:end]], rev= true)][1]
            write(f,str_prefix*"MAXIMUM NUM. ALLOC.--> $(results_dict[max_num_alloc_key][:name][1]) \n")
            pretty_table(f,results_dict[max_num_alloc_key][:num_alloc];header= (["ALLOCATIONS"],[""]),backend = Val(:text),tf = tf_markdown, alignment = :c)
            write(f,"\n")
            write(f,"\n")
            max_memory_key = result_keys[2:end][sortperm([results_dict[key][:memory][1] for key in result_keys[2:end]], rev= true)][1]
            write(f,str_prefix*"MAXIMUM MEMORY USAGE--> $(results_dict[max_memory_key][:name][1]) \n")
            pretty_table(f,results_dict[max_memory_key][:memory];header= (["MEMORY USAGE"],["bytes"]),backend = Val(:text),tf = tf_markdown, alignment = :c)
            write(f,"\n")
            write(f,"\n")
            write(f,"############################################ END OF RUN #######################################################\n")
        end
        str_buffer_array = []
        tmp_array = open(new_addition_tmp_file)
        new_lines = readlines(tmp_array)
        close(tmp_array)
        str_buffer_array = vcat(new_lines,str_buffer_array)
        summary_runner_array = open(summary_runner_file,"r")
        old_lines = readlines(summary_runner_array)
        close(summary_runner_array)
        str_buffer_array = vcat(str_buffer_array,old_lines)
        open(summary_runner_file, "w") do f
            for el in str_buffer_array
                write(f,el*"\n")
            end
        end
        rm(new_addition_tmp_file)
    end
    
    
end