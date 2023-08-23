
"""
        raghavan2021_lognormalization(x;scaling_factor=10_000,pseudocount=1.0,numi=nothing)
    This function takes in the data as a nested vector of vectors and performs the log-normalization procedure outline in Raghavan et al. 2021.
"""
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
    generate_modelStatisticsFromDataFrame(data)
Parses the model parrameters from ta Data Frame of the data.
"""  
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

"""
    generate_LabelInfoFromMetaDataFrame(metadata)
    Parses the label information from ta Data Frame of the metadata.

"""  
function generate_LabelInfoFromMetaDataFrame(metadata)
    unique_clusters = sort(unique(vec(Matrix(metadata[!,["cell_labels"]]))))
    total_cell_count = size(metadata)[1]; 
    KCalled = length(unique_clusters)
    KTrue = KCalled
    cluster_map = Dict( k => new_k for (k,new_k) in zip(sort(unique_clusters),collect(1:KTrue)) );
    remap_trueMembership = [cluster_map[vec(Matrix(metadata[!,["cell_labels"]]))[i]] for i in 1: total_cell_count];
    return KCalled,unique_clusters,remap_trueMembership,cluster_map
end

"""
    generate_modelInputsFromDataFrames(data,metadata;scaled_data = nothing,lognorm_data = nothing)
Combines outputs from generate_modelStatisticsFromDataFrame and generate_LabelInfoFromMetaDataFrame
"""  
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

"""
    getExperimentVariablesString(experiment_params_to_save,parametersDict)
Extracts the parameters from NCLUSION object to save in the results dictionary
"""  
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
    save_run(filepath,unique_time_id;outputs=nothing,parametersDict=nothing,initializationDict=nothing,posteriorSummarizationsDict=nothing,dataDescriptionDict=nothing)
Saves results in a JLD2 file
```math

```
"""  
function save_run(filepath,unique_time_id;outputs=nothing,parametersDict=nothing,initializationDict=nothing,posteriorSummarizationsDict=nothing,dataDescriptionDict=nothing)
    println("Saving Results")
    filename = filepath*unique_time_id*"/"*unique_time_id*".jld2"
    jldsave(filename; outputs, parametersDict, initializationDict,posteriorSummarizationsDict,dataDescriptionDict, compress=true)
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