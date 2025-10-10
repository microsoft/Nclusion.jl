
"""
        recursive_flatten(x::AbstractArray)
This function takes an arbitrarily nested set of vectors and recursively flattens them into one 1-D vector
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
This function recursively performs an inner melt of a vector input
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
    initialize_dict(;key_type=Any,val_type=Any)
This function initializes a dictionary with specified key and value types.
"""
function initialize_dict(;key_type=Any,val_type=Any)
    dict = Dict{key_type,val_type}()
    return dict
end

"""
    initialize_ordered_dict(;key_type=Any,val_type=Any)
This function initializes an ordered dictionary with specified key and value types.
"""
function initialize_ordered_dict(;key_type=Any,val_type=Any)
    od = OrderedDict{key_type,val_type}()
    return od
end

# """
#     addToOrderedDict!(ordered_dict,key_array,val_array)
# Adds a set of values to a previously initialized ordered dictionary
# """
# function addToOrderedDict!(ordered_dict,key_array,val_array)
#     num_var = length(key_array)
#     for i in 1:num_var
#         key = key_array[i]
#         val = val_array[i]
#         ordered_dict[key] = val
#     end
#     ordered_dict
# end
"""
    @add_variables_to_ordered_dict!(od, vars...)
This macro adds variables to an ordered dictionary with their names as keys and their values as values.
"""
macro add_variables_to_ordered_dict!(od, vars...)
    od_expr = esc(od)
    exprs = map(var -> :(push!($od_expr, $(string(var)) => $(esc(var)))), vars)
    return Expr(:block, exprs...)  # Combine the expressions into a single block
end

"""
    @add_variables_to_ordered_dict_as_string!(od, vars...)
This macro adds variables to an ordered dictionary with their names as keys and their string values as values.
"""
macro add_variables_to_ordered_dict_as_string!(od, vars...)
    od_expr = esc(od)
    exprs = map(var -> :(push!($od_expr, $(string(var)) => string($(esc(var))))), vars)
    return Expr(:block, exprs...)  # Combine the expressions into a single block
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
"""      
function setup_experiment_tag(experiment_filename)
    return "EXPERIMENT_$experiment_filename"
end


"""
    load_data(datafilename1,seed)
This function loads data from an HDF5 file.
"""
function load_data(datafilename1,seed)
    fid1 = h5open(datafilename1,"r")
    anndata_dict1 = read(fid1)
    return anndata_dict1
end


"""
    preparing_data(anndata_dict1;time_key=nothing,individuals_key=nothing)
This function prepares data from an HDF5 file.
"""
function preparing_data(anndata_dict1;time_key=nothing,individuals_key=nothing)
    gene_names = anndata_dict1["var"]["_index"]
    cell_ids = anndata_dict1["obs"]["_index"]
    cell_cluster_dict = nothing
    # time_ids = nothing
    time_vec = nothing
    # individuals_ids = nothing
    individuals_vec = nothing
    cell_cluster_labels = nothing
    if time_key != nothing
        if haskey(anndata_dict1["obs"],time_key)
            # time_ids = sort(unique(anndata_dict1["obs"][time_key]["codes"]))
            time_vec  = anndata_dict1["obs"][time_key]["codes"]
        end
    end
    if individuals_key != nothing
        if haskey(anndata_dict1["obs"],individuals_key)
            # individuals_ids = sort(unique(anndata_dict1["obs"][individuals_key]["codes"]))
            individuals_vec  = anndata_dict1["obs"][individuals_key]["codes"]
        end
    end
    if haskey(anndata_dict1["obs"],"cell_type")
        cell_cluster_labels = anndata_dict1["obs"]["cell_type"]["codes"]
        cell_cluster_labels = cell_cluster_labels .+ 1
        cell_cluster_dict = Dict(zip(anndata_dict1["obs"]["cell_type"]["categories"],1:length(unique(cell_cluster_labels))))
    end
    if isnothing(time_vec)
        # time_ids = Int8.([1])
        time_vec = ones(Int8,length(cell_ids))
    end
    if isnothing(individuals_vec)
        # individuals_ids = Int8.([1])
        individuals_vec = ones(Int8,length(cell_ids))
    end  
    # return gene_names,cell_ids, cell_cluster_labels, time_ids,time_vec,individuals_ids,individuals_vec
    return gene_names,cell_ids, cell_cluster_labels,time_vec,individuals_vec,cell_cluster_dict
end


"""
    make_nclusion_inputs(anndata_dict1;time_key=nothing,individuals_key=nothing, layer_name=nothing,layer_index=0,gene_set_file_path=nothing, min_genes_detected=10,standardization_of_used_representation=nothing,is_precomputed_latent_representation=false)
This function prepares inputs for the Nclusion model.
"""
function make_nclusion_inputs(anndata_dict1;time_key=nothing,individuals_key=nothing, layer_name=nothing,layer_index=0,gene_set_file_path=nothing, min_genes_detected=10,standardization_of_used_representation=nothing,is_precomputed_latent_representation=false)
    x_mat = anndata_dict1["X"]
    # gene_names, _, cell_cluster_labels, _ ,time_vec, _ ,individuals_vec,cell_cluster_dict=preparing_data(anndata_dict1;time_key=time_key,individuals_key=individuals_key)
    # gene_names, _, cell_cluster_labels, time_vec, individuals_vec,cell_cluster_dict=preparing_data(anndata_dict1;time_key=time_key,individuals_key=individuals_key)
    gene_names, cell_ids, cell_cluster_labels, time_vec, individuals_vec,cell_cluster_dict=preparing_data(anndata_dict1;time_key=time_key,individuals_key=individuals_key)
    new_order_features = sortperm(gene_names)
    # testing this approach 
    # time_vec_ = [ 1, 1 , 5, 2, 1, 2, 2, 4, 1, 1, 1, 4, 1, 3, 3, 1, 3, 2, 3]; individuals_vec_ = [1 ,1,3, 4, 1, 5, 1, 1, 2, 3, 3, 3, 4, 4, 3, 5, 5, 3, 1]; J_ = 5; N_ = length(individuals_vec_); x_mat_ = rand(J_, N_); sorting_keys_ = [(individuals_vec_[i], time_vec_[i]) for i in 1:N_];sorted_indices_ = sortperm(sorting_keys_); sorted_sorting_keys_ = sorting_keys_[sorted_indices_]; x_mat_sorted_ = x_mat_[:, sorted_indices_]
    sorting_keys = [(individuals_vec[i], time_vec[i], cell_ids[i]) for i in 1:length(cell_ids)]
    new_order_samples = sortperm(sorting_keys)
    gene_names = gene_names[new_order_features]
    gene_names_order_dict = OrderedDict(zip(gene_names,1:length(gene_names)))
    time_vec = time_vec[new_order_samples]
    individuals_vec = individuals_vec[new_order_samples]
    sorted_sorting_keys = [(el[1],el[2],el[3],n) for (n,el) in enumerate(sorting_keys[new_order_samples])]
    individuals_ids = sort(unique([el[1] for el in sorted_sorting_keys]))
    I = length(individuals_ids)
    time_ids = [sort(unique([el[2] for el in sorted_sorting_keys if el[1] == individuals_ids[i]])) for i in 1:I]
    individual_time_combinations_counts = countmap([(el[1],el[2]) for el in sorted_sorting_keys])
    total_unique_combinations = sort(collect(keys(individual_time_combinations_counts)))
    if !is_precomputed_latent_representation
        representation_key = "layers"
    else
        representation_key = "obsm"
    end
    if haskey(anndata_dict1,representation_key) && isnothing(layer_name)
        layer_names = collect(keys(anndata_dict1[representation_key]))
        if layer_index != 0
            layer_name = layer_names[layer_index]
        end
    end
    used_representation_feature_name= nothing
    if !isnothing(cell_cluster_labels)
        cell_cluster_labels = cell_cluster_labels[new_order_samples]
    end
    if Base.lowercase(layer_name) == "pca"
        maxpscs = 50
        M = fit(PCA, x_mat[new_order_features,new_order_samples]; maxoutdim=maxpscs)
        selected_representation = predict(M, x_mat[new_order_features,new_order_samples])
        gene_factor_loadings = DataFrame(permutedims(projection(M)),:auto);
        rename!(gene_factor_loadings,Symbol.(gene_names))
        used_representation_feature_name = ["PC$i" for i in 1:maxpscs]
        insertcols!(gene_factor_loadings,1,:feature_name =>used_representation_feature_name)
        alternative_representation = (gene_factor_loadings,x_mat[new_order_features,new_order_samples],new_order_features,new_order_samples)
    elseif Base.lowercase(layer_name) == "factor"
        if isnothing(gene_set_file_path)
            error("Please provide a gene set file path")
        end
        mask, gene_set_names = load_gmt_and_create_mask(gene_set_file_path, gene_names; min_genes_detected=min_genes_detected)
        new_order = sortperm(gene_set_names)
        mask = mask[:,new_order]
        gene_set_names = gene_set_names[new_order]
        detected_genes_bool = vec(sum(mask,dims=2) .!=0)
        X = deepcopy(x_mat[new_order_features,new_order_samples][detected_genes_bool,:]);
        mask = mask[detected_genes_bool,:]
        # gene_sets = load_gmt(gene_set_file_path)
        # mask, gene_set_names = create_binary_mask(gene_names, gene_sets);
        used_representation_feature_name = gene_set_names
        mask = mask';
        n_factors = length(gene_set_names);
        lam1,lam2,lr,n_iter,print_every =0,0.1,1e-1,150,1
        selected_representation, W_est = sparse_factor_model(X, mask, n_factors; lam1=lam1, lam2=lam2,lr=lr, n_iter=n_iter,print_every=print_every)
        selected_representation =  selected_representation';
        gene_factor_loadings = DataFrame(W_est,:auto);
        rename!(gene_factor_loadings,Symbol.(gene_names[detected_genes_bool]))
        insertcols!(gene_factor_loadings,1,:feature_name =>used_representation_feature_name)
        alternative_representation = (gene_factor_loadings,x_mat[new_order_features,new_order_samples],new_order_features,new_order_samples)
    elseif (Base.lowercase(layer_name) == Base.lowercase("scaledata")) || (Base.lowercase(layer_name) == Base.lowercase("raw")) || (Base.lowercase(layer_name) == Base.lowercase("logcounts"))
        selected_representation,selected_representation_feature_names,matching_projection_,anndata_dict1,layer_index,layer_name = select_data_representation(anndata_dict1;replace_x_keyvalue=false,is_precomputed_latent_representation=false,layer_index=layer_index,layer_name=layer_name,return_matching_projection_matrix=true)
        new_order_features=sortperm([gene_names_order_dict[gene] for gene in selected_representation_feature_names])
        selected_representation = selected_representation[new_order_features,new_order_samples]
        used_representation_feature_name = selected_representation_feature_names[sortperm([gene_names_order_dict[gene] for gene in selected_representation_feature_names])]
        alternative_representation = (nothing,nothing,new_order_features,new_order_samples)
    elseif  (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddings")) || (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddingsMeanCentered")) || (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddingsMedianCentered")) || (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddingsMeanThenMedianCentered")) || (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddingsMeanCenteredScaled")) || (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddingsMedianCenteredScaled")) || (Base.lowercase(layer_name) ==  Base.lowercase("ScanpyPcaEmbeddingsMedianCenteredScaled")) || (Base.lowercase(layer_name) ==  Base.lowercase("LdvaeLatentEmbeddingsZ10")) || (Base.lowercase(layer_name) ==  Base.lowercase("LdvaeLatentEmbeddingsZ50"))
        selected_representation,selected_representation_feature_names,matching_projection_,anndata_dict1,layer_index,layer_name = select_data_representation(anndata_dict1;replace_x_keyvalue=false,is_precomputed_latent_representation=is_precomputed_latent_representation,layer_index=layer_index,layer_name=layer_name,return_matching_projection_matrix=true);
        # print(selected_representation[:,1:5])
        selected_representation = selected_representation[:,new_order_samples]
        # selected_representation = selected_representation[:,new_order_samples]
        alternative_representation = (nothing,nothing,new_order_features,new_order_samples)
        used_representation_feature_name = selected_representation_feature_names
        if  typeof(matching_projection_) <: DataFrame 
            if eltype(matching_projection_[:,1]) <: String
                matching_projection_=matching_projection_[sortperm([gene_names_order_dict[el] for el in matching_projection_[:,1]]),:]
                matching_projection_genenames = matching_projection_[:,1]
                matching_projection_matrix = Matrix(matching_projection_[:,2:end])
                gene_factor_loadings = DataFrame(permutedims(matching_projection_matrix),:auto);
                rename!(gene_factor_loadings,Symbol.(matching_projection_genenames))
                insertcols!(gene_factor_loadings,1,:feature_name =>names(matching_projection_)[2:end])
                alternative_representation = (gene_factor_loadings,x_mat[new_order_features,new_order_samples],new_order_features,new_order_samples)
            end
        else 
            alternative_representation = (nothing,nothing,new_order_features,new_order_samples)
        end
    else
        selected_representation = deepcopy(x_mat[new_order_features,new_order_samples])
        alternative_representation = (nothing,nothing,new_order_features,new_order_samples)
        used_representation_feature_name = gene_names
        layer_name = "None"
    end
    if !isnothing(standardization_of_used_representation)
        if occursin("center",Base.lowercase(standardization_of_used_representation)) && occursin("scale",Base.lowercase(standardization_of_used_representation)) # && in(Base.lowercase(layer_name),["pca","factor"])
            selected_representation = center_and_scale_matrix_cols(selected_representation;center_cols = true,scale_cols = true);
        elseif occursin("center",Base.lowercase(standardization_of_used_representation)) && !occursin("scale",Base.lowercase(standardization_of_used_representation)) # && in(Base.lowercase(layer_name),["pca","factor"])
            selected_representation = center_and_scale_matrix_cols(selected_representation;center_cols = true,scale_cols = false);
        elseif !occursin("center",Base.lowercase(standardization_of_used_representation)) && occursin("scale",Base.lowercase(standardization_of_used_representation)) # && in(Base.lowercase(layer_name),["pca","factor"])
            selected_representation = center_and_scale_matrix_cols(selected_representation;center_cols = false,scale_cols = true);
        end
    end
    N = size(selected_representation)[2]
    G = size(selected_representation)[1]
    N_t = [[individual_time_combinations_counts[(individuals_ids[i],el)] for el in time_ids[i]] for i in 1:I]
    T = [length([individual_time_combinations_counts[(individuals_ids[i],el)] for el in time_ids[i]]) for i in 1:I]
    linear_sorted_sorting_keys = [(el[1],el[2],el[3],[0],el[4]) for el in sorted_sorting_keys]
    #[[individual_time_combinations_counts[ind] for ind in total_unique_combinations if ind[1] == individuals_ids[i]] for i in 1:I]
    z = nothing
    data_input = Vector{Vector{Vector{Vector{Float64}}}}(undef,I)
    # for i in 1:I
    #     data_input[i] = Vector{Vector{Vector{Float64}}}(undef,T[i])
    #     for t in 1:T[i]
    #         # N_t[i][t] = sum([(el[1] == i) && (el[2] == t) for el in sorted_sorting_keys])
    #         data_input[i][t] = [Float64.(collect(col)) for col in eachcol(selected_representation[:, (time_vec .== time_ids[i][t]) .&& (individuals_vec .== individuals_ids[i])])]
    #     end
    # end
    counter=0
    for i in 1:I
        data_input[i] = Vector{Vector{Vector{Float64}}}(undef,T[i])
        for t in 1:T[i]
            # N_t[i][t] = sum([(el[1] == i) && (el[2] == t) for el in sorted_sorting_keys])
            data_input[i][t] = Vector{Vector{Float64}}(undef,N_t[i][t])
            for n in 1:N_t[i][t]
                counter+=1
                data_input[i][t][n] = Vector{Float64}(undef,G)
                linear_sorted_sorting_keys[counter][4][1] = n
            end
        end
    end
    # curr_i = 1
    # data_input[curr_i] = Vector{Vector{Vector{Float64}}}(undef,T[curr_i])
    for (val,el) in enumerate(linear_sorted_sorting_keys)
        i = el[1]
        t = el[2]
        n = el[4][1]
        @test val == el[5]
        data_input[i][t][n] = Float64.(collect(selected_representation[:,val]))#[Float64.(collect(col)) for col in eachcol()]
    end
    # Check in N_t and data_input are consistent
    @test all([all([length(data_input[i][t]) for t in 1:T[i]] .== N_t[i]) for i in 1:I])
    # Check if timepoints are in order
    @test all(all.([[time_vec[(time_ids[i][t] .== time_vec) .&& (individuals_ids[i] .== individuals_vec) ][n-1]<=time_vec[(time_ids[i][t] .== time_vec) .&& (individuals_ids[i] .== individuals_vec) ][n] for t in 1:T[i] for n in 2:Int(N_t[i][t])] for i in 1:I]))  # Should be true
    # @test all([time_vec[i-1]<=time_vec[i] for i in 2:N])  # Should be true
    ####
    # N_t[1] = N
    # data_input[1] = [Float64.(collect(col)) for col in eachcol(x_mat)]
    # if !isnothing(cell_cluster_labels)
    #     z = Vector{Vector{Int}}(undef,1)
    #     z[1] = Int.(collect(cell_cluster_labels))
    # end
    return data_input,z,alternative_representation,used_representation_feature_name,layer_name
end

"""
    sort_projection_matrix_key_names(lst)
This function sorts projection matrix key names based on numerical order.
"""
function sort_projection_matrix_key_names(lst)
    sort(lst, by = s -> (
        occursin(r"\d", s), 
        occursin(r"\d", s) ? parse(Int, match(r"\d+", s).match) : -1
    ))
end

"""
    select_data_representation(anndata_dict1;replace_x_keyvalue=false,is_precomputed_latent_representation=false,layer_index=0,layer_name=nothing,return_matching_projection_matrix=true)
This function selects a data representation from an AnnData dictionary.
"""
function select_data_representation(anndata_dict1;replace_x_keyvalue=false,is_precomputed_latent_representation=false,layer_index=0,layer_name=nothing,return_matching_projection_matrix=true)
    if !is_precomputed_latent_representation
        representation_key = "layers"
    else
        representation_key = "obsm"
    end
    if haskey(anndata_dict1,representation_key)
        layer_names = collect(keys(anndata_dict1[representation_key]))
        if layer_index != 0
            layer_name = layer_names[layer_index]
            selected_representation = deepcopy(anndata_dict1[representation_key][layer_name])
        else
            if isnothing(layer_name)
                if representation_key == "layers"
                    layer_name = "None"
                    selected_representation = deepcopy(anndata_dict1["X"])
                else
                    error("Please provide a valid layer name for the representation: $representation_key")
                end
            else
                if !in(layer_name,layer_names)
                    error("Please provide a valid layer name for the representation: $representation_key")
                end
                selected_representation = deepcopy(anndata_dict1[representation_key][layer_name])
            end
        end
    end
    if representation_key == "obsm" && occursin("pca",lowercase(layer_name))
        num_features = size(selected_representation)[1]
        selected_representation_feature_names = ["PC_$(i)" for i in 1:num_features]
    elseif representation_key == "obsm" && occursin("ldvae",lowercase(layer_name))
        num_features = size(selected_representation)[1]
        selected_representation_feature_names = ["Z_$(i)" for i in 1:num_features]
    elseif representation_key == "layers"
        selected_representation_feature_names = anndata_dict1["var"]["_index"]
    else
        error("Please provide a valid layer name for the representation: $representation_key")
    end
    matching_projection_ = nothing
    if representation_key == "obsm" && haskey(anndata_dict1,"varm") && return_matching_projection_matrix
        if haskey(anndata_dict1["varm"],layer_name)
            if typeof(anndata_dict1["varm"][layer_name]) <: Dict
                cc = anndata_dict1["varm"][layer_name]
                matching_projection_ = DataFrame(OrderedDict(el => cc[el] for el in sort_projection_matrix_key_names(collect(keys(cc)))))
            else
                matching_projection_ = anndata_dict1["varm"][layer_name]
            end
        else
            matching_projection_ = nothing
        end
    else
        matching_projection_ = nothing
    end
    if replace_x_keyvalue
        anndata_dict1["X"] = selected_representation
        anndata_dict1["var"]["_index"] = selected_representation_feature_names
    end
    return selected_representation,selected_representation_feature_names,matching_projection_,anndata_dict1,layer_index,layer_name
end


"""
    select_cells_hvgs(x_mat,num_var_feat,num_cnts,gene_ids,cell_cluster_labels,scale_factor,N;chosen_cells=nothing)
This function selects highly variable genes from a given expression matrix.
"""
function select_cells_hvgs(x_mat,num_var_feat,num_cnts,gene_ids,cell_cluster_labels,scale_factor,N;chosen_cells=nothing)
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


"""
    verify_initialization_type(var_init,dimensions_tuple)
This function verifies the type of initialization for a given variable.
"""
function verify_initialization_type(var_init,dimensions_tuple)
    if typeof(var_init) <: Nothing
        var_init = fill(nothing, dimensions_tuple)#Vector{Nothing}(undef, prod(dimensions_tuple))
    elseif typeof(var_init) <: AbstractArray
        if length(var_init) != prod(dimensions_tuple)
            error("The length of the initialization vector does not match the dimensions of the variable")
        end
        var_init = reshape(var_init,dimensions_tuple)
    else
        error("The initialization type is not recognized")
    end
    return var_init
end


"""
    size_concentration_contractions(N;exp0=1.0)
This function computes the size concentration contraction based on the number of samples and an exponent.
"""
function size_concentration_contractions(N;exp0=1.0)
    return N^(-exp0)
end


"""
    initialize_model_parameters(data_input,KMax,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,seed;num_iter=500, size_concentration_contractions_exp0=1.0,rand_init = false,change_seeds = false, uniform_theta_init = true,g1_init=nothing,g2_init=nothing, m_mu_init=nothing, s_sq_mu_init=nothing,m_nu_init=nothing, s_sq_nu_init=nothing,y_init=nothing,u_init=nothing,v_init = nothing,a_init=nothing,b_init=nothing, h1_init=nothing, h2_init=nothing, w1_init=nothing, w2_init=nothing, d_init=nothing, c_init=nothing,r_init=nothing,condition_update_neighbors=nothing,condition_network_neighbors=nothing,update_clusterwise::Bool = false,samplebased_alpha0::Bool = false,eta_update_mode="Local",sigma_update_mode="Local", lambda_update_mode="Local",train_h::Bool = false,train_w::Bool = false,train_ab::Bool = false,train_uv::Bool = false)
This function initializes model parameters for a given dataset and configuration.
"""
function initialize_model_parameters(data_input,KMax,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,seed;num_iter=500, size_concentration_contractions_exp0=1.0,rand_init = false,change_seeds = false, uniform_theta_init = true,g1_init=nothing,g2_init=nothing, m_mu_init=nothing, s_sq_mu_init=nothing,m_nu_init=nothing, s_sq_nu_init=nothing,y_init=nothing,u_init=nothing,v_init = nothing,a_init=nothing,b_init=nothing, h1_init=nothing, h2_init=nothing, w1_init=nothing, w2_init=nothing, d_init=nothing, c_init=nothing,r_init=nothing,condition_update_neighbors=nothing,condition_network_neighbors=nothing,update_clusterwise::Bool = false,samplebased_alpha0::Bool = false,eta_update_mode="Local",sigma_update_mode="Local", lambda_update_mode="Local",train_h::Bool = false,train_w::Bool = false,train_ab::Bool = false,train_uv::Bool = false)
    # num_iter=500; size_concentration_contractions_exp0=1.0;rand_init = false;change_seeds = false; uniform_theta_init = true;g1_init=nothing;g2_init=nothing; m_mu_init=nothing; s_sq_mu_init=nothing;y_init=nothing;u_init=nothing;v_init = nothing;a_init=nothing;b_init=nothing; h1_init=nothing; h2_init=nothing; w1_init=nothing; w2_init=nothing; d_init=nothing; c_init=nothing;r_init=nothing; update_clusterwise=false; condition_update_neighbors=nothing;condition_network_neighbors=nothing; m_nu_init=nothing; s_sq_nu_init=nothing;
    
    if typeof(KMax) <: AbstractFloat
        KMax = Int(round(KMax))
    end
    Random.seed!(seed)
    K = KMax
    I = length(data_input)
    T = [length(data_input[i]) for i in 1:I]
    T_all = sum(T)
    J = length(data_input[1][1][1])
    float_type = eltype(data_input[1][1][1])
    N_t = [[length(data_input[i][t]) for t in  1:T[i]] for  i in 1:I]
    N = convert(typeof(I),sum([sum(el) for el in N_t]))
    linear_sample_index_as_ragged_array, linear_time_index_as_ragged_array = get_linear_index_as_ragged_array(N_t)
    if isnothing(condition_update_neighbors)
        condition_update_neighbors=get_linear_time_condition_update_neighbors(linear_time_index_as_ragged_array;get_ragged_array=false)
    end
    if isnothing(condition_network_neighbors)
        condition_network_neighbors=get_linear_time_condition_network_neighbors(linear_time_index_as_ragged_array;get_ragged_array=false)
    end
    #alpha0 = 1.0
    if !samplebased_alpha0 && typeof(alpha0) <: Number
        alpha = Vector{Vector{Float64}}(undef,I)
        for i in 1:I
            alpha[i] = Vector{Float64}(undef,T[i])
            for t in 1:T[i]
                alpha[i][t] = alpha0
            end
        end
        alpha0 = alpha
    elseif samplebased_alpha0
        alpha0 = [[size_concentration_contractions(N_t[i][t];exp0=size_concentration_contractions_exp0) for t in 1:T[i]] for i in 1:I]
    end
    m_mu_init = verify_initialization_type(m_mu_init,(J,K))
    s_sq_mu_init = verify_initialization_type(s_sq_mu_init,(J,K))
    m_nu_init = verify_initialization_type(m_nu_init,(J))
    s_sq_nu_init = verify_initialization_type(s_sq_nu_init,(J))
    y_init = verify_initialization_type(y_init,(J,K))
    g1_init = verify_initialization_type(g1_init,(K))
    g2_init = verify_initialization_type(g2_init,(K))
    # h1_init = verify_initialization_type(h1_init,(K))
    # h2_init = verify_initialization_type(h2_init,(K))
    w1_init = verify_initialization_type(w1_init,(T_all))
    w2_init = verify_initialization_type(w2_init,(T_all))
    d_init = verify_initialization_type(d_init,(K+1,T_all))
    c_init = verify_initialization_type(c_init,(T_all,N))
    r_init = verify_initialization_type(r_init,(K+1,N))
    if sigma_update_mode == "Clusterwise" && train_ab
        error("The sigma update for 'Clusterwise' is not implemented")
        # a_init = verify_initialization_type(a_init,(K))
        # b_init = verify_initialization_type(b_init,(K))
        # a_init = permutedims(hcat([a_init for j in 1:J]...))
        # b_init = permutedims(hcat([b_init for j in 1:J]...))
    elseif sigma_update_mode == "Genewise" && train_ab
        a_init = verify_initialization_type(a_init,(J))
        b_init = verify_initialization_type(b_init,(J))
        a_init = hcat([a_init for k in 1:K]...)
        b_init = hcat([b_init for k in 1:K]...)
    elseif sigma_update_mode == "Global" && train_ab
        error("The sigma update for 'Global' is not implemented")
        # a_init = verify_initialization_type(a_init,(1))
        # b_init = verify_initialization_type(b_init,(1))
        # if isnothing(a_init[1]) && rand_init
        #     a_init = [rand()]
        #     a_init =  a_init[1] .* hcat([ones(J) for k in 1:K]...)
        # elseif isnothing(a_init[1]) && !rand_init
        #     a_init = hcat([[a_init[1] for j in 1:J] for k in 1:K]...)
        # else
        #     a_init = a_init[1] .* hcat([ones(J) for k in 1:K]...)
        # end
        # if isnothing(b_init[1]) && rand_init
        #     b_init = [rand()]
        #     b_init = b_init[1] .* hcat([ones(J) for k in 1:K]...)
        # elseif isnothing(b_init[1]) && !rand_init
        #     b_init = hcat([[b_init[1] for j in 1:J] for k in 1:K]...)
        # else
        #     b_init = b_init[1] .* hcat([ones(J) for k in 1:K]...)
        # end
    elseif sigma_update_mode == "Local" && train_ab
        a_init = verify_initialization_type(a_init,(J,K))
        b_init = verify_initialization_type(b_init,(J,K))
    elseif !train_ab
        a_init = verify_initialization_type(a_init,(J,K))
        b_init = verify_initialization_type(b_init,(J,K))
    else
        error("The sigma update mode is not recognized")
    end
    # a_init = verify_initialization_type(a_init,(J))
    # b_init = verify_initialization_type(b_init,(J))
    if lambda_update_mode == "Clusterwise" && train_uv
        u_init = verify_initialization_type(u_init,(K))
        v_init = verify_initialization_type(v_init,(K))
        u_init = permutedims(hcat([u_init for j in 1:J]...))
        v_init = permutedims(hcat([v_init for j in 1:J]...))
    elseif lambda_update_mode == "Genewise" && train_uv
        u_init = verify_initialization_type(u_init,(J))
        v_init = verify_initialization_type(v_init,(J))
        u_init = hcat([u_init for k in 1:K]...)
        v_init = hcat([v_init for k in 1:K]...)
    elseif lambda_update_mode == "Global" && train_uv
        u_init = verify_initialization_type(u_init,(1))
        v_init = verify_initialization_type(v_init,(1))
        if isnothing(u_init[1]) && rand_init
            u_init = [rand()]
            u_init =  u_init[1] .* hcat([ones(J) for k in 1:K]...)
        elseif isnothing(u_init[1]) && !rand_init
            u_init = hcat([[u_init[1] for j in 1:J] for k in 1:K]...)
        else
            u_init = u_init[1] .* hcat([ones(J) for k in 1:K]...)
        end
        if isnothing(v_init[1]) && rand_init
            v_init = [rand()]
            v_init = v_init[1] .* hcat([ones(J) for k in 1:K]...)
        elseif isnothing(v_init[1]) && !rand_init
            v_init = hcat([[v_init[1] for j in 1:J] for k in 1:K]...)
        else
            v_init = v_init[1] .* hcat([ones(J) for k in 1:K]...)
        end
    elseif lambda_update_mode == "Local" && train_uv
        u_init = verify_initialization_type(u_init,(J,K))
        v_init = verify_initialization_type(v_init,(J,K))
    elseif !train_uv
        u_init = verify_initialization_type(u_init,(J,K))
        v_init = verify_initialization_type(v_init,(J,K))
    else
        error("The lambda update mode is not recognized")
    end
    if eta_update_mode == "Clusterwise" && train_h
        h1_init = verify_initialization_type(h1_init,(K))
        h2_init = verify_initialization_type(h2_init,(K))
        h1_init = permutedims(hcat([h1_init for j in 1:J]...))
        h2_init = permutedims(hcat([h2_init for j in 1:J]...))
    elseif eta_update_mode == "Genewise" && train_h
        h1_init = verify_initialization_type(h1_init,(J))
        h2_init = verify_initialization_type(h2_init,(J))
        h1_init = hcat([h1_init for k in 1:K]...)
        h2_init = hcat([h2_init for k in 1:K]...)
    elseif eta_update_mode == "Global" && train_h
        h1_init = verify_initialization_type(h1_init,(1))
        h2_init = verify_initialization_type(h2_init,(1))
        if isnothing(h1_init[1]) && rand_init
            h1_init = [rand()]
            h1_init =  h1_init[1] .* hcat([ones(J) for k in 1:K]...)
        elseif isnothing(h1_init[1]) && !rand_init
            h1_init = hcat([[h1_init[1] for j in 1:J] for k in 1:K]...)
        else
            h1_init = h1_init[1] .* hcat([ones(J) for k in 1:K]...)
        end
        if isnothing(h2_init[1]) && rand_init
            h2_init = [rand()]
            h2_init = h2_init[1] .* hcat([ones(J) for k in 1:K]...)
        elseif isnothing(h2_init[1]) && !rand_init
            h2_init = hcat([[h2_init[1] for j in 1:J] for k in 1:K]...)
        else
            h2_init = h2_init[1] .* hcat([ones(J) for k in 1:K]...)
        end
    elseif eta_update_mode == "Local" && train_h
        h1_init = verify_initialization_type(h1_init,(J,K))
        h2_init = verify_initialization_type(h2_init,(J,K))
    elseif !train_h
        h1_init = verify_initialization_type(h1_init,(J,K))
        h2_init = verify_initialization_type(h2_init,(J,K))
    else
        error("The eta update mode is not recognized")
    end
    cells = [CellFeature(i,t,n,KMax,T,data_input[i][t][n],rand_init =rand_init,c_init = c_init[:,linear_sample_index_as_ragged_array[i][t][n]],r_init = r_init[:,linear_sample_index_as_ragged_array[i][t][n]]) for i in 1:I for t in 1:T[i] for n in 1:N_t[i][t]];
    clusters = [ClusterFeature(k,J;float_type=float_type, m_mu_init = m_mu_init[:,k],s_sq_mu_init = s_sq_mu_init[:,k],m_nu_init = m_nu_init, s_sq_nu_init = s_sq_nu_init,y_init = y_init[:,k],g1_init = g1_init[k],g2_init = g2_init[k],h1_init = h1_init[:,k],h2_init = h2_init[:,k],a_init = a_init[:,k],b_init = b_init[:,k],u_init = u_init[:,k],v_init = v_init[:,k],rand_init =rand_init,update_clusterwise = update_clusterwise) for k in 1:KMax];
    dataparams = DataFeature(data_input);
    conditions = [ConditionFeature(i,t,KMax,T[i],condition_update_neighbors[i][t],condition_network_neighbors[i][t];float_type=float_type,d_init = d_init[:,linear_time_index_as_ragged_array[i][t]],w1_init = w1_init[linear_time_index_as_ragged_array[i][t]],w2_init = w2_init[linear_time_index_as_ragged_array[i][t]] ,rand_init =rand_init) for i in 1:I for t in 1:T[i]];
    matrixconditions = [MatrixConditionFeature(i,t,KMax,T[i],condition_update_neighbors[i][t],condition_network_neighbors[i][t];float_type=float_type) for i in 1:I for t in 1:T[i]];
    modelparams = ModelParameterFeature(data_input,K,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,num_iter,uniform_theta_init,rand_init,change_seeds,seed);
    training_logger = TrainFeature(1,T_all,K,J,Int64(num_iter+1));
    input_str_list = @name cells,clusters,conditions,matrixconditions,dataparams,modelparams,training_logger;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [cells,clusters,conditions,matrixconditions,dataparams,modelparams,training_logger];
    inputs = OrderedDict()
    addToDict!(inputs,input_key_list,input_var_list);
    return inputs
end

"""
    get_m_init_from_initialization_approach(m_initialization_approach,used_representation,cell_cluster_labels,K,seed;iseeds=nothing,m_mu_init=nothing)
This function initializes the cluster means based on a specified initialization approach.
"""
function get_m_init_from_initialization_approach(m_initialization_approach,used_representation,cell_cluster_labels,K,seed;iseeds=nothing,m_mu_init=nothing)
    J = size(used_representation)[1]
    N = size(used_representation)[2]
    Random.seed!(seed)
    if m_initialization_approach == "rand"
        m_mu_init = nothing
    elseif m_initialization_approach == "fcmeans"
        R = fuzzy_cmeans(used_representation, 100, 2, maxiter=200, display=:iter)
        m_mu_init = R.centers
    elseif m_initialization_approach == "kmeans"
        R = Clustering.kmeans(used_representation,K)
        m_mu_init = R.centers
    elseif m_initialization_approach == "kpp"
        iseeds = initseeds(:kmpp,used_representation,K)
        m_mu_init = used_representation[:,iseeds]
    elseif m_initialization_approach == "kmcen"
        iseeds = initseeds(:kmcen,used_representation,K)
        m_mu_init = used_representation[:,iseeds]
    elseif m_initialization_approach == "kpp+rand"
        iseeds = initseeds(:kmpp,used_representation,K)
        m_mu_init = used_representation[:,iseeds] .+ std(used_representation,dims=2)/sqrt(N) .* randn(size(used_representation[:,iseeds]))
    elseif m_initialization_approach == "cell_label"
        if !isnothing(cell_cluster_labels)
            m_mu_init = randn(J,K)
            unique_labels = unique(cell_cluster_labels)
            for k in 1:length(unique_labels)
                m_mu_init[:,k] .= mean(used_representation[:,cell_cluster_labels .== unique_labels[k]],dims=2)
            end
        else
            error("Cell Labels are not provided")
        end
    elseif m_initialization_approach == "cell_label+rand"
        if !isnothing(cell_cluster_labels)
            m_mu_init = randn(J,K)
            unique_labels = unique(cell_cluster_labels)
            for k in 1:length(unique_labels)
                nk = sum(cell_cluster_labels .== unique_labels[k])
                m_mu_init[:,k] .= mean(used_representation[:,cell_cluster_labels .== unique_labels[k]],dims=2)
                m_mu_init[:,k] .+= reshape(std(used_representation[:,cell_cluster_labels .== unique_labels[k]],dims=2),J)/sqrt(nk) .* randn(J)
            end
            # m_mu_init += std(anndata_dict1["X"])*randn(size(m_mu_init))
        else
            error("Cell Labels are not provided")
        end
    end
    return m_mu_init,iseeds
end

"""
    run_cavi(inputs;elbo_ep = 10^(-0),logger=nothing)
This function runs the Coordinate Ascent Variational Inference (CAVI) algorithm on the provided inputs.
"""
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

"""
    save_pips(pip,gene_names;unique_time_id="",filepath="")
This function saves the posterior inclusion probabilities (PIPs) to a CSV file.
"""
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

"""
    save_Nk(rtik_;unique_time_id="",filepath="")
This function saves the cluster sizes (Nk) to a CSV file.
"""
function save_Nk(rtik_;unique_time_id="",filepath="")
    KMax = length(rtik_[1][1])
    N = length(rtik_[1])
    T = 1
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


"""
    summarize_parameters(outputs_dict_vec,elapsed_time,final_elbo_vec,elbo_vec,rtik_vec,yjk_vec,perK_elbo_vec,delta_t_vec,nk_perL_vec,ηk_vec)
This function summarizes model parameters across multiple runs.
"""
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

"""
    _flushed_logger(msg;logger=nothing)
This function logs a message using the provided logger.
"""
function _flushed_logger(msg;logger=nothing)
    if !isnothing(logger)
        with_logger(logger) do
            @info msg
        end
    end
end


"""
    make_ids(dataset_name,G,N)
This function generates unique identifiers for the dataset, experiment, and time.
"""
function make_ids(dataset_name,G,N)
    unique_time_id = get_unique_time_id()
    dataset_used_id = "$(dataset_name)_$(G)HVGs-$(N)N"
    experiment_filename = "nclusion_$(dataset_name)"
    experiment_id = setup_experiment_tag(experiment_filename)

    return unique_time_id,dataset_used_id,experiment_id
end


"""
    mk_outputs_filepath(outdir,experiment_id,dataset_used,unique_time_id)
This function generates the output file path for the experiment.
"""
function mk_outputs_filepath(outdir,experiment_id,dataset_used,unique_time_id)
    filepath ="$outdir/outputs/$experiment_id/current/$(unique_time_id)_$(dataset_used)/"
    return filepath
end

"""
    mk_outputs_pathname(filepath)
This function creates the output directory if it does not exist.
"""
function mk_outputs_pathname(filepath)
    mkpath(filepath)
end

"""
    saving_summary_file(filepath;slurm_job_id="",script_name="",change_seeds = false,unique_time_id="",datafilename1="",KMax="", seed="",num_var_feat="",N="",elbo_ep="",notes_="",alpha0 = "",gamma0 = "" ,phi1 = "",phi2 = "",kappa1 = "",kappa2 = "",xi1 = "",xi2 = "",varphi1 = "",varphi2 = "",significance_prop="",min_number_cells="",min_percent_cells="",min_percent_of_genes="",max_percent_of_genes="", num_iter = "",dataset_name = "",outdir = "",time_key = "",individuals_key = "",m_initialization_approach="",update_clusterwise="", use_alt_representation="", check_cluster_interpretability_bool="", gene_set_file_path="", samplebased_alpha0="")
This function saves a summary of the experiment parameters to a text file.
"""
function saving_summary_file(filepath;slurm_job_id="",script_name="",change_seeds = false,unique_time_id="",datafilename1="",KMax="", seed="",num_var_feat="",N="",elbo_ep="",notes_="",alpha0 = "",gamma0 = "" ,phi1 = "",phi2 = "",kappa1 = "",kappa2 = "",xi1 = "",xi2 = "",varphi1 = "",varphi2 = "",significance_prop="",min_number_cells="",min_percent_cells="",min_percent_of_genes="",max_percent_of_genes="", num_iter = "",dataset_name = "",outdir = "",time_key = "",individuals_key = "",m_initialization_approach="",update_clusterwise="", use_alt_representation="", check_cluster_interpretability_bool="", gene_set_file_path="", samplebased_alpha0="")
    summary_file = filepath*"_QuickSummary_"*unique_time_id*".txt"
    vars = [datafilename1,dataset_name,slurm_job_id,script_name,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,significance_prop,min_percent_cells,min_number_cells,min_percent_of_genes,max_percent_of_genes,KMax,change_seeds, seed, m_initialization_approach, "Scanpy-Default", time_key,individuals_key,num_iter,num_var_feat,N,true,false,false,false,elbo_ep,outdir,update_clusterwise, use_alt_representation, check_cluster_interpretability_bool, gene_set_file_path, samplebased_alpha0,notes_]
    varnames = ["Filename","Dataset ID","Slurm Job ID","Script used to generate these results","alpha0","gamma0","phi1","phi2","kappa1","kappa2","xi1","xi2","varphi1","varphi2","significance_prop","min_percent_cells","min_number_cells","min_percent_of_genes","max_percent_of_genes","KMax","Seeds Changed During Training","Initial Seed used","Approach used to initalialize cluster means","Which Method did I use to select HVGs","Condition key in data","Individuals key in data","max number of iterations", "How many Genes were in this anaysis", "How many Cells were in this anaysis","(True/False) I used all cells in this anaysis",  "(True/False) I manually had to standardize the data and did not use the steps in the QC pipeline","(True/False) I started with raw counts","(True/False) I used all genes","Elbo intolerance threshold","Output directorry","Cluster-wise update of lambda","Alternative Representation used","Cluster interpretability checked","Gene set path","Sample - based alpha0 was used","Notes"]

     
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

"""
    save_embeddings(anndata_dict1,filepath;outputs_dict = nothing,logger = nothing,unique_time_id="",new_order_samples=nothing)
This function saves the embeddings (TSNE, PCA, UMAP) to CSV files.
"""
function save_embeddings(anndata_dict1,filepath;outputs_dict = nothing,logger = nothing,unique_time_id="",new_order_samples=nothing)
    G = size(anndata_dict1["X"])[1]
    N = size(anndata_dict1["X"])[2]
    tsne_data = nothing
    pca_data = nothing
    umap_data = nothing
    if isnothing(new_order_samples)
        new_order_samples = Int.(collect(1:N))
    end
    # @info "Getting TSNE Transform..."
    _flushed_logger("Getting TSNE Transform...";logger)
    
    if haskey(anndata_dict1["obsm"],"X_tsne")
        tsne_data =  permutedims(anndata_dict1["obsm"]["X_tsne"])
        tsne_data = tsne_data[new_order_samples,:];
        tsne_data_df  = DataFrame(tsne_data, :auto);
        ncols = size(tsne_data)[2]
        rename!(tsne_data_df,Symbol.(["TSNE_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-tsne_coordinates.csv",  tsne_data_df)
        if !isnothing(outputs_dict)
            outputs_dict[:X_tsne] = tsne_data
        end
    end

    
    _flushed_logger("Getting PCA Transform...";logger)
    if haskey(anndata_dict1["obsm"],"X_pca")
        pca_data =  permutedims(anndata_dict1["obsm"]["X_pca"])
        pca_data = pca_data[new_order_samples,:];
        pca_data_df  = DataFrame(pca_data, :auto);
        ncols = size(pca_data)[2]
        rename!(pca_data_df,Symbol.(["PC_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-pca_coordinates.csv",  pca_data_df)
        if !isnothing(outputs_dict)
            outputs_dict[:X_pca] = pca_data
        end
    end


    _flushed_logger("Getting UMAP Transform..."; logger)
    
    if haskey(anndata_dict1["obsm"],"X_umap")
        umap_data =  permutedims(anndata_dict1["obsm"]["X_umap"])
        umap_data = umap_data[new_order_samples,:];
        umap_data_df  = DataFrame(umap_data, :auto);
        ncols = size(umap_data)[2]
        rename!(umap_data_df,Symbol.(["UMAP_$i" for i in 1:ncols]));
        CSV.write(filepath*"$(G)G-"*unique_time_id*"-umap_coordinates.csv",  umap_data_df)
        if !isnothing(outputs_dict)
            outputs_dict[:X_umap] = umap_data
        end
    end
    return outputs_dict
end

"""
    make_labels(x_input,anndata_dict1,z_argmax)
This function generates a DataFrame containing cell labels and inferred cluster assignments.
"""
function make_labels(x_input,anndata_dict1,z_argmax)
    T = 1
    timepoint_map = [t * ones(Int,Int(length(x_input[t]))) for t in 1:T]
    timepoint_map = recursive_flatten(timepoint_map)
    # cluster_remap = Dict(v => k for (k,v) in cluster_map)
    cluster_results_df = nothing
    if haskey(anndata_dict1["obs"],"cell_type")
        cell_ids = anndata_dict1["obs"]["_index"]
        cell_cluster_labels = anndata_dict1["obs"]["cell_type"]["codes"]
        cell_cluster_labels = cell_cluster_labels .+ 1
        z = Vector{Vector{Int}}(undef,T)
        z[1] = Int.(collect(cell_cluster_labels))
        cluster_map = OrderedDict(k => v for (k,v) in enumerate(anndata_dict1["obs"]["cell_type"]["categories"]))
        cell_type_vec = [cluster_map[el] for el in recursive_flatten(z)]

        cluster_results_df = DataFrame(condtion=timepoint_map,cell_id = cell_ids,cell_type = cell_type_vec, called_label = recursive_flatten(z), inferred_label = recursive_flatten(z_argmax))
    else
        cell_type_vec = [ "NA" for el in recursive_flatten(z_argmax)]
        cluster_results_df = DataFrame(condtion=timepoint_map,cell_id = cell_ids,cell_type = cell_type_vec, called_label = cell_type_vec, inferred_label = recursive_flatten(z_argmax))
    end
    
    return cluster_results_df
    
end

"""
    save_labels(cluster_results_df;dataset_used="",G="",unique_time_id="",filepath="")
This function saves the cluster membership labels to a CSV file.
"""
function save_labels(cluster_results_df;dataset_used="",G="",unique_time_id="",filepath="")
    CSV.write(filepath*"$(dataset_used)_nclusion-"*unique_time_id*".csv",  cluster_results_df)
end

function load_constants()
    gene_set_file_path = ""
    significance_prop =  0.5 
    min_percent_of_genes = 	 0.05 
    max_percent_of_genes = 	 0.5 # 0.1
    time_key = nothing
    check_cluster_interpretability_bool = false
    individuals_key = nothing
    m_initialization_approach = "rand";
    slurm_job_id = "No SLURM Job ID Provided"
    change_seeds = true
    samplebased_alpha0 = false
    add_time = false
    add_individuals = false
    layer_index = 0
    layer_name = "xmat"
    is_precomputed_latent_representation =  occursin(lowercase("Embeddings"), lowercase(layer_name)) ?  true  : false
    update_clusterwise = false
    iseeds_init = nothing;
    m_mu_init_init = nothing;
    r_init = nothing;
    h1_init = nothing;
    h2_init = nothing;
    w1_init= nothing;
    w2_init= nothing;
    a_init = nothing;
    b_init = nothing;
    u_init = nothing;
    v_init = nothing;
    y_init = nothing;
    define_global_sparsity = true;
    define_global_autocorrelation = true;
    sparsity_based_on_used_representation_dim = false;
    train_w = true;
    train_ab = true;
    train_uv = true;
    train_h = true;
    train_m_nu = true;
    train_s_sq_nu = true;
    multiple_s_sq_updates = true;
    uv_scalar = 100.0;
    use_std_for_hvg = true;
    center_data_cols=true;
    scale_data_cols=true;
    min_number_cells = 	 1.0;
    min_percent_cells = 1/min_number_cells;#
    min_genes_detected_for_mask=10
    standardization_of_used_representation=nothing
    rand_init_inputs = true;
    uniform_theta_init = false;
    define_global_sparsity = true;
    eta_update_mode = "Global";
    sigma_update_mode = "Genewise";
    lambda_update_mode = "Local";
    center_data_cols=false;
    scale_data_cols=false;
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
    r_init_string = "r_init initialized to nothing";
    ab_init_string = "a_init and b_init initialized to nothing";
    uv_init_string = "u_init and v_init initialized to nothing";

    return gene_set_file_path,significance_prop,min_percent_of_genes,max_percent_of_genes,time_key,check_cluster_interpretability_bool,individuals_key,m_initialization_approach,slurm_job_id,change_seeds,samplebased_alpha0,add_time,add_individuals,layer_index,layer_name,is_precomputed_latent_representation,update_clusterwise,iseeds_init,m_mu_init_init,r_init,h1_init,h2_init,w1_init,w2_init,a_init,b_init,u_init,v_init,y_init,define_global_sparsity,define_global_autocorrelation,sparsity_based_on_used_representation_dim,train_w,train_ab,train_uv,train_h,train_m_nu,train_s_sq_nu,multiple_s_sq_updates,uv_scalar,use_std_for_hvg,center_data_cols,scale_data_cols,min_number_cells,min_percent_cells,min_genes_detected_for_mask,standardization_of_used_representation,rand_init_inputs,uniform_theta_init,define_global_sparsity,eta_update_mode,sigma_update_mode,lambda_update_mode,center_data_cols,scale_data_cols,data_structure_note, r_init_string, ab_init_string, uv_init_string
end

"""
    run_nclusion(datafilename1; KMax=25, seed=12345, dataset_name="foo", elbo_ep = 1.0e-6, outdir = "./", logger = nothing, num_iter=500, rand_init=false, alpha0 = 1.0, gamma0 = 1.0, phi1 = 1.0, phi2 = 1.0, kappa1 = 1.0, kappa2 = 1.0, xi1 = 1.0, xi2 = 1.0, varphi1 = 1.0, varphi2 = 1.0, nu0 = 0.0, sigma_sq_nu = 1e-12, delta_rsum_ep = 1.0, delta_rsum_lag = 2, delta_occymean_ep = 1.0, delta_occymean_lag = 2, elbo_sign_decrease_max_tolerance = 10, num_samples_posterior_samples = 1000, s_value_thresh=0.05, burnin = 10,
    )
This function runs the entire NCLUSION pipeline, from data loading to model training and result saving.
"""
function run_nclusion(datafilename1; 
    KMax=25,
    seed=12345,
    dataset_name="foo",
    elbo_ep = 1.0e-6,
    outdir = "./",
    logger = nothing,
    num_iter=500,
    rand_init=false,
    alpha0 = 1.0,
    gamma0 = 1.0,
    phi1 = 1.0,
    phi2 = 1.0,
    kappa1 = 1.0,
    kappa2 = 1.0,
    xi1 = 1.0,
    xi2 = 1.0,
    varphi1 = 1.0,
    varphi2 = 1.0,
    nu0 = 0.0,
    sigma_sq_nu = 1e-12,
    delta_rsum_ep = 1.0,
    delta_rsum_lag = 2,
    delta_occymean_ep = 1.0,
    delta_occymean_lag = 2,
    elbo_sign_decrease_max_tolerance = 10,
    num_samples_posterior_samples = 1000,
    s_value_thresh=0.05,
    burnin = 10,
    )
    Random.seed!(seed)
    gene_set_file_path,significance_prop,min_percent_of_genes,max_percent_of_genes,time_key,check_cluster_interpretability_bool,individuals_key,m_initialization_approach,slurm_job_id,change_seeds,samplebased_alpha0,add_time,add_individuals,layer_index,layer_name,is_precomputed_latent_representation,update_clusterwise,iseeds_init,m_mu_init_init,r_init,h1_init,h2_init,w1_init,w2_init,a_init,b_init,u_init,v_init,y_init,define_global_sparsity,define_global_autocorrelation,sparsity_based_on_used_representation_dim,train_w,train_ab,train_uv,train_h,train_m_nu,train_s_sq_nu,multiple_s_sq_updates,uv_scalar,use_std_for_hvg,center_data_cols,scale_data_cols,min_number_cells,min_percent_cells,min_genes_detected_for_mask,standardization_of_used_representation,rand_init_inputs,uniform_theta_init,define_global_sparsity,eta_update_mode,sigma_update_mode,lambda_update_mode,center_data_cols,scale_data_cols,data_structure_note, r_init_string, ab_init_string, uv_init_string = load_constants()

    parameter_settings_to_reproduce_run = initialize_ordered_dict(;key_type=String,val_type=String)
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,gene_set_file_path,significance_prop,min_percent_of_genes,max_percent_of_genes,time_key,check_cluster_interpretability_bool,individuals_key,m_initialization_approach,slurm_job_id,change_seeds,samplebased_alpha0,add_time,add_individuals,layer_index,layer_name,is_precomputed_latent_representation,update_clusterwise,iseeds_init,m_mu_init_init,r_init,h1_init,h2_init,w1_init,w2_init,a_init,b_init,u_init,v_init,y_init,define_global_sparsity,define_global_autocorrelation,sparsity_based_on_used_representation_dim,train_w,train_ab,train_uv,train_h,train_m_nu,train_s_sq_nu,multiple_s_sq_updates,uv_scalar,use_std_for_hvg,center_data_cols,scale_data_cols,min_number_cells,min_percent_cells,min_genes_detected_for_mask,standardization_of_used_representation,rand_init_inputs,uniform_theta_init,define_global_sparsity,eta_update_mode,sigma_update_mode,lambda_update_mode,center_data_cols,scale_data_cols,data_structure_note, r_init_string, ab_init_string, uv_init_string)
    K = KMax
    elbo_sign_change_max = elbo_sign_decrease_max_tolerance;
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,elbo_sign_change_max,num_iter,seed,KMax,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,delta_rsum_ep,delta_rsum_lag,delta_occymean_ep,delta_occymean_lag,num_samples_posterior_samples,s_value_thresh,burnin)
    _flushed_logger("Loading data and metadata...";logger)
    anndata_dict1= load_data(datafilename1,seed);
    N = size(anndata_dict1["X"])[2];
    J = size(anndata_dict1["X"])[1];
    _flushed_logger("Number of cells: $(N) ...";logger)
    _flushed_logger("Number of Genes: $(J) ...";logger)
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
    _flushed_logger("Preparing Dataset ...";logger)
    data_input,z,alternative_representation,used_representation_feature_name,layer_name = make_nclusion_inputs(anndata_dict1;time_key=time_key,individuals_key=individuals_key, layer_name=layer_name,gene_set_file_path=gene_set_file_path, min_genes_detected=min_genes_detected_for_mask,standardization_of_used_representation=standardization_of_used_representation,is_precomputed_latent_representation=is_precomputed_latent_representation);
    I_ = length(data_input);
    T = [ length(data_input[i]) for i in 1:I_]
    N_t = [[length(data_input[i][t]) for t in 1:T[i]] for i in 1:I_]
    T_all = sum(T);
    used_representation = hcat([data_input[i][t][n] for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]]...);
    used_representation = center_and_scale_matrix_cols(used_representation;center_cols = center_data_cols,scale_cols = scale_data_cols);
    original_mean = zeros(size(used_representation)[1])#vec(median(used_representation,dims=2))
    m_mu_init,iseeds = get_m_init_from_initialization_approach(m_initialization_approach,used_representation,cell_cluster_labels,K,seed;iseeds=iseeds_init,m_mu_init=m_mu_init_init);
    if !isnothing(m_mu_init)
        if size(m_mu_init)[2] != K
            m_mu_init = hcat(m_mu_init,randn(size(m_mu_init)[1],K-size(m_mu_init)[2]));
        end
    end
    nu0 = original_mean;
    inputs = initialize_model_parameters(data_input,K,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,seed;num_iter=num_iter,rand_init = rand_init_inputs,change_seeds = change_seeds,uniform_theta_init = uniform_theta_init, m_mu_init = m_mu_init,samplebased_alpha0=samplebased_alpha0,update_clusterwise=update_clusterwise,y_init=y_init,r_init = r_init,h1_init=h1_init,h2_init=h2_init,w1_init=w1_init,w2_init=w2_init,a_init=a_init,b_init=b_init,u_init=u_init,v_init=v_init,eta_update_mode=eta_update_mode,sigma_update_mode=sigma_update_mode,lambda_update_mode=lambda_update_mode,train_h=train_h,train_w=train_w,train_ab=train_ab,train_uv=train_uv,m_nu_init=[el for el in nu0]);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,iseeds,r_init_string,h1_init,h2_init,w1_init,w2_init,ab_init_string,uv_init_string,define_global_sparsity,define_global_autocorrelation);
    filepath = mk_outputs_filepath(outdir,experiment_id,dataset_used_id,unique_time_id);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,filepath);
    mk_outputs_pathname(filepath);
    notes_="$(data_structure_note)__*||*__$(layer_name)"
    summary_file = saving_summary_file(filepath;slurm_job_id=slurm_job_id,script_name="",change_seeds = change_seeds,unique_time_id=unique_time_id,datafilename1=datafilename1,KMax=KMax, seed=seed,num_var_feat=num_var_feat,N=N,elbo_ep=elbo_ep,notes_=notes_,alpha0 = alpha0,gamma0 =gamma0 ,phi1 = phi1,phi2 = phi2,kappa1 = kappa1,kappa2 = kappa2,xi1 =xi1,xi2 = xi2,varphi1 = varphi1,varphi2 = varphi2,significance_prop=significance_prop,min_percent_cells=min_percent_cells,min_number_cells=min_number_cells,min_percent_of_genes=min_percent_of_genes,max_percent_of_genes=max_percent_of_genes,num_iter = num_iter,dataset_name = dataset_name,outdir = outdir,time_key = time_key,individuals_key = individuals_key,m_initialization_approach = m_initialization_approach,update_clusterwise=update_clusterwise, use_alt_representation=layer_name, check_cluster_interpretability_bool=check_cluster_interpretability_bool, gene_set_file_path=gene_set_file_path, samplebased_alpha0=samplebased_alpha0,);
    @add_variables_to_ordered_dict_as_string!(parameter_settings_to_reproduce_run,inputs[:dataparams].I,inputs[:dataparams].T,inputs[:dataparams].J,inputs[:dataparams].N_t,inputs[:dataparams].N,inputs[:dataparams].Jlog,inputs[:dataparams].logpi,inputs[:modelparams].K,inputs[:modelparams].alpha0,inputs[:modelparams].gamma0,inputs[:modelparams].phi1,inputs[:modelparams].phi2,inputs[:modelparams].kappa1,inputs[:modelparams].kappa2,inputs[:modelparams].xi1,inputs[:modelparams].xi2,inputs[:modelparams].varphi1,inputs[:modelparams].varphi2,inputs[:modelparams].nu0,inputs[:modelparams].sigma_sq_nu,inputs[:modelparams].significance_prop,inputs[:modelparams].min_number_cells,inputs[:modelparams].min_percent_cells,inputs[:modelparams].min_percent_of_genes,inputs[:modelparams].max_percent_of_genes,inputs[:modelparams].num_iter,inputs[:modelparams].uniform_theta_init,inputs[:modelparams].rand_init,inputs[:modelparams].change_seeds,inputs[:modelparams].init_seed);
    training_features_dict = nothing;
    elapsed_time = @elapsed begin
            st = time();
            outputs_dict,training_features_dict = cavi(inputs;delta_rsum_ep =delta_rsum_ep,delta_rsum_lag = delta_rsum_lag,delta_occymean_ep = delta_occymean_ep,delta_occymean_lag = delta_occymean_lag,logger=logger,update_clusterwise=update_clusterwise,elbo_sign_change_max=elbo_sign_decrease_max_tolerance,check_cluster_interpretability_bool=check_cluster_interpretability_bool,train_h =train_h,train_w=train_w,train_ab=train_ab,train_uv=train_uv,multiple_s_sq_updates=multiple_s_sq_updates,eta_update_mode=eta_update_mode,sigma_update_mode=sigma_update_mode,lambda_update_mode=lambda_update_mode,burnin=burnin,remove_small_clusters=false,use_log=true,train_m_nu=train_m_nu,train_s_sq_nu=train_s_sq_nu);
    end
    dt = time() - st;
    _flushed_logger("\t \t Finished Training Model. Model took $dt seconds to run...";logger)
    _flushed_logger("\t \t Final ELBO $(outputs_dict[:LB][end])...";logger)
    _flushed_logger("Model Took a total of $elapsed_time seconds to run";logger)
    outputs_dict[:elapsed_time]=elapsed_time;
    inferred_labels = [argmax(outputs_dict[:r_][i][t][n]) for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]]
    num_clust = length(unique(inferred_labels))
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
   _flushed_logger("Calculating Posterior Summaries...";logger)
    temp_outputs_dict,temp_z_argmax,temp_posterior_gene_summaries_df,temp_cluster_summaries_df,temp_cell_summaries_df = posterior_summaries(outputs_dict, used_representation,used_representation_feature_name,seed;num_samples = num_samples_posterior_samples,s_value_thresh=s_value_thresh,return_data_frames = false,save_data_frames = true)
    if !isnothing(temp_outputs_dict)
        outputs_dict = temp_outputs_dict
    end
    if !isnothing(temp_z_argmax)
        z_argmax = temp_z_argmax
        outputs_dict[:z_argmax] = [z_argmax];
    end
    filename = "$filepath/output.jld2"
    _flushed_logger("Saving Outputs...";logger)
    jldsave(filename,true;outputs_dict=outputs_dict)
    trainingfilename = "$filepath/train_history.jld2"
    if !isempty(training_features_dict)
        _flushed_logger("Saving Training Features...";logger)
        jldsave(trainingfilename,true;outputs_dict=training_features_dict)
    end
    runParameterSettingsfilename = "$filepath/_RunParameterSettings_$unique_time_id.csv"
    runParameterSettingsdf = DataFrame(parameter_name=parameter_settings_to_reproduce_run.keys,parameter_value=[parameter_settings_to_reproduce_run[key] for key in parameter_settings_to_reproduce_run.keys]);
    CSV.write(runParameterSettingsfilename,runParameterSettingsdf);
    return outputs_dict
end

function _old_run_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = nothing,num_iter=500,save_metrics=false,rand_init=false)
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
    new_order = sortperm(gene_names)
    gene_names = gene_names[new_order]
    x_input,z = make_nclusion_inputs(anndata_dict1)
    

    filepath = mk_outputs_filepath(outdir,experiment_id,dataset_used_id,unique_time_id)
    _flushed_logger("Preparing saving Directory at $filepath ...";logger)
    mk_outputs_pathname(filepath)


    _flushed_logger("Saving Quick Summary...";logger)
    notes_=""
    summary_file = saving_summary_file(filepath;unique_time_id=unique_time_id,datafilename1=datafilename1,alpha1=alpha1,gamma1=gamma1,KMax=KMax, seed=seed,num_var_feat=num_var_feat,N=N,elbo_ep=elbo_ep,notes_=notes_)
    save_embeddings(anndata_dict1,filepath;logger = logger,unique_time_id=unique_time_id)

    

    _flushed_logger("Initializing Model parameters...";logger)
    inputs = initialize_model_parameters(x_input,KMax,alpha1,gamma1;num_iter=num_iter,rand_init=rand_init);

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
    
    if !isnothing(z)
        if save_metrics
            _flushed_logger("Calculating Extrinsic Metrics";logger)
            clustering_quality_metrics(cluster_results_df,filepath);
        end
    end

    return outputs_dict
end

"""
    append_summary(summary_file,vars,varnames)
This function appends additional information to the summary file.
"""
function append_summary(summary_file,vars,varnames)
    open(summary_file, "a") do f
        for i in eachindex(vars)
            write(f,"$(varnames[i]) = \t $(vars[i]) \n")
        end
    end
end

"""
    create_results_dict(run_,function_name)
This function creates a results dictionary from the benchmarking run.
"""
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

"""
    benchmark_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = nothing,num_iter=500)
This function benchmarks the NCLUSION model using the provided parameters.
"""
function benchmark_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = nothing,num_iter=500)
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
    new_order = sortperm(gene_names)
    gene_names = gene_names[new_order]
    x_input,_ = make_nclusion_inputs(anndata_dict1)
    

    filepath = mk_outputs_filepath(outdir,experiment_id,dataset_used_id,unique_time_id)
    _flushed_logger("Preparing saving Directory at $filepath ...";logger)
    mk_outputs_pathname(filepath)

    _flushed_logger("Saving Quick Summary...";logger)
    notes_=""
    summary_file = saving_summary_file(filepath;unique_time_id=unique_time_id,datafilename1=datafilename1,alpha1=alpha1,gamma1=gamma1,KMax=KMax, seed=seed,num_var_feat=num_var_feat,N=N,elbo_ep=elbo_ep,notes_=notes_)
    
    _flushed_logger("Initializing Model parameters...";logger)
    inputs = initialize_model_parameters(x_input,KMax,alpha1,gamma1;num_iter=num_iter);

    
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

"""
    get_highly_variable_genes_bool(anndata_dict1;n_hvgs=nothing, use_std=false)
This function determines the highly variable genes in the dataset.
"""
function get_highly_variable_genes_bool(anndata_dict1;n_hvgs=nothing, use_std=false)
    highly_variable_genes_bool = trues(size(anndata_dict1["X"])[1])
    if use_std
        gene_stds=vec(std(anndata_dict1["X"],dims=2))
        gene_names = anndata_dict1["var"]["_index"]
        sorted_gene_stds = sortperm(gene_stds,rev=true)
        if isnothing(n_hvgs)
            n_hvgs = min(2000,length(gene_stds)) 
        end
        highly_variable_genes_bool = [ in(name, gene_names[sorted_gene_stds[1:n_hvgs]]) ? true : false for name in gene_names ]
    else
        if isnothing(n_hvgs)
            if haskey(anndata_dict1["uns"],"hvg_meta_df")
                highly_variable_genes_col = "highly_variable-$(n_hvgs)"
                if haskey(anndata_dict1["uns"]["hvg_meta_df"],highly_variable_genes_col)
                    highly_variable_genes_bool .= anndata_dict1["uns"]["hvg_meta_df"][highly_variable_genes_col] .== 1
                end
            end
        end
    end
    return highly_variable_genes_bool
end

"""    subset_on_highly_variable_genes_bool(anndata_dict1,highly_variable_genes_bool,layer_index)
This function subsets the AnnData object to include only the highly variable genes.
"""
function subset_on_highly_variable_genes_bool(anndata_dict1,highly_variable_genes_bool,layer_index)
    anndata_dict1["X"] = anndata_dict1["X"][highly_variable_genes_bool,:]
    if haskey(anndata_dict1,"layers") && !isempty(anndata_dict1["layers"]) && layer_index != 0
        for layer in keys(anndata_dict1["layers"])
            anndata_dict1["layers"][layer] = anndata_dict1["layers"][layer][highly_variable_genes_bool,:]
        end
    end
    for key in keys(anndata_dict1["var"])
        anndata_dict1["var"][key] = anndata_dict1["var"][key][highly_variable_genes_bool]
    end
    return anndata_dict1
end

"""    center_and_scale_data_cols(anndata_dict1;center_cols = true,scale_cols = true)
This function centers and/or scales the data matrix in the AnnData object.
"""
function center_and_scale_data_cols(anndata_dict1;center_cols = true,scale_cols = true)
    xmat = anndata_dict1["X"]
    if center_cols
        xmat .= xmat .- mean(xmat,dims=2)
    end
    if scale_cols
        xmat .= xmat ./ std(xmat,dims=2)
    end
    anndata_dict1["X"] = xmat
    return anndata_dict1
end

"""    center_and_scale_matrix_cols(xmat;center_cols = true,scale_cols = true)
This function centers and/or scales the columns of the provided matrix.
"""
function center_and_scale_matrix_cols(xmat;center_cols = true,scale_cols = true)
    novariation = collect(1:size(xmat)[1])[[all(el .== el[1])  for el in eachrow(xmat)]]
    if length(novariation) > 0
       for i in novariation
        xmat[i,:] .+= 1e-32 .* randn(size(xmat)[2])
       end
    end
    if center_cols && !scale_cols
        xmat .= xmat .- nanmean(xmat,dims=2)
    elseif !center_cols && scale_cols
        xmat .= xmat ./ nanstd(xmat .+ 1e-8,dims=2)
    elseif center_cols && scale_cols
        xmat .= nanstandardize(xmat,dims=2)
    end
    if length(novariation) > 0
       for i in novariation
            xmat[i,:] .= 0.0
       end
    end
    return xmat
end


"""
    size_reorder_clusters(outputs_dict::OrderedDict{Symbol, Any})
This function reorders clusters based on their sizes in descending order.
"""
function size_reorder_clusters(outputs_dict::OrderedDict{Symbol, Any})
    #print(outputs_dict.keys)
    r = outputs_dict[:r_]
    I = length(r)
    T = [length(r[i]) for i in 1:I]
    N_t = [[length(r[i][t]) for t in 1:T[i]] for i in 1:I]
    Nkplus1 = sum([r[i][t][n] for i in 1:I for t in 1:T[i] for n in 1:N_t[i][t]])
    reindexing_Kplus1 = sortperm(Nkplus1,rev=true)
    reindexing_K = reindexing_Kplus1[reindexing_Kplus1 .!= maximum(reindexing_Kplus1)]
    for i in 1:I
        for t in 1:T[i]
            for n in 1:N_t[i][t]
                outputs_dict[:r_][i][t][n] = outputs_dict[:r_][i][t][n][reindexing_Kplus1]
            end
        end
    end
    for it in 1:sum(T)
        outputs_dict[:d_][it] = outputs_dict[:d_][it][reindexing_Kplus1]
    end
    # outputs_dict[:Nkplus1_] = outputs_dict[:Nkplus1_][reindexing_Kplus1]
    outputs_dict[:y_] = outputs_dict[:y_][reindexing_K]
    outputs_dict[:m_mu_] = outputs_dict[:m_mu_][reindexing_K]
    outputs_dict[:s_sq_mu_] = outputs_dict[:s_sq_mu_][reindexing_K]
    outputs_dict[:x_hat_] = outputs_dict[:x_hat_][reindexing_K]
    outputs_dict[:x_hat_sq_] = outputs_dict[:x_hat_sq_][reindexing_K]
    outputs_dict[:g1_] = outputs_dict[:g1_][reindexing_K]
    outputs_dict[:g2_] = outputs_dict[:g2_][reindexing_K]
    outputs_dict[:h1_] = outputs_dict[:h1_][reindexing_K]
    outputs_dict[:h2_] = outputs_dict[:h2_][reindexing_K]
    outputs_dict[:Nk_] = outputs_dict[:Nk_][reindexing_K]
    outputs_dict[:u_] = outputs_dict[:u_][reindexing_K]
    outputs_dict[:v_] = outputs_dict[:v_][reindexing_K]
    outputs_dict[:a_] = outputs_dict[:a_][reindexing_K]
    outputs_dict[:b_] = outputs_dict[:b_][reindexing_K]
    outputs_dict[:z_argmax] = [[argmax(outputs_dict[:r_][i][t][n]) for i in 1:I for t in 1:T[i] for n in 1:N_t[i][t]]]
    return outputs_dict
end

"""    calculate_s_values(sig_values;thresh=0.05)
This function calculates the Hoff S-value based on the provided significance values and threshold.
"""
function calculate_s_values(sig_values;thresh=0.05)
    hoff_s_value = mean(sig_values[sig_values .<= thresh])
    if isnan(hoff_s_value)
        hoff_s_value = 1.0
    end
    return hoff_s_value
end


"""
    return_occupied_clusters(r)
This function returns the indices of occupied clusters and their reindexing.
"""
function return_occupied_clusters(r)
    KMaxplus1 = length(r[1][1][1])
    I = length(r)
    T = [length(r[i]) for i in 1:I]
    N_t = [[length(r[i][t]) for t in 1:T[i]] for i in 1:I]
    z_argmax = [argmax(r[i][t][n]) for i in 1:I for t in 1:T[i] for n in 1:N_t[i][t]]
    cluster_occupancy_counts_dict = countmap(z_argmax)
    Nk = zeros(KMaxplus1)
    for k in 1:KMaxplus1
        if haskey(cluster_occupancy_counts_dict,k)
            Nk[k] = cluster_occupancy_counts_dict[k]
        end
    end
    N = Int(sum(Nk))
    occupied_cluster_indx = collect(1:KMaxplus1)[Nk .>=1]
    K_post = length(occupied_cluster_indx)
    remap_dict = Dict(occupied_cluster_indx[i] => i for i in 1:K_post)
    new_unique_cluster_indx = [i for i in 1:K_post]
    return KMaxplus1,K_post, N, Nk, occupied_cluster_indx,remap_dict,new_unique_cluster_indx
end

"""
    subset_on_occupied_clusters(outputs_dict::OrderedDict{Symbol, Any})
This function subsets the outputs dictionary to include only occupied clusters and reindexes them.
"""
function subset_on_occupied_clusters(outputs_dict::OrderedDict{Symbol, Any})
    r = outputs_dict[:r_];
    KMaxplus1,K_post, N, Nk, occupied_cluster_indx,remap_dict,new_unique_cluster_indx = return_occupied_clusters(r);
    outputs_dict[:Nk_] = Nk
    println("Reindexing clusters (Old Index => New Index)")
    for key in sort(collect(keys(remap_dict)))
        println("$(key) => $(remap_dict[key])")
    end
    println("___________________________________________________________")
    occupied_cluster_indx_kplus1 = [occupied_cluster_indx;KMaxplus1]
    I = length(r)
    T = [length(r[i]) for i in 1:I]
    N_t = [[length(r[i][t]) for t in 1:T[i]] for i in 1:I]
    for i in 1:I
        for t in 1:T[i]
            for n in 1:N_t[i][t]
                outputs_dict[:r_][i][t][n] = outputs_dict[:r_][i][t][n][occupied_cluster_indx_kplus1]
            end
        end
    end
    for it in 1:sum(T)
        outputs_dict[:d_][it] = outputs_dict[:d_][it][occupied_cluster_indx_kplus1]
    end
    outputs_dict[:y_] = outputs_dict[:y_][occupied_cluster_indx]
    outputs_dict[:m_mu_] = outputs_dict[:m_mu_][occupied_cluster_indx]
    outputs_dict[:s_sq_mu_] = outputs_dict[:s_sq_mu_][occupied_cluster_indx]
    outputs_dict[:x_hat_] = outputs_dict[:x_hat_][occupied_cluster_indx]
    outputs_dict[:x_hat_sq_] = outputs_dict[:x_hat_sq_][occupied_cluster_indx]
    outputs_dict[:g1_] = outputs_dict[:g1_][occupied_cluster_indx]
    outputs_dict[:g2_] = outputs_dict[:g2_][occupied_cluster_indx]
    outputs_dict[:h1_] = outputs_dict[:h1_][occupied_cluster_indx]
    outputs_dict[:h2_] = outputs_dict[:h2_][occupied_cluster_indx]
    outputs_dict[:Nk_] = outputs_dict[:Nk_][occupied_cluster_indx]
    outputs_dict[:u_] = outputs_dict[:u_][occupied_cluster_indx]
    outputs_dict[:v_] = outputs_dict[:v_][occupied_cluster_indx]
    outputs_dict[:a_] = outputs_dict[:a_][occupied_cluster_indx]
    outputs_dict[:b_] = outputs_dict[:b_][occupied_cluster_indx]
    outputs_dict[:z_argmax] = [remap_dict[outputs_dict[:z_argmax][1][i]] for i in 1:N]
    return outputs_dict
end


"""
    posterior_summaries(outputs_dict::OrderedDict{Symbol, Any}, used_representation::AbstractArray,used_representation_feature_name::AbstractArray,seed::Int;num_samples = 1000,s_value_thresh=0.05,return_data_frames = false,save_data_frames = true)
This function computes posterior summaries, including PIPs and S-values, and optionally saves them to CSV files.
"""
function posterior_summaries(outputs_dict::OrderedDict{Symbol, Any}, used_representation::AbstractArray,used_representation_feature_name::AbstractArray,seed::Int;num_samples = 1000,s_value_thresh=0.05,return_data_frames = false,save_data_frames = true)
    Random.seed!(seed)
    N = size(used_representation)[2]
    G = size(used_representation)[1]
    num_var_feat = G
    KMax_original = length(outputs_dict[:m_nu_])
    K_post = length(outputs_dict[:m_mu_])
    # used_representation_feature_name = anndata_dict1["var"]["_index"]
    # feature_median = median(used_representation,dims=2)
    feature_mean = mean(used_representation,dims=2)
    new_order = 1:length(used_representation_feature_name)#sortperm(used_representation_feature_name)
    used_representation_feature_name = used_representation_feature_name[new_order]
    x_mat = used_representation[new_order,:]
    # _prior_eta = 1/G
    # _prior_minus_eta = 1 - _prior_eta
    # prior_log_odds = log(_prior_eta / _prior_minus_eta)
    # _sigmasq = 1 ./ (outputs_dict[:b_] ./(outputs_dict[:a_]))
    # _slab_variance = 1 ./ hcat([((outputs_dict[:b_] .* outputs_dict[:v_][k]) ./ (outputs_dict[:a_] .*  outputs_dict[:u_][k] )  ) for k in 1:K]...)
    # var_dist_feature_included = _slab_variance .+  _slab_variance
    # _q_not_K_means = hcat([vec(mean(hcat([outputs_dict[:m_][k_prime] for k_prime in setdiff(eachindex(collect(1:K)), [k])]...), dims =2)) for k in 1:K]...)
    # _sq_feature_diff_q_mean = hcat([(outputs_dict[:m_][k] .-  _q_not_K_means[:,k] ).^2 for k in 1:K]...)
    # _sq_feature_nodiff_mean = hcat([(outputs_dict[:m_][k]).^2 for k in 1:K]...)
    # _sq_feature_diff_mean = hcat([(outputs_dict[:m_][k] .- feature_mean).^2 for k in 1:K]...)
    # # _sq_feature_diff_median = hcat([(outputs_dict[:m_][k] .- feature_median).^2 for k in 1:K]...)
    # # sigmoid((_slab_variance ./ (_sigmasq  .* var_dist_feature_included)) .* _sq_feature_diff_median + log.(_sigmasq ./var_dist_feature_included) .+ prior_log_odds)
    # hoff_pip = [ collect(col) for col in  eachcol(sigmoid((_slab_variance ./ (_sigmasq  .* var_dist_feature_included)) .* _sq_feature_diff_mean + log.(_sigmasq ./var_dist_feature_included) .+ prior_log_odds))]
    # sigmoid((_slab_variance ./ (_sigmasq  .* var_dist_feature_included)) .* _sq_feature_diff_q_mean + log.(_sigmasq ./var_dist_feature_included) .+ prior_log_odds)
    # sigmoid((_slab_variance ./ (_sigmasq  .* var_dist_feature_included)) .* _sq_feature_nodiff_mean + log.(_sigmasq ./var_dist_feature_included) .+ prior_log_odds)
    # hoff_pip = [vec(sigmoid((outputs_dict[:m_][k] .- feature_median).^2 .* _sigmasq  ./ (_sigmasq .* (outputs_dict[:v_][k] ./ outputs_dict[:u_][k] ))  .- log.( _sigmasq  ./ (_sigmasq .* (outputs_dict[:v_][k] ./ outputs_dict[:u_][k] ))) .+ log((1/G) / (1-1/G))) ) for k in 1:K]#outputs_dict[:y]#
    cluster_z_counts_dict = countmap(outputs_dict[:z_argmax][1])
    nk = zeros(length(outputs_dict[:m_mu_]))
    for k in 1:length(nk)
        if haskey(cluster_z_counts_dict,k)
            nk[k] = cluster_z_counts_dict[k]
        end
    end
    occupied_subset_bool = all(nk .>= 1)
    size_reordered_bool1 = all([outputs_dict[:Nk_][k-1] >= outputs_dict[:Nk_][k]  for k in collect(2:length(outputs_dict[:Nk_]))])
    size_reordered_bool2 = all([cluster_z_counts_dict[k-1] >= cluster_z_counts_dict[k]  for k in collect(2:length(outputs_dict[:Nk_]))])
    modified_outputs_dict_bool = false
    if !occupied_subset_bool
        println("Empty Clusters Detected! Subsetting on occupied clusters")
        println("___________________________________________________________")
        outputs_dict = subset_on_occupied_clusters(outputs_dict)
        modified_outputs_dict_bool = true
    end
    if (!size_reordered_bool1 ) || (!size_reordered_bool1 && !size_reordered_bool2)
        println("Clusters Not Ordered by Size! Reordering Clusters by Size")
        println("___________________________________________________________")
        outputs_dict = size_reorder_clusters(outputs_dict)
        modified_outputs_dict_bool = true
    end
    # m_ = outputs_dict[:m_]
    # s_sq_ = outputs_dict[:s_sq_]
    # y_pip_  = outputs_dict[:y_]
    r = outputs_dict[:r_]
    c = outputs_dict[:c]
    I_ = length(r)
    T = [length(r[i]) for i in 1:I_]
    N_t = [[length(r[i][t]) for t in 1:T[i]] for i in 1:I_]
    # cell_summaries_df = outputs_dict[:cell_summaries_df]
    filepath = outputs_dict[:filepath]
    unique_time_id = outputs_dict[:unique_time_id]
    # z_argmax_raw = outputs_dict[:z_argmax][1]
    G_post = length(used_representation_feature_name)
    # KMax = length(r[1][1])
    # N = length(r[1])
    # T = 1
    # Nk = sum(sum.(r))
    # occupied_cluster_indx = collect(1:KMax)[Nk .>=1]
    # K_post = length(occupied_cluster_indx)
    # remap_dict = Dict(occupied_cluster_indx[i] => i for i in 1:K_post)
    # 
    # KMax,K_post, _, Nk, occupied_cluster_indx,remap_dict,new_unique_cluster_indx = return_occupied_clusters(r)
    # z_argmax = [remap_dict[outputs_dict[:z_argmax][1][i]] for i in 1:N]
    z_argmax_post = outputs_dict[:z_argmax][1]
    a_post = hcat(outputs_dict[:a_]...)
    b_post = hcat(outputs_dict[:b_]...)
    ab_sigma_sq_post = b_post ./ (a_post .- 1)
    m_mu_post = hcat(outputs_dict[:m_mu_]...)
    s_sq_mu_post = hcat(outputs_dict[:s_sq_mu_]...)
    m_nu_post = hcat(outputs_dict[:m_nu_]...)[:,1]
    s_sq_nu_post = hcat(outputs_dict[:s_sq_nu_]...)[:,1]
    u_post = hcat(outputs_dict[:u_]...)
    v_post = hcat(outputs_dict[:v_]...)
    # hoff_pip_post = hcat(hoff_pip[occupied_cluster_indx]...)
    y_pip_post = hcat(outputs_dict[:y_]...)
    h1_post = hcat(outputs_dict[:h1_]...)
    h2_post = hcat(outputs_dict[:h2_]...)
    g1_post = outputs_dict[:g1_]
    g2_post = outputs_dict[:g2_]
    ragged_r = [r[i][t][n] for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]];
    ragged_c = [c[i][t][n] for i in 1:I_ for t in 1:T[i] for n in 1:N_t[i][t]];
    max_length_r = maximum(length.(ragged_r))
    max_length_c = maximum(length.(ragged_c))
    padded_r = [vcat(vec, fill(missing, max_length_r - length(vec))) for vec in ragged_r]
    padded_c = [vcat(vec, fill(missing, max_length_c - length(vec))) for vec in ragged_c]
    # y_post = hcat(y_[occupied_cluster_indx]...)
    d_post = [el for el in outputs_dict[:d_]]
    w1_post = outputs_dict[:w1_]
    w2_post = outputs_dict[:w2_]
    Nk_post = outputs_dict[:Nk_]
    K_post = length(Nk_post)
    new_unique_cluster_indx = [i for i in 1:K_post]
    mu_in = Matrix{Float64}(undef,G_post, K_post)
    mu_not_in = Matrix{Float64}(undef,G_post, K_post)
    sigma_sq_in = Matrix{Float64}(undef,G_post, K_post)
    sigma_sq_not_in = Matrix{Float64}(undef,G_post, K_post)
    posterior_es = Matrix{Float64}(undef,G_post, K_post)
    y_m_post = y_pip_post .* m_mu_post
    for k in new_unique_cluster_indx
        mu_in[:,k] = nanmean(x_mat[:,z_argmax_post .==k],dims=2)
        mu_not_in[:,k] = nanmean(x_mat[:,z_argmax_post .!=k],dims=2)
        if sum(z_argmax_post .==k) == 1
            sigma_sq_in[:,k] .= 0.0
        else
            sigma_sq_in[:,k] = nanvar(x_mat[:,z_argmax_post .==k], dims=2)
        end
        if sum(z_argmax_post .!=k) == 1
            sigma_sq_not_in[:,k] .= 0.0
        else
            sigma_sq_not_in[:,k] = nanvar(x_mat[:,z_argmax_post .!=k], dims=2)
        end
        # sigma_sq_in[:,k] = nanvar(x_mat[:,z_argmax_post .==k], dims=2)
        # sigma_sq_not_in[:,k] = nanvar(x_mat[:,z_argmax_post .!=k], dims=2)
        posterior_es[:,k] .= (y_m_post[:,k] .- nanmean(y_m_post[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)) ./sqrt.(0.5 .* (ab_sigma_sq_post[:,k]  .+ nanmean(ab_sigma_sq_post[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)))
    end
    # hoff_adj_weight = (1 .-sum(hoff_pip_post .>= 0.5,dims=2) ./K_post) ./(1 .- minimum(sum(hoff_pip_post .>= 0.5,dims=2)) ./K_post)
    y_adj_weight = (1 .-sum(y_pip_post .>= 0.5,dims=2) ./K_post) ./(1 .-  1.0 ./K_post) #minimum(sum(y_pip_post .>= 0.5,dims=2))
    # hoff_adj_pip_post = hoff_pip_post .*hoff_adj_weight
    y_adj_pip_post = y_pip_post .*y_adj_weight
    empirical_es = (mu_in .- mu_not_in) ./sqrt.(0.5 .* (sigma_sq_in .+ sigma_sq_not_in))
    empirical_ess = hcat([["0" for j in 1:G_post] for k in 1:K_post]...)
    posterior_ess = hcat([["0" for j in 1:G_post] for k in 1:K_post]...)
    empirical_ess[empirical_es .> 0] .= "+"
    empirical_ess[empirical_es .< 0] .= "-"
    posterior_ess[posterior_es .> 0] .= "+"
    posterior_ess[posterior_es .< 0] .= "-"
    # p_mu_lt_0 =  hoff_pip_post .* cdf.(Normal(0,1),- m_post ./sqrt.(s_sq_post)) 
    post_mu_dist = Normal.(m_mu_post,sqrt.(s_sq_mu_post))
    post_ab_sigma_sq_dist = InverseGamma.(a_post,b_post)
    # post_hoff_rho_dist = Bernoulli.(hoff_pip_post)
    # post_hoff_adj_rho_dist = Bernoulli.(hoff_adj_pip_post)
    post_y_rho_dist = Bernoulli.(clamp.(y_pip_post, 0.0, 1.0))
    post_y_adj_rho_dist = Bernoulli.(clamp.(y_adj_pip_post, 0.0, 1.0))
    post_mu_samples = [rand.(post_mu_dist) for s in 1:num_samples]
    post_ab_sigma_sq_samples = [rand.(post_ab_sigma_sq_dist) for s in 1:num_samples]
    # post_hoff_rho_samples = [rand.(post_hoff_rho_dist) for s in 1:num_samples]
    # post_hoff_adj_rho_samples = [rand.(post_hoff_adj_rho_dist) for s in 1:num_samples]
    post_y_rho_samples = [rand.(post_y_rho_dist) for s in 1:num_samples]
    post_y_adj_rho_samples = [rand.(post_y_adj_rho_dist) for s in 1:num_samples]
    # lsfr_mc_hoff_es_samples  = Vector{Matrix{Float64}}(undef,num_samples)
    # lsfr_mc_hoff_adj_es_samples  = Vector{Matrix{Float64}}(undef,num_samples)
    lsfr_mc_y_es_samples  = Vector{Matrix{Float64}}(undef,num_samples)
    lsfr_mc_y_adj_es_samples  = Vector{Matrix{Float64}}(undef,num_samples)
    for s in 1:num_samples
        # sampled_hoff_rho = post_hoff_rho_samples[s]
        # sampled_hoff_adj_rho = post_hoff_adj_rho_samples[s]
        sampled_y_rho = post_y_rho_samples[s]
        sampled_y_adj_rho = post_y_adj_rho_samples[s]
        sampled_mu = post_mu_samples[s]
        sampled_ab_sigma_sq = post_ab_sigma_sq_samples[s]
        # sampled_hoff_rho_mu = sampled_hoff_rho .* sampled_mu
        # sampled_hoff_adj_rho_mu = sampled_hoff_adj_rho .* sampled_mu
        sampled_y_rho_mu = sampled_y_rho .* sampled_mu
        sampled_y_adj_rho_mu = sampled_y_adj_rho .* sampled_mu
        # mc_hoff_es_sample =  Matrix{Float64}(undef,G_post, K_post)
        # mc_hoff_adj_es_sample =  Matrix{Float64}(undef,G_post, K_post)
        mc_es_y_sample =  Matrix{Float64}(undef,G_post, K_post)
        mc_y_adj_es_sample =  Matrix{Float64}(undef,G_post, K_post)
        for k in 1:K_post
            # mc_hoff_es_sample[:,k] .= sampled_hoff_rho_mu[:,k] .- mean(sampled_hoff_rho_mu[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)
            # mc_hoff_adj_es_sample[:,k] .= sampled_hoff_adj_rho_mu[:,k] .- mean(sampled_hoff_adj_rho_mu[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)
            mc_es_y_sample[:,k] .= sampled_y_rho_mu[:,k] .- mean(sampled_y_rho_mu[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2) ./ sqrt.(0.5 .* (sampled_ab_sigma_sq[:,k] .+ nanmean(sampled_ab_sigma_sq[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)))
            mc_y_adj_es_sample[:,k] .= sampled_y_adj_rho_mu[:,k] .- mean(sampled_y_adj_rho_mu[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2) ./ sqrt.(0.5 .* (sampled_ab_sigma_sq[:,k] .+ nanmean(sampled_ab_sigma_sq[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)))
        end
        # lsfr_mc_hoff_es_samples[s] = mc_hoff_es_sample
        # lsfr_mc_hoff_adj_es_samples[s] = mc_hoff_adj_es_sample
        lsfr_mc_y_es_samples[s] = mc_es_y_sample
        lsfr_mc_y_adj_es_samples[s] = mc_y_adj_es_sample
    end
    # lsfr_mc_hoff_es_samples = cat(lsfr_mc_hoff_es_samples...,dims=3)
    # lsfr_mc_hoff_adj_es_samples = cat(lsfr_mc_hoff_adj_es_samples...,dims=3)
    lsfr_mc_y_es_samples = cat(lsfr_mc_y_es_samples...,dims=3)
    lsfr_mc_y_adj_es_samples = cat(lsfr_mc_y_adj_es_samples...,dims=3)
    # adj_p_mu_lt_0 = hoff_adj_pip_post .* cdf.(Normal(0,1),- m_post ./sqrt.(lambda_sq_))
    # adj_p_mu_lt_0 = hoff_adj_pip_post .* cdf.(Normal(0,1),- m_post ./sqrt.(s_sq_post))
    # hoff_p_es_eq_0 = dropdims(mean(lsfr_mc_hoff_es_samples .== 0,dims=3),dims=3)
    # hoff_adj_p_es_eq_0 = dropdims(mean(lsfr_mc_hoff_adj_es_samples .== 0,dims=3),dims=3)
    y_p_es_eq_0 = dropdims(mean(lsfr_mc_y_es_samples .== 0,dims=3),dims=3)
    y_adj_p_es_eq_0 = dropdims(mean(lsfr_mc_y_adj_es_samples .== 0,dims=3),dims=3)

    # hoff_p_es_lt_0 = dropdims(mean(lsfr_mc_hoff_es_samples .< 0,dims=3),dims=3)
    # hoff_adj_p_es_lt_0 = dropdims(mean(lsfr_mc_hoff_adj_es_samples .< 0,dims=3),dims=3)
    y_p_es_lt_0 = dropdims(mean(lsfr_mc_y_es_samples .< 0,dims=3),dims=3)
    y_adj_p_es_lt_0 = dropdims(mean(lsfr_mc_y_adj_es_samples .< 0,dims=3),dims=3)


    # hoff_p_es_gt_0 = dropdims(mean(lsfr_mc_hoff_es_samples .> 0,dims=3),dims=3)
    # hoff_adj_p_es_gt_0 = dropdims(mean(lsfr_mc_hoff_adj_es_samples .> 0,dims=3),dims=3)
    y_p_es_gt_0 = dropdims(mean(lsfr_mc_y_es_samples .> 0,dims=3),dims=3)
    y_adj_p_es_gt_0 = dropdims(mean(lsfr_mc_y_adj_es_samples .> 0,dims=3),dims=3)
    # Make sure that the sum of the probabilities is 1
    # hoff_check_unity = vcat([hcat([isprobvec([hoff_p_es_eq_0[j,k],hoff_p_es_lt_0[j,k],hoff_p_es_gt_0[j,k] ]) for j in 1:G_post]...) for k in 1:K_post]...)
    # hoff_adj_check_unity = vcat([hcat([isprobvec([hoff_adj_p_es_eq_0[j,k],hoff_adj_p_es_lt_0[j,k],hoff_adj_p_es_gt_0[j,k] ]) for j in 1:G_post]...) for k in 1:K_post]...)
    y_check_unity = vcat([hcat([isprobvec([y_p_es_eq_0[j,k],y_p_es_lt_0[j,k],y_p_es_gt_0[j,k] ]) for j in 1:G_post]...) for k in 1:K_post]...)
    y_adj_check_unity = vcat([hcat([isprobvec([y_adj_p_es_eq_0[j,k],y_adj_p_es_lt_0[j,k],y_adj_p_es_gt_0[j,k] ]) for j in 1:G_post]...) for k in 1:K_post]...)
    # @test sum(hoff_check_unity .!= 1) == 0
    # @test sum(hoff_adj_check_unity .!= 1) == 0
    # @test sum(y_check_unity .!= 1) == 0
    # @test sum(y_adj_check_unity .!= 1) == 0

    # hoff_p_es_gteq_0 = hoff_p_es_gt_0 .+ hoff_p_es_eq_0
    # hoff_adj_p_es_gteq_0 = hoff_adj_p_es_gt_0 .+ hoff_adj_p_es_eq_0
    y_p_es_gteq_0 = y_p_es_gt_0 .+ y_p_es_eq_0
    y_adj_p_es_gteq_0 = y_adj_p_es_gt_0 .+ y_adj_p_es_eq_0
    # hoff_p_es_lteq_0 = hoff_p_es_lt_0 .+ hoff_p_es_eq_0
    # hoff_adj_p_es_lteq_0 = hoff_adj_p_es_lt_0 .+ hoff_adj_p_es_eq_0
    y_p_es_lteq_0 = y_p_es_lt_0 .+ y_p_es_eq_0
    y_adj_p_es_lteq_0 = y_adj_p_es_lt_0 .+ y_adj_p_es_eq_0
    # p_es_sum = hoff_p_es_eq_0 .+ hoff_p_es_lt_0 .+ hoff_p_es_gt_0
    # adj_p_es_sum = hoff_adj_p_es_eq_0 .+ hoff_adj_p_es_lt_0 .+ hoff_adj_p_es_gt_0
    # hoff_p_es_eq_0 = hoff_p_es_eq_0 ./p_es_sum
    # hoff_p_es_lt_0 = hoff_p_es_lt_0 ./p_es_sum
    # hoff_p_es_gt_0 = hoff_p_es_gt_0 ./p_es_sum
    # hoff_adj_p_es_eq_0 = hoff_adj_p_es_eq_0 ./adj_p_es_sum
    # hoff_adj_p_es_lt_0 = hoff_adj_p_es_lt_0 ./adj_p_es_sum
    # hoff_adj_p_es_gt_0 = hoff_adj_p_es_gt_0 ./adj_p_es_sum
    # # adj_p_mu_eq_0 = (1 .- hoff_adj_pip_post)
    # p_mu_gt_0 =  1 .- p_mu_lt_0 .- p_mu_eq_0
    # # adj_p_mu_gt_0 = 1 .- adj_p_mu_lt_0 .- adj_p_mu_eq_0
    # p_mu_lteq_0 = p_mu_lt_0 .+ p_mu_eq_0
    # # adj_p_mu_lteq_0 = adj_p_mu_lt_0 .+ adj_p_mu_eq_0
    # p_mu_gteq_0 = p_mu_gt_0 .+ p_mu_eq_0
    # # adj_p_mu_gteq_0 = adj_p_mu_gt_0 .+ adj_p_mu_eq_0
    # hoff_lfsr_kj = hcat([[minimum([p_mu_lteq_0[j,k],p_mu_gteq_0[j,k]]) for j in 1:G_post] for k in 1:K_post]...)
    # lfsr_j = minimum(hoff_lfsr_kj,dims=2)
    # hoff_lfsr_kj = hcat([[minimum([hoff_p_es_lteq_0[j,k],hoff_p_es_gteq_0[j,k]]) for j in 1:G_post] for k in 1:K_post]...)
    # hoff_adj_lfsr_kj = hcat([[minimum([hoff_adj_p_es_lteq_0[j,k],hoff_adj_p_es_gteq_0[j,k]]) for j in 1:G_post] for k in 1:K_post]...)
    y_lfsr_kj = hcat([[minimum([y_p_es_lteq_0[j,k],y_p_es_gteq_0[j,k]]) for j in 1:G_post] for k in 1:K_post]...)
    y_adj_lfsr_kj = hcat([[minimum([y_adj_p_es_lteq_0[j,k],y_adj_p_es_gteq_0[j,k]]) for j in 1:G_post] for k in 1:K_post]...)

    # hoff_min_lfsr_j = minimum(hoff_lfsr_kj,dims=2)
    # hoff_min_adj_lfsr_j = minimum(hoff_adj_lfsr_kj,dims=2)
    y_min_lfsr_j = minimum(y_lfsr_kj,dims=2)
    y_min_adj_lfsr_j = minimum(y_adj_lfsr_kj,dims=2)

    # hoff_avg_lfsr_j = mean(hoff_lfsr_kj,dims=2)
    # hoff_avg_adj_lfsr_j = mean(hoff_adj_lfsr_kj,dims=2)
    y_avg_lfsr_j = mean(y_lfsr_kj,dims=2)
    y_avg_adj_lfsr_j = mean(y_adj_lfsr_kj,dims=2)

    # hoff_adj_s_value = [calculate_s_values(hoff_adj_lfsr_kj[j,:];thresh=s_value_thresh) for j in 1:G_post]
    # hoff_s_value = [calculate_s_values(hoff_lfsr_kj[j,:];thresh=s_value_thresh) for j in 1:G_post]
    y_adj_s_value = [calculate_s_values(y_adj_lfsr_kj[j,:];thresh=s_value_thresh) for j in 1:G_post]
    y_s_value = [calculate_s_values(y_lfsr_kj[j,:];thresh=s_value_thresh) for j in 1:G_post]
    cell_summaries_df = DataFrame(cell_id = outputs_dict[:cell_ids], individuals_id = outputs_dict[:individuals_vec], time_id = outputs_dict[:time_vec], inferred_label = z_argmax_post)
    if haskey(outputs_dict,:true_labels)
        cell_summaries_df[:true_label] = outputs_dict[:true_labels]
    end
    if haskey(outputs_dict,:X_tsne)
        cell_summaries_df[!,:X_TSNE_1] = outputs_dict[:X_tsne][:,1]
        cell_summaries_df[!,:X_TSNE_2] = outputs_dict[:X_tsne][:,2]
    end
    if haskey(outputs_dict,:X_umap)
        cell_summaries_df[!,:X_UMAP_1] = outputs_dict[:X_umap][:,1]
        cell_summaries_df[!,:X_UMAP_2] = outputs_dict[:X_umap][:,2]
    end
    if haskey(outputs_dict,:X_pca)
        cell_summaries_df[!,:X_PCA_1] = outputs_dict[:X_pca][:,1]
        cell_summaries_df[!,:X_PCA_2] = outputs_dict[:X_pca][:,2]
    end
    if haskey(outputs_dict,:used_representation_tsne)
        cell_summaries_df[!,:Used_Representation_TSNE_1] = outputs_dict[:used_representation_tsne][:,1]
        cell_summaries_df[!,:Used_Representation_TSNE_2] = outputs_dict[:used_representation_tsne][:,2]
    end
    if haskey(outputs_dict,:used_representation_umap)
        cell_summaries_df[!,:Used_Representation_UMAP_1] = outputs_dict[:used_representation_umap][:,1]
        cell_summaries_df[!,:Used_Representation_UMAP_2] = outputs_dict[:used_representation_umap][:,2]
    end
    if haskey(outputs_dict,:used_representation_pca)
        cell_summaries_df[!,:Used_Representation_PCA_1] = outputs_dict[:used_representation_pca][:,1]
        cell_summaries_df[!,:Used_Representation_PCA_2] = outputs_dict[:used_representation_pca][:,2]
    end
    posterior_featurecluster_summaries_df = DataFrame(Gene = used_representation_feature_name)
    for k in 1:K_post
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_m_mu")] = m_mu_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_s_sq_mu")] = s_sq_mu_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_empirical_mu_in")] = mu_in[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_empirical_mu_not_in")] = mu_not_in[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_empirical_sigma_sq_in")] = sigma_sq_in[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_empirical_sigma_sq_not_in")] = sigma_sq_not_in[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_m_nu")] = m_nu_post
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_s_sq_nu")] = s_sq_nu_post
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_pip")] = hoff_pip_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_pip")] = y_pip_post[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_pip")] = hoff_adj_pip_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_pip")] = y_adj_pip_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_empirical_es")] = empirical_es[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_empirical_ess")] = empirical_ess[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_posterior_es")] = posterior_es[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_posterior_ess")] = posterior_ess[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_p_es_lt_0")] = hoff_p_es_lt_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_lt_0")] = y_p_es_lt_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_p_es_eq_0")] = hoff_p_es_eq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_eq_0")] = y_p_es_eq_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_p_es_gt_0")] = hoff_p_es_gt_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_gt_0")] = y_p_es_gt_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_p_es_lteq_0")] = hoff_p_es_lteq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_lteq_0")] = y_p_es_lteq_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_p_es_gteq_0")] = hoff_p_es_gteq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_gteq_0")] = y_p_es_gteq_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_p_es_lt_0")] = hoff_adj_p_es_lt_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_p_es_lt_0")] = y_adj_p_es_lt_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_p_es_eq_0")] = hoff_adj_p_es_eq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_p_es_eq_0")] = y_adj_p_es_eq_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_p_es_gt_0")] = hoff_adj_p_es_gt_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_p_es_gt_0")] = y_adj_p_es_gt_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_p_es_lteq_0")] = hoff_adj_p_es_lteq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_p_es_lteq_0")] = y_adj_p_es_lteq_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_p_es_gteq_0")] = hoff_adj_p_es_gteq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_p_es_gteq_0")] = y_adj_p_es_gteq_0[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_lfsr_kj")] = hoff_lfsr_kj[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_lfsr_kj")] = y_lfsr_kj[:,k]
        # posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_hoff_adj_lfsr_kj")] = hoff_adj_lfsr_kj[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_adj_lfsr_kj")] = y_adj_lfsr_kj[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_h1")] = h1_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_h2")] = h2_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_a")] = a_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_b")] = b_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_u")] = u_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_v")] = v_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_g1")] = g1_post[k] .* ones(G_post)
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_g2")] = g2_post[k] .* ones(G_post)
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_Nk")] =Nk_post[k] .* ones(G_post)
    end
    # cluster_summaries_df = DataFrame(cluster_id = ["Cluster-$k" for k in 1:K_post],Nk = Nk_post, h1 = h1_post, h2 = h2_post, g1 = g1_post, g2 = g2_post, u = u_post, v = v_post)
    conditions_summaries_df =DataFrame(individuals_id = [i for i in 1:I_ for t in 1:T[i]], time_id = [t for i in 1:I_ for t in 1:T[i]], w1 = w1_post, w2 = w2_post)
    d_post_df = DataFrame(permutedims(hcat(d_post...))[:,1:size(permutedims(hcat(d_post...)))[2]-1],:auto)
    # rename!(d_post_df, [ k < K_post+1 ? Symbol("Cluster-"*string(k)*"_d") : Symbol("Cluster-Kplus1_d") for k in 1:(K_post+1)])
    rename!(d_post_df, [Symbol("Cluster-"*string(k)*"_d") for k in 1:(K_post)])
    r_df = DataFrame(permutedims(hcat(padded_r...))[:,1:size(permutedims(hcat(padded_r...)))[2]-1],:auto)
    c_df = DataFrame(permutedims(hcat(padded_c...)),:auto)
    rename!(r_df, [Symbol("Cluster-"*string(k)*"_r") for k in 1:K_post])
    rename!(c_df, [Symbol("Condition-"*string(k)*"_c") for k in 1:max_length_c])
    cell_summaries_df = hcat(cell_summaries_df,r_df,c_df)
    conditions_summaries_df = hcat(conditions_summaries_df,d_post_df)
    # posterior_featurecluster_summaries_df[!,Symbol("hoff_min_lfsr_j")] = vec(hoff_min_lfsr_j)
    posterior_featurecluster_summaries_df[!,Symbol("y_min_lfsr_j")] = vec(y_min_lfsr_j)
    # posterior_featurecluster_summaries_df[!,Symbol("hoff_min_adj_lfsr_j")] = vec(hoff_min_adj_lfsr_j)
    posterior_featurecluster_summaries_df[!,Symbol("y_min_adj_lfsr_j")] = vec(y_min_adj_lfsr_j)
    # posterior_featurecluster_summaries_df[!,Symbol("hoff_avg_lfsr_j")] = vec(hoff_avg_lfsr_j)
    posterior_featurecluster_summaries_df[!,Symbol("y_avg_lfsr_j")] = vec(y_avg_lfsr_j)
    # posterior_featurecluster_summaries_df[!,Symbol("hoff_avg_adj_lfsr_j")] = vec(hoff_avg_adj_lfsr_j)
    posterior_featurecluster_summaries_df[!,Symbol("y_avg_adj_lfsr_j")] = vec(y_avg_adj_lfsr_j)
    # posterior_featurecluster_summaries_df[!,Symbol("hoff_s_$(join(split("$(s_value_thresh)","."),""))_value")] = hoff_s_value
    posterior_featurecluster_summaries_df[!,Symbol("y_s_$(join(split("$(s_value_thresh)","."),""))_value")] = y_s_value
    # posterior_featurecluster_summaries_df[!,Symbol("hoff_adj_s_$(join(split("$(s_value_thresh)","."),""))_value")] = hoff_adj_s_value
    posterior_featurecluster_summaries_df[!,Symbol("y_adj_s_$(join(split("$(s_value_thresh)","."),""))_value")] =y_adj_s_value
    filepath = outputs_dict[:filepath];
    unique_time_id = outputs_dict[:unique_time_id];
    if haskey(outputs_dict,:run_name)
        run_name = outputs_dict[:run_name]
    else
        run_name = ""
    end
    feature_summaries_filename = "$filepath/$(run_name)_FeaturesClusters_summaries-$(unique_time_id).csv";
    cell_summaries_filename = "$filepath/$(run_name)_Cells_summaries-$(unique_time_id).csv";
    # cluster_summaries_filename = "$filepath/$(run_name)_cluster_summaries-$(unique_time_id).csv";
    condition_summaries_filename = "$filepath/$(run_name)_Conditions_summaries-$(unique_time_id).csv";
    posterior_genecluster_summaries_df = nothing
    if save_data_frames
        CSV.write(feature_summaries_filename,posterior_featurecluster_summaries_df)
        CSV.write(cell_summaries_filename,cell_summaries_df)
        CSV.write(condition_summaries_filename,conditions_summaries_df)
        # CSV.write(cluster_summaries_filename,cluster_summaries_df)
    end
    if haskey(outputs_dict,:gene_factor_loadings)
        if !isnothing(outputs_dict[:gene_factor_loadings])
            projection_filename = "$filepath/$(run_name)_FactorLoading-$(unique_time_id).csv";
            GeneFactorSummaries_filename = "$filepath/$(run_name)_GeneFactorSummaries-$(unique_time_id).csv";
            if typeof(outputs_dict[:gene_factor_loadings]) <: DataFrame
                if eltype(outputs_dict[:gene_factor_loadings][:,1]) <: String
                    projection_matrix = Float64.(Matrix(outputs_dict[:gene_factor_loadings][:,2:end]))
                    posterior_genecluster_summaries_df = posterior_summaries_on_projection(y_pip_post, m_mu_post,s_sq_mu_post,a_post,b_post,projection_matrix,outputs_dict[:gene_factor_loadings][:,1],names(outputs_dict[:gene_factor_loadings][:,2:end]),seed,filepath;num_samples = num_samples,s_value_thresh=s_value_thresh,return_data_frames = true,save_data_frames = true,run_name = run_name,unique_time_id=unique_time_id)
                end 
            end
            if save_data_frames
                CSV.write(projection_filename,outputs_dict[:gene_factor_loadings])
                CSV.write(GeneFactorSummaries_filename,posterior_genecluster_summaries_df)
            end
        end
    end
    if return_data_frames && modified_outputs_dict_bool
        return outputs_dict,z_argmax_post,posterior_featurecluster_summaries_df, cell_summaries_df, conditions_summaries_df,posterior_genecluster_summaries_df
    elseif return_data_frames && !modified_outputs_dict_bool
        return nothing,nothing,posterior_featurecluster_summaries_df,cell_summaries_df,conditions_summaries_df,posterior_genecluster_summaries_df
    elseif !return_data_frames && modified_outputs_dict_bool
        return outputs_dict,z_argmax_post,nothing,nothing,nothing,nothing
    else
        return nothing,nothing,nothing,nothing,nothing,nothing  
    end
end


"""
    posterior_summaries_on_projection(y_pip_post::Matrix{Float64}, m_mu_post::Matrix{Float64},s_sq_mu_post::Matrix{Float64},a_post::Matrix{Float64},b_post::Matrix{Float64},projection_matrix::Matrix{Float64},used_representation_feature_name::AbstractArray,gene_names::AbstractArray,seed::Int,filepath;num_samples = 1000,s_value_thresh=0.05,return_data_frames = false,save_data_frames = true,run_name = "",unique_time_id="")
This function computes posterior summaries on a projection matrix.
""" 
function posterior_summaries_on_projection(y_pip_post::Matrix{Float64}, m_mu_post::Matrix{Float64},s_sq_mu_post::Matrix{Float64},a_post::Matrix{Float64},b_post::Matrix{Float64},projection_matrix::Matrix{Float64},used_representation_feature_name::AbstractArray,gene_names::AbstractArray,seed::Int,filepath;num_samples = 1000,s_value_thresh=0.05,return_data_frames = false,save_data_frames = true,run_name = "",unique_time_id="")
    Random.seed!(seed)
    G = size(projection_matrix)[1]
    J= size(projection_matrix)[2]
    num_latent_feat = G
    num_genes = J
    K_post = size(y_pip_post)[2]
    K_post = length(Nk_post)
    new_unique_cluster_indx = [i for i in 1:K_post]
    posterior_es = Matrix{Float64}(undef,J, K_post)
    ab_sigma_sq_post = b_post ./ (a_post .- 1)
    y_m_post = y_pip_post .* m_mu_post
    gene_mean_post = permutedims(permutedims(y_m_post) * projection_matrix)
    for k in new_unique_cluster_indx
        # sigma_sq_in[:,k] = nanvar(x_mat[:,z_argmax_post .==k], dims=2)
        # sigma_sq_not_in[:,k] = nanvar(x_mat[:,z_argmax_post .!=k], dims=2)
        posterior_es[:,k] .= (gene_mean_post[:,k] .- nanmean(gene_mean_post[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2)) ./ sqrt.(0.5 .* (diag(permutedims(projection_matrix) * diagm(ab_sigma_sq_post[:,k]) * projection_matrix) .+ diag(permutedims(projection_matrix) * diagm(vec(nanmean(ab_sigma_sq_post[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2))) * projection_matrix)))
    end
    # hoff_adj_weight = (1 .-sum(hoff_pip_post .>= 0.5,dims=2) ./K_post) ./(1 .- minimum(sum(hoff_pip_post .>= 0.5,dims=2)) ./K_post)
    # hoff_adj_pip_post = hoff_pip_post .*hoff_adj_weight
    posterior_ess = hcat([["0" for j in 1:J] for k in 1:K_post]...)
    posterior_ess[posterior_es .> 0] .= "+"
    posterior_ess[posterior_es .< 0] .= "-"
    # p_mu_lt_0 =  hoff_pip_post .* cdf.(Normal(0,1),- m_post ./sqrt.(s_sq_post)) 
    post_mu_dist = Normal.(m_mu_post,sqrt.(s_sq_mu_post))
    # post_hoff_rho_dist = Bernoulli.(hoff_pip_post)
    # post_hoff_adj_rho_dist = Bernoulli.(hoff_adj_pip_post)
    post_y_rho_dist = Bernoulli.(y_pip_post)
    post_mu_samples = [rand.(post_mu_dist) for s in 1:num_samples]
    post_ab_sigma_sq_dist = InverseGamma.(a_post,b_post)
    post_ab_sigma_sq_samples = [rand.(post_ab_sigma_sq_dist) for s in 1:num_samples]
    # post_hoff_rho_samples = [rand.(post_hoff_rho_dist) for s in 1:num_samples]
    # post_hoff_adj_rho_samples = [rand.(post_hoff_adj_rho_dist) for s in 1:num_samples]
    post_y_rho_samples = [rand.(post_y_rho_dist) for s in 1:num_samples]
    lsfr_mc_y_es_samples  = Vector{Matrix{Float64}}(undef,num_samples)
    for s in 1:num_samples
        sampled_y_rho = post_y_rho_samples[s]
        sampled_mu = post_mu_samples[s]
        sampled_ab_sigma_sq = post_ab_sigma_sq_samples[s]
        sampled_y_rho_mu = sampled_y_rho .* sampled_mu
        sampled_gene_mean = permutedims(permutedims(sampled_y_rho_mu) * projection_matrix)
        mc_es_y_sample =  Matrix{Float64}(undef,J, K_post)
        for k in 1:K_post
            mc_es_y_sample[:,k] .= sampled_gene_mean[:,k] .- mean(sampled_gene_mean[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2) ./ sqrt.(0.5 .* (diag(permutedims(projection_matrix) * diagm(sampled_ab_sigma_sq[:,k]) * projection_matrix) .+ diag(permutedims(projection_matrix) * diagm(vec(nanmean(sampled_ab_sigma_sq[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2))) * projection_matrix)))
        end
        lsfr_mc_y_es_samples[s] = mc_es_y_sample
    end
    lsfr_mc_y_es_samples = cat(lsfr_mc_y_es_samples...,dims=3)
    y_p_es_eq_0 = dropdims(mean(lsfr_mc_y_es_samples .== 0,dims=3),dims=3)

    y_p_es_lt_0 = dropdims(mean(lsfr_mc_y_es_samples .< 0,dims=3),dims=3)


    y_p_es_gt_0 = dropdims(mean(lsfr_mc_y_es_samples .> 0,dims=3),dims=3)
    y_check_unity = vcat([hcat([isprobvec([y_p_es_eq_0[j,k],y_p_es_lt_0[j,k],y_p_es_gt_0[j,k] ]) for j in 1:J]...) for k in 1:K_post]...)

    y_p_es_gteq_0 = y_p_es_gt_0 .+ y_p_es_eq_0
    y_p_es_lteq_0 = y_p_es_lt_0 .+ y_p_es_eq_0

    y_lfsr_kj = hcat([[minimum([y_p_es_lteq_0[j,k],y_p_es_gteq_0[j,k]]) for j in 1:J] for k in 1:K_post]...)

    y_min_lfsr_j = minimum(y_lfsr_kj,dims=2)

    y_avg_lfsr_j = mean(y_lfsr_kj,dims=2)


    y_s_value = [calculate_s_values(y_lfsr_kj[j,:];thresh=s_value_thresh) for j in 1:J]
    posterior_featurecluster_summaries_df = DataFrame(Gene = gene_names)
    for k in 1:K_post
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_gene_mean_in")] = gene_mean_post[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_ab_sigma_sq_in")] = diag(permutedims(projection_matrix) * diagm(ab_sigma_sq_post[:,k]) * projection_matrix) 
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_gene_mean_out")] = vec(nanmean(gene_mean_post[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2))
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_ab_sigma_sq_out")] = diag(permutedims(projection_matrix) * diagm(vec(nanmean(ab_sigma_sq_post[:,setdiff(eachindex(new_unique_cluster_indx), [k])],dims=2))) * projection_matrix)
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_posterior_es")] = posterior_es[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_posterior_ess")] = posterior_ess[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_lt_0")] = y_p_es_lt_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_eq_0")] = y_p_es_eq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_gt_0")] = y_p_es_gt_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_lteq_0")] = y_p_es_lteq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_p_es_gteq_0")] = y_p_es_gteq_0[:,k]
        posterior_featurecluster_summaries_df[!,Symbol("Cluster-"*string(k)*"_y_lfsr_kj")] = y_lfsr_kj[:,k]
    end
    posterior_featurecluster_summaries_df[!,Symbol("y_min_lfsr_j")] = vec(y_min_lfsr_j)
    posterior_featurecluster_summaries_df[!,Symbol("y_avg_lfsr_j")] = vec(y_avg_lfsr_j)
    posterior_featurecluster_summaries_df[!,Symbol("y_s_$(join(split("$(s_value_thresh)","."),""))_value")] = y_s_value
    feature_summaries_filename = "$filepath/$(run_name)_GeneProjectionsClusters_summaries-$(unique_time_id).csv";
    if save_data_frames
        CSV.write(feature_summaries_filename,posterior_featurecluster_summaries_df)
    end
    if return_data_frames
        return posterior_featurecluster_summaries_df
    end
end


"""
    make_new_embeddings(outputs_dict::Dict,used_representation::Matrix{Float64},used_representation_feature_name::AbstractArray,seed::Int; samples_as_rows=false)
This function generates new embeddings (PCA, t-SNE, UMAP) from the used representation matrix and adds them to the outputs dictionary.
"""
function make_new_embeddings(outputs_dict,used_representation,used_representation_feature_name,seed; samples_as_rows=false)
    Random.seed!(seed)
    new_order = sortperm(used_representation_feature_name)
    used_representation_feature_name = used_representation_feature_name[new_order]
    if !samples_as_rows
        x_mat = deepcopy(used_representation[new_order,:])
        x_mat = permutedims(x_mat)
    else
        x_mat = deepcopy(used_representation[:,new_order])
    end
    N = size(x_mat)[2]
    G = size(x_mat)[1]
    pca_model = fit(PCA, x_mat; maxoutdim=2)  # Reduce to 2 dimensions
    pca_transformed = MultivariateStats.transform(pca_model, x_mat)
    tsne_result = tsne(x_mat, 2, 0 ,1000, 30.0)
    umap_model = UMAP_(x_mat, 2;n_neighbors=15, min_dist=0.1)
    umap_result = UMAP.transform(umap_model, x_mat)
    outputs_dict[:used_representation_pca] = pca_transformed
    outputs_dict[:used_representation_tsne] = tsne_result
    outputs_dict[:used_representation_umap] = umap_result
    return outputs_dict
end

"""
    generate_report(script_path::String, arg1::String, arg2::String, arg3::String)::String
This function generates a report by executing an external Python script with the provided arguments.
"""
function generate_report(script_path::String, arg1::String, arg2::String, arg3::String)::String
    # Construct the shell command with the script path as an argument
    cmd = `python -u $script_path --resultsfiledir $arg1 --datadir $arg2 --make_abridged $arg3`
    
    try
        # Run the command and wait for it to complete
        run(cmd)
        # If it finishes successfully, return a success message
        return "Successfully created a report using '$script_path'."
    catch err
        # Capture the error message
        error_msg = String(err)
        return "Could not create a report. \nThe script's error is the following:\n\n$error_msg"
    end
end


"""
    submit_slurm_job_to_generate_report(script_path::String, arg1::String, arg2::String, arg3::String)::String
This function submits a SLURM job to generate a report by executing an external Python script with the provided arguments.
"""
function submit_slurm_job_to_generate_report(script_path::String, arg1::String, arg2::String, arg3::String)::String
    # Construct the shell command with the script path as an argument
    cmd = `sbatch -J report_generation -N 1 -c 1 -t 72:00:00 --mem=64GB -o /users/cnwizu/scratch/%x-log-%j.out -e /users/cnwizu/scratch/%x-log-%j.err --mail-type=END,FAIL --mail-user=chibuikem_nwizu@brown.edu --wrap="module load julia; module load llvm/16.0.2; module load r/4.4.0-yycctsj ; module load pcre2/10.42 cuda/12.1.1 texlive/20220321; module load cmake/3.26.3; module load libgit2/1.6.4; module load geos/3.11.2; module load libpng/1.6.39; module load gdal/3.7.0 proj/9.2.0; module load nlopt/2.7.1; source /users/cnwizu/data/cnwizu/cdHDPlmm/.venv/bin/activate ; export LD_PRELOAD=/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_def.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_avx2.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_core.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_lp64.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_thread.so:/gpfs/runtime/opt/intel/2020.2/lib/intel64_lin/libiomp5.so; python -u $script_path --resultsfiledir $arg1 --datadir $arg2 --make_abridged $arg3 "`
    # `python -u $script_path --resultsfiledir $arg1 --datadir $arg2 --make_abridged $arg3`
    
    try
        # Run the command and wait for it to complete
        run(cmd)
        # If it finishes successfully, return a success message
        return "Successfully submitted job created a report using '$script_path'."
    catch err
        # Capture the error message
        error_msg = String(err)
        return "Could not create a report. \nThe script's error is the following:\n\n$error_msg"
    end
end

# # Function to compute a cost matrix (e.g., based on overlap)
function compute_cost_matrix(labels1, labels2, K)
    cost_matrix = zeros(K, K)
    for i in 1:K
        for j in 1:K
            cost_matrix[i, j] = -sum((labels1 .== i) .& (labels2 .== j))  # Negative overlap (maximize match)
        end
    end
    return cost_matrix
end

# Function to relabel clusters using Hungarian algorithm
function relabel_clusters(labels, reference_labels, K)
    cost_matrix = compute_cost_matrix(reference_labels, labels, K)
    assignment, _ = hungarian(cost_matrix)  # Get optimal assignment

    label_map = Dict(j => i for (i, j) in enumerate(assignment))  # Map old labels to new
    return [label_map[l] for l in labels]  # Reassign labels
end
