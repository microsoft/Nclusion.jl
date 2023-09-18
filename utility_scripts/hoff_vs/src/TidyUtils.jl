module TidyUtils
    using Random
    using Distributions
    using Turing
    using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
    using Test
    import Debugger
    using CSV,DataFrames
    using RCall
    using Gnuplot

           #From synDataPreprocess.jl
    export recursive_flatten,
           outermelt, 
           innermelt,

           
           tidify_data,
           tidify_labels,
           tidify_mixture_weights,
           tidify_state_params,
           tidify_state_hyperprior_params,tidify_θ,
           tidify_τ_μ_expected_value,
           tidify_log_τ_kj_expected_value,
           tidify_log_π_expected_value,
           tidify_cttprime,
           tidify_rtik,
           tidify_clustering_params,
           tidify_time_params,
           tidify_scalar_hyperprior_params,
           tidify_ηtkj,
           tidify_error_state_params,
           tidy_get_Nt_from_xmat,
           tidy_get_Nt_from_rtikmat,
           tidy_get_next_chain_rows_indices,
           tidy_make_chain,
           tidy_truncate_chain,
           tidy_get_timeranges,

           IndexCounter,
           increment!,
           getCount


    function recursive_flatten(x::AbstractArray)
        if any(a->typeof(a)<:AbstractArray, x)#eltype(x) <: Vector
            recursive_flatten(vcat(x...))
        else
            return x
        end
    end
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
    function tidify_data(x;column_name="expression")
        T = length(x)
        N_t = [length(el) for el in x]
        G = length(x[1][1])
        gene_ids = collect(1:G)
        N = sum(N_t)
        cell_ids = collect(1:N)
        timepoints = collect(1:T)
        time = innermelt(timepoints,G .* N_t)
        gene = outermelt(gene_ids,N)
        cell = innermelt(cell_ids,G .* ones(Int,N))
        flatten_x = recursive_flatten(x)
        tidy_x = DataFrame(time=time,gene=gene,cell=cell)
        tidy_x[!,column_name] = flatten_x
        return tidy_x
    end
    function tidify_data(x)
        T = length(x)
        N_t = [length(el) for el in x]
        G = length(x[1][1])
        gene_ids = collect(1:G)
        N = sum(N_t)
        cell_ids = collect(1:N)
        timepoints = collect(1:T)
        time = innermelt(timepoints,G .* N_t)
        gene = outermelt(gene_ids,N)
        cell = innermelt(cell_ids,G .* ones(Int,N))
        flatten_x = recursive_flatten(x)
        nrows = length(flatten_x)
        ncols = 4
        tidy_x = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        tidy_x[:,1 ] = time 
        tidy_x[:,2 ] = cell
        tidy_x[:,3 ] = gene
        tidy_x[:,4 ] = flatten_x #hcat([time,cell,gene,flatten_x]...)
        # tidy_x[!,column_name] = flatten_x
        return tidy_x
    end
    function tidify_labels(z)
        T = length(z)
        N_t = [length(el) for el in z]
        K = length(unique(recursive_flatten(z)))
        state_ids = collect(1:K)
        N = sum(N_t)
        cell_ids = collect(1:N)
        timepoints = collect(1:T)
        time = innermelt(timepoints, N_t)
        # gene = outermelt(gene_ids,N)
        cell = innermelt(cell_ids, ones(Int,N))
        flatten_z = recursive_flatten(z)
        nrows = length(flatten_z)
        ncols = 3
        tidy_z = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        tidy_z[:,1 ] = time 
        tidy_z[:,2 ] = cell
        tidy_z[:,3 ] = flatten_z #hcat([time,cell,gene,flatten_x]...)
        # tidy_x[!,column_name] = flatten_x
        return tidy_z
    end
    function tidify_mixture_weights(π_)
        T = length(π_)
        K =  length(π_[1])
        state_ids = collect(1:K)
        timepoints = collect(1:T)
        time = innermelt(timepoints, K)
        states = outermelt(state_ids,T)
        
        flatten_π_ = recursive_flatten(π_)
        nrows = length(flatten_π_)
        ncols = 3
        tidy_π_ = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        tidy_π_[:,1 ] = time 
        tidy_π_[:,2 ] = states
        tidy_π_[:,3 ] = flatten_π_ #hcat([time,cell,gene,flatten_x]...)
        # tidy_x[!,column_name] = flatten_x
        return tidy_π_
    end
    function tidify_state_params(args...)
        K = length(args[1])
        G = length(args[1][1])
        num_params = length(args)
        state_ids = collect(1:K)
        gene_ids = collect(1:G)
        states = innermelt(state_ids,G)
        genes = outermelt(gene_ids,K)
        nrows = K*G
        ncols = num_params + 2 
        params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        for i in 1:num_params
            params_vec[i] = vcat(args[i]...)
        end
        params_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        params_mat[:,1] = states
        params_mat[:,2] = genes
        params_mat[:,3:end] = hcat(params_vec...)
        return params_mat
    end
    function tidify_state_hyperprior_params(args...)
        G = length(args[1])
        num_params = length(args)
        gene_ids = collect(1:G)
        nrows = G
        ncols = num_params + 1 
        params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        for i in 1:num_params
            params_vec[i] = vcat(args[i]...)
        end
        params_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        params_mat[:,1] = gene_ids
        params_mat[:,2:end] = hcat(params_vec...)
        return params_mat
    end
    function tidify_θ(θ_hat_vec)
        T = length(θ_hat_vec)
        Kplus = length(θ_hat_vec[1])
        timepoints = collect(1:T)
        state_ids = collect(1:Kplus)
        time = innermelt(timepoints,Kplus)
        states = outermelt(state_ids,T)
        nrows = Kplus*T
        ncols = 3 
        # params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        # for i in 1:num_params
        #     params_vec[i] = vcat(args[i]...)
        # end
        θ_hat_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        θ_hat_mat[:,1] = time
        θ_hat_mat[:,2] = states
        θ_hat_mat[:,end] = recursive_flatten(θ_hat_vec)
        return θ_hat_mat
    end
    function tidify_τ_μ_expected_value(e_τ_μ_tikj)
        T = length(e_τ_μ_tikj)
        N_t = [ length(el) for el in e_τ_μ_tikj]
        N = sum(N_t)
        G = length(e_τ_μ_tikj[1][1][1])
        K = length(e_τ_μ_tikj[1][1]) 
        gene_ids = collect(1:G)
        cell_ids = collect(1:N)
        timepoints = collect(1:T)
        state_ids= collect(1:K)
        timepoint_freq = countmap(Int.(xmat[1:G:end,1]))
        nrows = N*K*G
        ncols = 5
        e_τ_μ_tikj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* G .* N_t)
        cells = innermelt(cell_ids,K*G)
        states = outermelt(innermelt(state_ids,G),N)
        genes = outermelt(gene_ids,N*K)
        e_τ_μ_tikj_mat[:,1] = time
        e_τ_μ_tikj_mat[:,2] = cells
        e_τ_μ_tikj_mat[:,3] = states
        e_τ_μ_tikj_mat[:,4] = genes
        e_τ_μ_tikj_mat[:,5] = recursive_flatten(e_τ_μ_tikj)
        
        return e_τ_μ_tikj_mat
    end
    function tidify_log_τ_kj_expected_value(e_log_τ_kj_vec )
        K = length(e_log_τ_kj_vec)
        G = length(e_log_τ_kj_vec[1])
        state_ids = collect(1:K)
        gene_ids = collect(1:G)
        nrows = K*G
        ncols = 3
        states = innermelt(state_ids,G)
        genes = outermelt(gene_ids,K)
        e_log_τ_kj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        e_log_τ_kj_mat[:,1] = states
        e_log_τ_kj_mat[:,2] = genes
        e_log_τ_kj_mat[:,3] = recursive_flatten(e_log_τ_kj_vec)
        return e_log_τ_kj_mat
    end
    function tidify_log_π_expected_value(e_log_π)
        T = length(e_log_π)
        Kplus =length(e_log_π[1])
        timepoints = collect(1:T)
        state_ids = collect(1:Kplus)
        time = innermelt(timepoints,Kplus)
        states = outermelt(state_ids,T)
        nrows = Kplus*T
        ncols = 3 
        e_log_π_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        e_log_π_mat[:,1] = time
        e_log_π_mat[:,2] = states
        e_log_π_mat[:,3] = recursive_flatten(e_log_π)
        return e_log_π_mat
    end
    function tidify_cttprime(c_ttprime_vec)
        T = length(c_ttprime_vec)
        timepoints = collect(1:T)
        time1 = innermelt(timepoints,T)
        time2 = outermelt(timepoints,T)
        nrows = T*T
        ncols = 3 
        # params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        # for i in 1:num_params
        #     params_vec[i] = vcat(args[i]...)
        # end
        c_ttprime_hat_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        c_ttprime_hat_mat[:,1] = time1
        c_ttprime_hat_mat[:,2] = time2
        c_ttprime_hat_mat[:,end] = recursive_flatten(c_ttprime_vec)
        return c_ttprime_hat_mat
    end
    function tidify_rtik(rtik)
        T = length(rtik)
        N_t = [length(el) for el in rtik]
        K = length(rtik[1][1])
        N = sum(N_t)
        cell_ids = collect(1:N)
        timepoints = collect(1:T)
        state_ids = collect(1:K )
        
        nrows = N*K
        ncols = 4
        rtik_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* N_t)
        cells = innermelt(cell_ids,K)
        states = outermelt(state_ids,N)
        rtik_mat[:,1] = time
        rtik_mat[:,2] = cells
        rtik_mat[:,3] = states
        rtik_mat[:,4] =  recursive_flatten(rtik)
        #@test all([all(tidify_rtik(rtik)[:,col] .== rtik_mat[:,col]) for col in 1:size(rtik_mat)[2]])
        return rtik_mat
    end
    function tidify_clustering_params(args...)
        K = length(args[1])
        num_params = length(args)
        state_ids = collect(1:K)
        nrows = K
        ncols = num_params + 1 
        params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        for i in 1:num_params
            params_vec[i] = vcat(args[i]...)
        end
        params_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        params_mat[:,1] = state_ids
        params_mat[:,2:end] = hcat(params_vec...)
        return params_mat
    end
    function tidify_time_params(args...)
        T = length(args[1])
        num_params = length(args)
        timepoints = collect(1:T)
        nrows = T
        ncols = num_params + 1 
        params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        for i in 1:num_params
            params_vec[i] = vcat(args[i]...)
        end
        params_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        params_mat[:,1] = timepoints
        params_mat[:,2:end] = hcat(params_vec...)
        return params_mat
    end
    function tidify_scalar_hyperprior_params(args...)
        num_params = length(args)
        nrows = 1
        ncols = num_params + 1 
        params_vec = hcat(args...)
        params_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        params_mat[1,1] = 1
        params_mat[1,2:end] = params_vec
        return params_mat
    end
    function tidify_ηtkj(ηtkj)
        T = length(ηtkj)
        K = length(ηtkj[1])
        G = length(ηtkj[1][1])
        nrows = 2*G*K*T
        ncols = 5
        ηtkj_mat = Matrix{Union{Float64,Int64}}(undef,nrows,ncols)
        timepoints = collect(1:T)
        state_ids = collect(1:K)
        gene_ids = collect(1:G)
        choice_ids = collect(1:2)
        time = innermelt(timepoints,2*G*K)
        states = outermelt(innermelt(state_ids,2*G),T)
        genes = outermelt(outermelt(innermelt(gene_ids,2),K),T)
        importance_choice = outermelt(choice_ids,G*K*T)
        ηtkj_mat[:,1] = time
        ηtkj_mat[:,2] = states
        ηtkj_mat[:,3] = genes
        ηtkj_mat[:,4] = importance_choice
        ηtkj_mat[:,5] = recursive_flatten(ηtkj)
        return ηtkj_mat
    end
    function tidify_error_state_params(args...)
        K = 1
        G = length(args[1])
        num_params = length(args)
        state_ids = collect(1:K)
        gene_ids = collect(1:G)
        states = innermelt(state_ids,G)
        genes = outermelt(gene_ids,K)
        nrows = K*G
        ncols = num_params + 2 
        params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
        for i in 1:num_params
            params_vec[i] = vcat(args[i]...)
        end
        params_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        params_mat[:,1] = states
        params_mat[:,2] = genes
        params_mat[:,3:end] = hcat(params_vec...)
        return params_mat
    end
    function tidy_get_Nt_from_xmat(xmat)
        T = length(unique(xmat[:,1]))
        N = length(unique(xmat[:,2]))
        G = length(unique(xmat[:,3]))
        timepoint_freq = countmap(Int.(xmat[1:G:end,1]))
        N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
        return N_t
    end
    function tidy_get_Nt_from_rtikmat(rtik_mat)
        T = length(unique(rtik_mat[:,1]))
        N = length(unique(rtik_mat[:,2]))
        K = length(unique(rtik_mat[:,3]))
        timepoint_freq = countmap(Int.(rtik_mat[1:K:end,1]))
        N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
        return N_t
    end
    function tidy_get_timeranges(N_t)
        T = length(N_t)
        starts = Vector{Int}(undef,T)
        ends = Vector{Int}(undef,T)
        
        for t in 1:T
            if t == 1
                starts[1] = 0 + 1
                ends[1] = 0 + N_t[1]
                continue
            end
            starts[t] = ends[t-1] + 1
            ends[t] = ends[t-1] + N_t[t]
        end
        return zip(starts,ends)
    end
    function tidy_get_next_chain_rows_indices(chn_indx,value_mat)
        val_rows,_ = size(value_mat)
        return (chn_indx-1)*val_rows+1:chn_indx*val_rows
    end
    function tidy_make_chain(num_iter,value_mat;start = 0)
        value_mat_nrows,value_mat_ncols = size(value_mat)
        iteration_ids = collect(start:num_iter)
        nrows = value_mat_nrows*(num_iter+1-start)
        ncols = value_mat_ncols +1 
        iterations = innermelt(iteration_ids,value_mat_nrows)
        value_mat_chain_mat = Matrix{Union{Missing,Float64,Int}}(undef,nrows,ncols)
        value_mat_chain_mat[:,1] = iterations
    
        return value_mat_chain_mat
    end
    function tidy_truncate_chain(full_chain,value_mat,truncation_value)
        num_iter = length(unique(full_chain[:,1]))
        iteration_min = minimum(unique(full_chain[:,1]))
        value_mat_nrows,value_mat_ncols = size(value_mat)
        iteration_ids = collect(iteration_min:truncation_value)
        nonemptyiterations = length(iteration_ids)
        nrows = value_mat_nrows*nonemptyiterations
        trunc_chain = full_chain[1:nrows,:]
        return trunc_chain
    end

    mutable struct IndexCounter
        count
    end
    increment!(c::IndexCounter) = c.count = c.count + 1
    getCount(c::IndexCounter) = c.count
    
    
end