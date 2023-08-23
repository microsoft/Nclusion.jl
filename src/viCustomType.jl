"""
        Features
    This is an abstract type for all Features tracked by NCLUSION
"""
abstract type Features end 


"""
        CellFeatures
    This is type that allows NCLUSION to track all of the cell-specific features during inference
"""
struct CellFeatures{U <: AbstractFloat, W <:Int64,P <: Function, G} <: Features 
    t::W
    i::W
    # null_precision::U
    rtik::Vector{U}
    cache::Vector{U}
    x::NTuple{G, U}
    xsq::NTuple{G, U}
    _reset!::P
    BitType::DataType
    function CellFeatures(t,i,K,data)
        U = eltype(data)
        G = length(data)
        W = typeof(t)
        data_sq = data .^2
        function _reset!(val_vec::Vector{U},BitType::DataType) where U <: AbstractFloat
            for indx in eachindex(val_vec)
                val_vec[indx] = zero(BitType)
            end
            return val_vec
        end
        P = typeof(_reset!)
        new{U,W,P,G}(t,i,Vector{U}(undef,K),zeros(U,G),Tuple(data),Tuple(data_sq),_reset!,U)
        # new{U,W}
    end
end
"""
        ClusterFeatures
    This is type that allows NCLUSION to track all of the cluster-specific features during inference
"""
struct ClusterFeatures{U <: AbstractFloat, W <:Int64,P <: Function} <: Features #cluster_features
    k::W
    mk_hat::Vector{U}
    v_sq_k_hat::Vector{U}
    σ_sq_k_hat::Vector{U}
    var_muk::Vector{U}
    κk_hat::Vector{U}
    yjk_hat::Vector{U}
    Nk::Vector{U}#size 1 so we can mutate it
    x_hat::Vector{U}
    x_hat_sq::Vector{U}
    gk_hat::Vector{U}#size 1 so we can mutate it
    hk_hat::Vector{U}#size 1 so we can mutate it
    ak_hat::Vector{U}#size 1 so we can mutate it
    bk_hat::Vector{U}#size 1 so we can mutate it
    cache::Vector{U}# accumulation::Vector{U}
    _reset!::P
    BitType::DataType
    function ClusterFeatures(k,G;float_type=Float64)
        U = float_type
        W = typeof(k)
        function _reset!(val_vec::Vector{U},BitType::DataType) where U <: AbstractFloat
            for indx in eachindex(val_vec)
                val_vec[indx] = zero(BitType)
            end
            return val_vec
        end
        P = typeof(_reset!)
        new{U,W,P}(k,Vector{U}(undef,G),Vector{U}(undef,G),Vector{U}(undef,G),Vector{U}(undef,G),Vector{U}(undef,G),Vector{U}(undef,G),zeros(U,1),ones(U,G),ones(U,G),zeros(U,1),zeros(U,1),zeros(U,1),zeros(U,1),zeros(U,G),_reset!,U)#,zeros(U,1),zeros(U,1)
    end
end

"""
        GeneFeatures
    This is type that allows NCLUSION to track all of the gene specific features during inference
"""
struct GeneFeatures{U <: AbstractFloat, W <:Int64,P <: Function} <: Features #cluster_features
    j::W
    λ_sq::Vector{U}#size 1 so we can mutate it
    cache::Vector{U}# accumulation::Vector{U}
    _reset!::P
    BitType::DataType
    function GeneFeatures(j;float_type=Float64)
        U = float_type
        W = typeof(j)
        function _reset!(val_vec::Vector{U},BitType::DataType) where U <: AbstractFloat
            for indx in eachindex(val_vec)
                val_vec[indx] = zero(BitType)
            end
            return val_vec
        end
        P = typeof(_reset!)
        new{U,W,P}(j,Vector{U}(undef,1),zeros(U,1),_reset!,U)
    end
end

"""
        ConditionFeatures
    This is type that allows NCLUSION to track all of the condition specific features during inference
"""
struct ConditionFeatures{U <: AbstractFloat, W <:Int64} <: Features#cluster_features
    t::W
    d_hat_t_sum::Vector{U}#size 1 so we can mutate it
    d_hat_t::Vector{U}
    e_log_π_t_cache::Vector{U}
    Ntk::Vector{U}
    c_tt_prime::Vector{U}
    st_hat::Vector{U}#size 1 so we can mutate it
    BitType::DataType
    function ConditionFeatures(t,K,T;float_type=Float64)
        U = float_type
        W = typeof(t)
        Kplus = K+1
        new{U,W}(t,zeros(U,1),Vector{U}(undef,Kplus),Vector{U}(undef,Kplus),zeros(U,Kplus),Vector{U}(undef,T),zeros(U,1),U)
    end
end

"""
        DataFeatures
    This is type that allows NCLUSION to track all of other dataset-specific features during inference
"""
struct DataFeatures{U <: AbstractFloat,W <: Int64} <: Features#cluster_features
    T::W
    G::W
    N_t::Vector{W}
    N::W
    LinearAddress::Vector{Tuple{W,W}}
    TimeRanges::Vector{Tuple{W, W}}
    # GeneNames::Vector{String}
    # UniqueClusters::Vector{W}
    Glog::U
    logpi::U
    BitType::DataType
    function DataFeatures(dataset)
        U = eltype(dataset[1][1])
        T = length(dataset)
        N_t = [length(el) for el in dataset]
        N = convert(U,sum(N_t))
        G =length(dataset[1][1])
        Glog = convert(U,G*log(2π))
        logpi = convert(U,log(2π))
        W = typeof(T)
        LinearAddress = [(t,i) for t in 1:T for i in 1:N_t[t]]
        TimeRanges = get_timeranges(N_t)
        new{U,W}(T,G,N_t,N,LinearAddress,TimeRanges,Glog,logpi,U)
    end
end

"""
        ElboFeatures
    This is type that allows NCLUSION to track all the elbo during inference
"""
struct ElboFeatures{U <: AbstractFloat, W <:Int64} <: Features#cluster_features
    l::W
    elbo_::Vector{Union{Missing,U}}#size 1 so we can mutate it
    per_k_elbo::Matrix{Union{Missing,U}}
    BitType::DataType
    function ElboFeatures(l,K,num_iter;float_type=Float64)
        U = float_type
        W = typeof(l)
        new{U,W}(l,Vector{Union{Missing,Float64}}(undef,num_iter),Matrix{Union{Missing,Float64}}(undef,K,num_iter),U)
    end
end

"""
        ModelParameterFeatures
    This is type that allows NCLUSION to track all of the user defined model parameters during inference
"""
struct ModelParameterFeatures{U <: AbstractFloat,V <:AbstractFloat,W <: Int64} <: Features #cluster_features
    K::W
    ηk::Vector{U}#size 1 so we can mutate it
    α0::U
    γ0::U
    ϕ0::U
    num_iter::W
    num_local_iter::W
    uniform_theta_init::Bool
    rand_init::Bool
    BitType::DataType
    # ep::U
    # elbo_ep::U
    function ModelParameterFeatures(dataset,K,ηk,α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init)
        U = typeof(α0)
        V = eltype(dataset[1][1])
        W = typeof(K)
        new{U,V,W}(K,[ηk],α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init,V)
    end
end

"""
        get_timeranges(N_t)
    This function returns the linear indices that contain cells from the same condition.
    ```math

    ```
"""
function get_timeranges(N_t)
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
    return [(st,en) for (st,en) in zip(starts,ends)]
end

"""
        _reset!(val_vec::Vector{U},BitType::DataType)
    Performs an inplace setting of values in a vector to 0. Maintains the type of the variable prior to reset.
"""
function _reset!(val_vec::Vector{U},BitType::DataType) where U <: AbstractFloat
    for indx in eachindex(val_vec)
        val_vec[indx] = zero(BitType)
    end
    return val_vec
end