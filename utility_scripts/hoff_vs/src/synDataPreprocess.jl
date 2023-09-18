function getCustomersAtTimepoint(data_dict,tp)
    ka_ = collect(keys(data_dict))
    key_array_ = ka_[last.(ka_) .== tp]
    key_array_ = sort(key_array_, by = first)
    _tp = [(first(key),data_dict[key]) for key in key_array_]
    return _tp
end
# TODO ADD TESTS
# getCustomersAtTimepoint(data_dict,2)
# getCustomersAtTimepoint(data_dict,1)

function getTrueClusterMembershipDict(T,k,truth_dict)
    true_cluster_membership_dict = Dict(t => Dict(i => Set([]) for i in 1:k ) for t in 1:T)
    for t in 1:T
        cells = getCustomersAtTimepoint(truth_dict,t)
        for c in cells
            true_cluster_id = last(c)
            push!(true_cluster_membership_dict[t][true_cluster_id], first(c))
        end 
    end
    return true_cluster_membership_dict
end

function genData_vec2dict(time_vec)
    T = length(time_vec)
    C_t = [length(c) for c in time_vec]
    final_dict = Dict()
    for t in 1:T
        for c in 1:C_t[t]
            final_dict[(c,t)] = time_vec[t][c]
        end
    end
    # final_dict = Dict((c,t) => time_vec[t][c] for c in 1:C_t[t] for t in 1:T)
    return final_dict
end
# module syntheticDataPreprocessing

#     using Random
#     using Distributions
#     using Turing
#     using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     using Test
#     import Debugger
#     using CSV,DataFrames
    
#     export getCustomersAtTimepoint
#     function getCustomersAtTimepoint(data_dict,tp)
#         ka_ = collect(keys(data_dict))
#         key_array_ = ka_[last.(ka_) .== tp]
#         key_array_ = sort(key_array_, by = first)
#         _tp = [(first(key),data_dict[key]) for key in key_array_]
#         return _tp
#     end
#     # TODO ADD TESTS
#     # getCustomersAtTimepoint(data_dict,2)
#     # getCustomersAtTimepoint(data_dict,1)
#     export getTrueClusterMembershipDict, genData_vec2dict
#     function getTrueClusterMembershipDict(T,k,truth_dict)
#         true_cluster_membership_dict = Dict(t => Dict(i => Set([]) for i in 1:k ) for t in 1:T)
#         for t in 1:T
#             cells = getCustomersAtTimepoint(truth_dict,t)
#             for c in cells
#                 true_cluster_id = last(c)
#                 push!(true_cluster_membership_dict[t][true_cluster_id], first(c))
#             end 
#         end
#         return true_cluster_membership_dict
#     end

#     function genData_vec2dict(time_vec)
#         T = length(time_vec)
#         C_t = [length(c) for c in time_vec]
#         final_dict = Dict()
#         for t in 1:T
#             for c in 1:C_t[t]
#                 final_dict[(c,t)] = time_vec[t][c]
#             end
#         end
#         # final_dict = Dict((c,t) => time_vec[t][c] for c in 1:C_t[t] for t in 1:T)
#         return final_dict
#     end
# end

