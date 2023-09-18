



function getEmpricalCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
    empirical_cluster_membership_exclusive = Dict(t => Dict() for t in 1:T) #Exactly correct, no other cells included in the cluster
    empirical_cluster_membership_inclusive = Dict(t => Dict() for t in 1:T) #Correct, but other cells are in the cluster
    numCellsPerT = [length(getCustomersAtTimepoint(truth_dict,t)) for t in 1:T]
    true_cluster_membership_dict = getTrueClusterMembershipDict(T,k,truth_dict)
    z_df = getzDF(tchain,n_samples)
    for r in eachrow(z_df)
        for t in 1:T
            numCellsT = numCellsPerT[t]
            membership_dict_t = true_cluster_membership_dict[t]
            tp_col = [m.match for m in match.(time_re_func(t,"z"),names(r)) if !isnothing(m) ]
            unique_emp_cluster_assgin = unique(collect(r[tp_col]))
            
            temp_dict = Dict(clus => Set([]) for clus in unique_emp_cluster_assgin)
            
            tp_cells = r[tp_col]
            #Assume linear indexing of cell ids!!!!
            for c in 1:numCellsT
                cell_col = [m.match for m in match.(cellID_re_func(c,"z"),names(tp_cells)) if !isnothing(m) ]
                cell_cluster_assgn = collect(tp_cells[cell_col])[1]
                push!(temp_dict[cell_cluster_assgn],c)
            end
            all_empirical_clusters = [temp_dict[key] for key in keys(temp_dict)]
            for key in keys(membership_dict_t)
                true_cluster_membership = membership_dict_t[key]
                inclusive_check = [issubset(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
                exlusive_check = [issetequal(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
                if any(inclusive_check)
                    if haskey(empirical_cluster_membership_inclusive[t],true_cluster_membership)
                        empirical_cluster_membership_inclusive[t][true_cluster_membership] += 1
                    else
                        empirical_cluster_membership_inclusive[t][true_cluster_membership] = 1
                    end 
                end
                if any(exlusive_check)
                    if haskey(empirical_cluster_membership_exclusive[t],true_cluster_membership)
                        empirical_cluster_membership_exclusive[t][true_cluster_membership] +=1
                    else
                        empirical_cluster_membership_exclusive[t][true_cluster_membership] = 1
                    end
                end
            end
            
        end
    end
    emp_clus_mem_counts_inclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
    emp_clus_mem_counts_exclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
    for t in 1:T
        for clus in 1:k
            key = true_cluster_membership_dict[t][clus]
            if haskey(empirical_cluster_membership_exclusive[t],key)
                emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
            end
            # emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
            if haskey(empirical_cluster_membership_inclusive[t],key)
                emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
            end
            # emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
        end
    end
    exclusive_rate = [[emp_clus_mem_counts_exclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
    inclusive_rate = [[emp_clus_mem_counts_inclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
    return exclusive_rate, inclusive_rate
end

function getzIDs(tchain)
    z_ids = findall(map(name -> occursin("z", string(name)), names(tchain)));
    return z_ids
end
function getzDF(tchain,n_samples)
    z_ids = getzIDs(tchain)
    z_df = DataFrame(tchain[n_samples:end, z_ids, :])[:,3:end]
    return z_df
end


function getlambdaposteriorAvg(tchain,n_samples)
    λ_df = getlambdaposteriorDF(tchain,n_samples)
    cluster_means = [mean(col) for col in  eachcol(λ_df)]
    cluster_std = [std(col) for col in  eachcol(λ_df)]
    cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
    return cluster_mat
end
function getlambdaposteriorDF(tchain,n_samples)
    λ_ids = getlambdaIDs(tchain)
    λ_df = DataFrame(tchain[n_samples:end, λ_ids, :])[:,3:end]
    return λ_df
end
function getlambdaIDs(tchain)
    λ_ids = findall(map(name -> occursin("λ", string(name)), names(tchain)))
    return λ_ids
end


function getmuposteriorAvg(tchain,n_samples)
    m_df = getmuposteriorDF(tchain,n_samples)
    cluster_means = [mean(col) for col in  eachcol(m_df)]
    cluster_std = [std(col) for col in  eachcol(m_df)]
    cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
    return cluster_mat
end
function getmuposteriorDF(tchain,n_samples)
    m_ids = getmuIDs(tchain)
    m_df = DataFrame(tchain[n_samples:end, m_ids, :])[:,3:end]
    return m_df
end
function getmuIDs(tchain)
    m_ids = findall(map(name -> occursin("m", string(name)), names(tchain)))
    return m_ids
end

function getsigmaposteriorAvg(tchain,n_samples)
    s_df = getsigmaposteriorDF(tchain,n_samples)
    cluster_means = [mean(col) for col in  eachcol(s_df)]
    cluster_std = [std(col) for col in  eachcol(s_df)]
    cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
    return cluster_mat
end
function getsigmaposteriorDF(tchain,n_samples)
    s_ids = getsigmaIDs(tchain)
    s_df = DataFrame(tchain[n_samples:end, s_ids, :])[:,3:end]
    return s_df
end
function getsigmaIDs(tchain)
    s_ids = findall(map(name -> occursin("s²", string(name)), names(tchain)))
    return s_ids
end



 
function generateDataDF_from_TuringChain(chn)
    return DataFrame(chn)
end
function get_var_ids(chain, var_)
    # var_ids = findall(map(name -> occursin(var_, string(name)), names(chain)))
    var_ids = MCMCChains.namesingroup(chain, Symbol(var_))
    return var_ids
end
function get_chain_var(chain)
    return string.(collect(keys(get(chain; section=:parameters))))
end
function get_static_chn_param_univariate(var_, chn)
    _id = get_var_ids(chn, var_)
    _chn = chn[:, _id, :]
    _df = generateDataDF_from_TuringChain(_chn)
    col_names_ = [m.match for m in match.(general_re_func(var_),names(_df)) if !isnothing(m) ]
    _true =  [_df[:,c][1] for c in col_names_]
    return _true
end
function get_dynamic_chn_param_univariate(var_, chn, T)
    _id = get_var_ids(chn, var_)
    _chn = chn[:, _id, :]
    _df = generateDataDF_from_TuringChain(_chn)
    col_names_ = [[m.match for m in match.(time_re_func(t,var_),names(_df)) if !isnothing(m) ] for t in 1:T]
    _true =  [[ col for col in eachcol(_df[!,c])] for c in col_names_]
    return _true
end
function unpack_data_chn_param(_true)
    return [[c[1] for c in el ] for el in _true]
end






 
function get_πposterior(chain, num_burnin, T)
    v_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, :v), :])[!,3:end]
    v_t_vec = [[[r[j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(r)) if !isnothing(m) ]] for t in 1:T] for r in eachrow(v_df)]
    π_t_post = [stickbreak.(itr) for itr in v_t_vec]
    return π_t_post
end
function get_Poisson_cell_lhprob(chain,x, num_burnin, T, C_t)
    λdf =  getlambdaposteriorDF(chain,num_burnin)
    param_lh = [Poisson.(collect(r)) for r in eachrow(λdf)]
    cell_lhprob = [[[pdf.(itr_lh,x[t][c]) for c in 1:C_t[t]] for t in 1:T] for itr_lh in param_lh]
    return cell_lhprob
end
function get_Normal_cell_lhprob(chain,x, num_burnin, T, C_t)
    mdf =  getmuposteriorDF(chain,num_burnin)
    sdf =  getsigmaposteriorDF(chain,num_burnin)
    ms = [collect(r) for r in eachrow(mdf)]
    ss = [collect(r) for r in eachrow(sdf)]
    param_lh = [Normal.(ms[i],ss[i]) for i in 1:length(ms)]
    cell_lhprob = [[[pdf.(itr_lh,x[t][c]) for c in 1:C_t[t]] for t in 1:T] for itr_lh in param_lh]
    return cell_lhprob
end
function get_P_tensor(S,T,C_t,π_t_post,cell_lhprob; addPseudoCount = true, pseudoCount = 1)
    if  addPseudoCount
        pseudoCount = pseudoCount
    else
        pseudoCount = 0
    end
    P = Vector{Vector{Vector{Vector{Float64}}}}(undef, S)
    for s in 1:S
        P[s] = [[ (π_t_post[s][t] .* cell_lhprob[s][t][c] .+ pseudoCount) ./ sum((π_t_post[s][t] .* cell_lhprob[s][t][c] .+ pseudoCount)) for c in 1:C_t[t]] for t in 1:T ]
    end
    return P
end
function stephens_relabelling(P,S,T,C_t,KMax,num_itr;use_identity_perm = true,permute_rows = false)
    relabelling_results = Dict()
    relabelling_results["nu"] = Vector{Vector{Vector{Int64}}}(undef,num_itr)
    relabelling_results["cost"] = Vector{Vector{Float64}}(undef,num_itr)
    nu = init_nu(S,KMax;use_identity_perm = use_identity_perm)
    for i in 1:num_itr
        Q = calculate_Q_hat(P,S,T,C_t,KMax,nu)
        cost_per_t = calculate_Cost_per_t(P,Q,S,T,C_t,KMax,nu; permute_rows = permute_rows)#calculate_Cost_per_t(P,Q,S,T,C_t,KMax)
        cost_across_t = calculate_Cost_across_t(cost_per_t,S,T,KMax)
        CostMat_across_t = create_CostMat_across_t(cost_across_t,S,KMax)
        updatesDict = update_nu(CostMat_across_t,S)
        nu = updatesDict["nu"]
        cost = updatesDict["cost"]
        relabelling_results["nu"][i] = nu
        relabelling_results["cost"][i] = cost
    end
    return relabelling_results
    
end
function update_nu(CostMat_across_t,S)
    updatesDict = Dict()
    updatesDict["nu"] = Vector{Vector{Int64}}(undef,S)
    updatesDict["cost"] = Vector{Float64}(undef,S)
    for s in 1:S
        assignment, cost = hungarian(CostMat_across_t[s])
        updatesDict["nu"][s] = assignment
        updatesDict["cost"][s] = cost
    end
    return updatesDict
    
end
function init_nu(S,KMax;use_identity_perm = true)
    nu =  Vector{Vector{Int}}(undef,S)
    for s in 1:S
        if use_identity_perm
            _perm = collect(1:KMax)
        else
            _perm = randperm(KMax)
        end
        nu[s] = _perm#randperm(KMax)#
    end
    return nu
end
function calculate_Q_hat(P,S,T,C_t,KMax,nu)
    Q = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        num_cells = C_t[t]
        Q[t] = Vector{Vector{Float64}}(undef,num_cells)
        for c in 1:num_cells
            Q[t][c] = Vector{Float64}(undef,KMax)
            for k in 1:KMax
                # nu_k = nu[s][k]
                all_vals = [P[s][t][c][nu[s][k]] for s  in 1:S]
                Q[t][c][k] = mean(all_vals)
            end
        end
    end
    return Q
end
function calculate_Cost_per_t(P,Q,S,T,C_t,KMax,nu; permute_rows = false)
    cost_per_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
    for s in 1:S
        cost_per_t[s] = Vector{Vector{Vector{Float64}}}(undef,T)
        for t in 1:T
            cost_per_t[s][t] = Vector{Vector{Float64}}(undef,KMax)
            for j in 1:KMax
                cost_per_t[s][t][j] = Vector{Float64}(undef,KMax)
                for l in 1:KMax
                    _q_hat(c)  = Q[t][c][j]
                    # _p(c) = P[s][t][c][nu[s][l]]
                    _p(c) = P[s][t][c][l]
                    _kl(c) = _p(c)*(log(_p(c))-log(_q_hat(c)))
                    kl_cells = [_kl(c) for c in 1:C_t[t]]
                    cost_per_t[s][t][j][l] = sum( kl_cells)
                end
            end
            if permute_rows
                old_ref = deepcopy(cost_per_t[s][t])
                for j in 1:KMax
                    permuted_idx = nu[s][j]
                    cost_per_t[s][t][permuted_idx] = old_ref[j]
                end
            end
        end
    end
    return cost_per_t
end
function calculate_Cost_across_t(cost_per_t,S,T,KMax)
    cost_across_t = Vector{Vector{Vector{Float64}}}(undef,S)
    for s in 1:S
        cost_across_t[s] = Vector{Vector{Float64}}(undef,KMax)
        for j in 1:KMax
            cost_across_t[s][j] = Vector{Float64}(undef,KMax)
            for l in 1:KMax
                cost_across_t[s][j][l] = sum([cost_per_t[s][t][j][l] for t in 1:T])
            end
        end
    end
    return cost_across_t
end
function create_CostMat_across_t(cost_across_t,S,KMax)
    CostMat_across_t =Vector{Matrix{Float64}}(undef,S)
    for s in 1:S
        CostMat_across_t[s] = -188.0*ones(Float64,KMax,KMax)
        for l in 1:KMax
            for j in 1:KMax
                CostMat_across_t[s][j, l] = cost_across_t[s][j][l]
            end
        end
    end
    return CostMat_across_t
end
function get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    S = length(tchain)
    z_ids =  get_var_ids(tchain, "z")
    z_post_s = [unpack_data_chn_param(get_dynamic_chn_param_univariate("z", tchain[s,z_ids,:], T)) for s in 1:S]
    post_z_dict_s = genData_vec2dict.(z_post_s)
    true_z_dict = genData_vec2dict(true_z)
    counts_mat_s = zeros(Int,KMax,KTrue,T,S)
    for s in 1:S
        for t in 1:T
            for c in 1:C_t[t]
                row = Int(post_z_dict_s[s][c,t])
                col = Int(true_z_dict[c,t])
                counts_mat_s[row,col,t,s] += 1
            end
        end
    end
    avg_counts_mat = mean(counts_mat_s,dims=4)
    return avg_counts_mat
end
function get_average_posterior_cluster_frequency2(z_post_s,T,true_z,KMax,KTrue,C_t)
    S = length(z_post_s)
    post_z_dict_s = genData_vec2dict.(z_post_s)
    true_z_dict = genData_vec2dict(true_z)
    counts_mat_s = zeros(Int,KMax,KTrue,T,S)
    for s in 1:S
        for t in 1:T
            for c in 1:C_t[t]
                row = Int(post_z_dict_s[s][c,t])
                col = Int(true_z_dict[c,t])
                counts_mat_s[row,col,t,s] += 1
            end
        end
    end
    avg_counts_mat = mean(counts_mat_s,dims=4)
    return avg_counts_mat
end
function get_clus_w_maxPopAtT(avg_counts_mat,t)
    maxPopAtT = findmax(avg_counts_mat[:,:,t,1];dims=2)
    ClusIndx_w_maxPopAtT = maxPopAtT[2]
    clus_w_maxPopAtT = ["True Cluster "*string(k[2]) for k in ClusIndx_w_maxPopAtT]
    return vec(clus_w_maxPopAtT)
end
function FindNonunique_func(x)
    return  [k for (k, v) in countmap(x) if v > 1] # needs StatsBase
end





function get_chn_values_dict(chain,T,num_burnin,π_t_post,chn_var,dymamic_var_dict)
    
    chn_values_dict = Dict()
    for chr in chn_var
        if chr == "π_t"
            chn_values_dict[chr] = π_t_post
        else
            key = Symbol(chr)
            var_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, key), :])[!,3:end]
            if typeof(dymamic_var_dict[chr]) <: Bool
                if dymamic_var_dict[chr]
                    gg(t,chr,var_df) = [m.match for m in match.(_position1_re_func(t,chr),names(var_df)) if !isnothing(m) ]
                    var_vov = [[[r[c] for c in gg(t,chr,var_df)] for t in 1:T] for r in eachrow(var_df)]
                else
                    ff(chr) = [m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
                    var_vov = [[[r[c] for c in ff(chr)] for t in 1:T] for r in eachrow(var_df)]
                end
            else
                if dymamic_var_dict[chr] == "cluster"
                    hh(chr) = [m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
                    var_vov = [[[r[c] for c in hh(chr)]] for r in eachrow(var_df)]
                end
            end
            chn_values_dict[chr] = var_vov
        end
    end
    return chn_values_dict
end
function relabel_chain_variable_dict(relabel_chn_values_dict,S,T,KMax,chn_var,relabel_var_dict, relabel_perm)
    for chr in chn_var
        if typeof(relabel_var_dict[chr]) <: AbstractString
            for s in 1:S
                for t in 1:T
                    if relabel_var_dict[chr] == "cluster"
                        new_indx = relabel_perm[s]
                        relabel_chn_values_dict[chr][s][t] =  relabel_chn_values_dict[chr][s][t][new_indx]
                    elseif relabel_var_dict[chr] == "cluster-param"
                        new_indx = relabel_perm[s]
                        relabel_chn_values_dict[chr][s][1] =  relabel_chn_values_dict[chr][s][1][new_indx]
                    elseif relabel_var_dict[chr] == "assignment"
                        itr_relabel_dict = Dict(j => relabel_perm[s][j] for j in 1:KMax)
                        relabel_chn_values_dict[chr][s][t] = [itr_relabel_dict[el] for el in relabel_chn_values_dict[chr][s][t]]
                    end
                end
            end
        end
    end
    return relabel_chn_values_dict
end


function make_θ_t_var_names(T)
    var_ = "θ_t"
    return [var_*"[$(t)]" for t in 1:T]
end
function make_π_t_var_names(T,KMax)
    var_ = "π_t"
    return [var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
end
function make_v_var_names(T,KMax)
    var_ = "v"
    return [var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax-1 ]
end
function make_z_var_names(T,C_t)
    var_ = "z"
    return vcat([[var_*"[$(t)][$(c)]" for c in 1:C_t[t]] for t in 1:T]...)
end
function make_λ_var_names(KMax)
    var_ = "λ"
    return [var_*"[$(k)]" for k in 1:KMax]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
end
function make_m_var_names(KMax)
    var_ = "m"
    return [var_*"[$(k)]" for k in 1:KMax]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
end
function make_s_var_names(KMax)
    var_ = "s"
    return  [var_*"[$(k)]" for k in 1:KMax]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
end

 
function make_θ_t_vov2Mat(relabel_chn_values_dict)
    var_ = "θ_t"
    vov = relabel_chn_values_dict[var_]
    
    return permutedims(reduce(hcat,reduce.(vcat, vov) ))#reduce(vcat,reduce.(vcat,vov))
end
function make_π_t_vov2Mat(relabel_chn_values_dict)
    var_ = "π_t"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_v_vov2Mat(relabel_chn_values_dict)
    var_ = "v"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_z_vov2Mat(relabel_chn_values_dict)
    var_ = "z"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_λ_vov2Mat(relabel_chn_values_dict)
    var_ = "λ"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_m_vov2Mat(relabel_chn_values_dict)
    var_ = "m"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_s_vov2Mat(relabel_chn_values_dict)
    var_ = "s"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
##### GENERAL EXCEPT FOR θ_t ##########
function make_vectorParam_vov2Mat(param,relabel_chn_values_dict)
    var_ = param
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end






function partition_matrix_cluster_df_names(df,K)
    [ [m.match for m in match.(_position2_Matrix_re_func(t,"λ"),names(df)) if !isnothing(m) ] for t in 1:K]
end
function partition_var_matrix_cluster_df_names(df,K,var_)
    [ [m.match for m in match.(_position2_Matrix_re_func(t,var_),names(df)) if !isnothing(m) ] for t in 1:K]
end
function partition_matrix_cluster_df_names2(df,K)
    [ [m.match for m in match.(_position2_Matrix_re_func2(t,"λ"),names(df)) if !isnothing(m) ] for t in 1:K]
end
function _position1_Matrix_re_func(g,var)
    return Regex(string(var)*"\\["*string(g)*"(,\\[[0-9]+\\])*\\]")
end
function _position2_Matrix_re_func(g,var)
    return Regex(string(var)*"\\[[0-9]+,"*string(g)*"\\](\\[[0-9]+\\])*")
end
function _position2_Matrix_re_func2(g,var)
    return Regex(string(var)*"\\[([0-9]+,)*"*string(g)*"\\](\\[[0-9]+\\])*")
end
function π_t_post_vov2Mat(π_t_post)
    return permutedims(reduce(hcat,reduce.(vcat,π_t_post)))
end
function _get_πposterior(chain, num_burnin, T;islogProb= false)
    v_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, :v), :])[!,3:end]
    v_t_vec = [[[r[j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(r)) if !isnothing(m) ]] for t in 1:T] for r in eachrow(v_df)]
    π_t_post = [vec(map(el-> islogProb ? log.(stickbreak(el)) :  stickbreak(el) ,itr)) for itr in v_t_vec] #[stickbreak.(itr) for itr in v_t_vec]
    return π_t_post
end
function _get_HDP_πposterior(chain, num_burnin, T;islogProb= false)
    # DataFrame(tchain[thin_n_burin_samples:end, MCMCChains.namesingroup(tchain, :π_), :])[!,3:end]
    π_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, :π_), :])[!,3:end]
    π_t_vec = [[[r[j][1] for j in [m.match for m in match.(time_re_func(t,"π_"),names(r)) if !isnothing(m) ]] for t in 1:T] for r in eachrow(π_df)]
    π_t_post = [vec(map(el-> islogProb ? log.(el) :  el ,itr)) for itr in π_t_vec] #[stickbreak.(itr) for itr in v_t_vec]
    return π_t_post
end

function partition_cluster_df_names(df,K)
    [ [m.match for m in match.(_position1_re_func(t,"λ"),names(df)) if !isnothing(m) ] for t in 1:K]
end


######### I THINK SOMETHING IS UP HERE OR I NEED TO RUN IT LONGER(only ran for 200 itr for testing purposes)!!! #################
function get_cell_lhprob_mat(chain,x, lh_func, num_burnin, T, C_t,KMax;islog= false)
    λdf =  getlambdaposteriorDF(chain,num_burnin)
    S = size(λdf)[1]
    λdf_names = partition_matrix_cluster_df_names(λdf,KMax)
    get_row_parameterize(df,row , df_names) = collect(df[row,df_names])
    param_lh = Vector{Vector}(undef,S)
    if islog
        pdf_func = logpdf
    else
        pdf_func = pdf
    end
    for s in 1:S
        row_params = vec(map(el -> get_row_parameterize(λdf,s, el), λdf_names))
        param_lh[s] = map(el -> lh_func(el),row_params)
    end
    # param_lh = [Poisson.(collect(r)) for r in eachrow(λdf)]
    cell_lhprob = [[[ 
                    map(clus_lh -> pdf_func(clus_lh,x[t][c]),param_lh[s])  
                        for c in 1:C_t[t]]  
                            for t in 1:T]  
                                for s in 1:S] #[[[ vcat(map(el ->  pdf_func(itr_lh[el],x[t][c]), 1:KMax )...) for c in 1:C_t[t]] for t in 1:T] for itr_lh in param_lh]
    return cell_lhprob
end

function get_post_cluster_assgn(tchain,x,num_burnin,T, C_t)
    v_ids = MCMCChains.namesingroup(tchain, :v)
    vpost = Matrix(DataFrame(tchain[:,v_ids,:])[!,3:end])
    πpost = get_πposterior(tchain,num_burnin, T);
    S = size(πpost)[1]
    cell_lhprob = get_Poisson_cell_lhprob(tchain,x, num_burnin, T, C_t);
    post_ip = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
    zpost = Vector{Vector{Vector{Int}}}(undef,S)
    for s in 1:S 
        post_ip[s] = Vector{Vector{Vector{Float64}}}(undef,T)
        zpost[s] = Vector{Vector{Float64}}(undef,T)
        colend = 0
        for t in 1:T
            cells = C_t[t]
            post_ip[s][t] = Vector{Vector{Float64}}(undef,cells)
            # zpost[s][t] = zeros{Float64}(undef,cells)
            for c in 1:cells
                val_vec = πpost[s][t] .* cell_lhprob[s][t][c]
                val_sum = sum(val_vec)
                post_ip[s][t][c] = val_vec ./ val_sum
            end
            zpost[s][t] = rand.(Categorical.(post_ip[s][t]))
        end
    end
    
    return zpost,post_ip
    
end

####USE THIS ONE!!!!
function get_post_cluster_assgnMat2(tchain,x,num_burnin,T, C_t,KMax,lh_func;uselogProb= false)
    v_ids = MCMCChains.namesingroup(tchain, :v)
    vpost = Matrix(DataFrame(tchain[:,v_ids,:])[!,3:end])
    πpost = _get_πposterior(tchain, num_burnin, T;islogProb= uselogProb); #get_πposterior(tchain,num_burnin, T);
    S = size(πpost)[1]
    cell_lhprob = get_cell_lhprob_mat(tchain,x, lh_func, num_burnin, T, C_t,KMax;islog= uselogProb); #get_Poisson_cell_lhprob(tchain,x, num_burnin, T, C_t);
    post_ip = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
    zpost = Vector{Vector{Vector{Int}}}(undef,S)
    for s in 1:S 
        post_ip[s] = Vector{Vector{Vector{Float64}}}(undef,T)
        zpost[s] = Vector{Vector{Float64}}(undef,T)
        colend = 0
        for t in 1:T
            cells = C_t[t]
            post_ip[s][t] = Vector{Vector{Float64}}(undef,cells)
            # zpost[s][t] = zeros{Float64}(undef,cells)
            for c in 1:cells
                if uselogProb
                    val_vec = πpost[s][t] .+ cell_lhprob[s][t][c]
                    val_sum = StatsFuns.logsumexp(val_vec)
                    post_ip[s][t][c] = exp.(val_vec .- val_sum)
                else
                    val_vec = πpost[s][t] .* cell_lhprob[s][t][c]
                    val_sum = sum(val_vec)
                    post_ip[s][t][c] = val_vec ./ val_sum
                end
                
            end
            zpost[s][t] = rand.(Categorical.(post_ip[s][t]))
        end
    end
    
    return zpost,post_ip
    
end
function get_post_cluster_assgnMat3(tchain,x,num_burnin,T, C_t,KMax,lh_func;uselogProb= false,πpost = nothing)
    if isnothing(πpost)
        πpost = _get_HDP_πposterior(tchain, thin_n_burin_samples, T;islogProb= uselogProb); 
    end
    S = size(πpost)[1]
    cell_lhprob = get_cell_lhprob_mat(tchain,x, lh_func, num_burnin, T, C_t,KMax;islog= uselogProb); #get_Poisson_cell_lhprob(tchain,x, num_burnin, T, C_t);
    post_ip = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
    zpost = Vector{Vector{Vector{Int}}}(undef,S)
    for s in 1:S 
        post_ip[s] = Vector{Vector{Vector{Float64}}}(undef,T)
        zpost[s] = Vector{Vector{Float64}}(undef,T)
        colend = 0
        for t in 1:T
            cells = C_t[t]
            post_ip[s][t] = Vector{Vector{Float64}}(undef,cells)
            # zpost[s][t] = zeros{Float64}(undef,cells)
            for c in 1:cells
                if uselogProb
                    val_vec = πpost[s][t] .+ cell_lhprob[s][t][c]
                    val_sum = StatsFuns.logsumexp(val_vec)
                    post_ip[s][t][c] = exp.(val_vec .- val_sum)
                else
                    val_vec = πpost[s][t] .* cell_lhprob[s][t][c]
                    val_sum = sum(val_vec)
                    post_ip[s][t][c] = val_vec ./ val_sum
                end
                
            end
            zpost[s][t] = rand.(Categorical.(post_ip[s][t]))
        end
    end
    
    return zpost,post_ip
    
end

function _get_chnMatrix_values_dict(chain,T,KMax,num_burnin,π_t_post,post_ip,zpost,chn_var,dymamic_var_dict)
        
    chn_values_dict = Dict()
    for chr in chn_var
        if chr == "π_t"
            chn_values_dict[chr] = π_t_post
        elseif  chr == "PIP_t"
            chn_values_dict[chr] = post_ip
        elseif  chr == "z_t"
            chn_values_dict[chr] = zpost
        else
            key = Symbol(chr)
            var_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, key), :])[!,3:end]
            if typeof(dymamic_var_dict[chr]) <: Bool
                if dymamic_var_dict[chr]
                    gg(t,chr,var_df) = [m.match for m in match.(_position1_re_func(t,chr),names(var_df)) if !isnothing(m) ]
                    var_vov = [[[r[c] for c in gg(t,chr,var_df)] for t in 1:T] for r in eachrow(var_df)]
                else
                    ff(chr) = [m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
                    var_vov = [[[r[c] for c in ff(chr)] for t in 1:T] for r in eachrow(var_df)]
                end
            else
                if dymamic_var_dict[chr] == "cluster"
                    hh(chr) =  [ [m.match for m in match.(_position2_Matrix_re_func(t,chr),names(var_df)) if !isnothing(m) ] for t in 1:KMax]#[m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
                    var_vov = [[[collect(r[c])[1:end] for c in hh(chr)]] for r in eachrow(var_df)]
                end
            end
            chn_values_dict[chr] = var_vov
        end
    end
    return chn_values_dict
end
function _relabel_chain_variable_dict(relabel_chn_values_dict,S,T,KMax,C_t,chn_var,relabel_var_dict, relabel_perm)
    for chr in chn_var
        if typeof(relabel_var_dict[chr]) <: AbstractString
            for s in 1:S
                for t in 1:T
                    if relabel_var_dict[chr] == "cluster"
                        new_indx = relabel_perm[s]
                        relabel_chn_values_dict[chr][s][t] =  relabel_chn_values_dict[chr][s][t][new_indx]
                    # elseif relabel_var_dict[chr] == "cluster-param"
                    #     println("here son")
                    #     new_indx = relabel_perm[s]
                    #     relabel_chn_values_dict[chr][s][1] =  relabel_chn_values_dict[chr][s][1][new_indx]
                    elseif relabel_var_dict[chr] == "PIP"
                        new_indx = relabel_perm[s]
                        for c in 1:C_t[t]
                            relabel_chn_values_dict[chr][s][t][c] =  relabel_chn_values_dict[chr][s][t][c][new_indx] 
                        end
                    elseif relabel_var_dict[chr] == "assignment"
                        itr_relabel_dict = Dict(j => relabel_perm[s][j] for j in 1:KMax)
                        relabel_chn_values_dict[chr][s][t] = [itr_relabel_dict[el] for el in relabel_chn_values_dict[chr][s][t]]
                    end
                end
            if relabel_var_dict[chr] == "cluster-param"
                
                new_indx = relabel_perm[s]
                relabel_chn_values_dict[chr][s][1] =  relabel_chn_values_dict[chr][s][1][new_indx]
            end
            end
        end
    end
    return relabel_chn_values_dict
end

function make_PIP_t_var_names(T,C_t,KMax)
    var_ = "PIP_t"
    return [var_*"[$(t)][$(c)][$(k)]" for t in 1:T for c in 1:C_t[t] for k in 1:KMax]
end
function make_z_t_vov2Mat(relabel_chn_values_dict)
    var_ = "z_t"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_PIP_t_vov2Mat(relabel_chn_values_dict)
    var_ = "PIP_t"
    vov = relabel_chn_values_dict[var_]
    return permutedims(reduce(hcat , reduce.(vcat,reduce.(vcat,vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_λmat_vov2Mat(relabel_chn_values_dict)
    var_ = "λ"
    vov = relabel_chn_values_dict[var_]
    vov2 = reduce.(vcat,(reduce.(vcat, vov))) #permutedims(reduce.(vcat,(reduce.(vcat, vov))))
    return  permutedims(reduce(hcat, vov2))#permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
end
function make_λmat_var_names(KMax,G)
    var_ = "λ"
    return [var_*"[$(g),$(k)]"  for k in 1:KMax for g in 1:G]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
end

function getLambdaInferenceModel2_newNames(G,KMax)
    rename_dict = Dict()
    var = "λ"
    for g in 1:G
        for k in 1:KMax
            old_name = var*"[$k][$g]"
            new_name = var*"[$g,$k]"
            rename_dict[old_name] = new_name
        end
    end
    return rename_dict
end

function extract(tchain, var; burn=0, is_scalar=true)
    sym = Symbol(var)
    if is_scalar
        tail = tchain[sym].data[(burn + 1):end, 1]
    else
        tail = group(tchain, sym).value.data[(burn + 1):end, :, 1]
        # NOTE: `chain[:sym]` will not work in the future for parameter
        # vectors. Use `group` instead.
    end
    return tail
end


function re_func(t,var)
    return Regex(string(var)*"\\["*string(t)*"\\]\\[[0-9]+\\]")
end
function time_re_func(t,var)
    return Regex(string(var)*"\\["*string(t)*"\\]\\[[0-9]+\\]")
end
function cellID_re_func(id,var)
    return Regex(string(var)*"\\[[0-9]+\\]\\["*string(id)*"\\]")
end

function _position1_re_func(t,var)
    return Regex(string(var)*"\\["*string(t)*"\\](\\[[0-9]+\\])*")
end
function _position2_re_func(id,var)
    return Regex(string(var)*"\\[[0-9]+\\]\\["*string(id)*"\\](\\[[0-9]+\\])*")
end
function _position3_re_func(g,var)
    return Regex(string(var)*"\\[[0-9]+\\]\\[[0-9]+\\]\\["*string(g)*"\\](\\[[0-9]+\\])*")
end
function general_re_func(var)
    return  Regex(string(var)*"(\\[[0-9]+\\])+")
end

# module turingChainProcessing
#     include("synDataPreprocess.jl")
#     using .syntheticDataPreprocessing
#     using Random
#     using Distributions
#     using Turing
#     using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     using Test
#     import Debugger
#     using CSV,DataFrames
#     using Hungarian


#     export getEmpricalCorrectClusteringRates
#     function getEmpricalCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
#         empirical_cluster_membership_exclusive = Dict(t => Dict() for t in 1:T) #Exactly correct, no other cells included in the cluster
#         empirical_cluster_membership_inclusive = Dict(t => Dict() for t in 1:T) #Correct, but other cells are in the cluster
#         numCellsPerT = [length(getCustomersAtTimepoint(truth_dict,t)) for t in 1:T]
#         true_cluster_membership_dict = getTrueClusterMembershipDict(T,k,truth_dict)
#         z_df = getzDF(tchain,n_samples)
#         for r in eachrow(z_df)
#             for t in 1:T
#                 numCellsT = numCellsPerT[t]
#                 membership_dict_t = true_cluster_membership_dict[t]
#                 tp_col = [m.match for m in match.(time_re_func(t,"z"),names(r)) if !isnothing(m) ]
#                 unique_emp_cluster_assgin = unique(collect(r[tp_col]))
                
#                 temp_dict = Dict(clus => Set([]) for clus in unique_emp_cluster_assgin)
                
#                 tp_cells = r[tp_col]
#                 #Assume linear indexing of cell ids!!!!
#                 for c in 1:numCellsT
#                     cell_col = [m.match for m in match.(cellID_re_func(c,"z"),names(tp_cells)) if !isnothing(m) ]
#                     cell_cluster_assgn = collect(tp_cells[cell_col])[1]
#                     push!(temp_dict[cell_cluster_assgn],c)
#                 end
#                 all_empirical_clusters = [temp_dict[key] for key in keys(temp_dict)]
#                 for key in keys(membership_dict_t)
#                     true_cluster_membership = membership_dict_t[key]
#                     inclusive_check = [issubset(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
#                     exlusive_check = [issetequal(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
#                     if any(inclusive_check)
#                         if haskey(empirical_cluster_membership_inclusive[t],true_cluster_membership)
#                             empirical_cluster_membership_inclusive[t][true_cluster_membership] += 1
#                         else
#                             empirical_cluster_membership_inclusive[t][true_cluster_membership] = 1
#                         end 
#                     end
#                     if any(exlusive_check)
#                         if haskey(empirical_cluster_membership_exclusive[t],true_cluster_membership)
#                             empirical_cluster_membership_exclusive[t][true_cluster_membership] +=1
#                         else
#                             empirical_cluster_membership_exclusive[t][true_cluster_membership] = 1
#                         end
#                     end
#                 end
                
#             end
#         end
#         emp_clus_mem_counts_inclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
#         emp_clus_mem_counts_exclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
#         for t in 1:T
#             for clus in 1:k
#                 key = true_cluster_membership_dict[t][clus]
#                 if haskey(empirical_cluster_membership_exclusive[t],key)
#                     emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
#                 end
#                 # emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
#                 if haskey(empirical_cluster_membership_inclusive[t],key)
#                     emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
#                 end
#                 # emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
#             end
#         end
#         exclusive_rate = [[emp_clus_mem_counts_exclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
#         inclusive_rate = [[emp_clus_mem_counts_inclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
#         return exclusive_rate, inclusive_rate
#     end

#     export getzIDs, getzDF
#     function getzIDs(tchain)
#         z_ids = findall(map(name -> occursin("z", string(name)), names(tchain)));
#         return z_ids
#     end
#     function getzDF(tchain,n_samples)
#         z_ids = getzIDs(tchain)
#         z_df = DataFrame(tchain[n_samples:end, z_ids, :])[:,3:end]
#         return z_df
#     end

#     export getlambdaIDs, getlambdaposteriorAvg, getlambdaposteriorDF
#     function getlambdaposteriorAvg(tchain,n_samples)
#         λ_df = getlambdaposteriorDF(tchain,n_samples)
#         cluster_means = [mean(col) for col in  eachcol(λ_df)]
#         cluster_std = [std(col) for col in  eachcol(λ_df)]
#         cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
#         return cluster_mat
#     end
#     function getlambdaposteriorDF(tchain,n_samples)
#         λ_ids = getlambdaIDs(tchain)
#         λ_df = DataFrame(tchain[n_samples:end, λ_ids, :])[:,3:end]
#         return λ_df
#     end
#     function getlambdaIDs(tchain)
#         λ_ids = findall(map(name -> occursin("λ", string(name)), names(tchain)))
#         return λ_ids
#     end

#     export getmuposteriorAvg,getmuposteriorDF,getmuIDs
#     function getmuposteriorAvg(tchain,n_samples)
#         m_df = getmuposteriorDF(tchain,n_samples)
#         cluster_means = [mean(col) for col in  eachcol(m_df)]
#         cluster_std = [std(col) for col in  eachcol(m_df)]
#         cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
#         return cluster_mat
#     end
#     function getmuposteriorDF(tchain,n_samples)
#         m_ids = getmuIDs(tchain)
#         m_df = DataFrame(tchain[n_samples:end, m_ids, :])[:,3:end]
#         return m_df
#     end
#     function getmuIDs(tchain)
#         m_ids = findall(map(name -> occursin("m", string(name)), names(tchain)))
#         return m_ids
#     end

#     export getsigmaposteriorAvg,getsigmaposteriorDF,getsigmaIDs
#     function getsigmaposteriorAvg(tchain,n_samples)
#         s_df = getsigmaposteriorDF(tchain,n_samples)
#         cluster_means = [mean(col) for col in  eachcol(s_df)]
#         cluster_std = [std(col) for col in  eachcol(s_df)]
#         cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
#         return cluster_mat
#     end
#     function getsigmaposteriorDF(tchain,n_samples)
#         s_ids = getsigmaIDs(tchain)
#         s_df = DataFrame(tchain[n_samples:end, s_ids, :])[:,3:end]
#         return s_df
#     end
#     function getsigmaIDs(tchain)
#         s_ids = findall(map(name -> occursin("s²", string(name)), names(tchain)))
#         return s_ids
#     end


    
#     export generateDataDF_from_TuringChain,get_var_ids,get_chain_var,get_static_chn_param_univariate,get_dynamic_chn_param_univariate,unpack_data_chn_param
#     function generateDataDF_from_TuringChain(chn)
#         return DataFrame(chn)
#     end
#     function get_var_ids(chain, var_)
#         # var_ids = findall(map(name -> occursin(var_, string(name)), names(chain)))
#         var_ids = MCMCChains.namesingroup(chain, Symbol(var_))
#         return var_ids
#     end
#     function get_chain_var(chain)
#         return string.(collect(keys(get(chain; section=:parameters))))
#     end
#     function get_static_chn_param_univariate(var_, chn)
#         _id = get_var_ids(chn, var_)
#         _chn = chn[:, _id, :]
#         _df = generateDataDF_from_TuringChain(_chn)
#         col_names_ = [m.match for m in match.(general_re_func(var_),names(_df)) if !isnothing(m) ]
#         _true =  [_df[:,c][1] for c in col_names_]
#         return _true
#     end
#     function get_dynamic_chn_param_univariate(var_, chn, T)
#         _id = get_var_ids(chn, var_)
#         _chn = chn[:, _id, :]
#         _df = generateDataDF_from_TuringChain(_chn)
#         col_names_ = [[m.match for m in match.(time_re_func(t,var_),names(_df)) if !isnothing(m) ] for t in 1:T]
#         _true =  [[ col for col in eachcol(_df[!,c])] for c in col_names_]
#         return _true
#     end
#     function unpack_data_chn_param(_true)
#         return [[c[1] for c in el ] for el in _true]
#     end
    





#     export get_πposterior,get_Poisson_cell_lhprob,get_Normal_cell_lhprob,get_P_tensor,stephens_relabelling,update_nu, init_nu,calculate_Q_hat,calculate_Cost_per_t,calculate_Cost_across_t,create_CostMat_across_t,get_average_posterior_cluster_frequency,get_average_posterior_cluster_frequency2,get_clus_w_maxPopAtT,FindNonunique_func
#     function get_πposterior(chain, num_burnin, T)
#         v_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, :v), :])[!,3:end]
#         v_t_vec = [[[r[j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(r)) if !isnothing(m) ]] for t in 1:T] for r in eachrow(v_df)]
#         π_t_post = [stickbreak.(itr) for itr in v_t_vec]
#         return π_t_post
#     end
#     function get_Poisson_cell_lhprob(chain,x, num_burnin, T, C_t)
#         λdf =  getlambdaposteriorDF(chain,num_burnin)
#         param_lh = [Poisson.(collect(r)) for r in eachrow(λdf)]
#         cell_lhprob = [[[pdf.(itr_lh,x[t][c]) for c in 1:C_t[t]] for t in 1:T] for itr_lh in param_lh]
#         return cell_lhprob
#     end
#     function get_Normal_cell_lhprob(chain,x, num_burnin, T, C_t)
#         mdf =  getmuposteriorDF(chain,num_burnin)
#         sdf =  getsigmaposteriorDF(chain,num_burnin)
#         ms = [collect(r) for r in eachrow(mdf)]
#         ss = [collect(r) for r in eachrow(sdf)]
#         param_lh = [Normal.(ms[i],ss[i]) for i in 1:length(ms)]
#         cell_lhprob = [[[pdf.(itr_lh,x[t][c]) for c in 1:C_t[t]] for t in 1:T] for itr_lh in param_lh]
#         return cell_lhprob
#     end
#     function get_P_tensor(S,T,C_t,π_t_post,cell_lhprob; addPseudoCount = true, pseudoCount = 1)
#         if  addPseudoCount
#             pseudoCount = pseudoCount
#         else
#             pseudoCount = 0
#         end
#         P = Vector{Vector{Vector{Vector{Float64}}}}(undef, S)
#         for s in 1:S
#             P[s] = [[ (π_t_post[s][t] .* cell_lhprob[s][t][c] .+ pseudoCount) ./ sum((π_t_post[s][t] .* cell_lhprob[s][t][c] .+ pseudoCount)) for c in 1:C_t[t]] for t in 1:T ]
#         end
#         return P
#     end
#     function stephens_relabelling(P,S,T,C_t,KMax,num_itr;use_identity_perm = true,permute_rows = false)
#         relabelling_results = Dict()
#         relabelling_results["nu"] = Vector{Vector{Vector{Int64}}}(undef,num_itr)
#         relabelling_results["cost"] = Vector{Vector{Float64}}(undef,num_itr)
#         nu = init_nu(S,KMax;use_identity_perm = use_identity_perm)
#         for i in 1:num_itr
#             Q = calculate_Q_hat(P,S,T,C_t,KMax,nu)
#             cost_per_t = calculate_Cost_per_t(P,Q,S,T,C_t,KMax,nu; permute_rows = permute_rows)#calculate_Cost_per_t(P,Q,S,T,C_t,KMax)
#             cost_across_t = calculate_Cost_across_t(cost_per_t,S,T,KMax)
#             CostMat_across_t = create_CostMat_across_t(cost_across_t,S,KMax)
#             updatesDict = update_nu(CostMat_across_t,S)
#             nu = updatesDict["nu"]
#             cost = updatesDict["cost"]
#             relabelling_results["nu"][i] = nu
#             relabelling_results["cost"][i] = cost
#         end
#         return relabelling_results
        
#     end
#     function update_nu(CostMat_across_t,S)
#         updatesDict = Dict()
#         updatesDict["nu"] = Vector{Vector{Int64}}(undef,S)
#         updatesDict["cost"] = Vector{Float64}(undef,S)
#         for s in 1:S
#             assignment, cost = hungarian(CostMat_across_t[s])
#             updatesDict["nu"][s] = assignment
#             updatesDict["cost"][s] = cost
#         end
#         return updatesDict
        
#     end
#     function init_nu(S,KMax;use_identity_perm = true)
#         nu =  Vector{Vector{Int}}(undef,S)
#         for s in 1:S
#             if use_identity_perm
#                 _perm = collect(1:KMax)
#             else
#                 _perm = randperm(KMax)
#             end
#             nu[s] = _perm#randperm(KMax)#
#         end
#         return nu
#     end
#     function calculate_Q_hat(P,S,T,C_t,KMax,nu)
#         Q = Vector{Vector{Vector{Float64}}}(undef,T)
#         for t in 1:T
#             num_cells = C_t[t]
#             Q[t] = Vector{Vector{Float64}}(undef,num_cells)
#             for c in 1:num_cells
#                 Q[t][c] = Vector{Float64}(undef,KMax)
#                 for k in 1:KMax
#                     # nu_k = nu[s][k]
#                     all_vals = [P[s][t][c][nu[s][k]] for s  in 1:S]
#                     Q[t][c][k] = mean(all_vals)
#                 end
#             end
#         end
#         return Q
#     end
#     function calculate_Cost_per_t(P,Q,S,T,C_t,KMax,nu; permute_rows = false)
#         cost_per_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
#         for s in 1:S
#             cost_per_t[s] = Vector{Vector{Vector{Float64}}}(undef,T)
#             for t in 1:T
#                 cost_per_t[s][t] = Vector{Vector{Float64}}(undef,KMax)
#                 for j in 1:KMax
#                     cost_per_t[s][t][j] = Vector{Float64}(undef,KMax)
#                     for l in 1:KMax
#                         _q_hat(c)  = Q[t][c][j]
#                         # _p(c) = P[s][t][c][nu[s][l]]
#                         _p(c) = P[s][t][c][l]
#                         _kl(c) = _p(c)*(log(_p(c))-log(_q_hat(c)))
#                         kl_cells = [_kl(c) for c in 1:C_t[t]]
#                         cost_per_t[s][t][j][l] = sum( kl_cells)
#                     end
#                 end
#                 if permute_rows
#                     old_ref = deepcopy(cost_per_t[s][t])
#                     for j in 1:KMax
#                         permuted_idx = nu[s][j]
#                         cost_per_t[s][t][permuted_idx] = old_ref[j]
#                     end
#                 end
#             end
#         end
#         return cost_per_t
#     end
#     function calculate_Cost_across_t(cost_per_t,S,T,KMax)
#         cost_across_t = Vector{Vector{Vector{Float64}}}(undef,S)
#         for s in 1:S
#             cost_across_t[s] = Vector{Vector{Float64}}(undef,KMax)
#             for j in 1:KMax
#                 cost_across_t[s][j] = Vector{Float64}(undef,KMax)
#                 for l in 1:KMax
#                     cost_across_t[s][j][l] = sum([cost_per_t[s][t][j][l] for t in 1:T])
#                 end
#             end
#         end
#         return cost_across_t
#     end
#     function create_CostMat_across_t(cost_across_t,S,KMax)
#         CostMat_across_t =Vector{Matrix{Float64}}(undef,S)
#         for s in 1:S
#             CostMat_across_t[s] = -188.0*ones(Float64,KMax,KMax)
#             for l in 1:KMax
#                 for j in 1:KMax
#                     CostMat_across_t[s][j, l] = cost_across_t[s][j][l]
#                 end
#             end
#         end
#         return CostMat_across_t
#     end
#     function get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         S = length(tchain)
#         z_ids =  get_var_ids(tchain, "z")
#         z_post_s = [unpack_data_chn_param(get_dynamic_chn_param_univariate("z", tchain[s,z_ids,:], T)) for s in 1:S]
#         post_z_dict_s = genData_vec2dict.(z_post_s)
#         true_z_dict = genData_vec2dict(true_z)
#         counts_mat_s = zeros(Int,KMax,KTrue,T,S)
#         for s in 1:S
#             for t in 1:T
#                 for c in 1:C_t[t]
#                     row = Int(post_z_dict_s[s][c,t])
#                     col = Int(true_z_dict[c,t])
#                     counts_mat_s[row,col,t,s] += 1
#                 end
#             end
#         end
#         avg_counts_mat = mean(counts_mat_s,dims=4)
#         return avg_counts_mat
#     end
#     function get_average_posterior_cluster_frequency2(z_post_s,T,true_z,KMax,KTrue,C_t)
#         S = length(z_post_s)
#         post_z_dict_s = genData_vec2dict.(z_post_s)
#         true_z_dict = genData_vec2dict(true_z)
#         counts_mat_s = zeros(Int,KMax,KTrue,T,S)
#         for s in 1:S
#             for t in 1:T
#                 for c in 1:C_t[t]
#                     row = Int(post_z_dict_s[s][c,t])
#                     col = Int(true_z_dict[c,t])
#                     counts_mat_s[row,col,t,s] += 1
#                 end
#             end
#         end
#         avg_counts_mat = mean(counts_mat_s,dims=4)
#         return avg_counts_mat
#     end
#     function get_clus_w_maxPopAtT(avg_counts_mat,t)
#         maxPopAtT = findmax(avg_counts_mat[:,:,t,1];dims=2)
#         ClusIndx_w_maxPopAtT = maxPopAtT[2]
#         clus_w_maxPopAtT = ["True Cluster "*string(k[2]) for k in ClusIndx_w_maxPopAtT]
#         return vec(clus_w_maxPopAtT)
#     end
#     function FindNonunique_func(x)
#         return  [k for (k, v) in countmap(x) if v > 1] # needs StatsBase
#     end
    
    



#     export get_chn_values_dict,relabel_chain_variable_dict
#     function get_chn_values_dict(chain,T,num_burnin,π_t_post,chn_var,dymamic_var_dict)
        
#         chn_values_dict = Dict()
#         for chr in chn_var
#             if chr == "π_t"
#                 chn_values_dict[chr] = π_t_post
#             else
#                 key = Symbol(chr)
#                 var_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, key), :])[!,3:end]
#                 if typeof(dymamic_var_dict[chr]) <: Bool
#                     if dymamic_var_dict[chr]
#                         gg(t,chr,var_df) = [m.match for m in match.(_position1_re_func(t,chr),names(var_df)) if !isnothing(m) ]
#                         var_vov = [[[r[c] for c in gg(t,chr,var_df)] for t in 1:T] for r in eachrow(var_df)]
#                     else
#                         ff(chr) = [m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
#                         var_vov = [[[r[c] for c in ff(chr)] for t in 1:T] for r in eachrow(var_df)]
#                     end
#                 else
#                     if dymamic_var_dict[chr] == "cluster"
#                         hh(chr) = [m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
#                         var_vov = [[[r[c] for c in hh(chr)]] for r in eachrow(var_df)]
#                     end
#                 end
#                 chn_values_dict[chr] = var_vov
#             end
#         end
#         return chn_values_dict
#     end
#     function relabel_chain_variable_dict(relabel_chn_values_dict,S,T,KMax,chn_var,relabel_var_dict, relabel_perm)
#         for chr in chn_var
#             if typeof(relabel_var_dict[chr]) <: AbstractString
#                 for s in 1:S
#                     for t in 1:T
#                         if relabel_var_dict[chr] == "cluster"
#                             new_indx = relabel_perm[s]
#                             relabel_chn_values_dict[chr][s][t] =  relabel_chn_values_dict[chr][s][t][new_indx]
#                         elseif relabel_var_dict[chr] == "cluster-param"
#                             new_indx = relabel_perm[s]
#                             relabel_chn_values_dict[chr][s][1] =  relabel_chn_values_dict[chr][s][1][new_indx]
#                         elseif relabel_var_dict[chr] == "assignment"
#                             itr_relabel_dict = Dict(j => relabel_perm[s][j] for j in 1:KMax)
#                             relabel_chn_values_dict[chr][s][t] = [itr_relabel_dict[el] for el in relabel_chn_values_dict[chr][s][t]]
#                         end
#                     end
#                 end
#             end
#         end
#         return relabel_chn_values_dict
#     end

#     export make_θ_t_var_names,make_π_t_var_names,make_v_var_names,make_z_var_names,make_λ_var_names,make_m_var_names,make_s_var_names
#     function make_θ_t_var_names(T)
#         var_ = "θ_t"
#         return [var_*"[$(t)]" for t in 1:T]
#     end
#     function make_π_t_var_names(T,KMax)
#         var_ = "π_t"
#         return [var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
#     end
#     function make_v_var_names(T,KMax)
#         var_ = "v"
#         return [var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax-1 ]
#     end
#     function make_z_var_names(T,C_t)
#         var_ = "z"
#         return vcat([[var_*"[$(t)][$(c)]" for c in 1:C_t[t]] for t in 1:T]...)
#     end
#     function make_λ_var_names(KMax)
#         var_ = "λ"
#         return [var_*"[$(k)]" for k in 1:KMax]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
#     end
#     function make_m_var_names(KMax)
#         var_ = "m"
#         return [var_*"[$(k)]" for k in 1:KMax]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
#     end
#     function make_s_var_names(KMax)
#         var_ = "s"
#         return  [var_*"[$(k)]" for k in 1:KMax]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
#     end
    
#     export make_θ_t_vov2Mat,make_π_t_vov2Mat,make_v_vov2Mat,make_z_vov2Mat,make_λ_vov2Mat,make_m_vov2Mat,make_s_vov2Mat,make_vectorParam_vov2Mat
#     function make_θ_t_vov2Mat(relabel_chn_values_dict)
#         var_ = "θ_t"
#         vov = relabel_chn_values_dict[var_]
        
#         return permutedims(reduce(hcat,reduce.(vcat, vov) ))#reduce(vcat,reduce.(vcat,vov))
#     end
#     function make_π_t_vov2Mat(relabel_chn_values_dict)
#         var_ = "π_t"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_v_vov2Mat(relabel_chn_values_dict)
#         var_ = "v"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_z_vov2Mat(relabel_chn_values_dict)
#         var_ = "z"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_λ_vov2Mat(relabel_chn_values_dict)
#         var_ = "λ"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_m_vov2Mat(relabel_chn_values_dict)
#         var_ = "m"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_s_vov2Mat(relabel_chn_values_dict)
#         var_ = "s"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     ##### GENERAL EXCEPT FOR θ_t ##########
#     function make_vectorParam_vov2Mat(param,relabel_chn_values_dict)
#         var_ = param
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end



    
#     export partition_matrix_cluster_df_names,partition_matrix_cluster_df_names2,_position1_Matrix_re_func,_position2_Matrix_re_func,_position2_Matrix_re_func2,π_t_post_vov2Mat,_get_πposterior,partition_cluster_df_names

#     function partition_matrix_cluster_df_names(df,K)
#         [ [m.match for m in match.(_position2_Matrix_re_func(t,"λ"),names(df)) if !isnothing(m) ] for t in 1:K]
#     end
#     function partition_matrix_cluster_df_names2(df,K)
#         [ [m.match for m in match.(_position2_Matrix_re_func2(t,"λ"),names(df)) if !isnothing(m) ] for t in 1:K]
#     end
#     function _position1_Matrix_re_func(g,var)
#         return Regex(string(var)*"\\["*string(g)*"(,\\[[0-9]+\\])*\\]")
#     end
#     function _position2_Matrix_re_func(g,var)
#         return Regex(string(var)*"\\[[0-9]+,"*string(g)*"\\](\\[[0-9]+\\])*")
#     end
#     function _position2_Matrix_re_func2(g,var)
#         return Regex(string(var)*"\\[([0-9]+,)*"*string(g)*"\\](\\[[0-9]+\\])*")
#     end
#     function π_t_post_vov2Mat(π_t_post)
#         return permutedims(reduce(hcat,reduce.(vcat,π_t_post)))
#     end
#     function _get_πposterior(chain, num_burnin, T;islogProb= false)
#         v_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, :v), :])[!,3:end]
#         v_t_vec = [[[r[j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(r)) if !isnothing(m) ]] for t in 1:T] for r in eachrow(v_df)]
#         π_t_post = [vec(map(el-> islogProb ? log.(stickbreak(el)) :  stickbreak(el) ,itr)) for itr in v_t_vec] #[stickbreak.(itr) for itr in v_t_vec]
#         return π_t_post
#     end
#     function partition_cluster_df_names(df,K)
#         [ [m.match for m in match.(_position1_re_func(t,"λ"),names(df)) if !isnothing(m) ] for t in 1:K]
#     end

#     export get_cell_lhprob_mat,get_post_cluster_assgn,get_post_cluster_assgnMat2,_get_chnMatrix_values_dict,_relabel_chain_variable_dict,make_PIP_t_var_names,make_z_t_vov2Mat,make_PIP_t_vov2Mat,make_λmat_vov2Mat,make_λmat_var_names,getLambdaInferenceModel2_newNames,extract
#     ######### I THINK SOMETHING IS UP HERE OR I NEED TO RUN IT LONGER(only ran for 200 itr for testing purposes)!!! #################
#     function get_cell_lhprob_mat(chain,x, lh_func, num_burnin, T, C_t,KMax;islog= false)
#         λdf =  getlambdaposteriorDF(chain,num_burnin)
#         S = size(λdf)[1]
#         λdf_names = partition_matrix_cluster_df_names(λdf,KMax)
#         get_row_parameterize(df,row , df_names) = collect(df[row,df_names])
#         param_lh = Vector{Vector}(undef,S)
#         if islog
#             pdf_func = logpdf
#         else
#             pdf_func = pdf
#         end
#         for s in 1:S
#             row_params = vec(map(el -> get_row_parameterize(λdf,s, el), λdf_names))
#             param_lh[s] = map(el -> lh_func(el),row_params)
#         end
#         # param_lh = [Poisson.(collect(r)) for r in eachrow(λdf)]
#         cell_lhprob = [[[ 
#                         map(clus_lh -> pdf_func(clus_lh,x[t][c]),param_lh[s])  
#                             for c in 1:C_t[t]]  
#                                 for t in 1:T]  
#                                     for s in 1:S] #[[[ vcat(map(el ->  pdf_func(itr_lh[el],x[t][c]), 1:KMax )...) for c in 1:C_t[t]] for t in 1:T] for itr_lh in param_lh]
#         return cell_lhprob
#     end

#     function get_post_cluster_assgn(tchain,x,num_burnin,T, C_t)
#         v_ids = MCMCChains.namesingroup(tchain, :v)
#         vpost = Matrix(DataFrame(tchain[:,v_ids,:])[!,3:end])
#         πpost = get_πposterior(tchain,num_burnin, T);
#         S = size(πpost)[1]
#         cell_lhprob = get_Poisson_cell_lhprob(tchain,x, num_burnin, T, C_t);
#         post_ip = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
#         zpost = Vector{Vector{Vector{Int}}}(undef,S)
#         for s in 1:S 
#             post_ip[s] = Vector{Vector{Vector{Float64}}}(undef,T)
#             zpost[s] = Vector{Vector{Float64}}(undef,T)
#             colend = 0
#             for t in 1:T
#                 cells = C_t[t]
#                 post_ip[s][t] = Vector{Vector{Float64}}(undef,cells)
#                 # zpost[s][t] = zeros{Float64}(undef,cells)
#                 for c in 1:cells
#                     val_vec = πpost[s][t] .* cell_lhprob[s][t][c]
#                     val_sum = sum(val_vec)
#                     post_ip[s][t][c] = val_vec ./ val_sum
#                 end
#                 zpost[s][t] = rand.(Categorical.(post_ip[s][t]))
#             end
#         end
        
#         return zpost,post_ip
        
#     end

#     ####USE THIS ONE!!!!
#     function get_post_cluster_assgnMat2(tchain,x,num_burnin,T, C_t,KMax,lh_func;uselogProb= false)
#         v_ids = MCMCChains.namesingroup(tchain, :v)
#         vpost = Matrix(DataFrame(tchain[:,v_ids,:])[!,3:end])
#         πpost = _get_πposterior(tchain, num_burnin, T;islogProb= uselogProb); #get_πposterior(tchain,num_burnin, T);
#         S = size(πpost)[1]
#         cell_lhprob = get_cell_lhprob_mat(tchain,x, lh_func, num_burnin, T, C_t,KMax;islog= uselogProb); #get_Poisson_cell_lhprob(tchain,x, num_burnin, T, C_t);
#         post_ip = Vector{Vector{Vector{Vector{Float64}}}}(undef,S)
#         zpost = Vector{Vector{Vector{Int}}}(undef,S)
#         for s in 1:S 
#             post_ip[s] = Vector{Vector{Vector{Float64}}}(undef,T)
#             zpost[s] = Vector{Vector{Float64}}(undef,T)
#             colend = 0
#             for t in 1:T
#                 cells = C_t[t]
#                 post_ip[s][t] = Vector{Vector{Float64}}(undef,cells)
#                 # zpost[s][t] = zeros{Float64}(undef,cells)
#                 for c in 1:cells
#                     if uselogProb
#                         val_vec = πpost[s][t] .+ cell_lhprob[s][t][c]
#                         val_sum = StatsFuns.logsumexp(val_vec)
#                         post_ip[s][t][c] = exp.(val_vec .- val_sum)
#                     else
#                         val_vec = πpost[s][t] .* cell_lhprob[s][t][c]
#                         val_sum = sum(val_vec)
#                         post_ip[s][t][c] = val_vec ./ val_sum
#                     end
                    
#                 end
#                 zpost[s][t] = rand.(Categorical.(post_ip[s][t]))
#             end
#         end
        
#         return zpost,post_ip
        
#     end

#     function _get_chnMatrix_values_dict(chain,T,KMax,num_burnin,π_t_post,post_ip,zpost,chn_var,dymamic_var_dict)
            
#         chn_values_dict = Dict()
#         for chr in chn_var
#             if chr == "π_t"
#                 chn_values_dict[chr] = π_t_post
#             elseif  chr == "PIP_t"
#                 chn_values_dict[chr] = post_ip
#             elseif  chr == "z_t"
#                 chn_values_dict[chr] = zpost
#             else
#                 key = Symbol(chr)
#                 var_df = DataFrame(chain[num_burnin:end, MCMCChains.namesingroup(chain, key), :])[!,3:end]
#                 if typeof(dymamic_var_dict[chr]) <: Bool
#                     if dymamic_var_dict[chr]
#                         gg(t,chr,var_df) = [m.match for m in match.(_position1_re_func(t,chr),names(var_df)) if !isnothing(m) ]
#                         var_vov = [[[r[c] for c in gg(t,chr,var_df)] for t in 1:T] for r in eachrow(var_df)]
#                     else
#                         ff(chr) = [m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
#                         var_vov = [[[r[c] for c in ff(chr)] for t in 1:T] for r in eachrow(var_df)]
#                     end
#                 else
#                     if dymamic_var_dict[chr] == "cluster"
#                         hh(chr) =  [ [m.match for m in match.(_position2_Matrix_re_func(t,chr),names(var_df)) if !isnothing(m) ] for t in 1:KMax]#[m.match for m in match.(general_re_func(chr),names(var_df)) if !isnothing(m) ]
#                         var_vov = [[[collect(r[c])[1:end] for c in hh(chr)]] for r in eachrow(var_df)]
#                     end
#                 end
#                 chn_values_dict[chr] = var_vov
#             end
#         end
#         return chn_values_dict
#     end
#     function _relabel_chain_variable_dict(relabel_chn_values_dict,S,T,KMax,chn_var,relabel_var_dict, relabel_perm)
#         for chr in chn_var
#             if typeof(relabel_var_dict[chr]) <: AbstractString
#                 for s in 1:S
#                     for t in 1:T
#                         if relabel_var_dict[chr] == "cluster"
#                             new_indx = relabel_perm[s]
#                             relabel_chn_values_dict[chr][s][t] =  relabel_chn_values_dict[chr][s][t][new_indx]
#                         # elseif relabel_var_dict[chr] == "cluster-param"
#                         #     println("here son")
#                         #     new_indx = relabel_perm[s]
#                         #     relabel_chn_values_dict[chr][s][1] =  relabel_chn_values_dict[chr][s][1][new_indx]
#                         elseif relabel_var_dict[chr] == "PIP"
#                             new_indx = relabel_perm[s]
#                             for c in 1:C_t[t]
#                                 relabel_chn_values_dict[chr][s][t][c] =  relabel_chn_values_dict[chr][s][t][c][new_indx] 
#                             end
#                         elseif relabel_var_dict[chr] == "assignment"
#                             itr_relabel_dict = Dict(j => relabel_perm[s][j] for j in 1:KMax)
#                             relabel_chn_values_dict[chr][s][t] = [itr_relabel_dict[el] for el in relabel_chn_values_dict[chr][s][t]]
#                         end
#                     end
#                 if relabel_var_dict[chr] == "cluster-param"
                    
#                     new_indx = relabel_perm[s]
#                     relabel_chn_values_dict[chr][s][1] =  relabel_chn_values_dict[chr][s][1][new_indx]
#                 end
#                 end
#             end
#         end
#         return relabel_chn_values_dict
#     end

#     function make_PIP_t_var_names(T,C_t,KMax)
#         var_ = "PIP_t"
#         return [var_*"[$(t)][$(c)][$(k)]" for t in 1:T for c in 1:C_t[t] for k in 1:KMax]
#     end
#     function make_z_t_vov2Mat(relabel_chn_values_dict)
#         var_ = "z_t"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_PIP_t_vov2Mat(relabel_chn_values_dict)
#         var_ = "PIP_t"
#         vov = relabel_chn_values_dict[var_]
#         return permutedims(reduce(hcat , reduce.(vcat,reduce.(vcat,vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_λmat_vov2Mat(relabel_chn_values_dict)
#         var_ = "λ"
#         vov = relabel_chn_values_dict[var_]
#         vov2 = reduce.(vcat,(reduce.(vcat, vov))) #permutedims(reduce.(vcat,(reduce.(vcat, vov))))
#         return  permutedims(reduce(hcat, vov2))#permutedims(reduce(hcat,(reduce.(vcat, vov))))#permutedims(reduce(hcat,reduce(vcat,vov)))#permutedims(hcat(vcat(vov...)...))
#     end
#     function make_λmat_var_names(KMax,G)
#         var_ = "λ"
#         return [var_*"[$(g),$(k)]"  for k in 1:KMax for g in 1:G]#[var_*"[$(t)][$(k)]" for t in 1:T for k in 1:KMax]
#     end

#     function getLambdaInferenceModel2_newNames(G,KMax)
#         rename_dict = Dict()
#         var = "λ"
#         for g in 1:G
#             for k in 1:KMax
#                 old_name = var*"[$k][$g]"
#                 new_name = var*"[$g,$k]"
#                 rename_dict[old_name] = new_name
#             end
#         end
#         return rename_dict
#     end

#     function extract(tchain, var; burn=0, is_scalar=true)
#         sym = Symbol(var)
#         if is_scalar
#             tail = tchain[sym].data[(burn + 1):end, 1]
#         else
#             tail = group(tchain, sym).value.data[(burn + 1):end, :, 1]
#             # NOTE: `chain[:sym]` will not work in the future for parameter
#             # vectors. Use `group` instead.
#         end
#         return tail
#     end


    

#     export re_func, time_re_func,cellID_re_func, _position1_re_func,_position2_re_func,_position3_re_func,general_re_func
#     function re_func(t,var)
#         return Regex(string(var)*"\\["*string(t)*"\\]\\[[0-9]+\\]")
#     end
#     function time_re_func(t,var)
#         return Regex(string(var)*"\\["*string(t)*"\\]\\[[0-9]+\\]")
#     end
#     function cellID_re_func(id,var)
#         return Regex(string(var)*"\\[[0-9]+\\]\\["*string(id)*"\\]")
#     end

#     function _position1_re_func(t,var)
#         return Regex(string(var)*"\\["*string(t)*"\\](\\[[0-9]+\\])*")
#     end
#     function _position2_re_func(id,var)
#         return Regex(string(var)*"\\[[0-9]+\\]\\["*string(id)*"\\](\\[[0-9]+\\])*")
#     end
#     function _position3_re_func(g,var)
#         return Regex(string(var)*"\\[[0-9]+\\]\\[[0-9]+\\]\\["*string(g)*"\\](\\[[0-9]+\\])*")
#     end
#     function general_re_func(var)
#         return  Regex(string(var)*"(\\[[0-9]+\\])+")
#     end

# end