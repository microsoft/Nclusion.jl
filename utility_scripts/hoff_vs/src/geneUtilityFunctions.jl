
function create_vec_cluster_scores_mat(cell_ll_scores)
    T = length(cell_ll_scores);
    C_t = [length(el) for el in cell_ll_scores];
    KMax = length(cell_ll_scores[1][1]);
    G = length(cell_ll_scores[1][1][1]);
    
    total_number_cells = sum(C_t)
    cluster_scores = Vector{Matrix{Float64}}(undef,KMax)

    for k in 1:KMax
        k_scores = Matrix{Float64}(undef,total_number_cells,G)
        lin_indx = 0
        for t in 1:T
            num_cells = C_t[t]
            for i in 1:num_cells
                lin_indx +=1
                k_scores[lin_indx,:] = cell_ll_scores[t][i][k]

            end
        end
        cluster_scores[k] = k_scores

    end
    return cluster_scores
end
function create_vec_cluster_gene_weights_mat(cell_ll_gene_weights)
    T = length(cell_ll_gene_weights);
    C_t = [length(el) for el in cell_ll_gene_weights];
    KMax = length(cell_ll_gene_weights[1][1]);
    G = length(cell_ll_gene_weights[1][1][1]);
    
    total_number_cells = sum(C_t)
    cluster_gene_weights = Vector{Matrix{Float64}}(undef,KMax)

    for k in 1:KMax

        k_gene_weights = Matrix{Float64}(undef,total_number_cells,G)
        lin_indx = 0
        for t in 1:T
            num_cells = C_t[t]
            for i in 1:num_cells
                lin_indx +=1
                k_gene_weights[lin_indx,:] = cell_ll_gene_weights[t][i][k]
            end
        end

        cluster_gene_weights[k] = k_gene_weights
    end
    return cluster_gene_weights
end
function create_vec_cluster_gene_weightsbygene_mat(cell_ll_gene_weightsbygene)
    T = length(cell_ll_gene_weightsbygene);
    C_t = [length(el) for el in cell_ll_gene_weightsbygene];
    KMax = length(cell_ll_gene_weightsbygene[1][1][1]);
    G = length(cell_ll_gene_weightsbygene[1][1]);
    
    total_number_cells = sum(C_t)
    cluster_gene_weightsbygene = Array{Float64}(undef,KMax,G,total_number_cells)
    lin_indx = 0
    for t in 1:T
        num_cells = C_t[t]
        for i in 1:num_cells
            lin_indx +=1
            cluster_gene_weightsbygene[:,:,lin_indx] =  reduce(hcat,cell_ll_gene_weightsbygene[t][i])
        end
    end
    return cluster_gene_weightsbygene
end


function topGenesClusters_byCalledCluster(cluster_gene_weights,gene_names,z_post_s; topNGenes= 30)
    S = length(z_post_s);
    KMax = length(cluster_gene_weights);
    G = size(cluster_gene_weights[1])[2];

    ## Check if gene_names == G

    driver_genes = Dict()

    for k in 1:KMax
        cluster_membership_means = mean([vec(mean(cluster_gene_weights[k][vcat(z_post_s[s]...) .== k,:],dims=1)) for s in 1:S])
        if all(isnan.(cluster_membership_means))
            continue
        end
        driver_genes[k] = gene_names[first.(maxk!(collect(1:G), cluster_membership_means , topNGenes, initialized = true))];
    end
    return driver_genes
end
function topGenesClusters_byAllCells(cluster_gene_weights,gene_names; topNGenes= 30)
    KMax = length(cluster_gene_weights);
    G = size(cluster_gene_weights[1])[2];
    ## Check if gene_names == G


    driver_genes = Dict()
    for k in 1:KMax
        driver_genes[k] = gene_names[first.(maxk!(collect(1:G), vec(mean(cluster_gene_weights[k],dims=1)) , topNGenes, initialized = true))];
    end
    return driver_genes
end
function topGenesClustersbygene_byCalledCluster(cluster_gene_weightsbygene,gene_names,z_post_s; topNGenes= 30)
    S = length(z_post_s);
    KMax = size(cluster_gene_weightsbygene)[1];
    G = size(cluster_gene_weightsbygene)[2];
    ## Check if gene_names == G

    ## Check if gene_names == G

    driver_genes = Dict()

    for k in 1:KMax
        cluster_membership_means = mean([vec(mean(cluster_gene_weightsbygene[k,:,vcat(z_post_s[s]...) .== k],dims=2)) for s in 1:S])
        if all(isnan.(cluster_membership_means))
            continue
        end
        driver_genes[k] = gene_names[first.(maxk!(collect(1:G), cluster_membership_means , topNGenes, initialized = true))];
    end
    return driver_genes
end
function topGenesClustersbygene_byAllCells(cluster_gene_weightsbygene,gene_names; topNGenes= 30)
    KMax = size(cluster_gene_weightsbygene)[1];
    G = size(cluster_gene_weightsbygene)[2];
    ## Check if gene_names == G


    driver_genes = Dict()
    for k in 1:KMax
        driver_genes[k] = gene_names[first.(maxk!(collect(1:G), vec(mean(cluster_gene_weightsbygene[k,:,:],dims=2)) , topNGenes, initialized = true))];
    end
    return driver_genes
end



function cells_in_k_clus_index(z_post_s,s,k)
    return vcat(z_post_s[s]...) .== k
end
function cells_not_in_k_clus_index(z_post_s,s,k)
    return broadcast(!,vcat(z_post_s[s]...) .== k)
end
function gene_weights_clus(cluster_gene_weights,k,genename,s,gene_names,z_post_s)
    G = length(gene_names)
    cells_index = cells_in_k_clus_index(z_post_s,s,k)
    gene_weights =cluster_gene_weights[k][cells_index ,get_geneID(G,gene_names,genename)];
    return gene_weights
end
function gene_weights_not_clus(cluster_gene_weights,k,genename,s,gene_names,z_post_s)
    G = length(gene_names)
    cells_index = cells_not_in_k_clus_index(z_post_s,s,k)
    gene_weights =cluster_gene_weights[k][cells_index ,get_geneID(G,gene_names,genename)];
    return gene_weights
end
function gene_weights_clus_tensor(cluster_gene_weightsbygene,k,genename,s,gene_names,z_post_s)
    G = length(gene_names)
    cells_index = cells_in_k_clus_index(z_post_s,s,k)
    gene_weights =cluster_gene_weightsbygene[k,get_geneID(G,gene_names,genename),cells_index];
    return gene_weights
end
function gene_weights_not_clus_tensor(cluster_gene_weightsbygene,k,genename,s,gene_names,z_post_s)
    G = length(gene_names)
    cells_index = cells_not_in_k_clus_index(z_post_s,s,k)
    gene_weights =cluster_gene_weightsbygene[k,get_geneID(G,gene_names,genename),cells_index];
    return gene_weights
end
function gene_weights_all(cluster_gene_weights,k,genename,gene_names)
    G = length(gene_names)
    gene_weights = cluster_gene_weights[k][: ,get_geneID(G,gene_names,genename)];
    return gene_weights
end
function get_geneID(G,gene_names,genename)
    return collect(1:G)[gene_names .== genename][1]
end
function generate_topNGeneWeightValues_DF(cluster_gene_weights,k,driver_genes_clus,gene_names,z_post_s;post_sample =1)
    # m,s_sq,s,n,upb,lob
    topNGeneWeight_DF = DataFrame()
    for i in 1:length(driver_genes_clus[k])
        gene_weights_clus_array = gene_weights_clus(cluster_gene_weights,k,driver_genes_clus[k][i],post_sample,gene_names,z_post_s);
        # m,s_sq,s,n,upb,lob= calculate_gene_weight_summary_statistics(gene_weights_clus_array;pct =0.9);
        topNGeneWeight_DF[!,driver_genes_clus[k][i]] = gene_weights_clus_array
        # push!(,[driver_genes_clus[k][i] m s_sq s n upb lob])
    end
    return topNGeneWeight_DF
end
function generate_topNGeneWeightStats_DF(cluster_gene_weights,k,driver_genes_clus,gene_names,z_post_s;post_sample =1)
    # m,s_sq,s,n,upb,lob
    topNGeneWeightStats_DF = DataFrame(Gene = [], average = [],variance = [],stdev =[], n= [], upper_bound=[],lower_bound= [])
    for i in 1:length(driver_genes_clus[k])
        gene_weights_clus_array = gene_weights_clus(cluster_gene_weights,k,driver_genes_clus[k][i],post_sample,gene_names,z_post_s);
        m,s_sq,s,n,upb,lob= calculate_gene_weight_summary_statistics(gene_weights_clus_array;pct =0.9);
        push!(topNGeneWeightStats_DF,[driver_genes_clus[k][i] m s_sq s n upb lob])
    end
    return topNGeneWeightStats_DF
end
function generate_topNGeneWeightStats_DF_tensor(cluster_gene_weightsbygene,k,driver_genes_clus,gene_names,z_post_s;post_sample =1)
    # m,s_sq,s,n,upb,lob
    topNGeneWeightStats_DF = DataFrame(Gene = [], average = [],variance = [],stdev =[], n= [], upper_bound=[],lower_bound= [])
    for i in 1:length(driver_genes_clus[k])
        gene_weights_clus_array = gene_weights_clus_tensor(cluster_gene_weightsbygene,k,driver_genes_clus[k][i],post_sample,gene_names,z_post_s);
        m,s_sq,s,n,upb,lob= calculate_gene_weight_summary_statistics(gene_weights_clus_array;pct =0.9);
        push!(topNGeneWeightStats_DF,[driver_genes_clus[k][i] m s_sq s n upb lob])
    end
    return topNGeneWeightStats_DF
end

