function R_baseline_SOUP(data,gene_names,Kmax;Kmin = 2, scale_fator = 10000,num_cores =10,pseudocount=1)
    data.genenames= String.(collect(gene_names[!,1]))
    rename!(data,Symbol.(names(data)))
    data = permutedims(data,"genenames","")
    counts = data[!,2:end]
    row_names = data[!,1]
    for col in 1:length(names(counts))
        if all(typeof.(col) .<: AbstractString )
            counts[!,col] = parse.(Int64, counts[!,col])
        else
            continue
        end
    end
    
    
    @rput row_names;
    @rput counts;
    @rput Kmin;
    @rput Kmax;
    @rput scale_fator;
    @rput num_cores;
    @rput pseudocount;

    R"""
        library(devtools)
        library(SOUP)
        library(ggplot2)
        rownames(counts) <- row_names
        log_expr = log2(scaleRowSums(counts)*(scale_fator) + pseudocount)
        select_out = selectGenes(counts, type="count", n.cores = num_cores)
        select_genes = select_out$select.genes
        log_select_expr = log_expr[,colnames(log_expr) %in% select_genes]
        startTime <- Sys.time()
        soup_out = SOUP(log_select_expr, Ks=c(Kmin:Kmax),type="log")
        endTime <- Sys.time()
        diffTime =  endTime - startTime
        elapsedTime =  diffTime[[1]]
    """
    @rget soup_out;
    @rget elapsedTime;
    @rget log_select_expr;
    @rget select_genes;
    @rget select_out;
    @rget log_expr;

    return soup_out,elapsedTime,log_select_expr,select_genes,select_out,log_expr
end
function R_baseline_SOUP(data,gene_names,Kmax;Kmin = 2, scale_fator = 10000,num_cores =10,pseudocount=1, select_genes_pre = nothing)
    data.genenames= String.(collect(gene_names[!,1]))
    rename!(data,Symbol.(names(data)))
    data = permutedims(data,"genenames","")
    counts = data[!,2:end]
    row_names = data[!,1]
    for col in 1:length(names(counts))
        if all(typeof.(col) .<: AbstractString )
            counts[!,col] = parse.(Int64, counts[!,col])
        else
            continue
        end
    end
    
    
    @rput row_names;
    @rput counts;
    @rput Kmin;
    @rput Kmax;
    @rput scale_fator;
    @rput num_cores;
    @rput pseudocount;
    if !isnothing(select_genes_pre)
        select_genes_bool = true
        @rput select_genes_bool
        @rput select_genes_pre
    else
        select_genes_bool = false
        @rput select_genes_bool
    end

    R"""
        library(devtools)
        library(SOUP)
        library(ggplot2)
        library(profmem)
        rownames(counts) <- row_names
        log_expr = log2(scaleRowSums(counts)*(scale_fator) + pseudocount)
        if (select_genes_bool) {
            select_out = select_genes_pre
            select_genes = select_genes_pre
        } else {
            select_out = selectGenes(counts, type="count", n.cores = num_cores)
            select_genes = select_out$select.genes
        }
        
        log_select_expr = log_expr[,colnames(log_expr) %in% select_genes]
        startTime <- Sys.time()
        soup_out = SOUP(log_select_expr, Ks=c(Kmin:Kmax),type="log")
        endTime <- Sys.time()
        diffTime =  endTime - startTime
        elapsedTime =  diffTime[[1]]
        p <- profmem({
            profout<- SOUP(log_select_expr, Ks=c(Kmin:Kmax),type="log")
            })
        bytes_used <-  sum(p$bytes, na.rm=TRUE)
    """
    # println("heres")
    @rget soup_out;
    @rget elapsedTime;
    @rget bytes_used;
    @rget log_select_expr;
    @rget select_genes;
    @rget select_out;
    @rget log_expr;

    return soup_out,elapsedTime,bytes_used,log_select_expr,select_genes,select_out,log_expr
end


function all_unique(v::Vector)::Bool
    return length(unique(v)) == length(v)
end

function df_add_first_column(
    df::DataFrame,
    colname::Union{Symbol,String},
    col_data
)
    df1 = DataFrame([colname => col_data])
    hcat(df1, df)
end

function df_transpose(df::DataFrame, col::Union{Symbol, String})::DataFrame
    @assert all_unique(df[!, col]) "Column `col` contains non-unique elements"

    function foo(i)
        string(df[i, col]) => collect(df[i, Not(col)])
    end

    dft = DataFrame(map(foo, 1:nrow(df)))

    return df_add_first_column(dft, "Row", filter(x -> x != string(col), names(df)))
end


    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/cortex/all_genes/cortex_all_genes_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/18879-cortex_preprocessed-30-.csv" "" "" "" "" "" "" ""    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/cortex/5000hvg/cortex_5000hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/5000-cortex_preprocessed_5000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/cortex/2500hvg/cortex_2500hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/2500-cortex_preprocessed_2500-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/cortex/1000hvg/cortex_1000hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/1000-cortex_preprocessed_1000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/cortex/600hvg/cortex_600hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/600-cortex_preprocessed_600-30-.csv" "" "" "" "" "" "" ""
    
    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/cortex/all_genes/cortex_all_genes_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/18879-cortex_preprocessed-30-.csv" "" "" "" "" "" "" ""    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/cortex/5000_hvg/cortex_5000hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/5000-cortex_preprocessed_5000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/cortex/2500_hvg/cortex_2500hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/2500-cortex_preprocessed_2500-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/cortex/1000_hvg/cortex_1000hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/1000-cortex_preprocessed_1000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/cortex/600_hvg/cortex_600hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/cortex/600-cortex_preprocessed_600-30-.csv" "" "" "" "" "" "" ""


    ### ERROR WHEN RUN####
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/pbmc/all_genes/pbmc_all_genes_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/15720-pbmc_preprocessed-30-.csv" "" "" "" "" "" "" ""    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/pbmc/5000hvg/pbmc_5000hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/5000-pbmc_preprocessed_5000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/pbmc/2500hvg/pbmc_2500hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/2500-pbmc_preprocessed_2500-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/pbmc/1000hvg/pbmc_1000hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/1000-pbmc_preprocessed_1000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scLCA/pbmc/600hvg/pbmc_600hvg_scLCA.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/600-pbmc_preprocessed_600-30-.csv" "" "" "" "" "" "" ""
    
    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/pbmc/all_genes/pbmc_all_genes_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/15720-pbmc_preprocessed-30-.csv" "" "" "" "" "" "" ""    
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/pbmc/5000_hvg/pbmc_5000hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/5000-pbmc_preprocessed_5000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/pbmc/2500_hvg/pbmc_2500hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/2500-pbmc_preprocessed_2500-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/pbmc/1000_hvg/pbmc_1000hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/1000-pbmc_preprocessed_1000-30-.csv" "" "" "" "" "" "" ""
    # julia --project=. baselines/tsne_plotting_inferred_labels.jl "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/nclusion_experiments/quality_of_clustering/scCCESS_SIMLR/pbmc/600_hvg/pbmc_600hvg_scCCESS.csv" "/mnt/c/Users/cnwizu/Dropbox (Brown)/nclusion-datasets/pmbc/600-pbmc_preprocessed_600-30-.csv" "" "" "" "" "" "" ""