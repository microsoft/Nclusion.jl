function plotNMI(true_z,z_post_s;conf_level=0.95)
    T = length(true_z)
    nmi_vov = getNMI(true_z,z_post_s);
    conf_int = t_test.(nmi_vov; conf_level = conf_level);
    mean_vov = mean.(nmi_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([0.0 - 0.01 ,minimum(lob)  - 0.1])
    ylim_ub = minimum([1.0 + 0.01 , maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p=plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n Normalized Mutual Information\n over time series**", xlabel = "Time points", ylabel = "Normalized Mutual Information \n (NMI)", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5)))
     # "**ranges between [0,1]; Higher is Better" , annotation= (5,100,"**ranges between [0,1]; Higher is Better")
    return p
    
end

function plotNMI(nmi_vov;conf_level=0.95)
    T = length(nmi_vov)
    conf_int = t_test.(nmi_vov; conf_level = conf_level);
    mean_vov = mean.(nmi_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([0.0 - 0.01 ,minimum(lob)  - 0.1])
    ylim_ub = minimum([1.0 + 0.01 , maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p=plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n Normalized Mutual Information\n over time series**", xlabel = "Time points", ylabel = "Normalized Mutual Information \n (NMI)", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5)))
     # "**ranges between [0,1]; Higher is Better" , annotation= (5,100,"**ranges between [0,1]; Higher is Better")
    return p
    
end

function plotPosteriorClusteringNMIAndSave!(true_z,z_post_s,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext;conf_level=0.95)
    dataPosteriorClusteringNMI = plotNMI(true_z,z_post_s;conf_level=conf_level)
    static_plot_id = prefix *"PosteriorClusteringNMI"
    dataPosteriorClusteringNMI_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(dataPosteriorClusteringNMI, dataPosteriorClusteringNMI_name);
    push!(plot_name_vec, dataPosteriorClusteringNMI_name)
end

function plotAdjustedRandIndx(true_z,z_post_s;conf_level=0.95)
    T = length(true_z)
    ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov  = getRandIndices(true_z,z_post_s);
    conf_int = t_test.(ari_vov; conf_level = conf_level);
    mean_vov = mean.(ari_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([-1.0 - 0.01 ,minimum(lob)  - 0.1])
    ylim_ub = minimum([1.0 + 0.01, maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p = plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n Adjusted Rand Index (ARI)\n over time series**", xlabel = "Time points", ylabel = "Adjusted Rand Index \n (ARI)", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5)))
    # "**ranges between [-1,1]; Higher is Better" , annotation= (5,100,"**ranges between [-1,1]; Higher is Better")
    return p
    
end

function plotAdjustedRandIndx(ari_vov;conf_level=0.95)
    T = length(ari_vov)
    conf_int = t_test.(ari_vov; conf_level = conf_level);
    mean_vov = mean.(ari_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([-1.0 - 0.01 ,minimum(lob)  - 0.1])
    ylim_ub = minimum([1.0 + 0.01, maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p = plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n Adjusted Rand Index (ARI)\n over time series**", xlabel = "Time points", ylabel = "Adjusted Rand Index \n (ARI)", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5)))
    # "**ranges between [-1,1]; Higher is Better" , annotation= (5,100,"**ranges between [-1,1]; Higher is Better")
    return p
    
end


function plotPosteriorClusteringARIAndSave!(true_z,z_post_s,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext;conf_level=0.95)
    dataPosteriorClusteringARI = plotAdjustedRandIndx(true_z,z_post_s;conf_level=conf_level)
    static_plot_id = prefix *"PosteriorClusteringARI"
    dataPosteriorClusteringARI_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(dataPosteriorClusteringARI, dataPosteriorClusteringARI_name);
    push!(plot_name_vec, dataPosteriorClusteringARI_name)
end
function plotVMeasure(true_z,z_post_s;conf_level=0.95,beta=1.0)
    T = length(true_z)
    vmeasure_vov = getVmeasure(true_z,z_post_s;beta=beta);
    conf_int = t_test.(vmeasure_vov; conf_level = conf_level);
    mean_vov = mean.(vmeasure_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([0.0 - 0.01 ,minimum(lob)  - 0.1])
    ylim_ub = minimum([1.0 + 0.01, maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p  = plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n V-Measure \n over time series**", xlabel = "Time points", ylabel = "V-Measure", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5)))
    # "between 0,1 higher is better"
    return p
    
end

function plotVMeasure(vmeasure_vov;conf_level=0.95)
    T = length(vmeasure_vov)
    conf_int = t_test.(vmeasure_vov; conf_level = conf_level);
    mean_vov = mean.(vmeasure_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([0.0 - 0.01 ,minimum(lob)  - 0.1])
    ylim_ub = minimum([1.0 + 0.01, maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p  = plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n V-Measure \n over time series**", xlabel = "Time points", ylabel = "V-Measure", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5)))
    # "between 0,1 higher is better"
    return p
    
end

function plotPosteriorClusteringVMeasureAndSave!(true_z,z_post_s,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext;conf_level=0.95, beta = 1.0)
    dataPosteriorClusteringVMeasure = plotVMeasure(true_z,z_post_s;conf_level=conf_level,beta=beta)
    static_plot_id = prefix *"PosteriorClusteringVMeasure"
    dataPosteriorClusteringVMeasure_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(dataPosteriorClusteringVMeasure, dataPosteriorClusteringVMeasure_name);
    push!(plot_name_vec, dataPosteriorClusteringVMeasure_name)
end
function plotVarInfo(true_z,z_post_s;conf_level=0.95, KTrue = missing)
    T = length(true_z)
    varinfo_vov = getVarInfo(true_z,z_post_s);
    if ismissing(KTrue)
        KTrue =  maximum(unique(vcat(true_z...)))
    end
    conf_int = t_test.(varinfo_vov; conf_level = conf_level);
    mean_vov = mean.(varinfo_vov) ;
    lob = first.(conf_int) 
    hib = last.(conf_int)
    ylim_lb = maximum([0.0 - 0.01,minimum(lob)  - 0.1])
    ylim_ub = minimum([2*log(KTrue) + 0.01 , maximum(hib) + 0.1])
    # p = plot(collect(1:T),mean_vov)
    p = plot( collect(1:T),mean_vov,yerr = (mean_vov .- lob , hib .- mean_vov),grid =false,legend=false,title="Clustering Performance based on \n Variation of Information (VarInfo)\n over time series**", xlabel = "Time points", ylabel = "Variation of Information \n (VarInfo)", ylims = (ylim_lb,ylim_ub), seriestype = :scatter,color = :orange, m=(7.5,stroke(:black, 1.5))) #markersize = 5,
    # "Lower is better"
    return p
    
end

function plotPosteriorClusteringVarInfoAndSave!(true_z,z_post_s,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext;conf_level=0.95)
    dataPosteriorClusteringVarInfo = plotVarInfo(true_z,z_post_s;conf_level=conf_level)
    static_plot_id = prefix *"PosteriorClusteringVarInfo"
    dataPosteriorClusteringVarInfo_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(dataPosteriorClusteringVarInfo, dataPosteriorClusteringVarInfo_name);
    push!(plot_name_vec, dataPosteriorClusteringVarInfo_name)
end

function plotTimingsSubplot(plt, KMax_vec,y,series_name;sp = 1,label_font_size=7,legendfont_size=5,stroke_linewidth=2 )
    for i in 1:length(series_name)
        c = get(ColorSchemes.rainbow,i./length(series_name))
        plot!(plt,subplot=sp,KMax_vec,y[:,i],linewidth=stroke_linewidth,label=series_name[i],color=c, grid =false, yaxis=:log,legend=:outertop)   
    end
    
    plot!(plt,xlabel="Maximum Number of States (KMax)",xguidefont=font(label_font_size), legendfont=font(legendfont_size),subplot=sp)
    plot!(plt,ylabel="Log of Elapsed Time (sec)",yguidefont=font(label_font_size),subplot=sp) 

    return plt
end


function plotMetricsSubplot(results_dict,plt,sp_vec,KMax_vec,algo_str_vec,metric_str_vec,;conf_level = 0.95,label_font_size=7,legendfont_size=5,fill_alpha=0.95, marker_alpha=0.5,stroke_thickness=1)
    C_t = collect(results_dict["C_t"])
    tot_cells = sum(C_t)
    G  = results_dict["G"]
    T = results_dict["T"]
    mrk_shp_vec = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]
    for i in 1:length(sp_vec)
        sp = sp_vec[i]
        metric_str = metric_str_vec[i]
        ylim_lb_vec=[]
        ylim_ub_vec=[]
        for j in 1:length(algo_str_vec)
            algo_str = algo_str_vec[j]
            algo_metric_mat,algo_series_name,ylim_lb_vec,ylim_ub_vec = getMetrics_and_Stats(results_dict, algo_str,metric_str,KMax_vec,ylim_lb_vec,ylim_ub_vec;conf_level= 0.95)
            mrk_shp = mrk_shp_vec[j]
            for g in 1:length(algo_series_name)
                c = get(ColorSchemes.rainbow, g ./ length(algo_series_name))
                mean_vec = algo_metric_mat[:,1,g]
                lower_bound = algo_metric_mat[:,2,g]
                lower_err_bound = mean_vec .- lower_bound
                upper_bound = algo_metric_mat[:,3,g]
                upper_error_bound = upper_bound .- mean_vec
                label_name = algo_series_name[g]
                plot!(plt,subplot=sp,collect(1:T), mean_vec,yerr = (lower_err_bound ,upper_error_bound),alpha=fill_alpha, grid =false, linewidth=0.5,label=label_name,color=c, seriestype = :scatter,markershape=mrk_shp, m=(marker_alpha,stroke(:black, stroke_thickness)))
            end
        end
        plot_ylim_lb = minimum(ylim_lb_vec)
        plot_ylim_ub = maximum(ylim_ub_vec)
        plot!(plt,xlabel="Experimental Timepoints",xguidefont=font(label_font_size), legendfont=font(legendfont_size),legend=:outertop,subplot=sp)
        plot!(plt,ylabel="$(metric_str)",yguidefont=font(label_font_size) , ylims = (plot_ylim_lb,plot_ylim_ub),subplot=sp) 
       
    end

    return plt
end


