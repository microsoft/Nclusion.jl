function gnu_stackedhistogram_plot(tp,avg_counts_mat,KMax,KCalled,unique_clusters;filenamebase=nothing,fig_size=(800,600),color_order=[9,1,3,4,2,5,6,7,8,10],to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = colorschemes[:seaborn_colorblind][color_order]
    # Gnuplot.quitall()
    # @gp :- set output 'plot.jpg'


    pt_575_ll = "Patient 575 Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "set xrange [:$(KMax)]"
    # @gp :- "set ytics 10"
    @gp :- xlab = "Posterior Cluster Index" ylab = "Average Cluster Frequency"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    for t in tp:tp
        @gp :- subplt_indx title = "Time $t" "set size square" :-
        for l in 1:KCalled
            @gp :- avg_counts_mat[:,l,t,1] " using 1 t '$(unique_clusters[l])' lc rgb '#$(hex(linCol[l]))' fs transparent solid 0.75 " :- 
        end
        # @gp :- sum([avg_counts_mat[:,l,2,1] for l in 1:KCalled]) " using 1 t '$(pt_575_ll)' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :-
        # subplt_indx +=1
        @gp :- "set cbtics out nomirror" :-
        @gp :- "set key invert reverse Left outside"  :-
    end
    @gp

    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end

function gnu_allcellshistogram_plot(tp,avg_counts_mat,KMax,KCalled,pt_id;filenamebase=nothing,fig_size=(800,600), bar_color_hex_code="000000",to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = [bar_color_hex_code]
    # Gnuplot.quitall()
    # @gp :- set output 'plot.jpg'


    ll = "Patient $(pt_id) Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "set xrange [:$(KMax)]"
    # @gp :- "set ytics 10"
    @gp :- xlab = "Posterior Cluster Index" ylab = "Average Cluster Frequency"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    for t in tp:tp
        @gp :- subplt_indx title = "Time $t" "set size square" :-
        @gp :- sum([avg_counts_mat[:,l,2,1] for l in 1:KCalled]) " using 1 t '$(ll)' lc rgb '#$(linCol[1])' fs transparent solid 0.75 " :-
        subplt_indx +=1
        @gp :- "set cbtics out nomirror" :-
        @gp :- "set key invert reverse Left outside"  :-
    end
    @gp

    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end

function gnuplotcluster_geneimportanceweights(cluster_sumstats_df,k;ylim=10,filenamebase=nothing,fig_size=(800,600),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    num_topNGenes = length(cluster_sumstats_df[!,:Gene])
    # Gnuplot.quitall()
    # pos, vert = sort(abs.(randn(15))), LinRange(0.1,6,15)
    # poslow= sort(abs.(randn(15)))
    # poshigh= 5*poslow
    # poserr = abs.(randn(15)) #[[u,v] for (u,v) in zip(poslow,poshigh)]
    # lcval = 1:15

    yval_min=0.1
    ylim_max = ylim +0.5
    ylim_min = -0.05
    cluster_gene_names = string.(cluster_sumstats_df[!,:Gene])
    yval = sort(collect(LinRange(yval_min,ylim,num_topNGenes)),rev=true)
    tic_delta = mean([yval[i]-yval[i+1] for i in 1:length(yval)-1])
    pos = Float64.(cluster_sumstats_df[!,:average])
    poslow= Float64.(cluster_sumstats_df[!,:upper_bound])
    poshigh= Float64.(cluster_sumstats_df[!,:lower_bound])
    @gp "set ytics 0,$tic_delta,$(ylim+1)"
    @gp :- xlab = "Gene Importance Weight" ylab = "Genes in Cluster"
    @gp :- title = "Top $num_topNGenes By Importannce Weight in Cluster $k" 
    # @gp "set format x ''"
    @gp :- pos yval poslow poshigh cluster_gene_names "u 1:2:3:4:ytic(5) w xerrorbars t 'Gene Importance Weight (90% posterior interval)' lc 'black'" xrange= (0,maximum(poshigh)+1/2*maximum(poshigh)) yrange = (ylim_min,ylim_max)
    @gp :- "set key invert reverse Left bottom"
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end



function gnu_multiplot_1feature_type_binhistogram_plot(metadata,T, feature_column,KCalled,unique_clusters;filenamebase=nothing,fig_size=(800,600),color_order=[9,1,3,4,2,5,6,7,8,10],query_metadata=nothing,filter_label_function_dict=nothing,binsize = 0.1,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = colorschemes[:seaborn_colorblind][color_order]
    # Gnuplot.quitall()
    if isnothing(query_metadata)
        query_metadata = metadata
    end
    num_cols = ceil(sqrt(T));#round(sqrt(T))
    num_rows = round(sqrt(T));#ceil(sqrt(T))
    feature_name = replace(feature_column,"_"=>".")
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Count Histogram of $(feature_name) Across All Time Points' offset 0,0.1 rowsfirst scale 1.1,0.9"
    @gp :- "unset key"
    @gp :- xlab = "$(feature_name)" ylab = "Count"
    @gp :- "set xrange [$(minimum(query_metadata[!,feature_column])-1):$(maximum(query_metadata[!,feature_column])+1)]"
    subplt_indx = 1
    for t in 1:T
        metadata_to_use = get_metadata_at_timepoint(t,query_metadata, metadata)
        lowest_feature_val = minimum(metadata_to_use[!,feature_column])
        highest_feature_val = maximum(metadata_to_use[!,feature_column])

        for l in 1:KCalled
            label = unique_clusters[l]
            celltype_val = get_celltype_1feature_metadata(metadata_to_use, label,feature_column;filter_label_function_dict=filter_label_function_dict)
            c_celltype = hex(linCol[l])
            if length(celltype_val) == 1
                celltype_histbins = collect(lowest_feature_val:0.1:highest_feature_val)
                celltype_histcounts = (celltype_histbins .- 0.5*binsize) .< celltype_val[1] .<= (celltype_histbins .+ 0.5*binsize)
            elseif length(celltype_val) ==0
                celltype_histbins = collect(lowest_feature_val:0.1:highest_feature_val)
                celltype_histcounts = zeros(Float64,length(celltype_histbins))
            elseif length(unique(celltype_val)) ==1
                celltype_histbins = collect(lowest_feature_val:0.1:highest_feature_val)
                celltype_histcounts = (celltype_histbins .- 0.5*binsize) .< unique(celltype_val)[1] .<= (celltype_histbins .+ 0.5*binsize)
            else
                celltype_hist = hist(celltype_val, bs=binsize)
                celltype_histbins = celltype_hist.bins
                celltype_histcounts = celltype_hist.counts
            end
            @gp :- subplt_indx "set title 'Time Point $(t)'"
            @gp :- "set style fill transparent solid 0.5 noborder"
            if l == 1
                @gp :- celltype_histbins  celltype_histcounts  "with fillsteps tit '$(label)' lc '#$(c_celltype)' fs solid 1.0 noborder"
            else
                # @gp :- organoid_hist.bins  organoid_hist.counts  "with fillsteps tit 'Organoid' lc '$(c_organoid)' fs noborder"
                @gp :- celltype_histbins  celltype_histcounts  "with fillsteps tit '$(label)' lc '#$(c_celltype)' fs noborder"
            end
            @gp :- celltype_histbins  celltype_histcounts  "with steps tit '' lc '#$(c_celltype)' lw 2" "set grid"
        end
        subplt_indx +=1
    end
    @gp :- "set key vertical maxrows 1" :-
    @gp :- "set key at screen .85, screen .95"  :-
    @gp
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end

function gnu_multiplot_1feature_type_binhistogram_withmeans_plot(metadata,T, feature_column,KCalled,unique_clusters,avg_counts_mat,mean_μ_post,mk_hat_vec_;filenamebase=nothing,fig_size=(800,600),color_order=[9,1,3,4,2,5,6,7,8,10],query_metadata=nothing,filter_label_function_dict=nothing,binsize = 0.1,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = colorschemes[:seaborn_colorblind][color_order]
    Gnuplot.quitall()
    if isnothing(query_metadata)
        query_metadata = metadata
    end
    num_cols = ceil(sqrt(T));#round(sqrt(T))
    num_rows = round(sqrt(T));#ceil(sqrt(T))
    feature_name = replace(feature_column,"_"=>".")
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Count Histogram of $(feature_name) Across All Time Points' offset 0,0.1 rowsfirst scale 1.1,0.9"
    @gp :- "unset key"
    @gp :- xlab = "$(feature_name)" ylab = "Count"
    @gp :- "set xrange  [$(minimum(query_metadata[!,feature_column])-1):$(maximum(query_metadata[!,feature_column])+1)]"
    subplt_indx = 1
    for t in 1:T
        metadata_to_use = get_metadata_at_timepoint(t,query_metadata, metadata)
        celltype_histbins,celltype_histcounts,max_counts = extract_histbin_infomation(metadata_to_use, unique_clusters,feature_column;filter_label_function_dict=filter_label_function_dict,binsize = binsize)
        if !isnothing(avg_counts_mat) && !isnothing(mean_μ_post)
            non_zero_clusters_index = vec(sum(avg_counts_mat[:,:,t,1],dims=2) .>= 1.0);#vec(broadcast(!,iszero.(sum(avg_counts_mat[:,:,t,1],dims=2))))#
            non_zero_clusters = [el[1] for el in mean_μ_post[non_zero_clusters_index]]
            for clus in non_zero_clusters
                @gp :-  subplt_indx "set parametric"
                @gp :- "plot [t=0:$(max_counts+10)] $(clus),t  lc 'red' lw 2 linetype 2 dashtype 2"
                # @gp :- "plot $(clus) tit ''  with vectors nohead lc 'red' lw 2 linetype 2 dashtype 2"
            end
            # @gp :- "set tit 'Cluster Means'"
            # @gp :- non_zero_clusters "with vectors nohead tit 'Cluster Means' lc 'red' lw 2" 
        end
        for l in 1:KCalled
            label = unique_clusters[l]
            c_celltype = hex(linCol[l])
            # celltype_val = get_celltype_1feature_metadata(metadata_to_use, label,feature_column;filter_label_function_dict=filter_label_function_dict)
            
            @gp :- subplt_indx "set title 'Time Point $(t)'"
            @gp :- "set style fill transparent solid 0.5 noborder"
            # if l == 1
            #     @gp :- celltype_histbins[l]  celltype_histcounts[l]  "with fillsteps tit '$(label)' lc '#$(c_celltype)' fs solid 0.5 noborder"
            # else
            #     # @gp :- organoid_hist.bins  organoid_hist.counts  "with fillsteps tit 'Organoid' lc '$(c_organoid)' fs noborder"
            #     @gp :- celltype_histbins[l]  celltype_histcounts[l]  "with fillsteps tit '$(label)' lc '#$(c_celltype)' fs noborder"
            # end
            @gp :- celltype_histbins[l]  celltype_histcounts[l]  "with fillsteps tit '$(label)' lc '#$(c_celltype)' fs noborder"
            @gp :- celltype_histbins[l]  celltype_histcounts[l]  "with steps tit '' lc '#$(c_celltype)' lw 2" "set grid"
        end
        subplt_indx +=1
    end
    @gp :- "set key noautotitle"
    @gp :- "set key vertical maxrows 1" :- #
    @gp :- "set key at screen .85, screen .95"  :-
    @gp
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end
function gnu_plot_1gene_data_hist(data,genename;bs=0.5,fig_size=(800,600),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    @gp "set title 'Count Histogram of $genename'"
    @gp :- xlab = "$genename value" ylab = "Count"
    data_hist = hist(data,bs=bs)
    @gp :-  data_hist.bins  data_hist.counts "with fillsteps t '$genename value' lc 'black' fs solid 0.3 noborder"
    @gp :- data_hist.bins  data_hist.counts "with steps tit '' lc 'black' lw 2" "set grid"
end
function gnu_multiplot_1gene_data_hist(data,cell_labels,genename;bs=0.5, color_scheme = :tol_rainbow,color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256,fig_size=(800,600),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    cell_labels = join.(split.(cell_labels,"_")," ")
    color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
    linCol = colorschemes[color_scheme][color_order]
    # Gnuplot.options.term = "qt  size 1100,900"
    unique_labels = sort(unique(cell_labels))
    num_unique_labels = length(unique_labels)
    num_cols = ceil(sqrt(num_unique_labels));#round(sqrt(T))
    num_rows = round(sqrt(num_unique_labels));#ceil(sqrt(T))
    if num_cols == num_rows
        num_rows +=1
    end
    if num_cols == num_rows+1
        num_rows +=1
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Count Histogram of $genename' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- xlab = "$genename value" ylab = "Count"
    @gp :- "set xrange [$(minimum(data)+.1*minimum(data)):$(maximum(data)+.1*maximum(data))]"
    subplt_indx = 1
    for plt in 1:num_unique_labels+2
        # @gp :- "set title font ',16'"
        if plt == num_unique_labels+1
            sbplt_title = "Histogram for All Called Phenotypes"
            @gp :- plt title = sbplt_title  :- #"set size square"
            for lbl in 1:num_unique_labels
                data_hist = hist(data[cell_labels .== unique_labels[lbl]],bs=bs)
                @gp :-  data_hist.bins  data_hist.counts "with fillsteps t '$genename value for $(unique_labels[lbl]) ' lc rgb  '#$(hex(linCol[lbl]))' fs solid 0.3 noborder"
                @gp :- data_hist.bins  data_hist.counts "with steps tit '' lc rgb '#$(hex(linCol[lbl]))' lw 2" "set grid"
            end
        elseif plt == num_unique_labels+2
            sbplt_title = "Histogram for All Called Phenotypes"
            @gp :- plt title = sbplt_title  :- #"set size square
            data_hist = hist(data,bs=bs)
            @gp :-  data_hist.bins  data_hist.counts "with fillsteps t '$genename value for All Phenotypes' lc 'black' fs solid 0.3 noborder"
            @gp :- data_hist.bins  data_hist.counts "with steps tit '' lc 'black' lw 2" "set grid"
        else
            label = unique_labels[plt]
            sbplt_title = "Histogram for the Called Phenotype: \\n $label"
            @gp :- plt title = sbplt_title  :- #"set size square"
            data_hist = hist(data[cell_labels .== unique_labels[plt]],bs=bs)
            @gp :-  data_hist.bins  data_hist.counts "with fillsteps t '$genename value for $label ' lc rgb   '#$(hex(linCol[plt]))' fs solid 0.3 noborder"
            @gp :- data_hist.bins  data_hist.counts "with steps tit '' lc rgb   '#$(hex(linCol[plt]))' lw 2" "set grid"
        end

    end
    # @gp :- "set cbtics out nomirror" :-
    # @gp :- "set key font ',12'"
    # @gp :- "set key invert reverse Left outside"  :-
    @gp
end

function gnu_1feature_all_binhistogam_plot(all_to_use,feature_name)
    feature_name = replace(feature_name,"_"=>".")
    @gp "set title 'Count Histogram of $(feature_name) Across All Timepoints'"
    @gp :- xlab = "$(feature_name)" ylab = "Count"
    all_hist = hist(all_to_use)
    @gp :- all_hist.bins  all_hist.counts "with fillsteps t '$(feature_name) value' lc 'black' fs solid 0.3 noborder"
    @gp :- all_hist.bins  all_hist.counts "with steps tit '' lc 'black' lw 2" "set grid"
end
function gnu_multiplot_stackedhistogram_plot(avg_counts_mat,KMax,KCalled,unique_clusters;filenamebase=nothing,fig_size=(800,600),color_scheme = :tol_rainbow,color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29])
    # linCol = colorschemes[:seaborn_colorblind][color_order]
    # color_order=[9,1,3,4,2,5,6,7,8,10]
    
    linCol = nothing

    linCol_base= colorschemes[color_scheme][color_order]
    if KMax > length(linCol_base)
        linCol = linCol_base
        while KMax > length(linCol)
            append!(linCol,linCol_base)
        end
    else
        linCol = linCol_base
    end

    Gnuplot.quitall()
    # @gp :- set output 'plot.jpg'
    T = size(avg_counts_mat)[3]
    num_cols = ceil(sqrt(T));#round(sqrt(T))
    num_rows = round(sqrt(T));#ceil(sqrt(T))
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Average Posterior Cluster Count Across All Time Points' offset 0,0.1 rowsfirst scale 1.1,0.9"
    # @gp :- "set multiplot title offset 0,1 "
    pt_575_ll = "Patient 575 Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "unset key"
    @gp :- "set xrange [:$(KMax)]"
    
    # @gp :- "set ytics 10"
    @gp :- xlab = "Posterior Cluster Index" ylab = "Average \\n Cluster Count"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    posterior_cluster_labels = collect(1:KMax)#collect(1:2:KMax)
    for t in 1:T
        @gp :- subplt_indx title = "Time $t"  :- #"set size square"
        for l in 1:KCalled
            @gp :- avg_counts_mat[:,l,t,1] posterior_cluster_labels " using 1:xtic(2) t '$(unique_clusters[l])' lc rgb '#$(hex(linCol[l]))' fs transparent solid 0.75 " :- 
        end
        # @gp :- sum([avg_counts_mat[:,l,2,1] for l in 1:KCalled]) " using 1 t '$(pt_575_ll)' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :-
        subplt_indx +=1
        # @gp :- "set cbtics out nomirror" :-
        # @gp :- "set key invert reverse Left outside"  :-
    end
    @gp :- "set key vertical maxrows 1" :-
    @gp :- "set key at screen .85, screen .95"  :-
    @gp

    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end

function gnu_multiplot_allcellshistogram_plot(avg_counts_mat,KMax,KCalled,pt_id;filenamebase=nothing,fig_size=(800,600), bar_color_hex_code="000000",to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = [bar_color_hex_code]
    # Gnuplot.quitall()
    # @gp :- set output 'plot.jpg'

    T = size(avg_counts_mat)[3]
    num_cols = ceil(sqrt(T));#round(sqrt(T))
    num_rows = round(sqrt(T));#ceil(sqrt(T))
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Average Posterior Cluster Count Across All Time Points' offset 0,0.1 rowsfirst scale 1.1,0.9"
    ll = "Patient $(pt_id) Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "unset key"
    @gp :- "set xrange [:$(KMax)]"
    # @gp :- "set ytics 10"
    @gp :- xlab = "Posterior Cluster Index" ylab = "Average \\n Cluster Count"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    for t in 1:T
        @gp :- subplt_indx title = "Time $t"  :- #"set size square"
        @gp :- sum([avg_counts_mat[:,l,t,1] for l in 1:KCalled]) " using 1 t '$(ll)' lc rgb '#$(linCol[1])' fs transparent solid 0.75 " :-
        subplt_indx +=1
    end
    @gp :- "set key vertical maxrows 1" :-
    @gp :- "set key at screen .85, screen .95"  :-
    @gp

    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end


function gnu_plot_All_Time_Points_PCA(x,z,G,unique_clusters;unique_clusters_labels=nothing,filenamebase=nothing,fig_size=(1100,900),color_scheme = :tol_rainbow,color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],save_svg_copy=false,linCol=nothing,plt_name=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = length(unique_clusters)
    # linCol = nothing
    if isnothing(linCol)
        color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
        linCol_base = colorschemes[color_scheme][color_order]
        if KMax > length(linCol_base)
            linCol = linCol_base
            while KMax > length(linCol)
                append!(linCol,linCol_base)
            end
        else
            linCol = linCol_base
        end
    end
    # Gnuplot.quitall()
    
    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    data = hcat(Z,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    
    if isnothing(plt_name)
        plt_name = "PCA of Top $G Variably Expressed Genes in All Cells"
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 solid  noborder" #"set style fill  transparent solid 0.45 noborder" border lt .05 
    @gp :- "set style circle radius .20"
    @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+.1*minimum(permutedims(X_transformed)[:,1])):$(maximum(permutedims(X_transformed)[:,1])+.1*maximum(permutedims(X_transformed)[:,1]))]"
    @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+.1*minimum(permutedims(X_transformed)[:,2])):$(maximum(permutedims(X_transformed)[:,2])+.1*maximum(permutedims(X_transformed)[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- "set title font ',20'"
    @gp :- xlab = "PC1" ylab = "PC2"
    @gp :- title="$plt_name "
    if isnothing(unique_clusters_labels)
        unique_clusters_labels = unique_clusters
    end
    for l in 1:KMax
        clus = unique_clusters[l]
        @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " with circles t '$(unique_clusters_labels[l])' lc rgb '#$(hex(linCol[l]))' " :- 
        @gp :- "set cbtics out nomirror" :-
        @gp :- "set key font ',12'"
        @gp :- "set key invert reverse Left outside"  :-
    end
    @gp

    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
    return linCol
end

function gnu_plot_All_Time_Points_TSNE(x,z,G,unique_clusters;unique_clusters_labels=nothing,tsne_transform = nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),color_scheme = :tol_rainbow,color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],save_svg_copy=false,circle_rad=1,linCol = nothing,plt_name=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = length(unique_clusters)
    # linCol = nothing
    if isnothing(linCol)
        color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
        linCol_base = colorschemes[color_scheme][color_order]
        if KMax > length(linCol_base)
            linCol = linCol_base
            while KMax > length(linCol)
                append!(linCol,linCol_base)
            end
        else
            linCol = linCol_base
        end
    end
    # Gnuplot.quitall()

    rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())

    T = length(x)
    
    # vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T]))
    X = reduce(vcat,[permutedims(reduce(hcat,x[t])) for t in 1:T])#vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z = vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    K = Int(maximum(unique(Z)))
    clusterIDs =  sort(unique(Z))
    data = hcat(Z,X)
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end
    if isnothing(plt_name)
        plt_name = "TNSE of Top $G Variably Expressed Genes in All Cells"
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    plotTitle = "$plt_name "
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 solid  noborder" #"set style fill  transparent solid 0.45 noborder" border lt .05 
    @gp :- "set style circle radius $circle_rad"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- "set title font ',20'"
    @gp :- "set xrange [$(minimum(Y[:,1])+.1*minimum(Y[:,1])):$(maximum(Y[:,1])+.1*maximum(Y[:,1]))]"
    @gp :- "set yrange [$(minimum(Y[:,2])+.1*minimum(Y[:,2])):$(maximum(Y[:,2])+.1*maximum(Y[:,2]))]"
    @gp :- title = plotTitle
    @gp :- xlab = "TSNE1" ylab = "TSNE2"
    if isnothing(unique_clusters_labels)
        unique_clusters_labels = unique_clusters
    end
    for l in 1:KMax
        clus = unique_clusters[l]
        @gp :- Y[Z .==clus,:][:,1] Y[Z .==clus,:][:,2] "with circles t '$(unique_clusters_labels[l])' lc rgb '#$(hex(linCol[l]))' " :- 
        @gp :- "set cbtics out nomirror" :-
        @gp :- "set key font ',12'"
        @gp :- "set key invert reverse Left outside"  :-
        
        # with points t '$(unique_clusters_labels[l])' lc rgb '#$(hex(linCol[l]))' " :- 
        # @gp :- "set cbtics out nomirror" :-
        # @gp :- "set key invert reverse Left outside"  :-
    end
    @gp

    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
    return linCol
end

function gnu_multiplot_plot_All_Time_Points_TSNE(x,z_called,z_infer,G,unique_clusters_called,unique_clusters_infer;unique_clusters_labels_called=nothing,unique_clusters_labels_infer=nothing,tsne_transform = nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,500),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256,save_svg_copy=false,plt_name=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
    linCol_called = nothing
    linCol_infer = nothing

    color_order_called= vcat(color_order_called[1:4], shuffle(color_order_called[5:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end

    linCol_base_infer= colorschemes[color_scheme_infer]
    if KMax > length(linCol_base_infer)
        linCol_infer = linCol_base_infer
        while KMax > length(linCol_infer)
            append!(linCol_infer,linCol_base_infer)
        end
    else
        linCol_infer = linCol_base_infer
    end
    # Gnuplot.quitall()

    rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())

    T = length(x)
    
    # vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T]))
    X = reduce(vcat,[permutedims(reduce(hcat,x[t])) for t in 1:T])#vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
    # K = Int(maximum(unique(Z)))
    # clusterIDs =  sort(unique(Z))
    # data = hcat(Z,X)
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end
    if isnothing(plt_name)
        plt_name = "TSNE of Top $G Variably Expressed Genes in All Cells"
    end
    num_rows = 1
    num_cols= 2
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title '$plt_name' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set style fill  transparent solid 0.8 solid  noborder " #"set style fill  transparent solid 0.45 noborder border lt .05 " 
    @gp :- "set style circle radius 0.75"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set xrange [$(minimum(Y[:,1])+.1*minimum(Y[:,1])):$(maximum(Y[:,1])+.1*maximum(Y[:,1]))]"
    @gp :- "set yrange [$(minimum(Y[:,2])+.1*minimum(Y[:,2])):$(maximum(Y[:,2])+.1*maximum(Y[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "TSNE1" ylab = "TSNE2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    if isnothing(unique_clusters_labels_infer)
        unique_clusters_labels_infer = unique_clusters_infer
    end

    subplt_indx = 1
    for t in 1:2
        if t == 1
            sbplt_title = "Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        else
            sbplt_title = "Inferred Phenotypes" 
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = colorschemes[:glasbey_category10_n256]##linCol_infer
        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- Y[Z .==clus,:][:,1] Y[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',12'"
            @gp :- "set key invert reverse Left outside"  :-
            # with circles t '$(K_unique_clusters_labels_[l])' lc rgb '#$(hex(linCol[l]))' " :- 
            # @gp :- "set cbtics out nomirror" :-
            # @gp :- "set key invert reverse Left outside"  :-
        end
        @gp
        # @gp :- sum([avg_counts_mat[:,l,t,1] for l in 1:KCalled]) " using 1 t '$(ll)' lc rgb '#$(linCol[1])' fs transparent solid 0.75 " :-
        subplt_indx +=1
    end

    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
    

end
function gnu_multiplot_plot_All_Time_Points_PCA(x,z_called,z_infer,G,unique_clusters_called,unique_clusters_infer;unique_clusters_labels_called=nothing,unique_clusters_labels_infer=nothing,filenamebase=nothing,fig_size=(1100,500),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256,save_svg_copy=false,plt_name=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
    linCol_called = nothing
    linCol_infer = nothing

    color_order_called= vcat(color_order_called[1:4], shuffle(color_order_called[5:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end

    linCol_base_infer= colorschemes[color_scheme_infer]
    if KMax > length(linCol_base_infer)
        linCol_infer = linCol_base_infer
        while KMax > length(linCol_infer)
            append!(linCol_infer,linCol_base_infer)
        end
    else
        linCol_infer = linCol_base_infer
    end
    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
    data = hcat(Z_called,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(plt_name)
        plt_name = "PCA of Top $G Variably Expressed Genes in All Cells"
    end
    num_rows = 1
    num_cols= 2
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title '$plt_name' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set style fill  transparent solid 0.8 solid  noborder" #"set style fill  transparent solid 0.45 noborder" border lt .05  
    @gp :- "set style circle radius .20"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+.1*minimum(permutedims(X_transformed)[:,1])):$(maximum(permutedims(X_transformed)[:,1])+.1*maximum(permutedims(X_transformed)[:,1]))]"
    @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+.1*minimum(permutedims(X_transformed)[:,2])):$(maximum(permutedims(X_transformed)[:,2])+.1*maximum(permutedims(X_transformed)[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    if isnothing(unique_clusters_labels_infer)
        unique_clusters_labels_infer = unique_clusters_infer
    end

    subplt_indx = 1
    for t in 1:2
        if t == 1
            sbplt_title = "Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        else
            sbplt_title = "Inferred Phenotypes" 
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer
        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " with circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',12'"
            @gp :- "set key invert reverse Left outside"  :-
            
            # t '$(K_unique_clusters_labels_[l])' lc rgb '#$(hex(linCol[l]))' " :- 
            # @gp :- "set cbtics out nomirror" :-
            # @gp :- "set key invert reverse Left outside"  :-
        end
        @gp
        # @gp :- sum([avg_counts_mat[:,l,t,1] for l in 1:KCalled]) " using 1 t '$(ll)' lc rgb '#$(linCol[1])' fs transparent solid 0.75 " :-
        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end

end
function gnu_multiplot_plot_All_Time_Points_PCA_circles(x,z_called,z_infer,G,unique_clusters_called,unique_clusters_infer;unique_clusters_labels_called=nothing,unique_clusters_labels_infer=nothing,tsne_transform = nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256,annotate_biopsy=true,PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
    linCol_called = nothing
    linCol_infer = nothing

    color_order_called= vcat(color_order_called[1:4], shuffle(color_order_called[5:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end

    linCol_base_infer= colorschemes[color_scheme_infer]
    if KMax > length(linCol_base_infer)
        linCol_infer = linCol_base_infer
        while KMax > length(linCol_infer)
            append!(linCol_infer,linCol_base_infer)
        end
    else
        linCol_infer = linCol_base_infer
    end
    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
    data = hcat(Z_called,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    

    num_rows = 2
    num_cols= 1

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'PCA of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 0.2"
    @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+.1*minimum(permutedims(X_transformed)[:,1])):$(maximum(permutedims(X_transformed)[:,1])+.1*maximum(permutedims(X_transformed)[:,1]))]"
    @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+.1*minimum(permutedims(X_transformed)[:,2])):$(maximum(permutedims(X_transformed)[:,2])+.1*maximum(permutedims(X_transformed)[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    if isnothing(unique_clusters_labels_infer)
        unique_clusters_labels_infer = unique_clusters_infer
    end

    subplt_indx = 1
    for t in 1:2
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        else
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer
        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        # for l in 1:KCalled
        #     clus = K_unique_clusters_[l]
        #     @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " with p pt 7 '$(K_unique_clusters_labels_[l])' lc rgb '#$(hex(linCol[l]))' " :- 
        #     @gp :- "set cbtics out nomirror" :-
        #     @gp :- "set key invert reverse Left outside"  :-
        # end
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',12'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        # for l in 1:KCalled
        #     # clus = K_unique_clusters_[l]
        #     # 'Data1'
        #     @gp "  " :- 
        #     @gp :- "set cbtics out nomirror" :-
        #     @gp :- "set key invert reverse Left outside"  :-
        # end
        @gp
        # @gp :- sum([avg_counts_mat[:,l,t,1] for l in 1:KCalled]) " using 1 t '$(ll)' lc rgb '#$(linCol[1])' fs transparent solid 0.75 " :-
        subplt_indx +=1
    end
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end
function gnu_multiplot_plot_All_Time_Points_TSNE_circles(x,z_called,z_infer,G,unique_clusters_called,unique_clusters_infer;unique_clusters_labels_called=nothing,unique_clusters_labels_infer=nothing,tsne_transform = nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256,annotate_biopsy=true,PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
    linCol_called = nothing
    linCol_infer = nothing

    color_order_called= vcat(color_order_called[1:4], shuffle(color_order_called[5:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end

    linCol_base_infer= colorschemes[color_scheme_infer]
    if KMax > length(linCol_base_infer)
        linCol_infer = linCol_base_infer
        while KMax > length(linCol_infer)
            append!(linCol_infer,linCol_base_infer)
        end
    else
        linCol_infer = linCol_base_infer
    end
    # Gnuplot.quitall()

    T = length(x)
    
    # vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T]))
    X = reduce(vcat,[permutedims(reduce(hcat,x[t])) for t in 1:T])#vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
    # K = Int(maximum(unique(Z)))
    # clusterIDs =  sort(unique(Z))
    # data = hcat(Z,X)
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end
    

    num_rows = 2
    num_cols= 1

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'TSNE of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 0.2"
    @gp :- "set xrange [$(minimum(Y[:,1])+.1*minimum(Y[:,1])):$(maximum(Y[:,1])+.1*maximum(Y[:,1]))]"
    @gp :- "set yrange [$(minimum(Y[:,2])+.1*minimum(Y[:,2])):$(maximum(Y[:,2])+.1*maximum(Y[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    if isnothing(unique_clusters_labels_infer)
        unique_clusters_labels_infer = unique_clusters_infer
    end

    subplt_indx = 1
    for t in 1:2
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        else
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer
        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        # for l in 1:KCalled
        #     clus = K_unique_clusters_[l]
        #     @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " with p pt 7 '$(K_unique_clusters_labels_[l])' lc rgb '#$(hex(linCol[l]))' " :- 
        #     @gp :- "set cbtics out nomirror" :-
        #     @gp :- "set key invert reverse Left outside"  :-
        # end
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',12'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        # for l in 1:KCalled
        #     # clus = K_unique_clusters_[l]
        #     # 'Data1'
        #     @gp "  " :- 
        #     @gp :- "set cbtics out nomirror" :-
        #     @gp :- "set key invert reverse Left outside"  :-
        # end
        @gp
        # @gp :- sum([avg_counts_mat[:,l,t,1] for l in 1:KCalled]) " using 1 t '$(ll)' lc rgb '#$(linCol[1])' fs transparent solid 0.75 " :-
        subplt_indx +=1
    end
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end
function gnu_multiplot_plot_All_Time_Points_PCA_4Plots_circles(x,z_called,z_infer,G,unique_clusters_called;unique_clusters_labels_called=nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called]))
    linCol_called = nothing
    

    color_order_called= vcat(color_order_called[1:KMax], shuffle(color_order_called[KMax+1:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end
    # T = length(unique(xmat[:,1]));
    # N = length(unique(xmat[:,2]));
    # G = length(unique(xmat[:,3]));
    # timepoints = collect(1:T);
    # states_id = collect(1:K);
    # cell_ids = collect(1:N);
    # gene_ids = collect(1:G);
    # N_t = tidy_get_Nt_from_xmat(xmat);



    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
    data = hcat(Z_called,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')

    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    data = hcat(Z_called,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    

    num_rows = 2
    num_cols= 2

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes", "Biopsy Samples", "Timepoint Samples"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'PCA of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 0.3"
    @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+.1*minimum(permutedims(X_transformed)[:,1])):$(maximum(permutedims(X_transformed)[:,1])+.1*maximum(permutedims(X_transformed)[:,1]))]"
    @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+.1*minimum(permutedims(X_transformed)[:,2])):$(maximum(permutedims(X_transformed)[:,2])+.1*maximum(permutedims(X_transformed)[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    # unique_clusters_infer = nothing

    subplt_indx = 1
    for t in 1:4
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        elseif t == 2
            z_infer = z_infer
            color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Cluster $el" for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 3
            z_infer = get_biopsy_labels(z_called)
            color_scheme_infer = :ground_cover
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            color_order_=[1,14,2,3,4,5,6,7,8,9,10,11,12,13]
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end

            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [el ? " Biopsy" : "Organoid"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 4
            z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            color_scheme_infer = :viridis
            color_order_=[]
            push!(color_order_,1)
            tstart = 0
            for t in 2:T-1
                tmax = length(colorschemes[color_scheme_infer])-T-t
                if tstart == 2
                    tstart = 2
                end
                tindx = rand(collect(tstart:tmax))
                push!(color_order_,tindx)
                tstart = tindx
            end
            push!(color_order_,length(colorschemes[color_scheme_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Time $el"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',9'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        @gp

        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end


function gnu_multiplot_plot_All_Time_Points_TSNE_4Plots_circles(x,z_called,z_infer,G,unique_clusters_called;unique_clusters_labels_called=nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),tsne_transform = nothing, color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],PlotNames=nothing,to_display=true)
    # if to_display
    #     Gnuplot.quitall()
    #     Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    # else
    #     Gnuplot.quitall()
    #     Gnuplot.options.gpviewer = false
    # end
    KMax = maximum(length.([unique_clusters_called]))
    linCol_called = nothing
    

    color_order_called= vcat(color_order_called[1:KMax], shuffle(color_order_called[KMax+1:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end
    # T = length(unique(xmat[:,1]));
    # N = length(unique(xmat[:,2]));
    # G = length(unique(xmat[:,3]));
    # timepoints = collect(1:T);
    # states_id = collect(1:K);
    # cell_ids = collect(1:N);
    # gene_ids = collect(1:G);
    # N_t = tidy_get_Nt_from_xmat(xmat);


    rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())

    T = length(x)
    
    # vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T]))
    X = reduce(vcat,[permutedims(reduce(hcat,x[t])) for t in 1:T])#vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    # Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
    # K = Int(maximum(unique(Z)))
    # clusterIDs =  sort(unique(Z))
    # data = hcat(Z,X)
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end
    

    # data = hcat(Z_called,X)
    # M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    # X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')

    # Gnuplot.quitall()

    # T = length(x)

    # X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    # Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    # data = hcat(Z_called,X)
    # M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    # X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    

    num_rows = 2
    num_cols= 2

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes", "Biopsy Samples", "Timepoint Samples"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'TSNE Projection of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 1"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- "set xrange [$(minimum(Y[:,1])+.1*minimum(Y[:,1])):$(maximum(Y[:,1])+.1*maximum(Y[:,1]))]"
    @gp :- "set yrange [$(minimum(Y[:,2])+.1*minimum(Y[:,2])):$(maximum(Y[:,2])+.1*maximum(Y[:,2]))]"
    @gp :- xlab = "TSNE1" ylab = "TSNE2"


    # @gp "set size 1,1"  
    # @gp :- "set origin 0,0"
    # @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'PCA of $G Genes Expressed in Cells across All $T Time Points' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set grid y"
    # @gp :- "set grid x"
    # @gp :- "set style fill  transparent solid 0.45 noborder"
    # @gp :- "set style circle radius 0.2"
    # @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+1):$(maximum(permutedims(X_transformed)[:,1])+1)]"
    # @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+1):$(maximum(permutedims(X_transformed)[:,2])+1)]"
    # @gp :- "set xtics font ',11'"
    # @gp :- "set ytics font ',11'"
    # @gp :- "set xlabel font ',13'"
    # @gp :- "set ylabel font ',13'"
    # @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    # unique_clusters_infer = nothing

    subplt_indx = 1
    for t in 1:4
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        elseif t == 2
            z_infer = z_infer
            color_scheme_infer = :glasbey_category10_n256#:glasbey_bw_minc_20_minl_30_n256#:glasbey_bw_minc_20_n256#:glasbey_bw_n256#:dracula#:glasbey_hv_n256#
            unique_clusters_infer = sort(unique(vcat(z_infer...)), rev = true);
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Cluster $el" for el in sort(unique(vcat(z_infer...)), rev = true)]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 3
            z_infer = get_biopsy_labels(z_called)
            color_scheme_infer = :ground_cover
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            color_order_=[1,14,2,3,4,5,6,7,8,9,10,11,12,13]
            linCol_base_infer=[colorschemes[color_scheme_infer][11],colorschemes[:viridis][1]]#colorschemes[color_scheme_infer][color_order_]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end

            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [el ? " Biopsy" : "Organoid"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 4
            z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            color_scheme_infer = :roma
            color_order_=[]
            push!(color_order_,1)
            # tstart = 0
            # for t in 2:T-1
            #     tmax = length(colorschemes[color_scheme_infer])-T-t
            #     if tstart == 2
            #         tstart = 2
            #     end
            #     tindx = rand(collect(tstart:tmax))
            #     push!(color_order_,tindx)
            #     tstart = tindx
            # end
            # push!(color_order_,length(colorschemes[color_scheme_infer]))
            color_order_ = [(tt-1 )* floor(Int,(length(colorschemes[color_scheme_infer])-1)/(T-1)) + 1 for tt in 1:T]
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [isone(el) ? " Biopsy" : "Time $(el-1)"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            println("1")
            println("$KCalled")

            K_unique_clusters_ = unique_clusters_infer[[5,4,3,2,1]]
            println("2")
            K_unique_clusters_labels_ = unique_clusters_labels_infer[[5,4,3,2,1]]
            println("3")
            Z = Z_infer
            linCol = linCol_infer[vcat(collect(5:-1:1),collect(6:length(linCol_infer)))]

        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :-  Y[Z .==clus,:][:,1] Y[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',9'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        @gp

        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end


#:tol_rainbow
function gnu_tidy_multiplot_plot_All_Time_Points_PCA_ToyData_circles(xmat,zmat_called,zmat_infer,unique_clusters_called;unique_clusters_labels_called=nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called]))
    linCol_called = nothing
    

    color_order_called= vcat(color_order_called[1:KMax], shuffle(color_order_called[KMax+1:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    # timepoints = collect(1:T);
    # states_id = collect(1:K);
    # cell_ids = collect(1:N);
    # gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);


    x_vec = [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    x = [x_vec[st:en] for (st,en) in timeranges]
    z_1 = zmat_called[:,end]
    z_called = [z_1[st:en] for (st,en) in timeranges]
    z_2 = zmat_infer[:,end]
    z_infer = [z_2[st:en] for (st,en) in timeranges]

    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    data = hcat(Z_called,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    

    num_rows = 2
    num_cols= 2

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes", "Biopsy Samples", "Timepoint Samples"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'PCA of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 0.2"
    @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+.1*minimum(permutedims(X_transformed)[:,1])):$(maximum(permutedims(X_transformed)[:,1])+.1*maximum(permutedims(X_transformed)[:,1]))]"
    @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+.1*minimum(permutedims(X_transformed)[:,2])):$(maximum(permutedims(X_transformed)[:,2])+.1*maximum(permutedims(X_transformed)[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    # unique_clusters_infer = nothing

    subplt_indx = 1
    for t in 1:4
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        elseif t == 2
            z_infer = z_infer
            color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Cluster $el" for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 3
            z_infer = get_biopsy_labels(z_called)
            color_scheme_infer = :ground_cover
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            color_order_=[1,14,2,3,4,5,6,7,8,9,10,11,12,13]
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end

            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [el ? " Biopsy" : "Organoid"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 4
            # z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # unique_clusters_infer = sort(unique(vcat(z_infer...)));
            # sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            # linCol_infer = nothing
            # KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            
            z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            color_scheme_infer = :viridis
            color_order_=[]
            push!(color_order_,1)
            tstart = 0
            for t in 2:T-1
                tmax = length(colorschemes[color_scheme_infer])-T-t
                if tstart == 2
                    tstart = 2
                end
                tindx = rand(collect(tstart:tmax))
                push!(color_order_,tindx)
                tstart = tindx
            end
            push!(color_order_,length(colorschemes[color_scheme_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))


            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Time $el"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',9'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        @gp

        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end

function gnu_tidy_multiplot_plot_All_Time_Points_TSNE_ToyData_circles(xmat,zmat_called,zmat_infer,unique_clusters_called;unique_clusters_labels_called=nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),tsne_transform = nothing,color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called]))
    linCol_called = nothing
    

    color_order_called= vcat(color_order_called[1:KMax], shuffle(color_order_called[KMax+1:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end
    
    rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    # timepoints = collect(1:T);
    # states_id = collect(1:K);
    # cell_ids = collect(1:N);
    # gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);


    x_vec = [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    x = [x_vec[st:en] for (st,en) in timeranges]
    z_1 = zmat_called[:,end]
    z_called = [z_1[st:en] for (st,en) in timeranges]
    z_2 = zmat_infer[:,end]
    z_infer = [z_2[st:en] for (st,en) in timeranges]

    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    # data = hcat(Z_called,X)
    # M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    # X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end

    num_rows = 2
    num_cols= 2

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes", "Biopsy Samples", "Timepoint Samples"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'TSNE of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 1"
    @gp :- "set xrange [$(minimum(Y[:,1])+.1*minimum(Y[:,1])):$(maximum(Y[:,1])+.1*maximum(Y[:,1]))]"
    @gp :- "set yrange [$(minimum(Y[:,2])+.1*minimum(Y[:,2])):$(maximum(Y[:,2])+.1*maximum(Y[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "TSNE1" ylab = "TSNE2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    # unique_clusters_infer = nothing

    subplt_indx = 1
    for t in 1:4
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        elseif t == 2
            z_infer = z_infer
            color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Cluster $el" for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 3
            z_infer = get_biopsy_labels(z_called)
            color_scheme_infer = :ground_cover
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            color_order_=[1,14,2,3,4,5,6,7,8,9,10,11,12,13]
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end

            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [el ? " Biopsy" : "Organoid"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 4
            # z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # unique_clusters_infer = sort(unique(vcat(z_infer...)));
            # sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            # linCol_infer = nothing
            # KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            
            z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            color_scheme_infer = :viridis
            color_order_=[]
            push!(color_order_,1)
            tstart = 0
            for t in 2:T-1
                tmax = length(colorschemes[color_scheme_infer])-T-t
                if tstart == 2
                    tstart = 2
                end
                tindx = rand(collect(tstart:tmax))
                push!(color_order_,tindx)
                tstart = tindx
            end
            push!(color_order_,length(colorschemes[color_scheme_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))


            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Time $el"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :-  Y[Z .==clus,:][:,1] Y[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',9'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        @gp

        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end

function gnu_tidy_multiplot_plot_All_Time_Points_PCA_4Plots(xmat,zmat_called,zmat_infer,unique_clusters_called;unique_clusters_labels_called=nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called]))
    linCol_called = nothing
    

    color_order_called= vcat(color_order_called[1:KMax], shuffle(color_order_called[KMax+1:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    # timepoints = collect(1:T);
    # states_id = collect(1:K);
    # cell_ids = collect(1:N);
    # gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);


    x_vec = [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    x = [x_vec[st:en] for (st,en) in timeranges]
    z_1 = zmat_called[:,end]
    z_called = [z_1[st:en] for (st,en) in timeranges]
    z_2 = zmat_infer[:,end]
    z_infer = [z_2[st:en] for (st,en) in timeranges]

    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    data = hcat(Z_called,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    

    num_rows = 2
    num_cols= 2

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes", "Biopsy Samples", "Timepoint Samples"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'PCA of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 0.2"
    @gp :- "set xrange [$(minimum(permutedims(X_transformed)[:,1])+.1*minimum(permutedims(X_transformed)[:,1])):$(maximum(permutedims(X_transformed)[:,1])+.1*maximum(permutedims(X_transformed)[:,1]))]"
    @gp :- "set yrange [$(minimum(permutedims(X_transformed)[:,2])+.1*minimum(permutedims(X_transformed)[:,2])):$(maximum(permutedims(X_transformed)[:,2])+.1*maximum(permutedims(X_transformed)[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "PCA1" ylab = "PCA2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    # unique_clusters_infer = nothing

    subplt_indx = 1
    for t in 1:4
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        elseif t == 2
            z_infer = z_infer
            color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Cluster $el" for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 3
            z_infer = get_biopsy_labels(z_called)
            color_scheme_infer = :ground_cover
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            color_order_=[1,14,2,3,4,5,6,7,8,9,10,11,12,13]
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end

            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [el ? " Biopsy" : "Organoid"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 4
            # z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # unique_clusters_infer = sort(unique(vcat(z_infer...)));
            # sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            # linCol_infer = nothing
            # KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            
            z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            color_scheme_infer = :viridis
            color_order_=[]
            push!(color_order_,1)
            tstart = 0
            for t in 2:T-1
                tmax = length(colorschemes[color_scheme_infer])-T-t
                if tstart == 2
                    tstart = 2
                end
                tindx = rand(collect(tstart:tmax))
                push!(color_order_,tindx)
                tstart = tindx
            end
            push!(color_order_,length(colorschemes[color_scheme_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))


            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Time $el"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :- permutedims(X_transformed)[Z .==clus,:][:,1] permutedims(X_transformed)[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " :- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" :-
            @gp :- "set key font ',9'"
            @gp :- "set key invert reverse Left outside"  :-
        end
        @gp

        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end

function gnu_tidy_multiplot_plot_All_Time_Points_TSNE_4Plots(xmat,zmat_called,zmat_infer,unique_clusters_called;unique_clusters_labels_called=nothing, todisplay=true,filenamebase=nothing,fig_size=(1100,900),tsne_transform = nothing,color_scheme_called = :tol_rainbow,color_order_called=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],PlotNames=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    KMax = maximum(length.([unique_clusters_called]))
    linCol_called = nothing
    

    color_order_called= vcat(color_order_called[1:KMax], shuffle(color_order_called[KMax+1:end]))
    linCol_base_called = colorschemes[color_scheme_called][color_order_called]
    if KMax > length(linCol_base_called)
        linCol_called = linCol_base_called
        while KMax > length(linCol_called)
            append!(linCol_called,linCol_base_called)
        end
    else
        linCol_called = linCol_base_called
    end
    
    rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    # timepoints = collect(1:T);
    # states_id = collect(1:K);
    # cell_ids = collect(1:N);
    # gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);


    x_vec = [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    x = [x_vec[st:en] for (st,en) in timeranges]
    z_1 = zmat_called[:,end]
    z_called = [z_1[st:en] for (st,en) in timeranges]
    z_2 = zmat_infer[:,end]
    z_infer = [z_2[st:en] for (st,en) in timeranges]

    # Gnuplot.quitall()

    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    Z_called = vec(reduce(vcat,[permutedims(reduce(hcat,z_called[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    # data = hcat(Z_called,X)
    # M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    # X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end

    num_rows = 2
    num_cols= 2

    if isnothing(PlotNames)
        PlotNames = ["Called Phenotypes"  , "Inferred Phenotypes", "Biopsy Samples", "Timepoint Samples"] #[ "Plot $i" for i in 1:num_rows*num_cols]
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'TSNE of Top $G Variably Expressed Genes in All Cells' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill  transparent solid 0.8 noborder"#"set style fill  transparent solid 0.45 noborder border lt .05"
    @gp :- "set style circle radius 1"
    @gp :- "set xrange [$(minimum(Y[:,1])+.1*minimum(Y[:,1])):$(maximum(Y[:,1])+.1*maximum(Y[:,1]))]"
    @gp :- "set yrange [$(minimum(Y[:,2])+.1*minimum(Y[:,2])):$(maximum(Y[:,2])+.1*maximum(Y[:,2]))]"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- xlab = "TSNE1" ylab = "TSNE2"
    if isnothing(unique_clusters_labels_called)
        unique_clusters_labels_called = unique_clusters_called
    end
    # unique_clusters_infer = nothing

    subplt_indx = 1
    for t in 1:4
        if t == 1
            sbplt_title =  PlotNames[t] #"Called Phenotypes" 
            KCalled = length(unique_clusters_called)
            K_unique_clusters_ = unique_clusters_called
            K_unique_clusters_labels_ = unique_clusters_labels_called
            Z = Z_called
            linCol = linCol_called
        elseif t == 2
            z_infer = z_infer
            color_scheme_infer = :glasbey_bw_minc_20_maxl_70_n256
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = ["Cluster $el" for el in sort(unique(vcat(z_infer...)), rev = true)]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 3
            z_infer = get_biopsy_labels(z_called)
            color_scheme_infer = :ground_cover #11
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            color_order_=[1,14,2,3,4,5,6,7,8,9,10,11,12,13]
            linCol_base_infer= [colorschemes[:viridis][1], colorschemes[color_scheme_infer][11]]
            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end

            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [el ? " Biopsy" : "Organoid"  for el in sort(unique(vcat(z_infer...)))]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        elseif t == 4
            # z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # unique_clusters_infer = sort(unique(vcat(z_infer...)));
            # sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            # linCol_infer = nothing
            # KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            
            z_infer = get_timepoint_labels(z_called)
            # color_scheme_infer = :tol_rainbow
            # color_order_=vcat(collect(10:length(colorschemes[color_scheme_infer])), collect(9:-1:1))
            # linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            color_scheme_infer = :viridis
            color_order_=[]
            push!(color_order_,1)
            tstart = 0
            for t in 2:T-1
                tmax = length(colorschemes[color_scheme_infer])-T-t
                if tstart == 2
                    tstart = 2
                end
                tindx = rand(collect(tstart:tmax))
                push!(color_order_,tindx)
                tstart = tindx
            end
            push!(color_order_,length(colorschemes[color_scheme_infer]))
            linCol_base_infer= colorschemes[color_scheme_infer][color_order_]
            unique_clusters_infer = sort(unique(vcat(z_infer...)));
            sbplt_title = PlotNames[t] #"Inferred Phenotypes" 
            linCol_infer = nothing
            KMax = maximum(length.([unique_clusters_called,unique_clusters_infer]))


            if KMax > length(linCol_base_infer)
                linCol_infer = linCol_base_infer
                while KMax > length(linCol_infer)
                    append!(linCol_infer,linCol_base_infer)
                end
            else
                linCol_infer = linCol_base_infer
            end
            Z_infer = vec(reduce(vcat,[permutedims(reduce(hcat,z_infer[t])) for t in 1:T])) 
            unique_clusters_labels_infer = [isone(el) ? " Biopsy" : "Time $(el-1)"  for el in sort(unique(vcat(z_infer...)),rev=true)]
            KCalled = length(unique_clusters_infer)
            K_unique_clusters_ = unique_clusters_infer
            K_unique_clusters_labels_ = unique_clusters_labels_infer
            Z = Z_infer
            linCol = linCol_infer

        end
        @gp :- subplt_indx title = sbplt_title  :- #"set size square"
        @gp :- "set title font ',16'"
        for l in 1:KCalled
            clus = K_unique_clusters_[l]
            @gp :-  Y[Z .==clus,:][:,1] Y[Z .==clus,:][:,2] " w circles notitle lc rgb '#$(hex(linCol[l]))'  " #:- 
            @gp :- "plot keyentry w circles fs solid fc rgb '#$(hex(linCol[l]))'  title '$(K_unique_clusters_labels_[l])'  " #p pt 7        #circles fs solid fc 'dark-red' title "circles"
            @gp :- "set cbtics out nomirror" #:-
            @gp :- "set key font ',9'"
            @gp :- "set key invert reverse Left outside"  #:-
        end
        @gp

        subplt_indx +=1
    end
    
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end

end

function gnu_multiplot_stackedhistogram_plot2(avg_counts_mat,G,KMax,KCalled,unique_clusters;filenamebase=nothing,fig_size=(800,600),color_scheme = :tol_rainbow,color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = nothing

    linCol_base= colorschemes[color_scheme][color_order]
    if KMax > length(linCol_base)
        linCol = linCol_base
        while KMax > length(linCol)
            append!(linCol,linCol_base)
        end
    else
        linCol = linCol_base
    end

    # Gnuplot.quitall()
    # @gp :- set output 'plot.jpg'
    T = size(avg_counts_mat)[3]
    num_cols = ceil(sqrt(T));#round(sqrt(T))
    num_rows = round(sqrt(T));#ceil(sqrt(T))
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Average Posterior Cluster Count Across All Time Points ($G Genes)' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set multiplot title offset 0,1 "
    pt_575_ll = "Patient 575 Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "unset key"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set xrange [:$(KMax)]"
    
    # @gp :- "set ytics 10"
    @gp :- xlab = "Posterior Cluster Index" ylab = "Average \\n Cluster Count"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    posterior_cluster_labels = collect(1:KMax)#collect(1:2:KMax)
    for t in 1:T
        @gp :- subplt_indx title = "Time $t"  :- #"set size square"
        @gp :- "set title font ',15'"
        for l in 1:KCalled
            @gp :- avg_counts_mat[:,l,t,1] posterior_cluster_labels " using 1:xtic(2) t '$(unique_clusters[l])' lc rgb '#$(hex(linCol[l]))' fs transparent solid 0.75 " :- 
        end
        # @gp :- sum([avg_counts_mat[:,l,2,1] for l in 1:KCalled]) " using 1 t '$(pt_575_ll)' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :-
        subplt_indx +=1
        # @gp :- "set cbtics out nomirror" :-
        # @gp :- "set key invert reverse outside below"  :-
        # @gp :- "set key font ',12'"
        # @gp :- "set key vertical maxrows 1" :-
    end
    # @gp :- "set key vertical maxrows 1" :-
    # @gp :- "set key font ',12'" :-

    # @gp :- "set key at screen .85, screen 0.35"  :-
    @gp :- "set key noautotitle"
    @gp :- "set key vertical maxrows 1 font ',12'" :- #
    @gp :- "set key at screen .75, screen .97"  :-
    @gp

    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end

function gnu_multiplot_stackedhistogram_plot_proportions(avg_counts_mat,G,KMax,KCalled,unique_clusters;filenamebase=nothing,fig_size=(800,600),color_scheme = :tol_rainbow,color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29],to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    linCol = nothing

    linCol_base= colorschemes[color_scheme][color_order]
    if KMax > length(linCol_base)
        linCol = linCol_base
        while KMax > length(linCol)
            append!(linCol,linCol_base)
        end
    else
        linCol = linCol_base
    end

    # Gnuplot.quitall()
    # @gp :- set output 'plot.jpg'
    T = size(avg_counts_mat)[3]
    num_cols = ceil(sqrt(T));#round(sqrt(T))
    num_rows = round(sqrt(T));#ceil(sqrt(T))
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Average Posterior Cluster Count Across All Time Points ($G Genes)' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set multiplot title offset 0,1 "
    pt_575_ll = "Patient 575 Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "unset key"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set xrange [:$(KMax)]"
    KCalled_totals = vec(sum(sum(avg_counts_mat,dims = 3), dims = 1))
    # @gp :- "set ytics 10"
    @gp :- xlab = "Posterior Cluster Index" ylab = "Average \\n Cluster Count"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    posterior_cluster_labels = collect(1:KMax)#collect(1:2:KMax)
    for t in 1:T
        @gp :- subplt_indx title = "Time $t"  :- #"set size square"
        @gp :- "set title font ',15'"
        for l in 1:KCalled
            @gp :- avg_counts_mat[:,l,t,1] ./KCalled_totals[l]  posterior_cluster_labels " using 1:xtic(2) t '$(unique_clusters[l])' lc rgb '#$(hex(linCol[l]))' fs transparent solid 0.75 " :- 
        end
        # @gp :- sum([avg_counts_mat[:,l,2,1] for l in 1:KCalled]) " using 1 t '$(pt_575_ll)' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :-
        subplt_indx +=1
        # @gp :- "set cbtics out nomirror" :-
        # @gp :- "set key invert reverse outside below"  :-
        # @gp :- "set key font ',12'"
        # @gp :- "set key vertical maxrows 1" :-
    end
    # @gp :- "set key vertical maxrows 1" :-
    # @gp :- "set key font ',12'" :-

    # @gp :- "set key at screen .85, screen 0.35"  :-
    @gp :- "set key noautotitle"
    @gp :- "set key vertical maxrows 1 font ',12'" :- #
    @gp :- "set key at screen .75, screen .97"  :-
    @gp

    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end




function gnu_elbo_plot(elbo;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false, num_iterations=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    if isnothing(num_iterations)
        num_iterations = collect(1:length(elbo))
    end
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xrange [:$(length(elbo)+1)]"
    # @gp :- "set ytics 10"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set key font ',12'"
    @gp :- "set title font ',15'"
    @gp :- "set key invert reverse Left outside" 
    @gp :- num_iterations elbo  "w linespoints lw 1 lc 'blue' pt 7 t 'ELBO'"
    @gp :- xlab = "Elapsed Global Iterations" ylab = "ELBO"
    @gp :- title = "ELBO over Algorithm Iterations"

    @gp

    if !isnothing(filenamebase)
        plotID = "plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_elbovsK_plot(elbo_vec,Krange;filenamebase=nothing,fig_size=(1100,900),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    @gp "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xrange [:$(maximum(Krange)+1)]"
    # @gp :- "set ytics 10"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set key font ',12'"
    @gp :- "set title font ',15'"
    @gp :- "set key invert reverse Left outside" 
    @gp :- Krange elbo_vec  "w linespoints lw 1 lc 'blue' pt 7 t 'ELBO'"
    @gp :- xlab = "Maximum K Intialization" ylab = "Final ELBO"
    @gp :- title = "Final ELBO For different Ks "

    @gp

    if !isnothing(filenamebase)
        plotID = "plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end
function gnu_ARIvsK_plot(ari_vec,Krange;filenamebase=nothing,fig_size=(1100,900),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xrange [:$(maximum(Krange)+1)]"
    # @gp :- "set ytics 10"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set key font ',12'"
    @gp :- "set title font ',15'"
    @gp :- "set key invert reverse Left outside" 
    @gp :- Krange ari_vec  "w linespoints lw 1 lc 'red' pt 7 t 'ARI'"
    @gp :- xlab = "Maximum K Intialization" ylab = "Overall ARI"
    @gp :- title = "Overall ARI For different Ks "

    @gp

    if !isnothing(filenamebase)
        plotID = "plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
        end
    end
end

function gnu_multiplot_plot_inclusion_probabilities(gene_significance_weights;gene_labels=nothing,cluster_labels=nothing,color_order=nothing,linCol=nothing ,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,clusters_used_annotation=nothing,driver_genes_annotations=nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    # Gnuplot.options.gpviewer = true
    G = length(gene_significance_weights[1])
    K = length(gene_significance_weights)
    gene_significance_weights =[ [va <= 1/(10*G) ? 0.0 : va  for va in el] for el in gene_significance_weights]
    to_include_bool = [!(any(isnan.(el)) || all(iszero.(el))) for el in gene_significance_weights]
    gene_significance_weights2 = gene_significance_weights[to_include_bool]
    K_valid = sum(to_include_bool)
    K_valid_vec = collect(1:K)[to_include_bool]
    data = permutedims(reduce(hcat,gene_significance_weights2))
    # data =hcat([normToProb([data[k,j] for k in 1:K_valid]) for j in 1:G]...)
    clusters_used_annotation_valid = nothing
    if !isnothing(clusters_used_annotation)
        clusters_used_annotation_valid = clusters_used_annotation[to_include_bool]
    end
    if isnothing(cluster_labels)
        cluster_labels = ["Cluster $k" for k in 1:K]
    end
    new_cluster_labels = cluster_labels[to_include_bool]
    
    
    # cluster_labels = ["Cluster $k" for k in 1:K]
    if isnothing(linCol)
        # color_scheme = :tol_rainbow;
        # if isnothing(color_order)
        #     color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
        #     color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
        # end
        # linCol = colorschemes[color_scheme][color_order]
        color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
        linCol = colorschemes[color_scheme][color_order] 
    end
    num_cols = ceil(sqrt(K_valid));#round(sqrt(T))
    num_rows = round(sqrt(K_valid));#ceil(sqrt(T))

    # driver_gene_tic_linCol = colorschemes[:glasbey_bw_minc_20_hue_330_100_n256]
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title 'Posterior Inclusion Probabilities for Genes in Inferred Clusters' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set multiplot title offset 0,1 "
    pt_575_ll = "Patient 575 Organoid cells"
    @gp :- "set grid y"
    @gp :- "set style data histograms"
    # @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "unset key"
    @gp :- "set xtics font ',7'"
    @gp :- "set ytics font ',7'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',10'"
    @gp :- "set ylabel font ',10'"
    @gp :- "set xrange [:$(G)]"
    @gp :- "set yrange [0:$(maximum(maximum.(data)) + 1/G * maximum(maximum.(data)) )]"
    # @gp :- "set yrange [] writeback "
    # @gp :- "set ytics 10"
    @gp :- xlab = "Gene" ylab = "PIPs"#"Posterior Inclusion Probability (PIP^k_{j})"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    if isnothing(gene_labels)
        gene_labels = [["Gene $j" for j in 1:G] for k in 1:K]#collect(1:2:KMax)
    end

    for k in 1:K_valid
        @gp :- subplt_indx title = new_cluster_labels[k]  :- #"set size square"
        @gp :- "set title font ',15'"
        if !isnothing(clusters_used_annotation_valid)
            if clusters_used_annotation_valid[k]
                @gp :- "set title tc lt 7" :-
            
            else
                @gp :- "set title tc def" :-
            end
        end 

        # for j in 1:G
        #     # @gp :- data[k][j] gene_labels " using 1:xtic(2) t '$(gene_labels[j])' lc rgb '#$(hex(linCol[j]))' fs transparent solid 0.75 " :- 
            
        # end
        
        @gp :- data[k,:] gene_labels[k] " using 1:xtic(2) t '$(new_cluster_labels[k])' lc rgb '#$(hex(linCol[k]))' fs transparent solid 0.75 " :-
        
        if !isnothing(clusters_used_annotation_valid)
            if !isnothing(driver_genes_annotations) && clusters_used_annotation_valid[k]
                num_causal_clusters = length(driver_genes_annotations)
                for t in 1:num_causal_clusters
                    for j in 1:G
                        if driver_genes_annotations[t][j]
                            @gp :- "set xtics $j,$j textcolor lt 7" :-# ('Gene $j' $j)  textcolor lt 7
                            # println("here $t,$j")
                        else                        
                            @gp :- "set xtics $j,$j tc def" :- # @gp :- "set xtics ('Gene $j' $j)  tc def" # 
                        end
                    end
                end
            else
                @gp :- "set xtics tc def" :-
            end
        end 
    
        # @gp :- sum([avg_counts_mat[:,l,2,1] for l in 1:KCalled]) " using 1 t '$(pt_575_ll)' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :-
        subplt_indx +=1
        # @gp :- "set cbtics out nomirror" :-
        # @gp :- "set key invert reverse outside below"  :-
        # @gp :- "set key font ',12'"
        # @gp :- "set key vertical maxrows 1" :-
    end
    # @gp :- "set key vertical maxrows 1" :-
    # @gp :- "set key font ',12'" :-

    # @gp :- "set key at screen .85, screen 0.35"  :-
    # @gp :- "set key noautotitle"
    # @gp :- "set key vertical maxrows 1 font ',12'" :- #
    # @gp :- "set key at screen .75, screen .97"  :-
    
    if to_display
        @gp 
    end
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
    return linCol
end

function gnu_multiplot_plot_inclusion_probabilities(gene_significance_weights,topGgenes,gene_names;color_order=nothing,cluster_labels=nothing,linCol=nothing ,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,clusters_used_annotation=nothing,to_display=true,to_include_bool=nothing)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    G = length(gene_significance_weights[1])
    K = length(gene_significance_weights)
    gene_significance_weights =[ [va <= 1/(10*G) ? 0.0 : va  for va in el] for el in gene_significance_weights]
    if isnothing(to_include_bool)
        to_include_bool = [!(any(isnan.(el)) || all(iszero.(el))) for el in gene_significance_weights]
    end
    new_gene_significance_weights = gene_significance_weights[to_include_bool]
    # new_K = sum(to_include_bool)
    new_K_vec = collect(1:K)[to_include_bool]
    new_K = length(new_K_vec)
    data = permutedims(reduce(hcat,new_gene_significance_weights))
    # data =hcat([normToProb([data[k,j] for k in 1:K_valid]) for j in 1:G]...)
    clusters_used_annotation_valid = nothing
    if !isnothing(clusters_used_annotation)
        clusters_used_annotation_valid = clusters_used_annotation[to_include_bool]
    end
    if isnothing(cluster_labels)
        cluster_labels = ["Cluster $k" for k in 1:K]
    end
    new_cluster_labels = cluster_labels[to_include_bool]
    
   
    # new_gene_significance_weights = gene_significance_weights[to_include_bool]
    # new_K_vec = collect(1:K)[to_include_bool]
    

    

    if isnothing(linCol)
        # color_scheme = :tol_rainbow;
        # if isnothing(color_order)
        #     color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
        #     color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
        # end
        # linCol = colorschemes[color_scheme][color_order]
        color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
        linCol = colorschemes[color_scheme][color_order] 
    end
    num_cols = ceil(sqrt(new_K));#round(sqrt(T))
    num_rows = round(sqrt(new_K));#ceil(sqrt(T))
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $num_rows,$num_cols title 'Posterior Inclusion Probabilities for top $topGgenes Genes in Inferred Clusters' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set multiplot title offset 0,1 "

    @gp :- "set grid y"
    @gp :- "set style data histograms"
    # @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "unset key"
    @gp :- "set xtics font ',7'"
    @gp :- "set xtics rotate"
    @gp :- "set ytics font ',7'"
    @gp :- "set xlabel font ',7'"
    @gp :- "set ylabel font ',7'"
    @gp :- "set xrange [:$(topGgenes)]"
    @gp :- "set yrange [0:$(maximum(maximum.(data)) + 1/G * maximum(maximum.(data)) )]"
    # @gp :- "set ytics 10"
    @gp :- xlab = "Gene" ylab = "PIPs"#"Posterior Inclusion Probability (PIP^k_{j})"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    

    for k in 1:new_K
        @gp :- subplt_indx title = "Cluster $(new_K_vec[k])"  :- #"set size square"
        @gp :- "set title font ',15'"
        # for j in 1:G
        #     # @gp :- data[k][j] gene_labels " using 1:xtic(2) t '$(gene_labels[j])' lc rgb '#$(hex(linCol[j]))' fs transparent solid 0.75 " :- 
            
        # end
        top_Ggenes_bool = sortperm(new_gene_significance_weights[k],rev=true)[1:topGgenes]#[in(new_gene_significance_weights[k][g] ,sort(new_gene_significance_weights[k],rev=true)[1:topGgenes]) for g in 1:G]
        gene_labels = gene_names[top_Ggenes_bool]
        # gene_labels = collect(1:G)#collect(1:2:KMax)
        @gp :- data[k,top_Ggenes_bool] gene_labels " using 1:xtic(2) t '$(new_cluster_labels[k])' lc rgb '#$(hex(linCol[k]))' fs transparent solid 0.75 " :- 
        subplt_indx +=1
    end
    @gp
    # if to_display
         
    # end
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end

end

# function gnu_multiplot_plot_inclusion_probabilities2(gene_significance_weights,topGgenes,gene_names;color_order=nothing,cluster_labels=nothing,linCol=nothing ,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,clusters_used_annotation=nothing,to_display=true,to_include_bool=nothing)
#     if to_display
#         Gnuplot.quitall()
#         Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
#     else
#         Gnuplot.quitall()
#         Gnuplot.options.gpviewer = false
#     end
#     G = length(gene_significance_weights[1])
#     K = length(gene_significance_weights)
#     gene_significance_weights =[ [va <= 1/(10*G) ? 0.0 : va  for va in el] for el in gene_significance_weights]
#     if isnothing(to_include_bool)
#         to_include_bool = [!(any(isnan.(el)) || all(iszero.(el))) for el in gene_significance_weights]
#     end
#     new_gene_significance_weights = gene_significance_weights[to_include_bool]
#     # new_K = sum(to_include_bool)
#     new_K_vec = collect(1:K)[to_include_bool]
#     new_K = length(new_K_vec)
#     data = permutedims(reduce(hcat,new_gene_significance_weights))
#     # data =hcat([normToProb([data[k,j] for k in 1:K_valid]) for j in 1:G]...)
#     clusters_used_annotation_valid = nothing
#     if !isnothing(clusters_used_annotation)
#         clusters_used_annotation_valid = clusters_used_annotation[to_include_bool]
#     end
#     if isnothing(cluster_labels)
#         cluster_labels = ["Cluster $k" for k in 1:K]
#     end
#     new_cluster_labels = cluster_labels[to_include_bool]
    
   
#     # new_gene_significance_weights = gene_significance_weights[to_include_bool]
#     # new_K_vec = collect(1:K)[to_include_bool]
    

    

#     if isnothing(linCol)
#         # color_scheme = :tol_rainbow;
#         # if isnothing(color_order)
#         #     color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
#         #     color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
#         # end
#         # linCol = colorschemes[color_scheme][color_order]
#         color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
#         linCol = colorschemes[color_scheme][color_order] 
#     end
#     num_cols = ceil(sqrt(new_K));#round(sqrt(T))
#     num_rows = round(sqrt(new_K));#ceil(sqrt(T))
#     @gp "set size 1,1"  
#     @gp :- "set origin 0,0"
#     @gp :- "set multiplot layout $num_rows,$num_cols title 'Posterior Inclusion Probabilities for top $topGgenes Genes in Inferred Clusters' font ',20' offset 0,0.1 rowsfirst scale 1,1"
#     # @gp :- "set multiplot title offset 0,1 "

#     @gp :- "set grid y"
#     @gp :- "set style data histograms"
#     # @gp :- "set style histogram rowstacked"
#     # @gp :- "set boxwidth 0.5"
#     @gp :- "set style fill solid 1.0 border -1"
#     @gp :- "set ytics nomirror"
#     @gp :- "set xtics 1"
#     @gp :- "unset key"
#     @gp :- "set xtics font ',7'"
#     @gp :- "set xtics rotate"
#     @gp :- "set ytics font ',7'"
#     @gp :- "set xlabel font ',7'"
#     @gp :- "set ylabel font ',7'"
#     @gp :- "set xrange [:$(topGgenes)]"
#     @gp :- "set yrange [0:$(maximum(maximum.(data)) + 1/G * maximum(maximum.(data)) )]"
    
#     # @gp :- "set ytics 10"
#     @gp :- xlab = "Gene" ylab = "PIPs"#"Posterior Inclusion Probability (PIP^k_{j})"
#     # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
#     subplt_indx = 1
    

#     for k in 1:new_K
#         @gp :- subplt_indx title = "Cluster $(new_K_vec[k])"  :- #"set size square"
#         @gp :- "set title font ',15'"
#         # for j in 1:G
#         #     # @gp :- data[k][j] gene_labels " using 1:xtic(2) t '$(gene_labels[j])' lc rgb '#$(hex(linCol[j]))' fs transparent solid 0.75 " :- 
            
#         # end
#         top_Ggenes_bool = sortperm(new_gene_significance_weights[k],rev=true)[1:topGgenes]#[in(new_gene_significance_weights[k][g] ,sort(new_gene_significance_weights[k],rev=true)[1:topGgenes]) for g in 1:G]
#         gene_labels = gene_names[top_Ggenes_bool]
#         # gene_labels = collect(1:G)#collect(1:2:KMax)
#         @gp :- data[k,top_Ggenes_bool] gene_labels " using 1:xtic(2) t '$(new_cluster_labels[k])' lc rgb '#$(hex(linCol[k]))' fs transparent solid 0.75 " :- 
#         subplt_indx +=1
#     end
#     @gp
#     # if to_display
         
#     # end
#     if !isnothing(filenamebase)
#         plotID = "Multiplot"
#         filename_vec = split(filenamebase,".")
#         ext = filename_vec[end]
#         new_file =  filename_vec[end-1]*plotID  
#         filename_vec[end-1] =  new_file 
#         filename_vec[end]= "."*filename_vec[end]
#         new_filename = join(filename_vec)
#         if ext == "png"
#             Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
#             if save_svg_copy
#                 svg_filename_vec = deepcopy(filename_vec)
#                 svg_filename_vec[end]= ".svg"
#                 svg_filename = join(svg_filename_vec)
#                 Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
#             end
#         end
#         if ext == "svg"
#             Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
#         end
#         if ext == "pdf"
#             Gnuplot.options.term = "qt  size 100,100"
#             Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
#             if save_svg_copy
#                 svg_filename_vec = deepcopy(filename_vec)
#                 svg_filename_vec[end]= ".svg"
#                 svg_filename = join(svg_filename_vec)
#                 Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
#             end
#         end
#     end

# end


function gnu_plot_dotplot(gene_significance_weights,cluster_labels,gene_labels;color_order=nothing,fig_size=(1100,900),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    G = length(gene_significance_weights[1])
    gene_significance_weights =[ [va <= 1/(10*G) ? 0.0 : va  for va in el] for el in gene_significance_weights]
    to_include_bool = [!(any(isnan.(el)) || all(iszero.(el))) for el in gene_significance_weights]
    gene_significance_weights2 = gene_significance_weights[to_include_bool]
    K_valid = sum(to_include_bool)
    K_valid_vec = collect(1:K)[to_include_bool]
    data = permutedims(reduce(hcat,gene_significance_weights2))
    # data =hcat([normToProb([data[k,j] for k in 1:K_valid]) for j in 1:G]...)
    
    
    new_cluster_labels = cluster_labels[to_include_bool]
    new_gene_labels = gene_labels
    
    # cluster_labels = ["Cluster $k" for k in 1:K]
    color_scheme = :tol_rainbow;
    if isnothing(color_order)
        color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
        color_order= vcat(color_order[1:4], shuffle(color_order[5:end]))
    end
    linCol = colorschemes[color_scheme][color_order]

    
    y_vals = collect(1:K_valid)
    y_vec = innermelt(y_vals,G)
    x_vals = collect(LinRange(0.1,10,G))
    x_vec = outermelt(x_vals,K_valid)
    radii = recursive_flatten(gene_significance_weights2)
    cluster_labels_vec = repeat(new_cluster_labels,inner=G)
    gene_labels_vec = repeat(new_gene_labels,outer =K_valid)

    # Gnuplot.options.term = "qt  size 600,600"
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set style fill  transparent solid 0.8"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set xtics rotate"
    @gp :- "set xrange [-1:$(maximum(x_vals)+.1*maximum(x_vals))]"
    @gp :- "set yrange [0:$(maximum(y_vals)+.1*maximum(y_vals))]"
    @gp :- "set key invert top"  :-
    @gp :- "set key font ',9'"
    @gp :- x_vec y_vec radii cluster_labels_vec gene_labels_vec "using 1:2:3:ytic(4):xtic(5) with circles t 'Posterior Inclusion Probability (PIP^k_{j})' "
end

function gnu_plot_allocations_barchart(results_dict;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[22,10,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    data = [results_dict[key][:num_alloc][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Allocations for the $num_functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- xlab = "Function Name" ylab = "Allocations"
    @gp :- data function_labels " using 1:xtic(2) t 'Function Allocation(s)' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_plot_memory_barchart(results_dict;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[10,22,15,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    data = [results_dict[key][:memory][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    conversion_const= 1.0
    unit_type = "mem"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data ./conversion_const
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Memory Used for the $num_functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- xlab = "Function Name" ylab = "Memory $units_str"
    @gp :- data function_labels " using 1:xtic(2) t 'Memory Used' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_plot_avg_time_barchart(results_dict;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[15,10,22,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    data = [results_dict[key][:avg_time][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    
    unit_type = "time"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data ./conversion_const
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Avgerage Time spent on each of the $num_functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- xlab = "Function Name" ylab = "Elapsed Time $units_str"
    @gp :- data function_labels " using 1:xtic(2) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_plot_median_time_barchart(results_dict;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[26,10,22,15,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    data = [results_dict[key][:med_time][1] for key in collect(keys(results_dict))[2:end]];
    lb = [results_dict[key][:min_time][1] for key in collect(keys(results_dict))[2:end]];
    ub = [results_dict[key][:max_time][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");

    unit_type = "time"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data ./conversion_const
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Median Time spent on each of the $num_functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- xlab = "Function Name" ylab = "Elapsed Time $units_str"
    @gp :- data function_labels " using 1:xtic(2) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_plot_num_evals_barchart(results_dict;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[26,10,22,15,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    data = [results_dict[key][:num_evals][1] for key in collect(keys(results_dict))[2:end]];

    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Number of Evalations of each of the $num_functions with G=$G, K=$K,N=$N,T=$T (More is Better)'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- xlab = "Function Name" ylab = "Number of Evalations"
    @gp :- data function_labels " using 1:xtic(2) t 'Number of Evalations' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_plot_cumulative_avg_time(results_dict;filenamebase=nothing,fig_size=(1100,900),plot_prop=false,save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[15,10,22,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    
    
    unit_type = "time"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = cumsum([results_dict[key][:avg_time][1] for key in collect(keys(results_dict))[2:end]]);
    data = data ./conversion_const
    if plot_prop
        data = data ./ sum([results_dict[key][:avg_time][1] for key in collect(keys(results_dict))[2:end]])
        ylbl = "Cummulative Proportion of Total Time Spent"
    
    else
        ylbl = "Cummulative Total Time Spent $units_str"
    end
    function_range = collect(1:num_functions)
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Cumulative Time spent on each of the $num_functions functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- "set key left " 
    @gp :- xlab = "Function Name In Algorithmic Sequence" ylab = ylbl
    @gp :- function_range data function_labels "using 1:2:xtic(3)  w linespoints lw 1 lc 'red' pt 7 t 'Cummulative Time Spent'" :-
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_plot_cumulative_alloc(results_dict;filenamebase=nothing,fig_size=(1100,900),plot_prop=false,save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[15,10,22,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    data = cumsum([results_dict[key][:num_alloc][1] for key in collect(keys(results_dict))[2:end]]);
    if plot_prop
        data = data ./ sum([results_dict[key][:num_alloc][1] for key in collect(keys(results_dict))[2:end]])
        ylbl = "Cummulative Proportion of Total Allocations"
    
    else
        ylbl = "Cummulative Total Allocations"
    end
    function_range = collect(1:num_functions)
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Cumulative Allocations spent on each of the $num_functions functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- "set key left " 
    @gp :- xlab = "Function Name In Algorithmic Sequence" ylab = ylbl
    @gp :- function_range data function_labels "using 1:2:xtic(3)  w linespoints lw 1 lc 'red' pt 7 t 'Cummulative Allocations'" :-
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_plot_cumulative_memory(results_dict;filenamebase=nothing,fig_size=(1100,900),plot_prop=false,save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :tol_rainbow;
    color_order=[15,10,22,26,1,2,3,4,5,6,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
    linCol = colorschemes[color_scheme][color_order]
    G, K, N, T, C_t, rand_init, nothing_init= (; results_dict[collect(keys(results_dict))[1]]...)
    num_functions = length(collect(keys(results_dict))[2:end])
    
    unit_type = "mem"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    
    data = cumsum([results_dict[key][:memory][1] for key in collect(keys(results_dict))[2:end]]);
    data = data ./conversion_const
    if plot_prop
        data = data ./ sum([results_dict[key][:memory][1] for key in collect(keys(results_dict))[2:end]])
        ylbl = "Cummulative Proportion of Total Memory Used "
    
    else
        ylbl = "Cummulative Total Memory Used $units_str"
    end
    function_range = collect(1:num_functions)
    function_labels = [results_dict[key][:name][1] for key in collect(keys(results_dict))[2:end]];
    function_labels = replace.(function_labels,"_" => " ");
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Cumulative Memory Used spent on each of the $num_functions functions with G=$G, K=$K,N=$N,T=$T'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',10'"
    @gp :- "set key left " 
    @gp :- xlab = "Function Name In Algorithmic Sequence" ylab = ylbl
    @gp :- function_range data function_labels "using 1:2:xtic(3)  w linespoints lw 1 lc 'red' pt 7 t 'Cummulative Memory Used'" :-
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end


function gnu_plot_varying_allocations(summary_dict,varing_feature_levels, varing_features;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
    linCol = colorschemes[color_scheme]
    num_functions = length(collect(keys(summary_dict))[2:end])
    function_keys = collect(keys(summary_dict))[2:end]
    data = hcat([summary_dict[key][:num_alloc] for key in function_keys]...);

    function_range = collect(1:num_functions)
    function_labels = [summary_dict[key][:name][1] for key in function_keys];
    function_labels = replace.(function_labels,"_" => " ");
    
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Allocations spent on each of the $num_functions function as  $varing_features varies'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "set logscale x"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',7'"
    @gp :- "set key right outside " 
    @gp :- xlab = "$varing_features" ylab =  " Allocations"
    for i in 1:num_functions
        @gp :- varing_feature_levels data[:,i] varing_feature_levels "using 1:2:xtic(3)  w linespoints lw 1 lc rgb '#$(hex(linCol[i]))' pt $i t '$(function_labels[i])'" :-
    end
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_plot_varying_memory(summary_dict,varing_feature_levels, varing_features;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
    linCol = colorschemes[color_scheme]

    num_functions = length(collect(keys(summary_dict))[2:end])
    function_keys = collect(keys(summary_dict))[2:end]
    data = hcat([summary_dict[key][:memory] for key in function_keys]...);

    function_range = collect(1:num_functions)
    function_labels = [summary_dict[key][:name][1] for key in function_keys];
    function_labels = replace.(function_labels,"_" => " ");
    
    unit_type = "mem"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data ./conversion_const
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Memory Usage spent on each of the $num_functions function as  $varing_features varies'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    @gp :- "set logscale x"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',7'"
    @gp :- "set key right outside " 
    @gp :- xlab = "$varing_features" ylab =  " Memory $units_str"
    for i in 1:num_functions
        @gp :- varing_feature_levels data[:,i] varing_feature_levels "using 1:2:xtic(3)  w linespoints lw 1 lc rgb '#$(hex(linCol[i]))' pt $i t '$(function_labels[i])'" :-
    end
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end


function gnu_plot_varying_avg_time(summary_dict,varing_feature_levels, varing_features;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
    linCol = colorschemes[color_scheme]

    num_functions = length(collect(keys(summary_dict))[2:end])
    function_keys = collect(keys(summary_dict))[2:end]
    data = hcat([summary_dict[key][:avg_time] for key in function_keys]...);

    function_range = collect(1:num_functions)
    function_labels = [summary_dict[key][:name][1] for key in function_keys];
    function_labels = replace.(function_labels,"_" => " ");
    
    unit_type = "time"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data ./conversion_const
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Average Elapsed Time spent on each of the $num_functions function as  $varing_features varies'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set logscale x"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',7'"
    @gp :- "set key right outside " 
    @gp :- xlab = "$varing_features" ylab =  "Average Elapsed Time $units_str"
    for i in 1:num_functions
        @gp :- varing_feature_levels data[:,i] varing_feature_levels "using 1:2:xtic(3)  w linespoints lw 1 lc rgb '#$(hex(linCol[i]))' pt $i t '$(function_labels[i])'" :-
    end
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end
function gnu_plot_varying_median_time(summary_dict,varing_feature_levels, varing_features;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
    linCol = colorschemes[color_scheme]

    num_functions = length(collect(keys(summary_dict))[2:end])
    function_keys = collect(keys(summary_dict))[2:end]
    data = hcat([summary_dict[key][:med_time] for key in function_keys]...);

    function_range = collect(1:num_functions)
    function_labels = [summary_dict[key][:name][1] for key in function_keys];
    function_labels = replace.(function_labels,"_" => " ");
    unit_type = "time"
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data ./conversion_const
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Median Elapsed Time spent on each of the $num_functions function as  $varing_features varies'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set logscale x"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',7'"
    @gp :- "set key right outside " 
    @gp :- xlab = "$varing_features" ylab =  "Average Elapsed Time $units_str"
    for i in 1:num_functions
        @gp :- varing_feature_levels data[:,i] varing_feature_levels "using 1:2:xtic(3)  w linespoints lw 1 lc rgb '#$(hex(linCol[i]))' pt $i t '$(function_labels[i])'" :-
    end
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_plot_varying_num_evals(summary_dict,varing_feature_levels, varing_features;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
    linCol = colorschemes[color_scheme]

    num_functions = length(collect(keys(summary_dict))[2:end])
    function_keys = collect(keys(summary_dict))[2:end]
    data = hcat([summary_dict[key][:num_evals] for key in function_keys]...);

    function_range = collect(1:num_functions)
    function_labels = [summary_dict[key][:name][1] for key in function_keys];
    function_labels = replace.(function_labels,"_" => " ");
    
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Number Evalutions of each of the $num_functions function as  $varing_features varies'"
    @gp :- "set title font ',14'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xtics 1"
    # @gp :- "unset key"
    @gp :- "set xtics font ',10'"
    @gp :- "set ytics font ',10'"
    @gp :- "set xtics rotate"
    @gp :- "set logscale x"
    @gp :- "set xlabel font ',14'"
    @gp :- "set ylabel font ',14'"
    @gp :- "set key font ',7'"
    @gp :- "set key right outside " 
    @gp :- xlab = "$varing_features" ylab =  "Number of Evaluations"
    for i in 1:num_functions
        @gp :- varing_feature_levels data[:,i] varing_feature_levels "using 1:2:xtic(3)  w linespoints lw 1 lc rgb '#$(hex(linCol[i]))' pt $i t '$(function_labels[i])'" :-
    end
    #  data function_labels " using 1:2:xtic(3) t 'Elapsed Time Spent' lc rgb '#$(hex(linCol[1]))' fs transparent solid 0.75 " :- 
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_plot_inferred_signal_null_dist(x_to_use,goi,k,mean_μ_err_post,mean_τ_err_post,mean_μ_post,mean_τ_post;genename = nothing,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    # goi =5
    if isnothing(genename)
        label_ = "Gene $goi"
    else
        label_ = genename
    end
    # k =5
    fg(x,μ,τ) = √(τ) .* exp.(.-(x.-μ).^2 ./(2) .*τ)./(√(2π))
    
    data = permutedims(reduce(hcat, reduce(vcat,x_to_use)))[:,goi]
    xline_max= maximum(data) + 0.5*(maximum(data))
    xline_min = - xline_max
    xline = LinRange(xline_min,xline_max,200)
    @gp  "set title 'Count Histogram of $label_ in Cluster $k'"
    @gp :- xlab = "$label_ value" ylab = "Count"
    data_hist = hist(data,bs=0.05)
    #2 .* scale_factor .* data_hist.counts ./ sum(data_hist.counts)
    @gp :- "set y2tics"
    @gp :- "set y2range [$(minimum(data_hist.counts) ):$(maximum(data_hist.counts) )]"
    @gp :- "set ytics nomirror"
    @gp :- "set y2lab 'Count' "
    @gp :-  data_hist.bins  data_hist.counts  "with fillsteps t '# of Cells with value at $label_' lc 'black' fs solid 0.3 noborder axis x1y2"
    @gp :- data_hist.bins  data_hist.counts  "with steps tit '' lc 'black' lw 2 axis x1y2" "set grid"
    @gp :- "set style fill transparent solid 0.5"
    @gp :- "set key title 'Key' box 3" xlab = "$label_ value (x)" ylab = "P(x)"
    @gp :- xline fg(xline, mean_μ_err_post[1][goi], mean_τ_err_post[1][goi]) "w filledcurves lc '#E69F00' lw 2 t 'Inferred Noise Distribution' axis x1y1"
    @gp :- xline fg(xline, mean_μ_post[k][goi], mean_τ_post[k][goi]) "w filledcu lc rgb '#56B4E9' lw 2 t 'Inferred Signal Distribution'  fs transparent solid 0.25 axis x1y1"
    # @gp :- xline fg(xline, 5.0, 1) "w filledcu lc '#009E73' t '-1,2' fs transparent solid 0.1"
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_multiplot_all_inferred_signal_null_dist(x_to_use,mean_μ_err_post,mean_τ_err_post,mean_μ_post,mean_τ_post;genename = nothing,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    # goi =5
    G = length(x_to_use[1][1])
    K = length(mean_μ_post)
    if isnothing(genename)
        genename =["Gene $j" for j in 1:G]
    end
    # k =5
    fg(x,μ,τ) = √(τ) .* exp.(.-(x.-μ).^2 ./(2) .*τ)./(√(2π))
    
    
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(K),$(G) title 'Inferred Cluster-Specific and Null Gene Distribution' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set multiplot title offset 0,1 "

    @gp :- "set grid y"

    # @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"


    # @gp :- "set ytics 10"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1


    for j in 1:G
        for k in 1:K
            label_ = genename[j]
            @gp :- subplt_indx title = "Gene: $label_, Cluster: $k"  :- #"set size square"
            # @gp :- "set title font ',15'"
            # for j in 1:G
            #     # @gp :- data[k][j] gene_labels " using 1:xtic(2) t '$(gene_labels[j])' lc rgb '#$(hex(linCol[j]))' fs transparent solid 0.75 " :- 
                
            # end
            
            
            data = permutedims(reduce(hcat, reduce(vcat,x_to_use)))[:,j]
            xline_max= maximum(data) + 0.5*(maximum(data))
            xline_min = - xline_max
            xline = LinRange(xline_min,xline_max,200)
            @gp :- "set title 'Count Histogram of $label_ in Cluster $k'"
            @gp :- xlab = "$label_ value" ylab = "Count"
            data_hist = hist(data,bs=0.05)
            #2 .* scale_factor .* data_hist.counts ./ sum(data_hist.counts)
            @gp :- "set y2tics"
            @gp :- "set y2range [$(minimum(data_hist.counts) ):$(maximum(data_hist.counts) )]"
            @gp :- "set ytics nomirror"
            @gp :- "set y2lab 'Count' "
            @gp :-  data_hist.bins  data_hist.counts  "with fillsteps t '# of Cells with value at $label_' lc 'black' fs solid 0.3 noborder axis x1y2"
            @gp :- data_hist.bins  data_hist.counts  "with steps tit '' lc 'black' lw 2 axis x1y2" "set grid"
            @gp :- "set style fill transparent solid 0.5"
            @gp :- "set key title 'Key' box 3" xlab = "$label_ value (x)" ylab = "P(x)"
            @gp :- xline fg(xline, mean_μ_err_post[1][j], mean_τ_err_post[1][j]) "w filledcurves lc '#E69F00' lw 2 t 'Inferred Noise Distribution' axis x1y1"
            @gp :- xline fg(xline, mean_μ_post[k][j], mean_τ_post[k][j]) "w filledcu lc rgb '#56B4E9' lw 2 t 'Inferred Signal Distribution'  fs transparent solid 0.25 axis x1y1"
            subplt_indx +=1
            # @gp :- "set cbtics out nomirror" :-
            # @gp :- "set key invert reverse outside below"  :-
            # @gp :- "set key font ',12'"
            # @gp :- "set key vertical maxrows 1" :-
        end
    end

    # @gp :- xline fg(xline, 5.0, 1) "w filledcu lc '#009E73' t '-1,2' fs transparent solid 0.1"
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_plot_bayesfactors(bayes_factor;filenamebase=nothing,fig_size=(1100,900),linCol=nothing ,save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    G = length(bayes_factor[1])
    K = length(bayes_factor)
    # gene_significance_weights =[ [va <= 1/(10*G) ? 0.0 : va  for va in el] for el in gene_significance_weights]
    to_include_bool = [!(any(isnan.(el)) || all(iszero.(el))) for el in bayes_factor]
    bayes_factor2 = bayes_factor[to_include_bool]
    K_valid = sum(to_include_bool)
    K_valid_vec = collect(1:K)[to_include_bool]
    data = permutedims(reduce(hcat,bayes_factor2))
    # data =hcat([normToProb([data[k,j] for k in 1:K_valid]) for j in 1:G]...)
    
    
    new_cluster_labels = ["Cluster $k" for k in 1:K][to_include_bool]
    
    
    



    data = permutedims(reduce(hcat,bayes_factor))
    # cluster_labels = ["Cluster $k" for k in 1:K]
    if isnothing(linCol)
        color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
        linCol = colorschemes[color_scheme] 
    end
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    plotTitle = "Bayes Factor for Expressed Genes in All Cells"
    @gp :- "set grid y"
    @gp :- "set grid x"
    # @gp :- "set logscale y"
    @gp :- "set style fill  transparent solid 0.95 solid  border lt .1 " #"set style fill  transparent solid 0.45 noborder" border lt .05 
    @gp :- "set style circle radius .25"
    @gp :- "set xtics font ',11'"
    @gp :- "set ytics font ',11'"
    @gp :- "set xlabel font ',13'"
    @gp :- "set ylabel font ',13'"
    @gp :- "set title font ',15'"
    @gp :- "unset xtic"
    @gp :- "set xrange [0.5:$(K_valid*G+0.5)]"
    # @gp :- "set yrange [$(minimum(log.(data)) + 10*minimum(log.(data))):$(maximum(log.(data))+.1*maximum(log.(data)))]"
    @gp :- title = plotTitle
    @gp :- xlab = "Gene in Cluster" ylab = "log_{10}(BayesFactor)"
    for k in 1:K_valid
    
        @gp :- collect((k-1)*G + 1:k*G) log.(permutedims(reduce(hcat,bayes_factor2)))[k,:] "with circles t '$(new_cluster_labels[k])' lc rgb '#$(hex(linCol[k]))' " :- 
        @gp :- "set cbtics out nomirror" :-
        @gp :- "set key font ',12'"
        @gp :- "set key invert reverse Left outside"  :-
        
        # with points t '$(unique_clusters_labels[l])' lc rgb '#$(hex(linCol[l]))' " :- 
        # @gp :- "set cbtics out nomirror" :-
        # @gp :- "set key invert reverse Left outside"  :-
    end
    @gp :- "set arrow from graph 0,first 2 to graph 1,first 2 nohead lc rgb '#FF0000' front"
    @gp :- "set arrow from graph 0,first 1 to graph 1,first 1 nohead lc rgb '#FFA500' front"
    @gp :- "set arrow from graph 0,first 0 to graph 1,first 0 nohead lc rgb '#000000' front"
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_multiplot_grand_tour(x_to_use,z_labels;cluster_labels =nothing,genename = nothing,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,linCol=nothing ,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    # goi =5
    G = length(x_to_use[1][1])
    z_vec = vcat(z_labels...)
    K = length(unique(z_vec))
    N = length(z_vec)
    unique_clusters = unique(z_vec)
    # println(K)
    # println(unique_clusters)
    
    if isnothing(genename)
        genename =["Gene $j" for j in 1:G]
    end
    if isnothing(cluster_labels)
        cluster_labels = ["Cluster $k" for k in unique_clusters]
    end
    if isnothing(linCol)
        color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
        linCol = colorschemes[color_scheme] 
    end
    # k =5
    

    
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(G),$(G) title 'Grand Tour of Expression Data' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    # @gp :- "set multiplot title offset 0,1 "
    @gp :- "unset key"
    @gp :- "set grid y"

    # @gp :- "set style histogram rowstacked"
    # @gp :- "set boxwidth 0.5"
    @gp :- "set style fill solid 1.0 border -1"


    # @gp :- "set ytics 10"
    # @gp :- Matrix(avg_counts_mat[:,:,1,1]) "using 1 t 'Basal', '' using 2 t 'Classical', '' using 3 t 'Intermediate', '' using 4:xtic(1) t 'Organoid' "
    subplt_indx = 1
    xmat =  permutedims(reduce(hcat, reduce(vcat,x_to_use)))
    x1_jitter = rand(Uniform(-0.05,0.05),N)
    x2_jitter = rand(Uniform(-0.05,0.05),N)
    for j in 1:G
        for g in 1:G
            label1_ = genename[j]
            label2_ = genename[g]
            # @gp :- subplt_indx title = "x: $label1_,  y: $label2_"  :- #"set size square"
            # @gp :- "set title font ',15'"
            # for j in 1:G
            #     # @gp :- data[k][j] gene_labels " using 1:xtic(2) t '$(gene_labels[j])' lc rgb '#$(hex(linCol[j]))' fs transparent solid 0.75 " :- 
                
            # end
            
            if j == g
                data = xmat[:,j]
                @gp :- "set xrange [$(minimum(xmat)+.1*minimum(xmat)):$(maximum(xmat)+.1*maximum(xmat))]" #[$(minimum(data)+.1*minimum(data)):$(maximum(data)+.1*maximum(data))]
                xline_max= maximum(data) + 0.5*(maximum(data))
                xline_min = - xline_max
                xline = LinRange(xline_min,xline_max,200)
                # @gp :- "set title 'Grand Tour of Data'"
                # @gp :- xlab = "$(genename[j]) value" ylab = "$(genename[g]) value"
                data_hist = hist(data,bs=0.05)
                #2 .* scale_factor .* data_hist.counts ./ sum(data_hist.counts)
                # @gp :- "set y2tics"
                @gp :- subplt_indx title = "$label1_ counts" 
                @gp :- "set yrange [$(minimum(data_hist.counts) ):$(maximum(data_hist.counts) )]"
                @gp :- "set ytics nomirror"
                # @gp :- "set ylab 'Count' "
                @gp :-  data_hist.bins  data_hist.counts  "with fillsteps t '# of Cells with value at $label1_' lc 'black' fs solid 0.3"
                @gp :- data_hist.bins  data_hist.counts  "with steps tit '' lc 'black' lw 2 " "set grid"
                @gp :- "set style fill transparent solid 0.5"
            else
                x1 = xmat[:,j] .+ x1_jitter
                x2 = xmat[:,g] .+ x2_jitter
                xrange_x1= x1
                yrange_x2= x2
                @gp :- subplt_indx title = "x: $label1_,  y: $label2_"  
                @gp :- "set style fill transparent solid 0.85"
                @gp :- "set xrange [$(minimum(xrange_x1)+.1*minimum(xrange_x1)):$(maximum(xrange_x1)+.1*maximum(xrange_x1))]"
                @gp :- "set yrange [$(minimum(yrange_x2)+.1*minimum(yrange_x2)):$(maximum(yrange_x2)+.1*maximum(yrange_x2))]"
                # @gp :- title = plotTitle
                # @gp :- xlab = "$label1_" ylab = "$label2_"
                for l in 1:K
                    clus = unique_clusters[l]
                    @gp :- x1[z_vec .==clus] x2[z_vec .==clus] "with circles t '$(cluster_labels[l])' lc rgb '#$(hex(linCol[l]))' " :- 
                    @gp :- "set cbtics out nomirror" :-
                    # @gp :- "set key font ',12'"
                    # @gp :- "set key invert reverse Left outside"  :-
                    
                    # with points t '$(unique_clusters_labels[l])' lc rgb '#$(hex(linCol[l]))' " :- 
                    # @gp :- "set cbtics out nomirror" :-
                    # @gp :- "set key invert reverse Left outside"  :-
                end 
            end
            subplt_indx +=1
            # @gp :- "set cbtics out nomirror" :-
            # @gp :- "set key invert reverse outside below"  :-
            # @gp :- "set key font ',12'"
            # @gp :- "set key vertical maxrows 1" :-
        end
    end

    # @gp :- xline fg(xline, 5.0, 1) "w filledcu lc '#009E73' t '-1,2' fs transparent solid 0.1"
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
    return linCol
end

function gnu_ARIvsMaxIter_plot(ari_vec,max_iter_range;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xrange [:$(maximum(max_iter_range)+1)]"
    # @gp :- "set ytics 10"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set key font ',12'"
    @gp :- "set title font ',15'"
    @gp :- "set key invert reverse Left outside" 
    @gp :- max_iter_range ari_vec  "w linespoints lw 1 lc 'red' pt 7 t 'ARI'"
    @gp :- xlab = "Maximum Iterations Intialization" ylab = "Final ARI"
    @gp :- title = "Final ARI For Maximum Iterations "

    @gp

    if !isnothing(filenamebase)
        plotID = "plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_FinalElbovsMaxIter_plot(final_elbo_vec,max_iter_range;filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    @gp :- "set xrange [:$(maximum(max_iter_range)+1)]"
    # @gp :- "set ytics 10"
    @gp :- "set xtics font ',9'"
    @gp :- "set ytics font ',9'"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set key font ',12'"
    @gp :- "set title font ',15'"
    @gp :- "set key invert reverse Left outside" 
    @gp :- max_iter_range final_elbo_vec  "w linespoints lw 1 lc 'blue' pt 7 t 'ELBO'"
    @gp :- xlab = "Maximum Iterations Intialization" ylab = "Final ELBO"
    @gp :- title = "Final ELBO For Maximum Iterations "

    @gp

    if !isnothing(filenamebase)
        plotID = "plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_multiplot_elbo_vs_metrics_correlation(final_elbo_val_vec,summary_ari_vec,summary_nmi_vec,summary_vmeasure_vec;fig_size=(1100,900),to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    @gp "set multiplot layout 2,2 title 'Effect of Elbo on Metrics' offset 0,0.1 columnsfirst scale 1.1,0.9"
    # @gp :- "unset key"

    @gp :- 1 title = "Log(Elbo) Vs ARI" "set size square"
    # @gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
    @gp :- "set logscale x"
    @gp :- "set yrange [-1:1]"
    @gp :- xlab = "Log (Elbo)" ylab = "ARI"
    @gp :- final_elbo_val_vec summary_ari_vec "w p pt 7 ps 1 lc '#0072B2' lw 2 t 'ARI'"
    @gp :- "set key below"

    @gp :- 2 title = "Log(Elbo) Vs NMI" "set size square"
    # @gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
    @gp :- "set logscale x"
    @gp :- "set yrange [0:1]"
    @gp :- xlab = "Log (Elbo)" ylab = "NMI"
    @gp :- final_elbo_val_vec summary_nmi_vec "w p pt 7 ps 1 lc '#009E73' lw 2 t 'NMI'"
    @gp :- "set key below"

    @gp :- 3 title = "Log(Elbo) Vs VMeasure" "set size square"
    # @gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
    @gp :- "set logscale x"
    @gp :- "set yrange [0:1]"
    @gp :- xlab = "Log (Elbo)" ylab = "VMeasure"
    @gp :- final_elbo_val_vec summary_vmeasure_vec "w p pt 7 ps 1 lc '#E69F00' lw 2 t 'VMeasure'"
    @gp :- "set key below"
    # @gp :- "set key vertical maxrows 1" :-
    # @gp :- "set key at screen .85, screen .95"  :-
    @gp
end


function gnu_ParamVsMetric_scatterplot(metric_vec,params_range;unique_labels =nothing,filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false, lb= nothing,ub=nothing,linCol=nothing,point_labels=nothing,xlabel=nothing,ylabel=nothing,plt_title=nothing,metric_label=nothing,x_range_low=nothing,x_range_high=nothing,y_range_low=nothing,y_range_high=nothing,  conf_level = 0.95,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    num_points = length(metric_vec)
    if isnothing(linCol) 
        color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
        linCol = colorschemes[color_scheme] 
    end
    if isnothing(point_labels)
        point_labels = ones(Int,num_points)
    end
    if isnothing(xlabel)
        xlabel = "Varying Value"
    end
    if isnothing(ylabel)
        ylabel = "Metric of Interest"
    end
    if isnothing(plt_title)
        plt_title = "NCLUSION performance measured by the $ylabel v.s. $xlabel"
    end

    if isnothing(unique_labels)
        unique_labels = sort(unique(point_labels))
    end
    num_unique_labels = length(unique_labels)
    metrics_jitter = rand(Uniform(-0.005,0.005),num_points)
    ranges_jitter  = rand(Uniform(-1,1),num_points)
    params_range = params_range .+ ranges_jitter
    metric_vec = metric_vec .+ metrics_jitter

    if isnothing(metric_label)
        metric_label = ["Cluster Grouping $k" for k in unique_labels]
    end
    if isnothing(x_range_low)
        x_range_low = minimum(params_range) + 0.1 * minimum(params_range)
    end
    if isnothing(x_range_high)
        x_range_high = maximum(params_range) + 0.1 * maximum(params_range)
    end
    if isnothing(y_range_low)
        y_range_low = minimum(metric_vec) + 0.1 *  minimum(metric_vec)
    end
    if isnothing(y_range_high)
        y_range_high = maximum(metric_vec) + 0.1 * maximum(metric_vec)
    end
    

    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    # @gp :- "set xrange [:$(maximum(params_range)+1)]"
    # @gp :- "set ytics 10"
    # @gp :- "set logscale x"
    @gp :- "set xtics font ',7'"
    @gp :- "set ytics font ',7'"
    @gp :- "set xlabel font ',10'"
    @gp :- "set ylabel font ',10'"
    @gp :- "set key font ',7'"
    # @gp :- "unset key"
    @gp :- "set title font ',12'"
    # @gp :- "set key reverse right bmargin"
    @gp :- "set key reverse center bmargin"
    @gp :- "set yrange [$y_range_low:$y_range_high]" 
    @gp :- "set xrange [$x_range_low:$x_range_high]" 
    @gp :- "set pointsize 3"
    @gp :- "set style fill  transparent solid 0.8 border lt .05" # noborder
    @gp :- xlab = "$xlabel" ylab = "$ylabel"
    @gp :- title = "$plt_title "
    # metrics_jitter = rand(Uniform(-0.005,0.005),num_points)
    # ranges_jitter  = rand(Uniform(-1,1),num_points)
    # params_range = params_range .+ ranges_jitter
    # metric_vec = metric_vec .+metrics_jitter
    for k in 1:num_unique_labels
        point_color_indx = unique_labels[k]

        if !isnothing(lb) && !isnothing(ub)
            lo = lb .+metrics_jitter
            hi = ub .+metrics_jitter
            @gp :- "set bars 2"
            @gp :- params_range[k .== point_labels] metric_vec[k .== point_labels] lo[k .== point_labels] hi[k .== point_labels] "using 1:2:3:4 with yerr ps 1 lc rgb '#$(hex(linCol[point_color_indx]))' lw 2  t '$(conf_level * 100)% Confidence Interval'" :-
            #
        end
        @gp :- params_range[k .== point_labels] metric_vec[k .== point_labels]  "w p pt 7 ps 1 lc rgb '#$(hex(linCol[point_color_indx]))' lw 2 t '$(metric_label[point_color_indx])'" :-
        # if !isnothing(lb) && !isnothing(ub)
        #     lo = lb .+metrics_jitter
        #     hi = ub .+metrics_jitter
        #     @gp :- "set bars 2"
        #     @gp :- params_range[1 .== point_labels] metric_vec[1 .== point_labels] lo[1 .== point_labels] hi[1 .== point_labels] "using 1:2:3:4 with yerr ps 1 lc rgb '#$(hex(linCol[1]))' lw 2  t '$(conf_level * 100)% Confidence Interval'" :-
        #     #
        # end
        # @gp :- params_range[1 .== point_labels] metric_vec[1 .== point_labels]  "w p pt 7 ps 1 lc rgb '#$(hex(linCol[1]))' lw 2 t '$(metric_label[1])'" :-
        # # 
    
    end


    @gp

    if !isnothing(filenamebase)
        plotID = "plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_multiplot_ParamVsMetric_scatterplot(metric_vec,params_range;unique_labels =nothing, filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false, lb= nothing,ub=nothing,linCol=nothing,point_labels=nothing,xlabel=nothing,ylabel=nothing,plt_title=nothing,metric_label=nothing,x_range_low=nothing,x_range_high=nothing,y_range_low=nothing,y_range_high=nothing, conf_level = 0.95,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    num_points = length(metric_vec)
    if isnothing(linCol) 
        color_scheme = :glasbey_bw_minc_20_maxl_70_n256;
        linCol = colorschemes[color_scheme] 
    end
    if isnothing(point_labels)
        point_labels = ones(Int,num_points)
    end
    if isnothing(xlabel)
        xlabel = "Varying Value"
    end
    if isnothing(ylabel)
        ylabel = "Metric of Interest"
    end
    if isnothing(plt_title)
        plt_title = "NCLUSION performance measured by the $ylabel v.s. $xlabel"
    end
    
    if isnothing(unique_labels)
        unique_labels = sort(unique(point_labels))
    end
    num_unique_labels = length(unique_labels)
    num_cols = ceil(sqrt(num_unique_labels));#round(sqrt(T))
    num_rows = round(sqrt(num_unique_labels));#ceil(sqrt(T))

    metrics_jitter = rand(Uniform(-0.005,0.005),num_points)
    ranges_jitter  = rand(Uniform(-1,1),num_points)
    params_range = params_range .+ ranges_jitter
    metric_vec = metric_vec .+ metrics_jitter
    if isnothing(metric_label)
        metric_label = ["Cluster Grouping $k" for k in unique_labels]
    end
    if isnothing(x_range_low)
        x_range_low = minimum(params_range) + 0.1 * minimum(params_range)
    end
    if isnothing(x_range_high)
        x_range_high = maximum(params_range) + 0.1 * maximum(params_range)
    end
    if isnothing(y_range_low)
        y_range_low = minimum(metric_vec) + 0.1 *  minimum(metric_vec)
    end
    if isnothing(y_range_high)
        y_range_high = maximum(metric_vec) + 0.1 * maximum(metric_vec)
    end
    
    # driver_gene_tic_linCol = colorschemes[:glasbey_bw_minc_20_hue_330_100_n256]
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title '$plt_title' font ',20' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    # @gp :- "set xrange [:$(maximum(params_range)+1)]"
    # @gp :- "set ytics 10"
    # @gp :- "set logscale x"
    @gp :- "set xtics font ',7'"
    @gp :- "set ytics font ',7'"
    @gp :- "set xlabel font ',10'"
    @gp :- "set ylabel font ',10'"
    @gp :- "set key font ',7'"
    # @gp :- "unset key"
    @gp :- "set title font ',12'"
    # @gp :- "set key reverse right bmargin"
    @gp :- "set key reverse center bmargin"
    @gp :- "set yrange [$y_range_low:$y_range_high]" 
    @gp :- "set xrange [$x_range_low:$x_range_high]" 
    @gp :- "set pointsize 3"
    @gp :- "set style fill  transparent solid 0.8 noborder" # border lt .05
    @gp :- xlab = "$xlabel" ylab = "$ylabel"
    

    subplt_indx =1
    for k in 1:num_unique_labels
        point_color_indx = unique_labels[k]
        @gp :- subplt_indx title = "$(metric_label[point_color_indx]) "
        
        if k == 1
            if !isnothing(lb) && !isnothing(ub)
                lo = lb .+metrics_jitter
                hi = ub .+metrics_jitter
                @gp :- "set bars 2"
                @gp :- params_range[k .== point_labels] metric_vec[k .== point_labels] lo[k .== point_labels] hi[k .== point_labels] "using 1:2:3:4 with yerr ps 1 lc rgb '#$(hex(linCol[1]))' lw 2  t '$(conf_level * 100)%  Confidence Interval'" :-
                #
            end
            @gp :- params_range[k .== point_labels] metric_vec[k .== point_labels]  "w p pt 7 ps 1 lc rgb '#$(hex(linCol[k]))' lw 2 t '$(metric_label[k])'" :-
        
        else
            if !isnothing(lb) && !isnothing(ub)
                lo = lb .+metrics_jitter
                hi = ub .+metrics_jitter
                @gp :- "set bars 2"
                @gp :- params_range[k .== point_labels] metric_vec[k .== point_labels] lo[k .== point_labels] hi[k .== point_labels] "using 1:2:3:4 with yerr ps 1 lc rgb '#$(hex(linCol[point_color_indx]))' lw 2  t '$(conf_level * 100)% Confidence Interval'" :-
                #
            end
            @gp :- params_range[k .== point_labels] metric_vec[k .== point_labels]  "w p pt 7 ps 1 lc rgb '#$(hex(linCol[point_color_indx]))' lw 2 t '$(metric_label[point_color_indx])'" :-
            if !isnothing(lb) && !isnothing(ub)
                lo = lb .+metrics_jitter
                hi = ub .+metrics_jitter
                @gp :- "set bars 2"
                @gp :- params_range[1 .== point_labels] metric_vec[1 .== point_labels] lo[1 .== point_labels] hi[1 .== point_labels] "using 1:2:3:4 with yerr ps 1 lc rgb '#$(hex(linCol[1]))' lw 2  t '$(conf_level * 100)% Confidence Interval'" :-
                #
            end
            @gp :- params_range[1 .== point_labels] metric_vec[1 .== point_labels]  "w p pt 7 ps 1 lc rgb '#$(hex(linCol[1]))' lw 2 t '$(metric_label[1])'" :-
            # 
        end

        subplt_indx +=1
    end
    

    @gp

    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end


function gnu_plot_speedups(results_vec,function_labels,metric_name,function_name;K=5,T=10,G=720,N=10_000,rand_init=true,seed=12345, filenamebase=nothing,fig_size=(1100,900),save_svg_copy=false,linCol=nothing,units = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    if isnothing(linCol)
        color_scheme = :tol_rainbow;
        color_order=[10,22,15,26,6,1,2,3,4,5,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
        linCol = colorschemes[color_scheme][color_order]
    end
    num_func = length(results_vec)
    data = results_vec;
    num_func = length(results_vec)
    # function_labels = [el[:name][1] for key in collect(keys(results_dict))[2:end]];
    y_label_dict = Dict("num_alloc" => "Number of Allocations", "memory" => "Memory Used", "avg_time"=>"Average Time", "med_time"=>"Median Time","num_evals" => "Number of Evalations")
    y_label_units_dict = Dict("num_alloc" => "", "memory" => "mem", "avg_time"=>"time", "med_time"=>"time","num_evals" => "")
    unit_type = y_label_units_dict[metric_name]
    # conversion_const= 1.0
    units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
    data = data .* conversion_const

    y_label = y_label_dict[metric_name]
    function_labels = replace.(function_labels,"_" => "\\_");
    function_name = replace(function_name,"_" => "\\_");
    metric_name = replace(metric_name,"_" => " ");
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set title 'Comparing $(function_name) $y_label with G=$G, K=$K,N=$N,T=$T, rand\\_init=$rand_init, and seed = $seed'"
    @gp :- "set title font ',13'"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    # @gp :- "set xtics 1"
    # @gp :- "unset key"
    # @gp :- "set xtics font ',10'"
    @gp :- "unset xtics"
    @gp :- "set ytics font ',9'"
    # @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set key font ',10'"
    @gp :- "set key reverse center bmargin"
    @gp :- xlab = "$(replace(function_name,"_" => "\\_")) Variants" ylab = "$y_label $units_str"
    for i in 1:num_func
        @gp :- data[i] " t '$(function_labels[i])' lc rgb '#$(hex(linCol[i]))' fs transparent solid 0.75 " :- 
    end
    @gp
    if !isnothing(filenamebase)
        plotID = "Plot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_multiplot_metrics_speedups(results_vec,function_labels_vec,metric_name_vec,function_name;K=5,T=10,G=720,N=10_000,rand_init=true,seed=12345, filenamebase=nothing,fig_size=(1300,950),save_svg_copy=false,linCol=nothing,units_vector = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    if isnothing(linCol)
        color_scheme = :tol_rainbow;
        color_order=[10,22,15,26,6,1,2,3,4,5,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
        linCol = colorschemes[color_scheme][color_order]
    end
    num_features = length(results_vec)
    num_variants = length(results_vec[1])
    

    num_cols = ceil(sqrt(num_features));#round(sqrt(T))
    num_rows = round(sqrt(num_features));#ceil(sqrt(T))
    y_label_dict = Dict("num_alloc" => "Number of Allocations", "memory" => "Memory Used", "avg_time"=>"Average Time", "med_time"=>"Median Time","num_evals" => "Number of Evalations")
    y_label_units_dict = Dict("num_alloc" => "", "memory" => "mem", "avg_time"=>"time", "med_time"=>"time","num_evals" => "")
    
    # conversion_const= 1.0


    
    
    function_name = replace(function_name,"_" => "\\_");
    

    plt_title = "Comparing $(function_name) metrics with G=$G, K=$K,N=$N,T=$T, rand\\_init=$rand_init, and seed = $seed"
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title '$plt_title' font ',15' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    # @gp :- "set xtics 1"
    # @gp :- "unset key"
    # @gp :- "set xtics font ',10'"
    @gp :- "unset xtics"
    @gp :- "set ytics font ',9'"
    # @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set title font ',12'"
    @gp :- "set key font ',9'" # border lt .05
    @gp :- "set key reverse center bmargin"

    subplt_indx =1
    for k in 1:num_features
        metric_name = metric_name_vec[k]
        units = units_vector[k]
        function_labels = function_labels_vec[k]
        unit_type = y_label_units_dict[metric_name]
        units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)
        y_label = y_label_dict[metric_name]
        metric_name = replace(metric_name,"_" => " ");
        @gp :- subplt_indx title = "$(y_label) "
        
        @gp :- xlab = "$(replace(function_name,"_" => "\\_")) Variants" ylab = "$y_label $units_str"
        data = results_vec[k];
        data = data .* conversion_const
        function_labels = replace.(function_labels,"_" => "\\_");
        for i in 1:num_variants
            @gp :- data[i] " t '$(function_labels[i])' lc rgb '#$(hex(linCol[i]))' fs transparent solid 0.75 " :- 
        end

        subplt_indx +=1
    end
    

    @gp


    # function_labels = [el[:name][1] for key in collect(keys(results_dict))[2:end]];
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

function gnu_multiplot_functions_speedups(results_vec,function_labels_vec,metric_name,function_name_vec;K=5,T=10,G=720,N=10_000,rand_init=true,seed=12345, filenamebase=nothing,fig_size=(1300,950),save_svg_copy=false,linCol=nothing,units_vector = nothing,to_display=true)
    if to_display
        Gnuplot.quitall()
        Gnuplot.options.term = "qt  size $(fig_size[1]),$(fig_size[2])"
    else
        Gnuplot.quitall()
        Gnuplot.options.gpviewer = false
    end
    if isnothing(linCol)
        color_scheme = :tol_rainbow;
        color_order=[10,22,15,26,6,1,2,3,4,5,7,8,9,11,12,13,14,16,17,18,19,20,21,23,24,25,27,28,29]
        linCol = colorschemes[color_scheme][color_order]
    end
    num_features = length(results_vec)
    num_variants = length(results_vec[1])
    

    num_cols = ceil(sqrt(num_features));#round(sqrt(T))
    num_rows = round(sqrt(num_features));#ceil(sqrt(T))
    y_label_dict = Dict("num_alloc" => "Number of Allocations", "memory" => "Memory Used", "avg_time"=>"Average Time", "med_time"=>"Median Time","num_evals" => "Number of Evalations")
    y_label_units_dict = Dict("num_alloc" => "", "memory" => "mem", "avg_time"=>"time", "med_time"=>"time","num_evals" => "")
    
    # conversion_const= 1.0

    unit_type = y_label_units_dict[metric_name]
    y_label = y_label_dict[metric_name]
    metric_name = replace(metric_name,"_" => " ");
    
    
    # function_name = replace(function_name,"_" => "\\_");
    

    plt_title = "Comparing the $(y_label) across functions with G=$G, K=$K,N=$N,T=$T, rand\\_init=$rand_init, and seed = $seed"
    @gp "set size 1,1"  
    @gp :- "set origin 0,0"
    @gp :- "set multiplot layout $(num_rows),$(num_cols) title '$plt_title' font ',15' offset 0,0.1 rowsfirst scale 1,1"
    @gp :- "set grid y"
    @gp :- "set grid x"
    @gp :- "set style data histograms"
    @gp :- "set style fill solid 1.0 border -1"
    @gp :- "set ytics nomirror"
    # @gp :- "set xtics 1"
    # @gp :- "unset key"
    # @gp :- "set xtics font ',10'"
    @gp :- "unset xtics"
    @gp :- "set ytics font ',9'"
    # @gp :- "set xtics rotate"
    @gp :- "set xlabel font ',12'"
    @gp :- "set ylabel font ',12'"
    @gp :- "set title font ',12'"
    @gp :- "set key font ',9'" # border lt .05
    @gp :- "set key reverse center bmargin"

    subplt_indx =1
    for k in 1:num_features
        # metric_name = metric_name_vec[k]
        units = units_vector[k]
        function_labels = function_labels_vec[k]
        # function_name = replace(function_name_vec[k],"_" => "\\_")
        units_str,conversion_const = get_conversion_units(units;unit_type=unit_type)

        @gp :- subplt_indx title = "$(replace(function_name_vec[k],"_" => " "))"
        
        @gp :- xlab = "$(replace(function_name_vec[k],"_" => " ")) Variants" ylab = "$y_label $units_str"
        data = results_vec[k];
        data = data .* conversion_const
        function_labels = replace.(function_labels,"_" => "\\_");
        for i in 1:num_variants
            @gp :- data[i] " t '$(function_labels[i])' lc rgb '#$(hex(linCol[i]))' fs transparent solid 0.75 " :- 
        end

        subplt_indx +=1
    end
    

    @gp


    # function_labels = [el[:name][1] for key in collect(keys(results_dict))[2:end]];
    if !isnothing(filenamebase)
        plotID = "Multiplot"
        filename_vec = split(filenamebase,".")
        ext = filename_vec[end]
        new_file =  filename_vec[end-1]*plotID  
        filename_vec[end-1] =  new_file 
        filename_vec[end]= "."*filename_vec[end]
        new_filename = join(filename_vec)
        if ext == "png"
            Gnuplot.save(term = "pngcairo size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$svg_filename")
            end
        end
        if ext == "svg"
            Gnuplot.save(term = "svg size $(fig_size[1]),$(fig_size[2])", output = "$new_filename")
        end
        if ext == "pdf"
            Gnuplot.options.term = "qt  size 100,100"
            Gnuplot.save(term = "pdfcairo size $(100),$(100)", output = "$new_filename")
            if save_svg_copy
                svg_filename_vec = deepcopy(filename_vec)
                svg_filename_vec[end]= ".svg"
                svg_filename = join(svg_filename_vec)
                Gnuplot.save(term = "svg size $(100),$(100)", output = "$svg_filename")
            end
        end
    end
end

