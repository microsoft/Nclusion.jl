r_library_path <- NULL
organism = "org.Hs.eg.db"
path_to_data <- NULL
path_to_pips<- NULL
path_to_labels <- NULL
outfile_base <- NULL
thresh = 0.5
dataset_name <- ''

# args <- commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]
args <- commandArgs(trailingOnly=TRUE)
# keys <- attachLocally(args)
# print(keys)
# str(mget(keys, envir=globalenv()))
path_to_data <- args[1]
path_to_pips <- args[2]
path_to_labels <- args[3]
outfile_base <- args[4]

if (length(args) > 4){
dataset_name <- args[5]
}

if (length(args) == 6){
r_library_path <- args[6]
}

if (!is.null(r_library_path)){
.libPaths(r_library_path)
}

seed = 12345
library(devtools)
library(ggplot2)
library(reticulate)
library(R.utils)
library(stringr)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(optparse)
set.seed(seed)

library(ComplexHeatmap)
library(circlize)

library(clusterProfiler)
library(enrichplot)
library(geneset)
library(genekitr)
library(patchwork)
library(organism, character.only = TRUE)

# path_to_data <- "/home/v-mahughes/nclusion_preprocessed_data/galenAML2019/2000hvgs_galenAML_preprocessed.h5ad"
# path_to_pips<- "/home/v-mahughes/RESULTS_81523/galenAML/NCLUSION/2000hvgs/2000G-2023-08-10T073505-pips.csv"
# path_to_labels <- "/home/v-mahughes/RESULTS_81523/galenAML/NCLUSION/2000hvgs/galen-AML_scanpy-2000G_nclusion_2000genes-2023-08-10T073505.csv"
# outfile_base <- "/home/v-mahughes/nclusion_figures/galen_fig4_81623/"

data <- read.csv(path_to_pips, row.names = 1, header= TRUE)
num_col = dim(data)[2]
num_row = dim(data)[1]
mat = data.matrix(data)
pips = mat
K = dim(pips)[1]
G = dim(pips)[2]

# Adjust PIPs to reflect significance more accurately
genes_names <- colnames(pips)
sig_genes_df = pips >= thresh
genes_sig_sum = colSums(sig_genes_df)
Kmax_occupied = max(genes_sig_sum)
gene_cluster_utilization_occ = (1 - (genes_sig_sum / Kmax_occupied))/((1 - (1 /Kmax_occupied)))
adj_pips_mat_occ =gene_cluster_utilization_occ * t(pips)

labels <- read.csv(file = path_to_labels)


# Filtering for genes above adjusted PIP value threshold specified
ro_tc <- unique(labels[,dim(labels)[2]])
row_bool = apply(adj_pips_mat_occ[,ro_tc] >= thresh, 1, any)

# this reorders the columns, but initially you wont know that order do id set it to the identity permutation vec
num_clus = length(ro_tc)
permutation_vec = c(1:num_clus)
new_clust <- ro_tc

# Generate Adjusted PIP Heatmap (Only shows genes above specified PIP threshold)
only_sig_genes_names = genes_names[row_bool]
only_sig_mat = t(adj_pips_mat_occ[row_bool,ro_tc])
rownames(only_sig_mat) <- new_clust
data_set<-""
suffix <- ""
file_name1 = paste0(outfile_base,data_set,"heatmap_AdjustedPips_KbyG_defaultclustering.pdf")
pdf(file=file_name1, width = 8.5, height = 11)

anno_labels<-sapply(ro_tc,function(x)paste0(x))[permutation_vec] 

anno_fill <- c("#D70000", "#8C3CFF", "#028800", "#00ACC7", "#E7A500", "#FF7FD1", "#6C004F", "#583B00", "#005759", "#15E18D", "#0000DD", "#A2766A", "#BCB7FF", "#C004B9", "#645473", "#790000", "#0774D8", "#739B7D", "#FF7852", "#004B00", "#8F7B01", "#F3007B", "#8FBA00", "#A67BB8", "#5A02A3", "#E3AFAF", "#A03A52", "#A2C8C8", "#9E4B00", "#546745", "#BBC389", "#5F7B88", "#60383C", "#8388FF", "#390000", "#E353FF", "#305382", "#7FCAFF", "#C5668F", "#00816A", "#929EB7", "#CC7407", "#7F2B8E", "#00BEA4", "#2DB152", "#4E33FF", "#00E500", "#FF00CE", "#C85848", "#E59CFF", "#1DA1FF", "#6E70AB", "#C89A69", "#78573B", "#04DAE6", "#C1A3C4", "#FF6A8A", "#BB00FE", "#925380", "#9F0274", "#94A150", "#374425", "#AF6DFF", "#596D00", "#FF3147", "#838057", "#006D2E", "#8956AF", "#5A4AA3", "#773516", "#86C39A", "#5F1123", "#D58581", "#A42918", "#0088B1", "#CB0044", "#FFA056", "#EB4E00", "#6C9700", "#538649", "#755A00", "#C8C440", "#92D370", "#4B9894", "#4D230D", "#61345C", "#8400CF", "#8B0031", "#9F6E32", "#AC8499", "#C63189", "#025438", "#086B84", "#87A8EC", "#6466EF", "#C45DBA", "#019F70", "#815159", "#836F8C", "#B3C0DA", "#B99129", "#FF97B2", "#A793E1", "#698DBE", "#4C5001", "#4802CC", "#61006E", "#456A66", "#9D5743", "#7BACB5", "#CD84BD", "#0054C1", "#7B2F4F", "#FB7C00", "#34C000", "#FF9C88", "#E1B769", "#536177", "#5C3A7C", "#EDA5DA", "#F053A3", "#5D7E69", "#C47750", "#D14868", "#6E00EB", "#1F3400", "#C14104", "#6DD5C2", "#46709F", "#A201C4", "#0A8289", "#AFA601", "#A65C6B", "#FE77FF", "#8B85AE", "#C77FE9", "#9AAB85", "#876CD9", "#01BAF7", "#AF5ED2", "#59512B", "#B6005F", "#7CB66A", "#4985FF", "#00C282", "#D295AB", "#A34BA8", "#E306E3", "#16A300", "#392E00", "#843033", "#5E95AA", "#5A1000", "#7B4600", "#6F6F31", "#335826", "#4D60B6", "#A29564", "#624028", "#45D458", "#70AAD0", "#2E6B4E", "#73AF9E", "#FD1500", "#D8B492", "#7A893B", "#7DC6D9", "#DC9137", "#EC615E", "#EC5FD4", "#E57BA7", "#A66C98", "#009744", "#BA5F22", "#BCAD53", "#88D830", "#873573", "#AEA8D2", "#E38C63", "#D1B1EC", "#37429F", "#3ABEC2", "#669D4D", "#9E0399", "#4E4E7A", "#7B4C86", "#C33531", "#8D6677", "#AA002D", "#7F0175", "#01824D", "#734A67", "#727791", "#6E0099", "#A0BA52", "#E16E31", "#C56A71", "#6D5B96", "#A33C74", "#326200", "#880050", "#335869", "#BA8D7C", "#1959FF", "#919202", "#2C8BD5", "#1726FF", "#21D3FF", "#A490AF", "#8B6D4F", "#5E213E", "#DC03B3", "#6F57CA", "#652821", "#AD7700", "#A3BFF7", "#B58446", "#9738DC", "#B25194", "#7242A3", "#878FD1", "#8A70B1", "#6BAF36", "#5A7AC9", "#C79FFF", "#56841A", "#00D6A7", "#824739", "#11431D", "#5AAB75", "#915B01", "#F64570", "#FF9703", "#E14231", "#BA92CF", "#34584D", "#F8807D", "#913400", "#B3CD00", "#2E9FD3", "#798B9F", "#51817D", "#C136D7", "#EC0553", "#B9AC7E", "#487032", "#849565", "#D99D89", "#0064A3", "#4C9078", "#8F6198", "#FF5338", "#A7423B", "#006E70", "#98843E", "#DCB0C8")[permutation_vec]

#### MAKE ADJUSTED PIP HEATMAP
ht = Heatmap(t(only_sig_mat)[,permutation_vec], name = "adj. pips", col = colorRamp2(c(0, 0.9), c("#e4e1e1", "#b40000")), column_title = "Adj. PIPS", cluster_rows = TRUE,  row_dend_side = "right",cluster_columns = FALSE, column_names_rot = 90, use_raster=FALSE ,show_column_names = FALSE,row_names_side = "left",column_split = c(1:length(new_clust)),bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =anno_fill),labels = anno_labels, labels_gp = gpar(col = "white", fontsize = 10))))

ht = draw(ht)
dev.off()
col_order_ht <- column_order(ht)
ro_order_ht <- row_order(ht)

# Load the adata 
scale_fator <- 10000
sc <- import("scanpy")
np <- import("numpy")
pd <- import("pandas")

adata <- sc$read_h5ad(path_to_data)
labels <- read.csv(file = path_to_labels)
X <- as.data.frame(t(adata$raw$X)) # adata$X
rownames(X) <- rownames(adata$var)
colnames(X) <- adata$obs_names$values
genes_names <- colnames(pips)
num_genes = length(rownames(X))
inclusion_bool <- rep(FALSE, num_genes)
for (i in 1:num_genes){
  x = rownames(X)[i]
    if(str_replace_all(x, "[^[:alnum:]]", "")  %in% str_replace_all(only_sig_genes_names, "[^[:alnum:]]", "")){
      inclusion_bool[i]<- TRUE
      cat("Element present in vector at pos",i)
      cat("\n")
    }
}
X <-  X[inclusion_bool,]

seurat_dataset <- CreateSeuratObject(counts = X, min.cells = 0, min.features = 0)
Idents(object = seurat_dataset)<- labels[,dim(labels)[2]]
seurat_dataset <- NormalizeData(seurat_dataset, normalization.method = "LogNormalize", scale.factor = scale_fator)
seurat_dataset <- ScaleData(seurat_dataset)
seurat_dataset[["cell.types"]] <- Idents(object = seurat_dataset)
seurat_dataset[["new.cell.types"]] <- seurat_dataset[["cell.types"]]

levels(seurat_dataset$new.cell.types) <- new_clust
scaled_mat <- seurat_dataset[["RNA"]]@scale.data
counts_mat <- seurat_dataset@assays$RNA@counts
num_cells = length(labels$cell_id)
inclusion_bool_cells <- rep(FALSE, num_cells)
for (i in 1:num_cells){
x = labels$cell_id[i]
  if(x  %in%  colnames(scaled_mat)){
      inclusion_bool_cells[i]<- TRUE
    next
  }
cat("Element not present in vector at pos",i)
cat("\n")
}

clus_id = ro_tc
num_clus = length(clus_id)
num_g = dim(scaled_mat)[1]

mean_mat <- matrix(, nrow = num_clus, ncol = num_g) 
var_mat <- matrix(, nrow = num_clus, ncol = num_g) 

not_mean_mat <- matrix(, nrow = num_clus, ncol = num_g) 
not_var_mat <- matrix(, nrow = num_clus, ncol = num_g) 

for (i in 1:num_clus){
  clus <- clus_id[i]
  new_bool = labels$inferred_label == clus
  new_not_bool = labels$inferred_label != clus

  mean_mat[i,] <- apply(counts_mat[,new_bool],1, mean)
  var_mat[i,] <- apply(counts_mat[,new_bool],1, var)
  not_mean_mat[i,] <- apply(counts_mat[,new_not_bool],1, mean)
  not_var_mat[i,] <- apply(counts_mat[,new_not_bool],1, var)

}


colnames(mean_mat) <- colnames(only_sig_mat) 
rownames(mean_mat) <- rownames(only_sig_mat) 

colnames(var_mat) <- colnames(only_sig_mat) 
rownames(var_mat) <- rownames(only_sig_mat) 

colnames(not_mean_mat) <- colnames(only_sig_mat) 
rownames(not_mean_mat) <- rownames(only_sig_mat) 

colnames(not_var_mat) <- colnames(only_sig_mat) 
rownames(not_var_mat) <- rownames(only_sig_mat) 

effect_size_mat <- matrix(, nrow = num_clus, ncol = num_g)
effect_size_calc <- function(pop_mean,notpop_mean,pop_var,notpop_var){
  (pop_mean-notpop_mean) / sqrt(0.5 * (pop_var + notpop_var))
}

for (i in 1:num_clus){
  effect_size_mat[i,] <- effect_size_calc(mean_mat[i,],not_mean_mat[i,],var_mat[i,],not_var_mat[i,])
}

colnames(effect_size_mat) <- colnames(only_sig_mat) 
rownames(effect_size_mat) <- rownames(only_sig_mat) 
signed_effect_size_mat <- matrix(, nrow = num_clus, ncol = num_g)
for (i in 1:num_clus){
  signed_effect_size_mat[i,effect_size_mat[i,]==0] <- "0"#3
  signed_effect_size_mat[i,effect_size_mat[i,]>0] <- "+"#1
  signed_effect_size_mat[i,effect_size_mat[i,]<0] <- "-" #2
}

colnames(signed_effect_size_mat) <- colnames(only_sig_mat) 
rownames(signed_effect_size_mat) <- rownames(only_sig_mat) 
c("+","-","0") 


#### MAKE EFFECT SIZE HEATMAP

file_name_ = paste0(outfile_base,data_set,"heatmap_EffectSizeSignMat_KbyG_basedonAdjPips",suffix,".pdf")
pdf(file=file_name_, width = 8.5, height = 11)
he =Heatmap(t(signed_effect_size_mat)[,permutation_vec], name = "Effect Size Sign", col = structure( c("#ffdb3b","#8f0388", "#000000"), names = c("+","-","0") ),cluster_columns = FALSE, row_order = ro_order_ht,  column_title = "Sign of Effect Size", column_names_rot = 90, use_raster=FALSE , show_column_names = FALSE,row_names_side = "left",column_split = c(1:Kmax_occupied),bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =anno_fill),labels = anno_labels, labels_gp = gpar(col = "white", fontsize = 10))))
he = draw(he)
dev.off()


#### MAKE VIOLIN PLOTS TO SHOW MODULE EXPRESSION 
sdata <- seurat_dataset
cluster_idx <- ro_tc
new_cluster_idx <- ro_tc


for (i in c(1:num_clus)){ 
  clus <- cluster_idx[i]
  new_clus <- new_cluster_idx[i]
  pltname= paste0("Cluster_",clus,"Genes")

  genes_names2 = colnames(signed_effect_size_mat)
  print(colnames(adj_pips_mat_occ))
  row_bool1 = adj_pips_mat_occ[,clus] >= thresh
  row_bool2 = signed_effect_size_mat[rownames(signed_effect_size_mat)==paste0(clus),] == "+"

  only_sig_genes_names = genes_names2[row_bool2]
 
  flist <- sub('[.]', '-', only_sig_genes_names)
  flist2 <- list(flist)

  sdata <- AddModuleScore(object = sdata,features =flist2,ctrl = 2,name=pltname)

  norm_pltname <- paste0("NormCluster_",clus,"Genes1")

  sdata[[norm_pltname]] <- (sdata[[paste0(pltname,"1")]] - min(sdata[[paste0(pltname,"1")]]))/(max(sdata[[paste0(pltname,"1")]]) - min(sdata[[paste0(pltname,"1")]]))

 
  sdata_reorder = sdata
  levels(seurat_dataset$new.cell.types) <- new_clust[permutation_vec]

   
  file_name_feature = paste0(outfile_base,"violinplot_normedmodulescore_colored_nolegend_",pltname,"-",suffix,".pdf")
  pdf(file=file_name_feature, width = 8.5, height = 7)
  abc <- ggplot(sdata_reorder[[c('new.cell.types',norm_pltname)]], aes(x = new.cell.types, y= .data[[norm_pltname]],fill=new.cell.types)) + geom_violin(trim=FALSE, scale = 'width')+ scale_fill_manual(values=anno_fill) + theme_classic()+ geom_boxplot(width=0.1, fill="white", outlier.colour=rgb(.5,.5,.5, 0.5)) + ylim(-0.1, 1.1)
  print(abc)
  dev.off()

  file_name_feature = paste0(outfile_base,"violinplot_modulescore_colored_nolegend_",pltname,"-",suffix,".pdf")
  pdf(file=file_name_feature, width = 8.5, height = 7)
  abc <- ggplot(sdata_reorder[[c('new.cell.types',paste0(pltname,"1"))]], aes(x = new.cell.types, y= .data[[paste0(pltname,"1")]],fill=new.cell.types)) + geom_violin(trim=FALSE, scale = 'width')+ scale_fill_manual(values=anno_fill) + theme_classic()+ geom_boxplot(width=0.1, fill="white", outlier.colour=rgb(.5,.5,.5, 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab(element_blank()) + ylab(element_blank()) + theme(legend.position="none")
  print(abc)
  dev.off()


}



#### FIG 3F

for (i in c(1:num_clus)){ 
c = i
pltname= paste0("Cluster_",cluster_idx[c],"Genes")
clus <- paste0(cluster_idx[c])
row_bool2 = only_sig_mat[clus,] >= thresh & signed_effect_size_mat[clus,] == "+"
only_sig_genes_names2 = colnames(only_sig_mat)[row_bool2]
flist <- sub('[.]', '-', only_sig_genes_names2)

genekitr::transId(id=flist , transTo = "entrez", org = "human")
sym2ent <- genekitr::transId(id=flist , transTo = "entrez", org = "human")

gs <- geneset::getGO(org = "human",ont = "bp")

ego <- genORA(sym2ent$entrezid,
              geneset = gs,
              p_cutoff = 0.05,
              q_cutoff = 0.05
)
file_name_feature_lolli = paste0(outfile_base,"GO_bp_enrichmentplots_",pltname,"-",suffix,".pdf")
pdf(file=file_name_feature_lolli, width = 5, height = 11)
efg <-ego %>% arrange(desc(-log2(qvalue))) %>% head(10)  %>%  ggplot( aes(x=-log2(qvalue) , y= reorder(Description, -log2(qvalue))))  +  geom_segment( aes( xend=0, yend=Description)) + geom_point( size=4, color=anno_fill[new_cluster_idx[c]]) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  + xlab("-Log(FDR q value)") + ylab(element_blank()) + ggtitle(paste0("Cluster ",cluster_idx[c],"\n Biological Pathway Enrichment \n n = ",length(only_sig_genes_names2)," genes")) +theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,20)+ scale_y_discrete(labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))
print(efg)
dev.off()
}