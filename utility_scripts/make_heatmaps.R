option_list = list(
  make_option(c("-f", "--path_to_pct_table"), type="character", default=NULL, 
              help="path to the cell distribution table", metavar="character"),
  make_option(c("-c", "--column_title"), type="character", default="Inferred Cell Type", 
              help="heatmap column title", metavar="character"),
  make_option(c("-r", "--row_title"), type="character", default="Annotated Cell Type", 
              help="heatmap row title", metavar="character"),
    make_option(c("-s", "--save_path"), type="character", default="complex_heatmap.pdf", 
              help="save path of heatmap pdf", metavar="character"),
  make_option(c("-p", "--purity_score"), type="logical", default=FALSE, 
              help="if included, then it makes the purity score heatmap", metavar="logical")
); 

path_to_pct_table <- NULL
coltitle <- 'Inferred Cell Type'
rowtitle <- 'Annotated Cell Type'
save_path <- 'cell_type_distribution_heatmap.pdf'
r_library_path <- NULL

args <- commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]
keys <- attachLocally(args)

.libPaths(r_library_path)
library(circlize)
library(reticulate)
library(devtools)
library(ClueR)
library(R.utils)
library(ComplexHeatmap)
pd <- import("pandas")

data <- read.csv(path_to_pct_table, row.names=1)

colnames(data) <- seq(from=1,to=length(colnames(data)))
col_fun = colorRamp2(c(0, 0.1, 1), c("white", "pink", "purple"))
col_fun(seq(0, 1))

mat <- as.matrix(data)
mat_name <- ' ' 

ht <- Heatmap(mat, name = mat_name, col = col_fun, width = ncol(mat)*unit(5, "mm"), height = nrow(mat)*unit(5, "mm"), cluster_rows=FALSE, cluster_columns=FALSE, column_names_side="bottom", row_names_side="left", column_title = coltitle, row_title = rowtitle, column_title_side = "bottom", column_names_rot = 90, column_names_gp = gpar(fontsize=8), row_names_centered = TRUE, border_gp = gpar(col='black', fontsize=2), column_names_centered = TRUE)


ht_draw <- draw(ht)
w <- ComplexHeatmap:::width(ht_draw)
w <- convertX(w, "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(ht_draw)
h <- convertY(h, "inch", valueOnly = TRUE)

pdf(file=save_path, width=w, height=h)
draw(ht)
dev.off()


