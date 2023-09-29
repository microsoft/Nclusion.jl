coltitle <- 'Inferred Cell Type'
rowtitle <- 'Annotated Cell Type'

args <- commandArgs(trailingOnly=TRUE)

path_to_pct_table <- args[1]
save_path <- args[2]
r_library_path <- args[3]

if (r_library_path != 'NULL'){
.libPaths(r_library_path)
}

library(circlize)
library(devtools)
library(R.utils)
library(ComplexHeatmap)

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


