r_library_path <- NULL
path_to_data <- "tutorial_data/vanGalen2019/Seurat_AML.rds"

args <- commandArgs(trailingOnly=TRUE)

r_library_path <- args[1]

.libPaths(r_library_path)

library(devtools)
library(R.utils)
library(Seurat)
# library(SeuratDisk)
library(reticulate)

data <- readRDS(path_to_data)

all(data@meta.data$nCount_RNA == colSums(data@assays$RNA@layers$counts))

countsdf <- as.data.frame(data@assays$RNA@layers$counts)
metadf <- as.data.frame(data@meta.data)

rownames(countsdf) <- data@assays$RNA@features[["counts"]]
colnames(countsdf) <- data@assays$RNA@cells[["counts"]]

write.csv(countsdf, 'tutorial_data/vanGalen2019/galenAML_countsdf.csv')
write.csv(metadf, 'tutorial_data/vanGalen2019/galenAML_metadf.csv')

