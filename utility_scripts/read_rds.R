r_library_path <- NULL
path_to_data <- NULL 

args <- commandArgs(trailingOnly=TRUE)

r_library_path <- args[1]
path_to_data <- args[2]

.libPaths(r_library_path)

library(devtools)
library(R.utils)
library(Seurat)
library(reticulate)

data <- readRDS(path_to_data)

all(data@meta.data$nCount_RNA == colSums(data@assays$RNA@layers$counts))

countsdf <- as.data.frame(data@assays$RNA@layers$counts)
metadf <- as.data.frame(data@meta.data)

rownames(countsdf) <- data@assays$RNA@features[["counts"]]
colnames(countsdf) <- data@assays$RNA@cells[["counts"]]

write.csv(countsdf, 'tutorial_data/vanGalen2019/galenAML_countsdf.csv')
write.csv(metadf, 'tutorial_data/vanGalen2019/galenAML_metadf.csv')

