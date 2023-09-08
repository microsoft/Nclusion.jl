.libPaths("/home/v-mahughes/r_packages/")
# install.packages('hdf5r')

# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")

library(devtools)
library(R.utils)
library(Seurat)
library(SeuratDisk)
library(reticulate)

# sessionInfo()


data <- readRDS('/home/v-mahughes/nclusion-datasets/galenAML/Seurat_AML.rds')
head(data)
print(typeof(data))
# print(data@meta.data)
# print(dim(data@raw.data))
# print(dim(data@meta.data))
# c <- GetAssayData(object = data, slot = "counts")
# print(dim(data@assays$RNA@layers$counts))
# print(head(data@meta.data))

# print('scaled')
# print(head(data@assays$RNA@layers$scale.data))
# print('log')
# print(head(data@assays$RNA@layers$data))
# print('raw')
# print(head(data@assays$RNA@layers$counts))

# data.seurat <- as.Seurat(data)

# adata <- Convert(from=data, to="anndata", filename="galen_aml_preprocessed.h5ad")

# SaveH5Seurat(data, filename = "aml_data.h5Seurat")

# Convert("~/SCoOP-sc/preprocess/aml_data.h5Seurat", dest = "h5ad")

# seuratobject_ad <- Convert(from=data, to="anndata", filename="seuratobject.h5ad")

all(data@meta.data$nCount_RNA == colSums(data@assays$RNA@layers$counts))

countsdf <- as.data.frame(data@assays$RNA@layers$counts)
metadf <- as.data.frame(data@meta.data)

rownames(countsdf) <- data@assays$RNA@features[["counts"]]
colnames(countsdf) <- data@assays$RNA@cells[["counts"]]
print(head(countsdf))
print(head(metadf))

write.csv(countsdf, 'galenAML_countsdf.csv')
write.csv(metadf, 'galenAML_metadf.csv')

