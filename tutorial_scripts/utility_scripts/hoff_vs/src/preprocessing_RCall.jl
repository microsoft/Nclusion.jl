function R_lognorm_FeatSelect_Scale_noMetaDataGenenames(data_prep;num_var_feat = 1070, scale_fator = 10000,only_var_feat=false,seed =12345)
    @rput data_prep;
    @rput num_var_feat; 
    @rput scale_fator; 
    @rput seed;
    if only_var_feat
        R"""
            library(stringr)
            library(Seurat)
            library(dplyr)
            # library(scater)
            library(SingleCellExperiment)
            library(tidyverse)
            library(Matrix)
            library(scales)
            library(cowplot)
            library(RCurl)
            set.seed(seed)
            rownames(data_prep) <- data_prep[,1]
            data_prep <- data_prep[, -1]
            dat <- CreateSeuratObject(counts = data_prep, min.cells = 1, min.features = 1)
            dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = scale_fator)
            dat_pca <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = num_var_feat)
            dat_pca <- ScaleData(dat_pca, features = dat_pca@assays$RNA@var.features) #
            x_int_scaled_postp <- dat_pca@assays$RNA@scale.data
            x_int_lognorm_postp <- as.matrix(dat_pca@assays$RNA@data)
            gene_names_lognorm_postp <- rownames(dat_pca@assays$RNA@data)
            gene_names_scaled_postp <- rownames(dat_pca@assays$RNA@scale.data)
        """
    else
        R"""
            library(stringr)
            library(Seurat)
            library(dplyr)
            # library(scater)
            library(SingleCellExperiment)
            library(tidyverse)
            library(Matrix)
            library(scales)
            library(cowplot)
            library(RCurl)
            set.seed(seed)
            rownames(data_prep) <- data_prep[,1]
            data_prep <- data_prep[, -1]
            dat <- CreateSeuratObject(counts = data_prep, min.cells = 1, min.features = 1)
            dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = scale_fator)
            dat_pca <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = num_var_feat)
            dat_pca <- ScaleData(dat_pca, features = rownames(dat_pca)) #dat_pca@assays$RNA@var.features
            x_int_scaled_postp <- dat_pca@assays$RNA@scale.data
            x_int_lognorm_postp <- as.matrix(dat_pca@assays$RNA@data)
            gene_names_lognorm_postp <- rownames(dat_pca@assays$RNA@data)
            gene_names_scaled_postp <- rownames(dat_pca@assays$RNA@scale.data)
        """
    end

    @rget x_int_lognorm_postp;
    @rget x_int_scaled_postp;
    @rget gene_names_lognorm_postp;
    @rget gene_names_scaled_postp;

    return x_int_lognorm_postp,x_int_scaled_postp,gene_names_lognorm_postp,gene_names_scaled_postp
end

function R_lognorm_FeatSelect_Scale(data_prep,metadata_prep,gene_names_prep;num_var_feat = 1070, scale_fator = 10000,only_var_feat=false,seed =12345)
    @rput data_prep;
    @rput metadata_prep;
    @rput gene_names_prep;
    @rput num_var_feat; 
    @rput scale_fator; 
    @rput seed;
    if only_var_feat
        R"""
            library(stringr)
            library(Seurat)
            library(dplyr)
            # library(scater)
            library(SingleCellExperiment)
            library(tidyverse)
            library(Matrix)
            library(scales)
            library(cowplot)
            library(RCurl)
            set.seed(seed)
            # rownames(data_prep) <- data_prep[,1]
            # data_prep <- data_prep[, -1]
            rownames(metadata_prep) <- metadata_prep[,1]
            rownames(data_prep) <- gene_names_prep[,1]
            dat <- CreateSeuratObject(counts = data_prep, meta.data=metadata_prep, min.cells = 1, min.features = 1)
            dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = scale_fator)
            dat_pca <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = num_var_feat)
            dat_pca <- ScaleData(dat_pca, features = dat_pca@assays$RNA@var.features) #
            x_int_scaled_postp <- dat_pca@assays$RNA@scale.data
            x_int_lognorm_postp <- as.matrix(dat_pca@assays$RNA@data)
            gene_names_lognorm_postp <- rownames(dat_pca@assays$RNA@data)
            gene_names_scaled_postp <- rownames(dat_pca@assays$RNA@scale.data)
            metadata_postp <- dat_pca@meta.data
        """
    else
        R"""
            library(stringr)
            library(Seurat)
            library(dplyr)
            # library(scater)
            library(SingleCellExperiment)
            library(tidyverse)
            library(Matrix)
            library(scales)
            library(cowplot)
            library(RCurl)
            set.seed(seed)
            # rownames(data_prep) <- data_prep[,1]
            # data_prep <- data_prep[, -1]
            # dat <- CreateSeuratObject(counts = data_prep, min.cells = 1, min.features = 1)
            # dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
            # dat_pca <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 3000)
            # dat_pca <- ScaleData(dat_pca, features = rownames(dat_pca)) #dat_pca@assays$RNA@var.features
            # x_int_scaled_postp <- dat_pca@assays$RNA@scale.data
            # x_int_lognorm_postp <- as.matrix(dat_pca@assays$RNA@data)
            # gene_names_lognorm_postp <- rownames(dat_pca@assays$RNA@data)
            # gene_names_scaled_postp <- rownames(dat_pca@assays$RNA@scale.data)
            # metadata_postp <- dat_pca@meta.data

            rownames(metadata_prep) <- metadata_prep[,1]
            rownames(data_prep) <- gene_names_prep[,1]
            dat <- CreateSeuratObject(counts = data_prep, meta.data=metadata_prep, min.cells = 1, min.features = 1)
            dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = scale_fator)
            dat_pca <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = num_var_feat)
            dat_pca <- ScaleData(dat_pca, features = rownames(dat_pca)) #
            x_int_scaled_postp <- dat_pca@assays$RNA@scale.data
            x_int_lognorm_postp <- as.matrix(dat_pca@assays$RNA@data)
            gene_names_lognorm_postp <- rownames(dat_pca@assays$RNA@data)
            gene_names_scaled_postp <- rownames(dat_pca@assays$RNA@scale.data)
            metadata_postp <- dat_pca@meta.data
        """
    end

    @rget x_int_lognorm_postp;
    @rget x_int_scaled_postp;
    @rget metadata_postp;
    @rget gene_names_lognorm_postp;
    @rget gene_names_scaled_postp;

    return x_int_lognorm_postp,x_int_scaled_postp,metadata_postp,gene_names_lognorm_postp,gene_names_scaled_postp
end

function R_VariableFeatSelect(data_prep, num_var_feat;gene_names=nothing, scale_fator = 10000,seed =12345)
    if isnothing(gene_names)
        gene_names = data_prep[!,1];
        data_prep = data_prep[!,2:end]; 
    elseif typeof(gene_names) <: Vector
        gene_names = DataFrame(gene_id = gene_names)
    end
    @rput data_prep;
    @rput gene_names;
    @rput num_var_feat; 
    @rput scale_fator; 
    @rput seed;
    R"""
        library(stringr)
        library(Seurat)
        library(dplyr)
        library(scater)
        library(SingleCellExperiment)
        library(tidyverse)
        library(Matrix)
        library(scales)
        library(cowplot)
        library(RCurl)
        set.seed(seed)
        rownames(data_prep) <- gene_names[,1]
        dat <- CreateSeuratObject(counts = data_prep, min.cells = 1, min.features = 1)
        dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = scale_fator)
        dat_pca <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = num_var_feat)
        var_genes <- dat_pca@assays$RNA@var.features #
    """

    @rget var_genes;
    
    return var_genes
end
