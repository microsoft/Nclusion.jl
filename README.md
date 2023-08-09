# Project

> This repo has been populated by an initial template to help get you started. Please
> make sure to update the content to build a great experience for community-building.

As the maintainer of this project, please make a few updates:

- Improving this README.MD file to provide a great experience
- Updating SUPPORT.MD with content about this project's support experience
- Understanding the security reporting process in SECURITY.MD
- Remove this section from the README

## Contributing

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## Trademarks

This project may contain trademarks or logos for projects, products, or services. Authorized use of Microsoft 
trademarks or logos is subject to and must follow 
[Microsoft's Trademark & Brand Guidelines](https://www.microsoft.com/en-us/legal/intellectualproperty/trademarks/usage/general).
Use of Microsoft trademarks or logos in modified versions of this project must not cause confusion or imply Microsoft sponsorship.
Any use of third-party trademarks or logos are subject to those third-party's policies.

## Introduction
Profiling phenotypic heterogenity in cellular populations has recieved a renewed interested recently due to the increased accuracy of single cell sequencing and the decrease in sequencing costs. As a result, new statistical methods have been developed to characterize the phenotypic heterogenity that result from single cell sequencing experiments such as single cell RNA-sequencing (scRNA-Seq) experiments. One popular way of carrying out this task is through clustering followed by performing additional statistical tests post-clustering to select salient genes that contribute to cluster identity. Current methods are limited in that the number of clusters in the data must be known prior to clustering and post-clustering variable selection does not leverage information about the cluster formation in order to elect salient genes.

We present a method for the Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION) - a sparse Bayesian Nonparametric model based on a Hierarchical Dirichlet Process mixture model. By combining a sparsity-promoting spike-and-slab prior on the cluster means and a nonparametric Hierarchical Dirichlet Process prior on the number of clusters, NCLUSION jointly learns the salient genes defining cluster identity and an optimal number of clusters needed to partition the data. Together with reasonable model approximations using varaiational inference, our proposed approach is scalable to large single cell RNA sequencing experiments.
## The Method
The Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION) is a sparse Bayesian Nonparametric model based on a Hierarchical Dirichlet Process mixture model which aims to identify single cell RNA sequencing clusters and the genes defining cluster identity. It does through the use of a sparsity-promoting spike-and-slab prior on the cluster means and a nonparametric Hierarchical Dirichlet Process prior on the number of clusters. The key idea behind the concept NCLUSION is the treatment of cluster formation as a shift in expression of particular genes away from a global mean. This mean shift is modeled by the spike-and-slab prior. The number of clusters is modeled via a stick-breaking proccess that allocates more clusters as the amount of variation in the data grows. The learning of the model parameters proceeds via a variational-EM algorithm that allows the method to scale with the number of cells in the data. As an overview of NCLUSION and its corresponding software implementation, we will assume that we have access to an expression matrix from a scRNA-seq study on N cells and G genes denoted as X where X is an N x G matrix of expression counts with G denoting the number of genes.

The goal of NCLUSION is to identify the phenotypic clusters in the data while jointly identifying the genes that define these clusters. To accomplish this, we asses the quality of clustering using both labels and label-free metrics. We also perform gene-pathway enrichment to asses the quality of elected genes.
## Installation

This package requires that Julia 1.6+ and its dependencies are installed. Instructions for installation can be found <a href="https://github.com/JuliaLang/julia"> here. </a>

To run the evaluation scripts we provide in the tutorials, the following R and python packages are required. We recommend that you install these R packages into their own library directory and the python packages into a new conda environment to avoid configuration issues. 

### R 
- <a href="https://cran.r-project.org/web/packages/devtools/index.html"> devtools </a>
- <a href="https://cran.r-project.org/web/packages/ggplot2/index.html"> ggplot2 </a>
- <a href="https://cran.r-project.org/web/packages/reticulate/index.html"> reticulate </a>
- <a href="https://cran.r-project.org/web/packages/R.utils/index.html"> R.utils </a>
- <a href="https://cran.r-project.org/web/packages/stringr/index.html"> stringr </a>
- <a href="https://cran.r-project.org/web/packages/Matrix/index.html"> Matrix </a>
- <a href="https://cran.r-project.org/web/packages/Seurat/index.html"> Seurat </a>
- <a href="https://cran.r-project.org/web/packages/dplyr/index.html"> dplyr </a>
- <a href="https://cran.r-project.org/web/packages/RColorBrewer/index.html"> RColorBrewer </a>
- <a href="https://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html"> SingleCellExperiment </a>
- <a href="https://cran.r-project.org/web/packages/tidyverse/index.html"> tidyverse </a>
- <a href="https://cran.r-project.org/web/packages/scales/index.html"> scales </a>
- <a href="https://cran.r-project.org/web/packages/cowplot/index.html"> cowplot </a>
- <a href="https://cran.r-project.org/web/packages/RCurl/index.html"> RCurl </a>
- <a href="https://cran.r-project.org/web/packages/optparse/index.html"> optparse </a>
- <a href="https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html"> ComplexHeatmap </a>
- <a href="https://cran.r-project.org/web/packages/circlize/index.html"> circilize </a>
- <a href="https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html"> clusterProfiler </a>
- <a href="https://bioconductor.org/packages/release/bioc/html/enrichplot.html"> enrichplot </a>
- <a href="https://cran.r-project.org/web/packages/geneset/index.html"> geneset </a>
- <a href="https://cran.r-project.org/web/packages/genekitr/index.html"> genekitr </a>
- <a href="https://cran.r-project.org/web/packages/patchwork/index.html"> patchwork </a>
- <a href="https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html"> org.Hs.eg.db </a>

These packages can be installed by entering the following commands into an R shell:
```
install.packages(c('devtools', 'ggplot2', 'reticulate','R.utils', 'stringr', 'Matrix', 'Seurat', 'dplyr', 'RcolorBrewer', 'tidyverse', 'scales', 'cowplot', 'RCurl', 'optparse', 'circilize', 'geneset', 'genekitr', 'patchwork'))
```
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "ComplexHeatmap", "clusterProfiler", "enrichplot", "org.Hs.eg.db"))
```

### python 

- <a href="https://scanpy.readthedocs.io/en/stable/"> scanpy </a>
- <a href="https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html"> pandas </a>
- <a href="https://matplotlib.org/stable/users/installing/index.html"> matplotlib </a> 

These packages can be installed with the following commands in an activated conda environment:

```
conda install -c conda-forge scanpy
```
```
conda install -c anaconda pandas
```
```
conda install -c conda-forge matplotlib
```
## Questions and Feedback

If you have any questions or concerns regarding NCLUSION or the tutorials, please contact <a href="mailto:chibuikem_nwizu@brown.edu"> Chibuikem Nwizu</a>, <a href="mailto:lcrawford@microsoft.com"> Lorin Crawford</a>, or <a href="mailto:v-mahughes@microsoft.com"> Madeline Hughes</a>. All feedback on the repository, manuscript, and tutorials is appreciated.

## References


