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

## The Method

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


