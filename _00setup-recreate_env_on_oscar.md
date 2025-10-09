0. Deactivate any active virtual envionments and enter the directory where you will create the environments. 
1. Ensure the following files are in the directory:
   ```sh
        pyproject.toml
        Project.toml
   ```
2. Ensure that the following file(s)/director(y/ies) are NOT present in the directory
    ```sh
        pyproject.toml
        Manifest.toml
        poetry.lock
        .venv/  #may need to change .venv to actual environment name if thats the environments name
        pyjuliapkg/
   ```
   If they are, please delete. 
3. Create the following variable
   ```sh
   CURRDIR=$(pwd)
   ```
4. In the ``pyproject.toml`` file, ensure that the following values are set to the following in the ``[tool.poetry]`` section
   ```toml
        name = "{directory_name}"
        package-mode = false
   ```
   where ``"{directory_name}"`` is the name of the directory that contains the files
5. Load the required modules on oscar
   ```sh
        module load python/3.11.0s-ixrhc3q; module load r/4.4.0-yycctsj; module load llvm/16.0.2; module load pcre2/10.42 texlive/20220321; module load cmake/3.26.3; module load libgit2/1.6.4; module load geos/3.11.2; module load libpng/1.6.39; module load gdal/3.7.0 proj/9.2.0; module load cuda/12.2.0 cudnn/8.9.6.50 openssl libarchive/3.6.2; module load graphviz inkscape; module load hdf5; module load gsl;module load julia;module load jags;export LD_PRELOAD=/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_def.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_avx2.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_core.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_lp64.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_thread.so:/gpfs/runtime/opt/intel/2020.2/lib/intel64_lin/libiomp5.so:/lib64/libssl.so.3:/lib64/libssl.so:/lib64/libcrypto.so.3;
   ```
6. To use the ``pyproject.toml`` file and need to install packages in a new virtual environment, here are the different ways you can do it:
   1.  Using poetry (Recommended since was used the poetry manages the project)
       If your ``pyproject.toml`` is managed by Poetry, create and activate a new virtual environment, then install dependencies:
       ```sh
            # Ensure Poetry is installed (https://python-poetry.org/docs/#installing-with-pipx)
            pipx install poetry
            pipx inject poetry poetry-plugin-shell
            # Create a new virtual environment and install dependencies
            poetry install
            #Usually creates an environment with the name .venv. may need to change actual environment if thats the environments name
       ``` 
   2.  If the ``pyproject.toml`` specifies dependencies via \[build-system\] but does not use poetry:
       ```sh
            # Create a virtual environment
            python -m venv .venv
            source .venv/bin/activate  # On Windows use: .venv\Scripts\activate
            # Use pip-tools to extract dependencies and install them
            pip install pip-tools
            pip-compile pyproject.toml
            pip install -r requirements.txt
       ```
       This approach extracts dependencies from ``pyproject.toml`` into a ``requirements.txt`` file and installs them with pip.
   3.  For others ways see [here](https://chatgpt.com/share/67b0e9fa-75fc-800a-a76e-15c8d60c7710)
7. Next, start the new environment
   ```sh
        source $CURRDIR/.venv/bin/activate
   ```
8. Next, we must unset the library path to allow us to create the Julia environment:
   ```sh
        unset LD_LIBRARY_PATH
        julia --project=$CURRDIR
   ```
9.  In Julia, we load enviroment via the ``Project.toml`` file:
   ```julia
        using Pkg
        Pkg.instantiate()
        Pkg.precompile()
        using PythonCall
        Pkg.resolve();
        Pkg.pin(name="OpenSSL_jll", version="3.0.15");
        Pkg.pin(name="PythonCall", version="0.9.24"); 
        Pkg.resolve();
        Pkg.instantiate();
        Pkg.precompile();
   ```
10. Next, we must connect python and julia. We only have to do this one time:
    1.  In the command line in OSCAR:
        ```sh
            module purge; module load python/3.11.0s-ixrhc3q; module load r/4.4.0-yycctsj; module load llvm/16.0.2; module load pcre2/10.42 texlive/20220321; module load cmake/3.26.3; module load libgit2/1.6.4; module load geos/3.11.2; module load libpng/1.6.39; module load gdal/3.7.0 proj/9.2.0; module load cuda/12.2.0 cudnn/8.9.6.50 openssl libarchive/3.6.2; module load graphviz inkscape; module load hdf5; module load gsl;module load julia;module load jags;export LD_PRELOAD=/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_def.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_avx2.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_core.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_lp64.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_thread.so:/gpfs/runtime/opt/intel/2020.2/lib/intel64_lin/libiomp5.so:/lib64/libssl.so.3:/lib64/libssl.so:/lib64/libcrypto.so.3;
        ```
    2.  Start python as follows:
        ```sh
            $CURRDIR/.venv/bin/python -u -X juliapkg-project=$CURRDIR -X juliacall-threads=auto -X juliacall-handle-signals=yes -X juliapkg-offline=yes  #may need to change actual environment if thats the environments name
        ```
        where ``$CURRDIR`` is the full path to the current directory
    3.  In python:
        ```python
            import juliapkg
            from juliacall import Main as jl
        ```
11. Reset everything:
    ```sh
        deactivate
        module purge
        module load python/3.11.0s-ixrhc3q; module load r/4.4.0-yycctsj; module load llvm/16.0.2; module load gnuplot/5.4.3; module load pcre2/10.42 texlive/20220321; module load cmake/3.26.3; module load libgit2/1.6.4; module load geos/3.11.2; module load libpng/1.6.39; module load gdal/3.7.0 proj/9.2.0; module load cuda/12.2.0 cudnn/8.9.6.50 openssl libarchive/3.6.2; module load graphviz inkscape; module load hdf5; module load gsl;module load julia;module load jags;export LD_PRELOAD=/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_def.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_avx2.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_core.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_lp64.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_thread.so:/gpfs/runtime/opt/intel/2020.2/lib/intel64_lin/libiomp5.so;
        source $CURRDIR/.venv/bin/activate
    ```
12. (Optional) Create an alias for loading the environment faster
    ```bash
        echo "alias start_nclusion_manuscript='CURRDIR=/absolute/path/to/curr/directory; cd $CURRDIR;module load python/3.11.0s-ixrhc3q; module load r/4.4.0-yycctsj; module load llvm/16.0.2; module load gnuplot/5.4.3; module load pcre2/10.42 texlive/20220321; module load cmake/3.26.3; module load libgit2/1.6.4; module load geos/3.11.2; module load libpng/1.6.39; module load gdal/3.7.0 proj/9.2.0; module load cuda/12.2.0 cudnn/8.9.6.50 openssl libarchive/3.6.2; module load graphviz inkscape; module load hdf5; module load gsl;module load julia;module load jags;export LD_PRELOAD=/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_def.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_avx2.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_core.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_lp64.so:/gpfs/runtime/opt/intel/2020.2/mkl/lib/intel64/libmkl_intel_thread.so:/gpfs/runtime/opt/intel/2020.2/lib/intel64_lin/libiomp5.so;source $CURRDIR/.venv/bin/activate; '" >> ~/.bash_aliases #may need to change actual environment if thats the environments name
        source ~/.bash_aliases
    ```
    where ``/absolute/path/to/curr/directory`` is the full path to the current directory (``/users/cnwizu/data/cnwizu/nclusion_manuscript_figure_reproducibility``)
13. Make sure you have created an alias that loads all of the dependencies on oscar (including R and its dependences)
    ```bash
        start_nclusion_manuscript
        R
    ```
14. Then in R
    ```R
        install.packages("renv")
        #Warning in install.packages("renv") :
        #  'lib = "/oscar/rt/9.2/software/0.20-generic/0.20.1/opt/spack/linux-rhel9-x86_64_v3/gcc-11.3.1/r-4.4.0-yycctsjvszuj5o2q4gfbaehsq7rkl4bz/rlib/R/library"' is not writable
        #Would you like to use a personal library instead? (yes/No/cancel) 'yes' [ENTER]
        # Would you like to create a personal library
        # ‘/oscar/home/cnwizu/R/x86_64-pc-linux-gnu-library/4.4’
        # to install packages into? (yes/No/cancel) 'yes' [ENTER]
        # --- Please select a CRAN mirror for use in this session ---
        # Secure CRAN mirrors 
        # Selection: '69' [ENTER]
        library(renv)
        renv::init(bare = TRUE)
        q()
    ```
15. Make sure and file called ".Rprofile" is created. In it (or if you have to make one) paste the following lines
    ```R
        # options(renv.download.override = utils::download.file)
        options(repos = c(CRAN = "https://cloud.r-project.org/", RForge = "https://r-forge.r-project.org", BioCsoft = "https://bioconductor.org/packages/3.12/bioc", BioCann = "https://bioconductor.org/packages/3.12/data/annotation", BioCexp = "https://bioconductor.org/packages/3.12/data/experiment", BioCworkflows = "https://bioconductor.org/packages/3.12/workflows",RCran="http://cran.us.r-project.org"))
        Sys.setenv(CXXSTD = "CXX11")
        Sys.setenv(CXXFLAGS = "-std=c++11")
        Sys.setenv(CXX14 = "g++")  # Instructs R to use g++ for C++14, but we downgrade to C++11
        Sys.setenv(CXX14FLAGS = "-std=c++11")  # Override the C++14 standard to C++11
        Sys.setenv(RENV_DOWNLOAD_FILE_METHOD = "libcurl")
        source("renv/activate.R")
    ```
16. In bash make sure you are in the current directory and that the ``/absolute/path/to/curr/directory/renv/`` directory exists. Then
    ```bash
        cd  renv/
        mkdir cellar
        cd cellar
        wget http://download.r-forge.r-project.org/src/contrib/lcmix_0.3.tar.gz
        wget https://github.com/XiDsLab/Festem/releases/download/v1.2.1/Festem_1.2.1.tar.gz
        wget https://github.com/XiDsLab/Festem_paper/raw/refs/heads/main/EMDE_V0.tar.gz
        mv EMDE_V0.tar.gz EMDE_0.0.0.tar.gz
        wget https://github.com/prabhakarlab/DUBStepR/archive/refs/tags/v1.1.3.tar.gz
        mv v1.1.3.tar.gz DUBStepR_1.1.3.tar.gz
        cd ../../ #make sure this gets you back to the /absolute/path/to/curr/directory/ directory
    ```
18. Restart the R session and type (Need make sure your github token is up to date so download packages to help with this)
    ```R
        library(renv)
        renv::install(c("gitcreds","usethis"))
        library(usethis)
        library(gitcreds)
        #0.
        usethis::create_github_token() # (OPTIONAL if not done before or token has expired) follow instruction in browser
        #0.5 copy and save this token somewhere (in my case in my Joplin Notebook)
        #1. When done copy token and type
        gitcreds::gitcreds_set()
        # Paste token "ghp_XXXXXXXXXXXXXXXXXXX..."
        # Verify status
        gitcreds::gitcreds_set()
        # -> Your current credentials for 'https://github.com':

        #   protocol: https
        #   host    : github.com
        #   username: PersonalAccessToken
        #   password: <-- hidden -->

        # -> What would you like to do? 

        # 1: Abort update with error, and keep the existing credentials
        # 2: Replace these credentials
        # 3: See the password / token

        # Selection: '1' [ENTER]
    ```
19. Finally load all of the needed R packages in R
    ```R
        library(renv)
        library(usethis)
        library(gitcreds)
        gitcreds::gitcreds_set()
        # Paste token "ghp_XXXXXXXXXXXXXXXXXXX..."
        # Verify status
        gitcreds::gitcreds_set()
        # -> Your current credentials for 'https://github.com':

        #   protocol: https
        #   host    : github.com
        #   username: PersonalAccessToken
        #   password: <-- hidden -->

        # -> What would you like to do? 

        # 1: Abort update with error, and keep the existing credentials
        # 2: Replace these credentials
        # 3: See the password / token

        # Selection: '1' [ENTER]
        Sys.setenv(GITHUB_PAT = gitcreds::gitcreds_get()$password)
        renv::install(c("ggplot2","dplyr","tidyverse","devtools", "BiocManager","reticulate","R.utils","stringr","Seurat","RColorBrewer","Matrix","scales","cowplot","RCurl","optparse","cluster","seriation","circlize","geneset","patchwork","xlsx", "ClueR","gitcreds","usethis","BH@1.72.0-3","devtools"))
        renv::init(bioconductor = TRUE)
        renv::install(c("bioc::singleCellHaystack","bioc::SingleCellExperiment","bioc::ComplexHeatmap","bioc::VennDetail","bioc::splatter","bioc::VariantAnnotation", "bioc::scater", "bioc::scDesign3", "bioc::org.Hs.eg.db","stats","bioc::DESeq2","ROSeq","scry","scran","TruncatedNormal","peakRAM","pracma","ggpubr","truncnorm","bioc::fgsea","bioc::genekitr","bioc::clusterProfiler","bioc::enrichplot","bioc::harmony","bioc::MAST","bioc::DuoClustering2018","bioc::zellkonverter","bioc::infercnv","bioc::SC3"))
        renv::install(c("mojaveazure/seurat-disk","VCCRI/CIDR","PYangLab/scCCESS","lingxuez/SOUP","lingxuez/SOUP","yulijia/SIMLR","fanyue322/TDEseq"))
        renv::install("bitbucket::scLCA/single_cell_lca")
        renv::install(c("biclust","clusterSim","fpc","ggalluvial","pROC","Rfast2","satijalab/seurat-data","survminer"))
        renv::install(c("Festem","DUBStepR","lcmix","EMDE","SIMLR"))#
        renv::install("immunogenomics/presto")
        renv::snapshot()
    ```
20. 
21. 
22. 
23. 
24. 
25. 
26. 
27. 
28. 
29. 
30. 