# Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION)

NCLUSION package documentation, examples, and tutorials can be found here: <a href="https://microsoft.github.io/nclusion"> NCLUSION Documentation and Tutorials </a>

## Introduction

Profiling phenotypic heterogenity in cellular populations has recieved a renewed interested recently due to the increased accuracy of single cell sequencing and the decrease in sequencing costs. As a result, new statistical methods have been developed to characterize the phenotypic heterogenity that result from single cell sequencing experiments such as single cell RNA-sequencing (scRNA-Seq) experiments. One popular way of carrying out this task is through clustering followed by performing additional statistical tests post-clustering to select salient genes that contribute to cluster identity. Current methods are limited in that the number of clusters in the data must be known prior to clustering and post-clustering variable selection does not leverage information about the cluster formation in order to elect salient genes.

We present a method for the Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION) - a sparse Bayesian Nonparametric model based on a Hierarchical Dirichlet Process mixture model. By combining a sparsity-promoting spike-and-slab prior on the cluster means and a nonparametric Hierarchical Dirichlet Process prior on the number of clusters, NCLUSION jointly learns the salient genes defining cluster identity and an optimal number of clusters needed to partition the data. Together with reasonable model approximations using varaiational inference, our proposed approach is scalable to large single cell RNA sequencing experiments.

## The Method

The Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION) is a sparse Bayesian Nonparametric model based on a Hierarchical Dirichlet Process mixture model which aims to identify single cell RNA sequencing clusters and the genes defining cluster identity. It does through the use of a sparsity-promoting spike-and-slab prior on the cluster means and a nonparametric Hierarchical Dirichlet Process prior on the number of clusters. The key idea behind the concept NCLUSION is the treatment of cluster formation as a shift in expression of particular genes away from a global mean. This mean shift is modeled by the spike-and-slab prior. The number of clusters is modeled via a stick-breaking proccess that allocates more clusters as the amount of variation in the data grows. The learning of the model parameters proceeds via a variational-EM algorithm that allows the method to scale with the number of cells in the data. As an overview of NCLUSION and its corresponding software implementation, we will assume that we have access to an expression matrix from a scRNA-seq study on N cells and G genes denoted as X where X is an N x G matrix of expression counts with G denoting the number of genes.

The goal of NCLUSION is to identify the phenotypic clusters in the data while jointly identifying the genes that define these clusters. To accomplish this, we asses the quality of clustering using both labels and label-free metrics. We also perform gene-pathway enrichment to asses the quality of elected genes.

## Installation

This package requires that Julia 1.6+ and its dependencies are installed. Instructions for installation can be found <a href="https://github.com/JuliaLang/julia"> here. </a>

To create a new Julia project environment and install the NCLUSION package, follow these steps:

<ol><li> Start a Julia session by entering <code>julia</code> in a bash terminal while in
your working directory.</li> <li>  Enter a Julia Pkg repl by pressing the right square bracket <code>]</code>.
  Create a new Julia project environment by typing <code>generate</code> followed by the name of your project (i.e. <code>generate nclusion_demo</code>). This creates a project
  directory named after the environment you created, and two files: <code>Project.toml</code>
  and <code>Manifest.toml</code>.
  These files, which can be found in the project
directory, specify the packages that are downloaded in the environment.</li>
<li>Activate the project environment in the Julia Pkg repl by typing <code>activate</code>
followed by the name of the Julia environment (i.e. <code>activate
nclusion_demo</code>).</li> <li>Install NCLUSION in the project environment, by running the following command in the
  Julia Pkg repl:<pre><code>using Pkg;Pkg.add("nclusion")</code></pre></li></ol>

Alternatively, to install NCLUSION into an existing Julia project environment, follow these steps:

<ol><li>Activate the project environment and start a Julia repl by entering <code>julia --project=/path/to/julia/env --thread=auto</code> into a bash terminal, where <code>--project</code> is set equal to the path to your existiing Julia project environment.</li>
  <li>Install NCLUSION in the activated project environment by entering the
  following command in the Julia repl: <pre><code>using Pkg;Pkg.add("nclusion")</code></pre></li></ol>

You can exit the Julia Pkg repl by typing <code>Ctrl + C</code>, and the Julia repl by entering <code>exit()</code>.

## Questions and Feedback

If you have any questions or concerns regarding NCLUSION or the tutorials, please contact <a href="mailto:chibuikem_nwizu@brown.edu"> Chibuikem Nwizu</a>, <a href="mailto:lcrawford@microsoft.com"> Lorin Crawford</a>, <a href="mailto:ava.amini@microsoft.com"> Ava Amini</a>, or <a href="mailto:v-mahughes@microsoft.com"> Madeline Hughes</a>. All feedback on the repository, manuscript, and tutorials is appreciated.

