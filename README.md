# Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION)

NCLUSION package documentation, examples, and tutorials can be found here: <a href="https://microsoft.github.io/nclusion"> NCLUSION Documentation and Tutorials </a>

## Introduction

<p>Profiling phenotypic heterogenity in cellular populations has recieved a
renewed interested due to the increased accuracy of single-cell
sequencing and the decrease in sequencing costs. As a result, new statistical
methods have been developed to characterize the phenotypic heterogenity derived
from single-cell sequencing experiments such as single-cell RNA-sequencing
(scRNA-Seq). Current methods often utilize clustering followed
by additional statistical tests post-clustering to select marker
genes that contribute to cluster identity. These methods are limited in that
the number of clusters in the data must be known prior to clustering and
post-clustering variable selection does not leverage information about the
cluster formation in order to elect marker genes.</p>
<p>We present a method for the Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION) - a sparse Bayesian Nonparametric model based on a Hierarchical Dirichlet Process mixture model. By combining a sparsity-promoting spike-and-slab prior on the cluster means and a nonparametric Hierarchical Dirichlet Process prior on the number of clusters, NCLUSION jointly learns the marker genes defining cluster identity and an optimal number of clusters needed to partition the data. Together with reasonable model approximations using varaiational inference, our proposed approach is scalable to large single-cell RNA sequencing experiments.</p>

## The Method

<p>The Nonparametric CLUstering of SIngle cell PopulatiONs (NCLUSION) is a
sparse Bayesian Nonparametric model based on a Hierarchical Dirichlet Process
mixture model, which aims to cluster cells from single-cell RNA sequencing
experiments and extract genes defining cluster
identity. It does this through the use of a
sparsity-promoting spike-and-slab prior on the cluster means and a nonparametric
Hierarchical Dirichlet Process prior on the number of clusters. The key idea
behind the concept NCLUSION is the treatment of cluster formation as a shift in
expression of particular genes away from a global mean. This mean shift is
modeled by the spike-and-slab prior. The number of clusters is modeled via a
stick-breaking proccess that allocates more clusters as the amount of variation
in the data grows. The learning of the model parameters proceeds via a
variational-EM algorithm that allows the method to scale with the number of
cells in the data.</p>

## Installation

This package requires that Julia 1.6+ and its dependencies are installed. Instructions for installation can be found <a href="https://github.com/JuliaLang/julia"> here. </a>

To create a new Julia project environment and install the NCLUSION package, follow these steps:

<ol><li> Start a Julia session by entering <code>julia</code> in a bash terminal while in
your working directory.</li> <li>  Enter a Julia Pkg REPL by pressing the right square bracket <code>]</code>.
  Create a new Julia project environment by typing <code>generate</code> followed by the name of your project (i.e. <code>generate nclusion_demo</code>). This creates a project
  directory named after the environment you created, and two files: <code>Project.toml</code>
  and <code>Manifest.toml</code>.
  These files, which can be found in the project
directory, specify the packages that are downloaded in the environment.</li>
<li>Activate the project environment in the Julia Pkg REPL by typing <code>activate</code>
followed by the name of the Julia environment (i.e. <code>activate
nclusion_demo</code>).</li> <li>Install NCLUSION in the project environment, by running the following command in the
  Julia Pkg REPL:<pre><code>add nclusion</code></pre></li></ol>

Alternatively, to install NCLUSION into an existing Julia project environment, follow these steps:

<ol><li>Activate the project environment and start a Julia REPL by entering <code>julia --project=/path/to/julia/env --thread=auto</code> into a bash terminal, where <code>--project</code> is set equal to the path to your existing Julia project environment.</li>
  <li>Install NCLUSION in the activated project environment by entering the
  following command in the Julia REPL: <pre><code>using Pkg;Pkg.add("nclusion")</code></pre></li></ol>

You can exit the Julia Pkg REPL by typing <code>Ctrl + C</code>, and the Julia REPL by entering <code>exit()</code>.

## Questions and Feedback

If you have any questions or concerns regarding NCLUSION or the tutorials, please contact <a href="mailto:chibuikem_nwizu@brown.edu"> Chibuikem Nwizu</a>, <a href="mailto:lcrawford@microsoft.com"> Lorin Crawford</a>, <a href="mailto:ava.amini@microsoft.com"> Ava Amini</a>, or <a href="mailto:v-mahughes@microsoft.com"> Madeline Hughes</a>. All feedback on the repository, manuscript, and tutorials is appreciated.

