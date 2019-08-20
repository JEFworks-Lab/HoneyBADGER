# [![HoneyBADGER](http://jef.works/HoneyBADGER/assets/img/logo.png)](http://jef.works/HoneyBADGER)

[![Build Status](http://travis-ci.org/JEFworks/HoneyBADGER.svg?branch=master)](https://travis-ci.org/JEFworks/HoneyBADGER)

`HoneyBADGER` (**h**idden Markov model integrated **B**ayesian **a**pproach for **d**etecting CNV and LOH events from sin**g**le-c**e**ll **R**NA-seq data) identifies and infers the presence of CNV and LOH events in single cells and reconstructs subclonal architecture using allele and expression information from single-cell RNA-sequencing data. 

The overall approach is detailed in the following publication:  
[Fan J*, Lee HO*, Lee S, et al. Linking transcriptional and genetic tumor heterogeneity through allele analysis of single-cell RNA-seq data. Genome Res. 2018;](https://genome.cshlp.org/content/early/2018/06/13/gr.228080.117)

---

## Benefits and Capabilities

### (1) Iterative HMM approach detects CNVs
![](http://jef.works/HoneyBADGER/assets/img/approach.png)

### (2) Bayesian hierarchical model uses allele and expression data to infer probability of CNVs in single cells
![](http://jef.works/HoneyBADGER/assets/img/example.png)

### (3) CNV inference from transcriptional data enables transcriptional characterization of subclones and other integrative analyses
![](http://jef.works/HoneyBADGER/assets/img/integration.png)


## Please see the `HoneyBADGER` [website](http://jef.works/HoneyBADGER/) for more tutorials and information.

---

## Installation

To install `HoneyBADGER`, we recommend using `devtools`:

```
require(devtools)
devtools::install_github('JEFworks/HoneyBADGER')
```

`HoneyBADGER` uses [`JAGS` (Just Another Gibbs Sampler)](http://mcmc-jags.sourceforge.net/) through `rjags`. Therefore, `JAGS` must be installed per your operating system requirements. Please see this [R-bloggers tutorial](https://www.r-bloggers.com/getting-started-with-jags-rjags-and-bayesian-modelling/) for additional tips for installing `JAGS` and `rjags`.

Additional dependencies may need to be installed from [`Bioconductor`](https://www.bioconductor.org/install/) such as `GenomicRanges` and others:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```

---

## Contributing

We welcome any bug reports, enhancement requests, and other contributions. To submit a bug report or enhancement request, please use the `HoneyBADGER` <a href="http://github.com/JEFworks/HoneyBADGER/issues">GitHub issues tracker</a>. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.

