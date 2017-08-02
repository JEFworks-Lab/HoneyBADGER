[![Build Status](http://travis-ci.org/JEFworks/HoneyBADGER.svg?branch=master)](https://travis-ci.org/JEFworks/HoneyBADGER)

[![HoneyBADGER](http://jef.works/HoneyBADGER/assets/img/logo.png)](http://jef.works/HoneyBADGER)

# Overview of HoneyBADGER

`HoneyBADGER` (**h**idden Markov model integrated **B**ayesian **a**pproach for **d**etecting CNVs and LOHs from sin**g**le-c**e**ll **R**NA-seq data) identifies and infers the presence of CNV and LOH events in single cells and reconstruct subclonal architecture using scRNA-seq data. 

The overall approach is detailed in the following publication: COMING SOON!

---

## Sample analysis and images

### Iterative HMM approach for CNV detection
![](http://jef.works/HoneyBADGER/assets/img/approach.png)

### Bayesian model to infer CNVs in single cells using allele and expression data
![](http://jef.works/HoneyBADGER/assets/img/allele.png)  
![](http://jef.works/HoneyBADGER/assets/img/expression.png)

---

## Installation

To install `HoneyBADGER`, we recommend using `devtools`:

```
require(devtools)
devtools::install_github('JEFworks/HoneyBADGER')
```

---

## Tutorials
- [Preparing data](docs/Preparing_Data.md)
- [Getting Started](docs/Getting_Started.md)
- Integration with `scde` and `liger`

---

## Help

For questions, suggests, and other comments, please submit through <a href="http://github.com/JEFworks/HoneyBADGER/issues">Github Issues</a>.
