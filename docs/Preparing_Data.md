---
title: "Preparing Data"
author: "Jean Fan"
date: '2018-01-01'
output:
  pdf_document: default
md_document:
  variant: markdown_github
vignette: |
  %\VignetteIndexEntry{Vignette Title} %\VignetteEngine{knitr::rmarkdown} \usepackage[utf8]{inputenc}
---





## Allele Data

For the allele-based model, you will need a set of heterozygous SNP positions. This can be ideally obtained from previous WES data from the same sample or estimated from the [ExAC database](http://exac.broadinstitute.org/).

Example using your own VCF file:


```r
## Use your own vcf file with heterozygous variants
vcfFile <- "hets.vcf.gz"
## For testing purposes, restrict to set of SNPs on region on chromosome 1
require(GenomicRanges)
testRanges <- GRanges('1', IRanges(start=1e5, width=1e6))
require(VariantAnnotation)
param <- ScanVcfParam(which=testRanges)
vcf <- readVcf(vcfFile, "hg19", param=param)
## limit to common snps by MAF
info <- info(vcf)
maf <- info[, 'AF'] # AF is Integer allele frequency for each Alt allele
maft <- 0.1
vi <- sapply(maf, function(x) any(x > maft))
print(table(vi))
snps <- rowData(vcf)
snps <- snps[vi,]
## get rid of non single nucleotide changes
vi <- width(snps@elementMetadata$REF) == 1
snps <- snps[vi,]
## look at final SNPs
print(snps)
```

Example using ExAC:

```r
## available for all autosomes (Chr1 to Chr22)
load(system.file("ExAC", "ExAC_chr1.RData", package = "HoneyBADGER"))
print(head(snps))
```

```
## GRanges object with 6 ranges and 5 metadata columns:
##               seqnames         ranges strand | paramRangeID            REF
##                  <Rle>      <IRanges>  <Rle> |     <factor> <DNAStringSet>
##   1:17365_C/G        1 [17365, 17365]      * |         <NA>              C
##   1:17385_G/A        1 [17385, 17385]      * |         <NA>              G
##   1:69270_A/G        1 [69270, 69270]      * |         <NA>              A
##    rs75062661        1 [69511, 69511]      * |         <NA>              A
##   1:69761_A/T        1 [69761, 69761]      * |         <NA>              A
##   1:69897_T/C        1 [69897, 69897]      * |         <NA>              T
##                              ALT        QUAL                 FILTER
##               <DNAStringSetList>   <numeric>            <character>
##   1:17365_C/G                  G    826621.1 InbreedingCoeff_Filter
##   1:17385_G/A                  A    592354.5 InbreedingCoeff_Filter
##   1:69270_A/G                  G   1758695.3                   PASS
##    rs75062661                  G 120729371.2                   PASS
##   1:69761_A/T                  T   2041395.6                   PASS
##   1:69897_T/C                  C   1171733.0                   PASS
##   -------
##   seqinfo: 85 sequences from hg19 genome
```

Now we can get the number of reads corresponding to each SNP site for each cell using their `.bam` files. Here, we have placed all `.bam` and corresponding `.bai` index files in the `data-raw/` folder. There is one `.bam` and `.bai` for each cell. 


```r
library(HoneyBADGER)
path <- "data-raw/"
files <- list.files(path = path)
bamFiles <- files[grepl('.bam$', files)]
bamFiles <- paste0(path, bamFiles) ## list of paths to bam files
indexFiles <- files[grepl('.bai$', files)] 
indexFiles <- paste0(path, indexFiles) ## list of paths to index files
results <- getSnpMats(snps, bamFiles, indexFiles)
```

Now we have a matrix of SNP coverage as well as reference and allele count for use in our `HoneyBADGER` allele model. 


```r
ref <- results$refCount
alt <- results$altCount
cov <- results$cov
```

## Gene expression data

For gene expression data, we recommend quantification by counts transformed to log CPM. The same processing pipeline and transformation is highly recommended for the normal reference. Normal references can be ideally obtained from matched normal but can also be estimated using [GTeX](https://www.gtexportal.org/home/). 

## Accomodating 10X Data

For 10X data, you can use the output of `CellRanger`. For example, the `Gene / cell matrix (filtered)` can be normalized to CPMs and log transformmed to serve as the gene expression matrix. For the allele matrix, `Genome-aligned BAM` and `Genome-aligned BAM index` will be used as `bamFile` and `indexFile` respectively. However, as all cells will be contained in the same bam, we will use a different function to get the allele counts for each cell `getCellAlleleCount`. The column names of the expression matrix will be your cell barcodes `cellBarcodes`.  


```r
results <- getCellAlleleCount(snps, bamFile, indexFile, cellBarcodes)
```



