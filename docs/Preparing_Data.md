---
title: "Preparing Data"
author: "Jean Fan"
date: '2017-05-19'
output: html_document
md_document:
  variant: markdown_github
vignette: |
  %\VignetteIndexEntry{Vignette Title} %\VignetteEngine{knitr::rmarkdown} \usepackage[utf8]{inputenc}
---


Starting with aligned `.bam` files, there are a number of ways to prepare your data for `HoneyBADGER`. 

## Allele Data

For the allele-based model, you will need a set of heterozygous SNP positions. This can be ideally obtained from previous WES data from the same sample or estimated from the [ExAC database](http://exac.broadinstitute.org/).


```r
## Download ExAC vcf or use your own vcf
vcfFile <- "../data-raw/ExAC.r0.3.sites.vep.vcf.gz"
## For testing purposes, restrict to set of SNPs on region on chromosome 3
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
```

```
## vi
## FALSE  TRUE 
##  7960    76
```

```r
snps <- rowData(vcf)
snps <- snps[vi,]
## get rid of non single nucleotide changes
vi <- width(snps@elementMetadata$REF) == 1
snps <- snps[vi,]
## look at final SNPs
print(snps)
```

```
## GRanges object with 72 ranges and 5 metadata columns:
##                seqnames             ranges strand   | paramRangeID
##                   <Rle>          <IRanges>  <Rle>   |     <factor>
##   1:139213_A/G        1   [139213, 139213]      *   |         <NA>
##   1:139233_C/A        1   [139233, 139233]      *   |         <NA>
##   1:325075_G/C        1   [325075, 325075]      *   |         <NA>
##   1:325127_T/C        1   [325127, 325127]      *   |         <NA>
##   1:325155_C/A        1   [325155, 325155]      *   |         <NA>
##            ...      ...                ...    ... ...          ...
##      rs4333796        1 [1007432, 1007432]      *   |         <NA>
##     rs10907177        1 [1021346, 1021346]      *   |         <NA>
##      rs3737728        1 [1021415, 1021415]      *   |         <NA>
##      rs4074137        1 [1026707, 1026707]      *   |         <NA>
##      rs4562563        1 [1026801, 1026801]      *   |         <NA>
##                           REF                ALT       QUAL      FILTER
##                <DNAStringSet> <DNAStringSetList>  <numeric> <character>
##   1:139213_A/G              A                  G 1634387.96        PASS
##   1:139233_C/A              C                  A 1489200.18        PASS
##   1:325075_G/C              G                  C  343566.82        PASS
##   1:325127_T/C              T                  C  136558.82        PASS
##   1:325155_C/A              C                  A   82321.84        PASS
##            ...            ...                ...        ...         ...
##      rs4333796              G                  A    4289637        PASS
##     rs10907177              A                  G   26083658        PASS
##      rs3737728              A                  G   72706171        PASS
##      rs4074137              C                  A    3815651        PASS
##      rs4562563              T                  A    8608083        PASS
##   -------
##   seqinfo: 85 sequences from hg19 genome
```

Now we can get the number of reads corresponding to each SNP site for each cell using their `.bam` files. 


```r
library(HoneyBADGER)
path <- "../data-raw/"
files <- list.files(path = path)
bamFiles <- files[grepl('.bam$', files)]
bamFiles <- paste0(path, bamFiles)
indexFiles <- files[grepl('.bai$', files)]
indexFiles <- paste0(path, indexFiles)
results <- getSnpMats(snps, bamFiles, indexFiles)
print(names(results))
```

```
[1] "refCount" "altCount" "cov" 
```

Now we have a matrix of SNP coverage as well as reference and allele count for use in our `HoneyBADGER` allele model. 


```r
print(head(results$refCount))
```

```
                data-raw/MM16.MM16_15_CGTACTAG-TAGATCGC.recalibrated.bam
1:567783-567783                                                        0
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                        4
1:949654-949654                                                        0
                data-raw/MM16.MM16_16_CGTACTAG-CTCTCTAT.recalibrated.bam
1:567783-567783                                                        0
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                       51
1:949654-949654                                                        0
                data-raw/MM16.MM16_17_CGTACTAG-TATCCTCT.recalibrated.bam
1:567783-567783                                                        0
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                        0
1:949654-949654                                                        0
```

```r
print(head(results$altCount))
```

```
                data-raw/MM16.MM16_15_CGTACTAG-TAGATCGC.recalibrated.bam
1:567783-567783                                                       45
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                        4
1:949654-949654                                                        0
                data-raw/MM16.MM16_16_CGTACTAG-CTCTCTAT.recalibrated.bam
1:567783-567783                                                      107
1:888639-888639                                                        2
1:888659-888659                                                        2
1:948921-948921                                                        8
1:949608-949608                                                       51
1:949654-949654                                                        0
                data-raw/MM16.MM16_17_CGTACTAG-TATCCTCT.recalibrated.bam
1:567783-567783                                                        5
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                        0
1:949654-949654                                                        0
```

```r
print(head(results$cov))
```

```
                data-raw/MM16.MM16_15_CGTACTAG-TAGATCGC.recalibrated.bam
1:567783-567783                                                       45
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                        4
1:949654-949654                                                        3
                data-raw/MM16.MM16_16_CGTACTAG-CTCTCTAT.recalibrated.bam
1:567783-567783                                                      107
1:888639-888639                                                        2
1:888659-888659                                                        2
1:948921-948921                                                        8
1:949608-949608                                                       51
1:949654-949654                                                       38
                data-raw/MM16.MM16_17_CGTACTAG-TATCCTCT.recalibrated.bam
1:567783-567783                                                        5
1:888639-888639                                                        0
1:888659-888659                                                        0
1:948921-948921                                                        0
1:949608-949608                                                        0
1:949654-949654                                                        0
```

## Gene expression data

For gene expression data, we recommend quantification by TPM or FPM transformed to log scale. The same processing pipeline and transformation is highly recommended for the normal reference. Normal references can be ideally obtained from matched normal but can also be estimated using [GTeX](https://www.gtexportal.org/home/). 
