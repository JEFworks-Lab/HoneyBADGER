`HoneyBADGER` enables detection of megabase-level copy number alterations such as deletions, amplifications, and copy-neutral loss-of-heterozygosity in single cells from single-cell RNA-seq data. `HoneyBADGER` relies on allele-information and gene-expression information derived from single-cell RNA-seq data. In this tutorial, we walk you through how to prepare your dataset for `HoneyBADGER`.

``` r
library(HoneyBADGER)
```

Preparing allele data
=====================

`HoneyBADGER` provides an HMM-integrated Bayesian hierarchical allele-based model to identify and infer the probability of copy number alterations in each single cell on the basis of persistent allelic imbalance.

To run the allele-based model, you will need matrices of heterozygous SNP counts. Specifically, we will need the counts of the reference allele, alternate allele, and total coverage for each SNP in each cell. Functions within the `HoneyBADGER` package `getSnpMats`, `getAlleleCount`, and `getCellAlleleCount` will help you create these matrices given a list of indexed bams where each cell corresponds to one bam (common for single-cell datasets generated from plate-based approaches), or a single bam with multiple cell barcodes (common for droplet-based single-cell datasets).

Heterozygous SNP positions will also need to be provided. These will ideally obtained from previous WES data from the same sample. When WES data from the same sample is not available, common heterozygous SNPs can be derived from databses such as [ExAC database](http://exac.broadinstitute.org/).

In this example, we will create a list of heterozygous SNPs as GRanges from a VCF file:

``` r
# Use your own vcf file with heterozygous variants
vcfFile <- "hets.vcf.gz"
# For testing purposes, restrict to set of SNPs on region on chromosome 1
require(GenomicRanges)
testRanges <- GRanges('1', IRanges(start=1e5, width=1e6))
require(VariantAnnotation)
param <- ScanVcfParam(which=testRanges)
# Be sure to use the correct genome species/version
vcf <- readVcf(vcfFile, "hg19", param=param)

snps <- rowRanges(vcf) 
# AF is the allele frequency for each alternate allele
info <- info(vcf)
maf <- info[, 'AF'] 
# limit to common snps by MAF (ie. > 10% in population)
maft <- 0.1
vi <- sapply(maf, function(x) any(x > maft))
snps <- snps[vi,]
# get rid of non single nucleotide changes
vi <- width(snps@elementMetadata$REF) == 1
# result should be GRanges object
snps <- snps[vi,]
```

This process has already been done for common heterozygous SNPs from ExAC (hg19) and can be loaded directly from `HoneyBADGER`:

``` r
# available for all autosomes (Chr1 to Chr22) for hg19 only
load(system.file("ExAC", "ExAC_chr1.RData", package = "HoneyBADGER"))
print(head(snps))
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

Now, given this list of potential heterozygous SNPs, we can get the number of reads corresponding to each SNP site for each cell using their `.bam` files. Here, we have placed all `.bam` and corresponding `.bai` index files in the `data-raw/` folder. There is one `.bam` and `.bai` for each cell.

``` r
library(HoneyBADGER)

path <- "data-raw/"
files <- list.files(path = path)
# list of paths to bam files
bamFiles <- files[grepl('.bam$', files)]
bamFiles <- paste0(path, bamFiles) 
# list of paths to index files
indexFiles <- files[grepl('.bai$', files)] 
indexFiles <- paste0(path, indexFiles) 

results <- getSnpMats(snps, bamFiles, indexFiles)
```

`getSnpMats` creates a matrix of SNP coverage as well as reference and allele count for use in our `HoneyBADGER` allele model.

``` r
r <- results$refCount
#cov <- results$cov
cov <- results$refCount + results$altCount ## alternatively
```

Preparing gene expression data
==============================

`HoneyBADGER` provides an HMM-integrated Bayesian hierarchical expression-based model to identify and infer the probability of copy number alterations in each single cell on the basis of persistent deviations in gene expression from a normal expression reference. Normal references can be ideally obtained from matched normal cells or sorted samples from the same patient but can also be estimated using [GTeX](https://www.gtexportal.org/home/).

To run the expression-based model, we recommend quantification by counts transformed to log CPM. The same processing pipeline and transformation is highly recommended for the normal reference.

Accomodating 10X Data
=====================

For 10X data, you can use the output of `CellRanger`. For example, the `Gene / cell matrix (filtered)` can be normalized to CPMs and log transformmed to serve as the gene expression matrix. For the allele matrix, `Genome-aligned BAM` and `Genome-aligned BAM index` will be used as `bamFile` and `indexFile` respectively. However, as all cells will be contained in the same bam, we will use a different function to get the allele counts for each cell `getCellAlleleCount`. The column names of the expression matrix will be your cell barcodes `cellBarcodes`.

``` r
results <- getSnpMats10X(snps, bamFile, indexFile, cellBarcodes)
r <- results$refCount
cov <- results$refCount + results$altCount
```

An alternative and much faster way of obtaining these allele-specific count tables is with the [scAlleleCount package](https://github.com/barkasn/scAlleleCount). Once installed you can issue the following command to obtain the tables:

``` r
results <- getFastCellAlleleCount(snps, bamFile, cellBarcodes)
```

Additional trouble shooting
=====================

If you are using the set of ExAC SNPs that comes with the `HoneyBADGER` package, you may need to remove certain chromosomes from the `seqlevels` depending on your hg19 build.

``` r
library(HoneyBADGER)
load(system.file("ExAC", "ExAC_chr1.RData", package = "HoneyBADGER"))
## ignore alt contigs (restrict to autosomes here)
vi <- seqlevels(snps) %in% as.character(c(1:22))
table(vi)
seqlevels(snps) <- seqlevels(snps)[vi]
```

The chromosome names used in your `snps` parameter needs to match the chromosome names used in your aligned `bams`. In this example, 10X used non-canonical chromosome names (other than UCSC RefSeq or Ensemble annotations) in order to distinguish human from mouse so we need to modify the `seqlevels` to match.

``` r
# need to match snp name with seqlevel in bam
seqlevels(snps) <- paste0('hg19_', seqlevels(snps))
```

To obtain a set of valid cell barcodes from 10X, you can use 10X Cell Ranger's whitelisted set of barcodes. However, most of these will not be present in your bam (because none of the cells sequenced have the barcode). So alternatively, you can restrict to the set of filtered barcodes with corresponding expression counts. 

``` r
## whitelist
cellBarcodes <- readLines('cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt')
## from filtered expression matrix
cellBarcodes <- readLines(gzfile('filtered_feature_bc_matrix/barcodes.tsv.gz')) ## barcodes with expression counts (filtered)
head(cellBarcodes)
```

---

# Additional tutorials
- [Getting started](Getting_Started.md)
- [Integrating with other analyses](Integrating.md)
- [Interactive visualization](https://jef.works/blog/2018/04/15/interactive-honeybadger-laf-profiles/)

