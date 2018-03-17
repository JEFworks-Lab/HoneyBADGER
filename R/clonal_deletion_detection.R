#' @author Catherine Xue
#' Functions for assessing for clonal deletions


#' Helper function to get total allele counts and coverage across all cells at each 
#' putative heterozygous SNP, given the allele counts and coverage from each cell
#'
#' @param snpMats List of data frames with allele counts at putative heterozygous
#'   SNPs for each cell
#' @return
#'   snpSummary data frame of allele counts and coverage at putative heterozygous
#'     SNPs aggregated across all cells in sample
#' 
#' @examples
#' \dontrun{
#' seqlevels(snps) <- paste0('chr', seqlevels(snps)) #convert to UCSC format
#' snpMat <- getSnpMats(snps, bams_pt, bais_pt)
#' snpSummary <- summarizeMatrix(snpMat)
#' head(snpSummary)        
#' }
#'
#' @export
#' 
summarizeMatrix <- function(snpMats) {
  snpSummary <- data.frame(rowSums(snpMats$ref), rowSums(snpMats$alt), rowSums(snpMats$cov))
  colnames(snpSummary) <- c('ref', 'alt', 'cov')
  filter <- apply(snpSummary, 1, function(x) {x[1] + x[2] == x[3]})
  snpSummary <- snpSummary[filter,]

  return(snpSummary)
}

#' Helper function to determine heterzogosity at a position from allele counts 
#' Returns TRUE for heterozygous SNPs, FALSE otherwise based on coverage
#' If reads are recorded from only one of reference or alternate allele, then it is
#' automakically classified as homozygous
#' If reads from both alleles are recorded, then we adopt a null hypothesis that
#' the number of minor allele reads sequenced is distributed  according to
#' Binomial(p = 0.5, n), where n is the total read depth at that position. SNPs at
#' which the null hypothesis p = 0.5 was rejected were classified as homozygous,
#' while the remainder were classified as heterozygous. 
#'
#' @param cov vector of allele counts and coverage for a putative heterozygous SNP
#' @param t threshold for hypothesis testing for heterozygosity, default 1e-8
#'   (Illumina sequencing error rate)
#' @return
#'   Boolean - TRUE for heterozygous, FALSE for homozygous
#' 
#' @examples
#' \dontrun{
#' }
#'
#' @export
#' 
isHet <- function (cov, t=1e-8) {
  x <- min(as.numeric(cov[1]), as.numeric(cov[2]))
  if (x == 0)  return(FALSE)
  return(pbinom(x, as.numeric(cov[3]), 0.5) > t)
}

#' Alternate helper function to determine heterzogosity at a position from allele  
#' counts. Returns TRUE for heterozygous SNPs, FALSE otherwise based on coverage.
#' If only reads are recorded from reference allele, then SNP is automatically
#' classified as homozygous.
#' If reads from both alleles or just the alternate allele are recorded, then the
#' hypothesis test is applied
#'
#' @param cov vector of allele counts and coverage for a putative heterozygous SNP
#' @param t threshold for hypothesis testing for heterozygosity, default 1e-8
#'   (Illumina sequencing error rate)
#' @return
#'   Boolean - TRUE for heterozygous, FALSE for homozygous
#' 
#' @examples
#' \dontrun{
#' }
#'
#' @export
#' 
isHetAlt <- function (cov, t=1e-8) {
  a <- as.numeric(cov[1])
  b <- as.numeric(cov[2])
  c <- as.numeric(cov[3])

  if (b == 0)  return(FALSE)
  x <- min(a, b)
  return(pbinom(x, c , 0.5) > t)
}

#' Helper function to get het rate from allele counts at a list of common variant
#' SNPs
#'
#' @param snpSummary data frame of common variant SNPs and their allele counts and
#'   coverage
#' @param t threshold for hypothesis testing for heterozygosity, default 1e-8
#'   (Illumina sequencing error rate)
#' @param alt Boolean of whether or not to use alternate het filter
#' @return
#'   hetRate percentage of heterozygous SNPs in the input data frame
#' 
#' @examples
#' \dontrun{
#' }
#'
#' @export
#' 
getHetRate <- function (snpSummary, t=1e-8, alt=FALSE) {
  if (alt == TRUE)  hetFilter <- apply(snpSummary, 1, isHetAlt)
  else  hetFilter <- apply(snpSummary, 1, isHet)
  hetRate <- mean(hets)
  return(hetRate)
}

#' Function to filter out SNPs from ExAC database
#'
#' @param vcf readVcf object from ExAC database read over a certain range
#' @param maft minor allele frequency cutoff for filtering
#' @return
#'   snps GenomicRanges object with positions of putative heterozygous SNPs
#'   info data fram with metadata for putative heterozygous SNPs
#' 
#' @examples
#' \dontrun{
#' vcfFile <- "../data-raw/ExAC.r0.3.sites.vep.vcf.gz"
#' # Filter for SNPs over all of chr1	
#' testRanges <- GRanges(seqnames='1', IRanges(start = 0, width=249250621))
#' param = ScanVcfParam(which=testRanges)
#' vcf <- readVcf(vcfFile, "hg19", param=param)
#' temp <- trimSnps(vcf)
#' snps <- temp$snps
#' info <- temp$info
#' head(info)
#' head(snps)
#' }
#'
#' @export
#'
trimSnps <- function (vcf, maft=0.1) {
    snps <- rowRanges(vcf)
    info <- info(vcf)
    if (length(snps) == 0) {return(list('snps'=snps, 'info'=info))}
    maf <- info[, 'AF']
    vi <- sapply(maf, function(x) any(x > maft))
    info <- info[vi,]
    snps <- snps[vi,]
    vi <- width(snps@elementMetadata$REF) == 1
    info <- info[vi,]
    snps <- snps[vi,]
    vi <- width(snps@elementMetadata$ALT@partitioning) == 1
    info <- info[vi,]
    snps <- snps[vi,]

    mat <- info[, c('Het_AFR', 'Hom_AFR', 'Het_AMR', 'Hom_AFR', 'Het_EAS', 'Hom_AMR', 'Het_FIN', 'Hom_EAS', 'Het_NFE', 'Hom_NFE', 'Het_OTH', 'Hom_OTH', 'Het_SAS', 'Hom_SAS')]
    ## convert to regular integer matrix
    mat <- do.call(cbind, lapply(mat, unlist))
    vi <- rowSums(mat) > 0
    info <- info[vi,]
    snps <- snps[vi,]
    return(list('snps'=snps, 'info'=info))
}
