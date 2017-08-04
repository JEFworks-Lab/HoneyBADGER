## Helper functions to prepare allele data for HoneyBADGER

#' Get alternative allele count for positions of interest
#'
#' @param gr GenomicRanges object for positions of interest
#' @param bamFile bam file
#' @param indexFile bai index file
#' @param verbose Boolean of whether or not to print progress and info
#' @return
#'   refCount reference allele count information for each position of interest
#'   altCount alternative allele count information for each position of interest
#'
#' @examples
#' \dontrun{
#' # Sites of interest (chr1:4600000, chr2:2000)
#' gr <- GRanges(c('chr1', 'chr2'), IRanges(start=c(4600000, 2000), width=1))
#' # we can get the coverage at these SNP sites from our bams
#' path <- '../data-raw/bams/'
#' files <- list.files(path = path)
#' files <- files[grepl('.bam$', files)]
#' alleleCounts <- lapply(files, function(f) {
#'    bamFile <- paste0(path, f)
#'    indexFile <- paste0(path, paste0(f, '.bai'))
#'    getAlleleCount(gr, bamFile, indexFile)
#' })
#' altCounts <- do.call(cbind, lapply(1:length(gr), function(i) alleleCounts[[i]][[1]]))
#' refCounts <- do.call(cbind, lapply(1:length(gr), function(i) alleleCounts[[i]][[2]]))
#' colnames(altCounts) <- colnames(refCounts) <- files
#' }
#'
#' @export
#'
getAlleleCount <- function (gr, bamFile, indexFile, verbose = FALSE) {
  df <- data.frame(seqnames(gr), ranges(gr))
  names <- paste(df$seqnames.gr., paste(df$start, df$end, sep='-'), sep=':')
  if (verbose) {
    print("Getting allele counts for...")
    print(names)
  }
  
  pp <- PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = TRUE, max_depth = 1e+07, min_base_quality = 20, min_mapq = 10)
  if (verbose) {
    print("Getting pileup...")
  }
  pu <- pileup(file = bamFile, index = indexFile, scanBamParam = ScanBamParam(which = gr), pileupParam = pp)
  
  if (verbose) {
    print("Getting allele read counts...")
  }
  refCount <- unlist(lapply(seq_along(names), function(i) {
    b = as.character(pu[pu$which_label==names[i],]$nucleotide) == as.character(data.frame(gr$REF)$value[i])
    if(length(b)==0) { return(0) } # neither allele observed
    else if(sum(b)==0) { return(0) } # alt allele observed only
    else { return(pu[pu$which_label==names[i],]$count[b]) }
  }))
  altCount <- unlist(lapply(seq_along(names), function(i) {
    b = as.character(pu[pu$which_label==names[i],]$nucleotide) == as.character(data.frame(gr$ALT)$value[i])
    if(length(b)==0) { return(0) } # neither allele observed
    else if(sum(b)==0) { return(0) } # ref allele observed only
    else { return(pu[pu$which_label==names[i],]$count[b]) }
  }))
  names(refCount) <- names(altCount) <- names
  
  if (verbose) {
    print("Done!")
  }
  return(list('ref'=refCount, 'alt'=altCount))
}


#' Get coverage count for positions of interest
#'
#' @param gr GenomicRanges object for positions of interest 
#' @param bamFile bam file
#' @param indexFile bai index file
#' @param verbose Boolean of whether or not to print progress and info
#' @return totCount Total coverage count information for each position of interest
#'
#' @examples
#' \dontrun{
#' # Sites of interest (chr1:4600000, chr2:2000)
#' gr <- GRanges(c('chr1', 'chr2'), IRanges(start=c(4600000, 2000), width=1))
#' # we can get the coverage at these SNP sites from our bams
#' path <- '../data-raw/bams/'
#' files <- list.files(path = path)
#' files <- files[grepl('.bam$', files)]
#' cov <- do.call(cbind, lapply(files, function(f) {
#'     bamFile <- paste0(path, f)
#'     indexFile <- paste0(path, paste0(f, '.bai'))
#'     getCoverage(gr, bamFile, indexFile)
#' }))
#' colnames(cov) <- files
#' }
#'
#' @export
#'
getCoverage <- function (gr, bamFile, indexFile, verbose = FALSE) {
  df <- data.frame(seqnames(gr), ranges(gr))
  names <- paste(df$seqnames.gr., paste(df$start, df$end, sep='-'), sep=':')
  if (verbose) {
    print("Getting coverage for...")
    print(names)
  }
  
  pp <- PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = FALSE, max_depth = 1e+07, min_base_quality = 20, min_mapq = 10)
  if (verbose) {
    print("Getting pileup...")
  }
  pu <- pileup(file = bamFile, index = indexFile, scanBamParam = ScanBamParam(which = gr), pileupParam = pp)
  rownames(pu) <- pu$which_label
  if (verbose) {
    print("Getting coverage counts...")
  }
  totCount <- pu[names, ]$count
  totCount[is.na(totCount)] <- 0
  names(totCount) <- names
  if (verbose) {
    print("Done!")
  }
  return(totCount)
}



#' Helper function to get coverage and allele count matrices given a set of putative heterozygous SNP positions
#'
#' @param snps GenomicRanges object for positions of interest
#' @param bamFiles list of bam file
#' @param indexFiles list of bai index file
#' @param n.cores number of cores
#' @param verbose Boolean of whether or not to print progress and info
#' @return
#'   refCount reference allele count matrix for each cell and each position of interest
#'   altCount alternative allele count matrix for each cell and each position of interest
#'   cov total coverage count matrix for each cell and each position of interest
#' 
#' @examples
#' \dontrun{
#' # Get putative hets from ExAC
#' vcfFile <- "../data-raw/ExAC.r0.3.sites.vep.vcf.gz"
#' testRanges <- GRanges(chr, IRanges(start = 1, width=1000))
#' param = ScanVcfParam(which=testRanges)
#' vcf <- readVcf(vcfFile, "hg19", param=param)
#' ## common snps by MAF
#' info <- info(vcf)
#' if(nrow(info)==0) {
#'     if(verbose) {
#'         print("ERROR no row in vcf")
#'     }
#'     return(NA)
#' }
#' maf <- info[, 'AF'] # AF is Integer allele frequency for each Alt allele
#' if(verbose) {
#'     print(paste0("Filtering to snps with maf > ", maft))
#' }
#' vi <- sapply(maf, function(x) any(x > maft))
#' if(verbose) {
#'     print(table(vi))
#' }
#' snps <- rowRanges(vcf)
#' snps <- snps[vi,]
#' ## get rid of non single nucleotide changes
#' vi <- width(snps@elementMetadata$REF) == 1
#' snps <- snps[vi,]
#' ## also gets rid of sites with multiple alt alleles though...hard to know which is in our patient
#' vi <- width(snps@elementMetadata$ALT@partitioning) == 1
#' snps <- snps[vi,]
#' ## Get bams
#' files <- list.files(path = "../data-raw")
#' bamFiles <- files[grepl('.bam$', files)]
#' bamFiles <- paste0(path, bamFiles)
#' indexFiles <- files[grepl('.bai$', files)]
#' indexFiles <- paste0(path, indexFiles)
#' results <- getSnpMats(snps, bamFiles, indexFiles)
#' }
#'
#' @export
#' 
getSnpMats <- function(snps, bamFiles, indexFiles, n.cores=10, verbose=FALSE) {
  if (is.list(snps))  snps <- snps[[1]]
  if (length(snps) == 0)  return(list('ref'=NA, 'alt'=NA, 'cov'=NA))

  ## loop
  cov <- do.call(cbind, mclapply(seq_along(bamFiles), function(i) {
    bamFile <- bamFiles[i]
    indexFile <- indexFiles[i]
    getCoverage(snps, bamFile, indexFile, verbose)
  }, mc.cores=n.cores))
  colnames(cov) <- bamFiles
  
  ## any coverage?
  if(verbose) {
    print("Snps with coverage:")
    print(table(rowSums(cov)>0))
  }
  vi <- rowSums(cov)>0; table(vi)
  cov <- cov[vi,]
  snps <- snps[vi,]
  
  if (length(snps) == 0) {
    results <- list('ref'=NA, 'alt'=NA, 'cov'=NA)
    return(results)
  }

  if(verbose) {
    print("Getting allele counts...")
  }
  alleleCount <- mclapply(seq_along(bamFiles), function(i) {
    bamFile <- bamFiles[i]
    indexFile <- indexFiles[i]
    getAlleleCount(snps, bamFile, indexFile, verbose)
  }, mc.cores=n.cores)
  refCount <- do.call(cbind, lapply(alleleCount, function(x) x[[1]]))
  altCount <- do.call(cbind, lapply(alleleCount, function(x) x[[2]]))
  colnames(refCount) <- colnames(altCount) <- bamFiles
  
  ## check correspondence
  if(verbose) {
    print("altCount + refCount == cov:")
    print(table(altCount + refCount == cov))
    print("altCount + refCount < cov: sequencing errors")
    print(table(altCount + refCount < cov))
    ##vi <- which(altCount + refCount != cov, arr.ind=TRUE)
    ## some sequencing errors evident
    ##altCount[vi]
    ##refCount[vi]
    ##cov[vi]
  }
  
  results <- list('ref'=refCount, 'alt'=altCount, 'cov'=cov)
  return(results)
}

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
