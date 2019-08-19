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
    b = as.character(pu[pu$which_label==names[i],]$nucleotide) == as.character(as.data.frame(gr$REF)$x[i])
    if(length(b)==0) { return(0) } # neither allele observed
    else if(sum(b)==0) { return(0) } # alt allele observed only
    else { return(pu[pu$which_label==names[i],]$count[b]) }
  }))
  altCount <- unlist(lapply(seq_along(names), function(i) {
    b = as.character(pu[pu$which_label==names[i],]$nucleotide) == as.character(as.data.frame(gr$ALT)$value[i])
    if(length(b)==0) { return(0) } # neither allele observed
    else if(sum(b)==0) { return(0) } # ref allele observed only
    else { return(pu[pu$which_label==names[i],]$count[b]) }
  }))
  names(refCount) <- names(altCount) <- names
  
  if (verbose) {
    print("Done!")
  }
  return(list(refCount, altCount))
}


#' Get alternative allele count for positions of interest for multiple cells from one bam
#'
#' @param gr GenomicRanges object for positions of interest
#' @param bamFile bam file containing all cells indexed by barcode (such as 10X)
#' @param cellBarcodes vector of valid cell barcodes  
#' @param tag BAM tag used to identify cells (e.g. 
#'            'CB' = Chromium cellular barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.
#'            'CR' = Chromium cellular barcode sequence as reported by the sequencer)
#' @param indexFile bai index file
#' @param verbose Boolean of whether or not to print progress and info
#' @param n.cores Number of cores; Can parallelize across cells
#' @return
#'   refCount reference allele count information for each position of interest
#'   altCount alternative allele count information for each position of interest
#'
#' @examples
#' \dontrun{
#' ## Get putative hets from ExAC
#' load(system.file("ExAC", "ExAC_chr1.RData", package = "HoneyBADGER"))
#' head(snps)
#' cov <- getCoverage(snps, bamFile, indexFile, verbose=TRUE)
#' ## for the sites with coverage, get all cells
#' vi <- cov>0
#' snps.cov <- snps[vi,]
#' barcodes <- c('AAACATACAAAACG-1', 'AAACATACAAAAGC-1', 'AAACATACAAACAG-1')
#' mats.chr1 <- getCellAlleleCount(snps.cov, bamFile, indexFile, barcodes, verbose=TRUE, n.cores=10)
#' }
#'
#' @export
#' 
getCellAlleleCount <- function (gr, bamFile, indexFile, cellBarcodes, tag = 'CB', verbose = FALSE, n.cores = 1) {
  df <- data.frame(seqnames(gr), ranges(gr))
  names <- paste(df$seqnames.gr., paste(df$start, df$end, sep = "-"), 
                 sep = ":")
  if (verbose) {
    print("Getting allele counts for...")
    print(names)
  }
  pp <- PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = TRUE, 
                    max_depth = 1e+07, min_base_quality = 20, min_mapq = 10)
  if (verbose) {
    print("Getting pileup...")
  }
  counts <- mclapply(cellBarcodes, function(cell) {
    tf <- list(cell); names(tf) <- tag
    sbp <- ScanBamParam(tagFilter = tf, which = gr)
    pu <- pileup(file = bamFile, index = indexFile, scanBamParam = sbp, 
                 pileupParam = pp)
    if (verbose) {
      print(paste0("Getting allele read counts for ", cell, 
                   "..."))
    }
    refCount <- unlist(lapply(seq_along(names), function(i) {
      b = as.character(pu[pu$which_label == names[i], ]$nucleotide) == 
        as.character(as.data.frame(gr$REF)$x[i])
      if (length(b) == 0) {
        return(0)
      }
      else if (sum(b) == 0) {
        return(0)
      }
      else {
        return(pu[pu$which_label == names[i], ]$count[b])
      }
    }))
    altCount <- unlist(lapply(seq_along(names), function(i) {
      b = as.character(pu[pu$which_label == names[i], ]$nucleotide) == 
        as.character(as.data.frame(gr$ALT)$value[i])
      if (length(b) == 0) {
        return(0)
      }
      else if (sum(b) == 0) {
        return(0)
      }
      else {
        return(pu[pu$which_label == names[i], ]$count[b])
      }
    }))
    names(refCount) <- names(altCount) <- names
    return(list(ref = refCount, alt = altCount))
  }, mc.cores = n.cores)
  refCount <- do.call(cbind, lapply(counts, function(x) x$ref))
  altCount <- do.call(cbind, lapply(counts, function(x) x$alt))
  colnames(refCount) <- colnames(altCount) <- cellBarcodes
  return(list(ref = refCount, alt = altCount))
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
#' snps <- rowData(vcf)
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
getSnpMats <- function(snps, bamFiles, indexFiles, n.cores=1, verbose=FALSE) {
  
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
  vi <- rowSums(cov)>0
  if(sum(vi)==0) {
    print('ERROR: NO SNPS WITH COVERAGE. Check if snps and bams are using the same alignment reference.')
    return(NULL)
  }
  cov <- cov[vi,]
  snps <- snps[vi,]
  
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
  
  results <- list(refCount=refCount, altCount=altCount, cov=cov)
  return(results)
}


#' Helper function to get coverage and allele count matrices given a set of putative heterozygous SNP positions specific for 10X data
#'
#' @param snps GenomicRanges object for positions of interest
#' @param bamFiles list of bam file
#' @param indexFiles list of bai index file
#' @param barcodes Vector of possible cell barcodes 
#' @param tag BAM tag used to identify cells (e.g. 
#'            'CB' = Chromium cellular barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.
#'            'CR' = Chromium cellular barcode sequence as reported by the sequencer)
#' @param n.cores number of cores
#' @param verbose Boolean of whether or not to print progress and info
#' 
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
#' snps <- rowData(vcf)
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
#' barcodes <- read.table("filtered_gene_bc_matrices/hg19/barcodes.tsv", header=FALSE)
#' results <- getSnpMats10X(snps, bamFiles, indexFiles, barcodes)
#' }
#'
#' @export
#' 
getSnpMats10X <- function (snps, bamFile, indexFile, barcodes, tag = 'CB', verbose = FALSE, n.cores = 1) {
  cov <- getCoverage(snps, bamFile, indexFile, verbose)
  if (verbose) {
    print("Snps with coverage:")
    print(table(cov > 0))
  }
  vi <- cov > 0
  if (sum(vi) == 0) {
    print("ERROR: NO SNPS WITH COVERAGE. Check if snps and bams are using the same alignment reference.")
    return(NULL)
  }
  cov <- cov[vi]
  snps.cov <- snps[vi, ]
  if (verbose) {
    print("Getting allele counts...")
  }
  alleleCount <- getCellAlleleCount(snps.cov, bamFile, indexFile, 
                                    barcodes, verbose = verbose, n.cores = n.cores, tag = tag)
  refCount <- alleleCount[[1]]
  altCount <- alleleCount[[2]]
  if (verbose) {
    print("altCount + refCount == cov:")
    print(table(altCount + refCount == cov))
    print("altCount + refCount < cov: sequencing errors")
    print(table(altCount + refCount < cov))
  }
  results <- list(refCount = refCount, altCount = altCount, 
                  cov = cov)
  return(results)
}
