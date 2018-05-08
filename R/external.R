#' Functions ported over from ChIPseeker
#' so that it can be removed as a dependency
#' since its installation is proving troublesome 
#' 
#' All functions originally by
#' @author G Yu

##' Annotate peaks
##'
##'
##' @title annotatePeak
##' @param peak peak file or GRanges object
##' @param tssRegion Region Range of TSS
##' @param TxDb TxDb object
##' @param level one of transcript and gene
##' @param assignGenomicAnnotation logical, assign peak genomic annotation or not
##' @param genomicAnnotationPriority genomic annotation priority
##' @param annoDb annotation package
##' @param addFlankGeneInfo logical, add flanking gene information from the peaks
##' @param flankDistance distance of flanking sequence
##' @param sameStrand logical, whether find nearest/overlap gene in the same strand
##' @param ignoreOverlap logical, whether ignore overlap of TSS with peak
##' @param ignoreUpstream logical, if True only annotate gene at the 3' of the peak.
##' @param ignoreDownstream logical, if True only annotate gene at the 5' of the peak.
##' @param overlap one of 'TSS' or 'all', if overlap="all", then gene overlap with peak will be reported as nearest gene, no matter the overlap is at TSS region or not.
##' @param verbose print message or not
##' @return data.frame or GRanges object with columns of:
##'
##' all columns provided by input.
##'
##' annotation: genomic feature of the peak, for instance if the peak is
##' located in 5'UTR, it will annotated by 5'UTR. Possible annotation is
##' Promoter-TSS, Exon, 5' UTR, 3' UTR, Intron, and Intergenic.
##'
##' geneChr: Chromosome of the nearest gene
##'
##' geneStart: gene start
##'
##' geneEnd: gene end
##'
##' geneLength: gene length
##'
##' geneStrand: gene strand
##'
##' geneId: entrezgene ID
##'
##' distanceToTSS: distance from peak to gene TSS
##'
##' if annoDb is provided, extra column will be included:
##'
##' ENSEMBL: ensembl ID of the nearest gene
##'
##' SYMBOL: gene symbol
##'
##' GENENAME: full gene name
##' @import BiocGenerics S4Vectors GenomeInfoDb
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
##' peakAnno
##' }
##' @author G Yu
annotatePeak <- function(peak,
                         tssRegion=c(-3000, 3000),
                         TxDb=NULL,
                         level = "transcript",
                         assignGenomicAnnotation=TRUE,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         annoDb=NULL,
                         addFlankGeneInfo=FALSE,
                         flankDistance=5000,
                         sameStrand = FALSE,
                         ignoreOverlap=FALSE,
                         ignoreUpstream=FALSE,
                         ignoreDownstream=FALSE,
                         overlap = "TSS",
                         verbose=TRUE) {
  
  is_GRanges_of_TxDb <- FALSE
  if (is(TxDb, "GRanges")) {
    is_GRanges_of_TxDb <- TRUE
    assignGenomicAnnotation <- FALSE
    annoDb <- NULL
    addFlankGeneInfo <- FALSE
    message("#\n#.. 'TxDb' is a self-defined 'GRanges' object...\n#")
    message("#.. Some parameters of 'annotatePeak' will be disable,")
    message("#.. including:")
    message("#..\tlevel, assignGenomicAnnotation, genomicAnnotationPriority,")
    message("#..\tannoDb, addFlankGeneInfo and flankDistance.")
    message("#\n#.. Some plotting functions are designed for visualizing genomic annotation")
    message("#.. and will not be available for the output object.\n#")
  }
  
  if (is_GRanges_of_TxDb) {
    level <- "USER_DEFINED"
  } else {
    level <- match.arg(level, c("transcript", "gene"))
  }
  
  if (assignGenomicAnnotation && all(genomicAnnotationPriority %in% c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")) == FALSE) {
    stop('genomicAnnotationPriority should be any order of c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")')
  }
  
  if ( is(peak, "GRanges") ){
    ## this test will be TRUE
    ## when peak is an instance of class/subclass of "GRanges"
    input <- "gr"
    peak.gr <- peak
  } else {
    input <- "file"
    peak.gr <- loadPeak(peak, verbose)
  }
  
  peakNum <- length(peak.gr)
  
  if (verbose)
    cat(">> preparing features information...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  
  if (is_GRanges_of_TxDb) {
    features <- TxDb
  } else {
    TxDb <- loadTxDb(TxDb)
    
    if (level=="transcript") {
      features <- getGene(TxDb, by="transcript")
    } else {
      features <- getGene(TxDb, by="gene")
    }
  }
  if (verbose)
    cat(">> identifying nearest features...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  
  ## nearest features
  idx.dist <- getNearestFeatureIndicesAndDistances(peak.gr, features,
                                                   sameStrand, ignoreOverlap,
                                                   ignoreUpstream,ignoreDownstream,
                                                   overlap=overlap)
  
  if (verbose)
    cat(">> calculating distance from peak to TSS...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  ## distance
  distance <- idx.dist$distance
  
  ## update peak, remove un-map peak if exists.
  peak.gr <- idx.dist$peak
  
  ## annotation
  if (assignGenomicAnnotation == TRUE) {
    if (verbose)
      cat(">> assigning genomic annotation...\t\t",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    
    anno <- getGenomicAnnotation(peak.gr, distance, tssRegion, TxDb, level, genomicAnnotationPriority, sameStrand=sameStrand)
    annotation <- anno[["annotation"]]
    detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
  } else {
    annotation <- NULL
    detailGenomicAnnotation <- NULL
  }
  
  ## append annotation to peak.gr
  if (!is.null(annotation))
    mcols(peak.gr)[["annotation"]] <- annotation
  
  
  has_nearest_idx <- which(idx.dist$index <= length(features))
  nearestFeatures <- features[idx.dist$index[has_nearest_idx]]
  
  ## duplicated names since more than 1 peak may annotated by only 1 gene
  names(nearestFeatures) <- NULL
  nearestFeatures.df <- as.data.frame(nearestFeatures)
  if (is_GRanges_of_TxDb) {
    colnames(nearestFeatures.df)[1:5] <- c("geneChr", "geneStart", "geneEnd",
                                           "geneLength", "geneStrand")
  } else if (level == "transcript") {
    colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                      "geneLength", "geneStrand", "geneId", "transcriptId")
    nearestFeatures.df$geneId <- TXID2EG(as.character(nearestFeatures.df$geneId), geneIdOnly=TRUE)
  } else {
    colnames(nearestFeatures.df) <- c("geneChr", "geneStart", "geneEnd",
                                      "geneLength", "geneStrand", "geneId")
  }
  
  for(cn in colnames(nearestFeatures.df)) {
    mcols(peak.gr)[[cn]][has_nearest_idx] <- unlist(nearestFeatures.df[, cn])
  }
  
  mcols(peak.gr)[["distanceToTSS"]] <- distance
  
  if (!is.null(annoDb)) {
    if (verbose)
      cat(">> adding gene annotation...\t\t\t",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    .idtype <- IDType(TxDb)
    if (length(.idtype) == 0 || is.na(.idtype) || is.null(.idtype)) {
      n <- length(peak.gr)
      if (n > 100)
        n <- 100
      sampleID <- peak.gr$geneId[1:n]
      
      if (all(grepl('^ENS', sampleID))) {
        .idtype <- "Ensembl Gene ID"
      } else if (all(grepl('^\\d+$', sampleID))) {
        .idtype <- "Entrez Gene ID"
      } else {
        warning("Unknown ID type, gene annotation will not be added...")
        .idtype <- NA
      }
    }
    
    if (!is.na(.idtype)) {
      peak.gr %<>% addGeneAnno(annoDb, .idtype)
    }
  }
  
  if (addFlankGeneInfo == TRUE) {
    if (verbose)
      cat(">> adding flank feature information from peaks...\t",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    
    flankInfo <- getAllFlankingGene(peak.gr, features, level, flankDistance)
    
    if (level == "transcript") {
      mcols(peak.gr)[["flank_txIds"]] <- NA
      mcols(peak.gr)[["flank_txIds"]][flankInfo$peakIdx] <- flankInfo$flank_txIds
    }
    
    mcols(peak.gr)[["flank_geneIds"]] <- NA
    mcols(peak.gr)[["flank_gene_distances"]] <- NA
    
    mcols(peak.gr)[["flank_geneIds"]][flankInfo$peakIdx] <- flankInfo$flank_geneIds
    mcols(peak.gr)[["flank_gene_distances"]][flankInfo$peakIdx] <- flankInfo$flank_gene_distances
    
  }
  
  if (!is_GRanges_of_TxDb) {
    if(verbose)
      cat(">> assigning chromosome lengths\t\t\t",
          format(Sys.time(), "%Y-%m-%d %X"), "\n")
    
    peak.gr@seqinfo <- seqinfo(TxDb)[names(seqlengths(peak.gr))]
  }
  
  if(verbose)
    cat(">> done...\t\t\t\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  
  if (assignGenomicAnnotation) {
    res <- new("csAnno",
               anno = peak.gr,
               tssRegion = tssRegion,
               level=level,
               hasGenomicAnnotation = TRUE,
               detailGenomicAnnotation=detailGenomicAnnotation,
               annoStat=getGenomicAnnoStat(peak.gr),
               peakNum=peakNum
    )
  } else {
    res <- new("csAnno",
               anno = peak.gr,
               tssRegion = tssRegion,
               level=level,
               hasGenomicAnnotation = FALSE,
               peakNum=peakNum
    )
  }
  return(res)
}


##' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
loadTxDb <- function(TxDb) {
  if ( is.null(TxDb) ) {
    warning(">> TxDb is not specified, use 'TxDb.Hsapiens.UCSC.hg19.knownGene' by default...")
    TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  return(TxDb)
}

##' @importFrom AnnotationDbi get
##' @importFrom GenomicFeatures genes
##' @importFrom GenomicFeatures transcriptsBy
getGene <- function(TxDb, by="gene") {
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  by <- match.arg(by, c("gene", "transcript"))
  
  if (by == "gene") {
    if ( exists("Genes", envir=ChIPseekerEnv, inherits=FALSE) ) {
      features <- get("Genes", envir=ChIPseekerEnv)
    } else {
      features <- genes(TxDb)
      assign("Genes", features, envir=ChIPseekerEnv)
    }
  } else {
    if ( exists("Transcripts", envir=ChIPseekerEnv, inherits=FALSE) ) {
      features <- get("Transcripts", envir=ChIPseekerEnv)
    } else {
      features <- transcriptsBy(TxDb)
      features <- unlist(features)
      assign("Transcripts", features, envir=ChIPseekerEnv)
    }
  }
  
  return(features)
}

##' get index of features that closest to peak and calculate distance
##'
##'
##' @title getNearestFeatureIndicesAndDistances
##' @param peaks peak in GRanges
##' @param features features in GRanges
##' @param sameStrand logical, whether find nearest gene in the same strand
##' @param ignoreOverlap logical, whether ignore overlap of TSS with peak
##' @param ignoreUpstream logical, if True only annotate gene at the 3' of the peak.
##' @param ignoreDownstream logical, if True only annotate gene at the 5' of the peak.
##' @param overlap one of "TSS" or "all"
##' @return list
##' @import BiocGenerics IRanges GenomicRanges
##' @author G Yu
getNearestFeatureIndicesAndDistances <- function(peaks, features,
                                                 sameStrand = FALSE,
                                                 ignoreOverlap=FALSE,
                                                 ignoreUpstream=FALSE,
                                                 ignoreDownstream=FALSE,
                                                 overlap = "TSS") {
  
  overlap <- match.arg(overlap, c("TSS", "all"))
  
  if (!ignoreOverlap && overlap == "all") {
    overlap_hit <- findOverlaps(peaks, unstrand(features))
  }
  
  ## peaks only conatin all peak records, in GRanges object
  ## feature is the annotation in GRanges object
  
  ## only keep start position based on strand
  ## start(features) <- end(features) <- ifelse(strand(features) == "+", start(features), end(features))
  features <- resize(features, width=1) # faster
  
  ## add dummy NA feature for peaks that are at the last or first feature
  ## suggested by Michael Kluge
  features.bak <- features
  seqlevels(features) <- c(seqlevels(features), "chrNA")
  dummy <- GRanges("chrNA", IRanges(1,1))
  
  ## dummy$tx_id <- -1
  ## dummy$tx_name <- "NA"
  
  cns <- names(mcols(features))
  for (cn in cns) {
    if (grepl('id', cn)) {
      mcols(dummy)[[cn]] <- -1
    } else {
      mcols(dummy)[[cn]] <- NA
    }
  }
  
  features <- append(features, dummy)
  dummyID <- length(features)
  
  if (sameStrand) {
    ## nearest from peak start
    ps.idx <- follow(peaks, features)
    
    ## nearest from peak end
    pe.idx <- precede(peaks, features)
  } else {
    ps.idx <- follow(peaks, unstrand(features))
    pe.idx <- precede(peaks, unstrand(features))
  }
  
  na.idx <- is.na(ps.idx) & is.na(pe.idx)
  if (sum(na.idx) > 0) { ## suggested by Thomas Schwarzl
    ps.idx <- ps.idx[!na.idx]
    pe.idx <- pe.idx[!na.idx]
    ##peaks <- peaks[!na.idx]
  }
  
  # set NA values to dummy value if only one entry is affected
  ps.idx[is.na(ps.idx)] <- dummyID
  pe.idx[is.na(pe.idx)] <- dummyID
  
  ## features from nearest peak start
  psF <- features[ps.idx]
  
  ## feature distances from peak start
  psD <- ifelse(strand(psF) == "+", 1, -1) *
    (start(peaks[!na.idx]) - start(psF))
  psD[ps.idx == dummyID] <- Inf # ensure that there is even no match if a seq with name "chrNA" exists
  
  ## features from nearest peak end
  peF <- features[pe.idx]
  ## feature distances from peak end
  peD <- ifelse(strand(peF) == "+", 1, -1) *
    (end(peaks[!na.idx]) - start(peF))
  peD[pe.idx == dummyID] <- Inf # ensure that there is even no match if a seq with name "chrNA" exists
  
  ## restore the old feature object
  features <- features.bak
  
  pse <- data.frame(ps=psD, pe=peD)
  if (ignoreUpstream) {
    j <- rep(2, nrow(pse))
  } else if (ignoreDownstream) {
    j <- rep(1, nrow(pse))
  } else {
    j <- apply(pse, 1, function(i) which.min(abs(i)))
  }
  
  ## index
  idx <- ps.idx
  idx[j==2] <- pe.idx[j==2]
  
  ## distance
  dd <- psD
  dd[j==2] <- peD[j==2]
  
  index <- distanceToTSS <- rep(NA, length(peaks))
  distanceToTSS[!na.idx] <- dd
  index[!na.idx] <- idx
  
  if (!ignoreOverlap) {
    ## hit <- findOverlaps(peaks, unstrand(features))
    
    if (overlap == "all") {
      hit <- overlap_hit
      if ( length(hit) != 0 ) {
        qh <- queryHits(hit)
        hit.idx <- getFirstHitIndex(qh)
        hit <- hit[hit.idx]
        peakIdx <- queryHits(hit)
        featureIdx <- subjectHits(hit)
        
        index[peakIdx] <- featureIdx
        distance_both_end <- data.frame(start=start(peaks[peakIdx]) - start(features[featureIdx]),
                                        end = end(peaks[peakIdx]) - start(features[featureIdx]))
        distance_idx <- apply(distance_both_end, 1, function(i) which.min(abs(i)))
        distance_minimal <- distance_both_end[,1]
        distance_minimal[distance_idx == 2] <- distance_both_end[distance_idx==2, 2]
        
        distanceToTSS[peakIdx] <- distance_minimal * ifelse(strand(features[featureIdx]) == "+", 1, -1)
        
      }
    }
    
    hit <- findOverlaps(peaks, unstrand(features))
    
    if ( length(hit) != 0 ) {
      qh <- queryHits(hit)
      hit.idx <- getFirstHitIndex(qh)
      hit <- hit[hit.idx]
      peakIdx <- queryHits(hit)
      featureIdx <- subjectHits(hit)
      
      index[peakIdx] <- featureIdx
      distanceToTSS[peakIdx] <- 0
    }
    
  }
  
  j <- is.na(distanceToTSS) | is.na(index)
  
  res <- list(index=index[!j],
              distance=distanceToTSS[!j],
              peak=peaks[!j])
  
  return(res)
}

isPeakFeatureOverlap <- function(peak, feature) {
  peakRange <- ranges(peak)
  featureRange <- ranges(feature)
  x <- intersect(peakRange, featureRange)
  return(length(x) != 0)
}

##' get Genomic Annotation of peaks
##'
##'
##' @title getGenomicAnnotation
##' @param peaks peaks in GRanges object
##' @param distance distance of peak to TSS
##' @param tssRegion tssRegion, default is -3kb to +3kb
##' @param TxDb TxDb object
##' @param level one of gene or transcript
##' @param genomicAnnotationPriority genomic Annotation Priority
##' @param sameStrand whether annotate gene in same strand
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @return character vector
##' @author G Yu
getGenomicAnnotation <- function(peaks,
                                 distance,
                                 tssRegion=c(-3000, 3000),
                                 TxDb,
                                 level,
                                 genomicAnnotationPriority,
                                 sameStrand = FALSE
) {
  
  ##
  ## since some annotation overlap,
  ## a priority is assign based on *genomicAnnotationPriority*
  ## use the following priority by default:
  ##
  ## 1. Promoter
  ## 2. 5' UTR
  ## 3. 3' UTR
  ## 4. Exon
  ## 5. Intron
  ## 6. Downstream
  ## 7. Intergenic
  ##
  
  .ChIPseekerEnv(TxDb)
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  annotation <- rep(NA, length(distance))
  
  flag <- rep(FALSE, length(distance))
  detailGenomicAnnotation <- data.frame(
    genic=flag,
    Intergenic=flag,
    Promoter=flag,
    fiveUTR=flag,
    threeUTR=flag,
    Exon=flag,
    Intron=flag,
    downstream=flag,
    distal_intergenic=flag)
  
  anno <- list(annotation=annotation,
               detailGenomicAnnotation=detailGenomicAnnotation)
  
  genomicAnnotationPriority <- rev(genomicAnnotationPriority)
  for (AP in genomicAnnotationPriority) {
    if (AP == "Intergenic") {
      ## Intergenic
      annotation[is.na(annotation)] <- "Intergenic"
      anno[["annotation"]] <- annotation
    } else if (AP == "Intron") {
      ## Introns
      intronList <- get_intronList(ChIPseekerEnv)
      anno <- updateGenomicAnnotation(peaks, intronList, "Intron", anno, sameStrand=sameStrand)
    } else if (AP == "Exon") {
      ## Exons
      exonList <- get_exonList(ChIPseekerEnv)
      anno <- updateGenomicAnnotation(peaks, exonList, "Exon", anno, sameStrand=sameStrand)
    } else if (AP == "3UTR") {
      ## 3' UTR Exons
      if ( exists("threeUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        threeUTRList <- get("threeUTRList", envir=ChIPseekerEnv)
      } else {
        threeUTRList <- threeUTRsByTranscript(TxDb)
        assign("threeUTRList", threeUTRList, envir=ChIPseekerEnv)
      }
      anno <- updateGenomicAnnotation(peaks, threeUTRList, "threeUTR", anno, sameStrand=sameStrand)
    } else if (AP == "5UTR") {
      ## 5' UTR Exons
      if ( exists("fiveUTRList", envir=ChIPseekerEnv, inherits=FALSE) ) {
        fiveUTRList <- get("fiveUTRList", envir=ChIPseekerEnv)
      } else {
        fiveUTRList <- fiveUTRsByTranscript(TxDb)
        assign("fiveUTRList", fiveUTRList, envir=ChIPseekerEnv)
      }
      anno <- updateGenomicAnnotation(peaks, fiveUTRList, "fiveUTR", anno, sameStrand=sameStrand)
    }
    
    annotation <- anno[["annotation"]]
    detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
    
    if (AP == "Promoter") {
      ## TSS
      tssIndex <- distance >= tssRegion[1] & distance <= tssRegion[2]
      annotation[tssIndex] <- "Promoter"
      detailGenomicAnnotation[tssIndex, "Promoter"] <- TRUE
      
      pm <- max(abs(tssRegion))
      if (pm/1000 >= 2) {
        dd <- seq(1:ceiling(pm/1000))*1000
        for (i in 1:length(dd)) {
          if (i == 1) {
            lbs <- paste("Promoter", " (<=", dd[i]/1000, "kb)", sep="")
            annotation[abs(distance) <= dd[i] &
                         annotation == "Promoter"] <- lbs
          } else {
            lbs <- paste("Promoter", " (", dd[i-1]/1000, "-", dd[i]/1000, "kb)", sep="")
            annotation[abs(distance) <= dd[i] &
                         abs(distance) > dd[i-1] &
                         annotation == "Promoter"] <- lbs
          }
        }
      }
    }
    
  }
  
  
  genicIndex <- which(apply(detailGenomicAnnotation[, c("Exon", "Intron")], 1, any))
  detailGenomicAnnotation[-genicIndex, "Intergenic"] <- TRUE
  detailGenomicAnnotation[genicIndex, "genic"] <- TRUE
  
  ## intergenicIndex <- anno[["annotation"]] == "Intergenic"
  ## anno[["detailGenomicAnnotation"]][intergenicIndex, "Intergenic"] <- TRUE
  ## anno[["detailGenomicAnnotation"]][!intergenicIndex, "genic"] <- TRUE
  
  
  features <- getGene(TxDb, by=level)
  
  ## nearest from gene end
  if (sameStrand) {
    idx <- precede(peaks, features)
  } else {
    idx <- precede(peaks, unstrand(features))
  }
  na.idx <- which(is.na(idx))
  if (length(na.idx)) {
    idx <- idx[-na.idx]
    peaks <- peaks[-na.idx]
  }
  peF <- features[idx]
  dd <- ifelse(strand(peF) == "+",
               start(peaks) - end(peF),
               end(peaks) - start(peF))
  if (length(na.idx)) {
    dd2 <- numeric(length(idx) + length(na.idx))
    dd2[-na.idx] <- dd
  } else {
    dd2 <- dd
  }
  
  for (i in 1:3) { ## downstream within 3k
    j <- which(annotation == "Intergenic" & abs(dd2) <= i*1000 & dd2 != 0)
    if (length(j) > 0) {
      if (i == 1) {
        lbs <- "Downstream (<1kb)"
      } else {
        lbs <- paste("Downstream (", i-1, "-", i, "kb)", sep="")
      }
      annotation[j] <- lbs
    }
  }
  annotation[which(annotation == "Intergenic")] = "Distal Intergenic"
  
  dsd <- getOption("ChIPseeker.downstreamDistance")
  if (is.null(dsd))
    dsd <- 3000 ## downstream 3k by default
  
  downstreamIndex <- dd2 > 0 & dd2 < dsd
  detailGenomicAnnotation[downstreamIndex, "downstream"] <- TRUE
  detailGenomicAnnotation[which(annotation == "Distal Intergenic"), "distal_intergenic"] <- TRUE
  return(list(annotation=annotation, detailGenomicAnnotation=detailGenomicAnnotation))
}

##' @importFrom GenomicFeatures intronsByTranscript
get_intronList <- function(ChIPseekerEnv) {
  TxDb <- get("TXDB", envir=ChIPseekerEnv)
  if ( exists("intronList", envir=ChIPseekerEnv, inherits=FALSE) ) {
    intronList <- get("intronList", envir=ChIPseekerEnv)
  } else {
    intronList <- intronsByTranscript(TxDb)
    assign("intronList", intronList, envir=ChIPseekerEnv)
  }
  return(intronList)
}

updateGenomicAnnotation <- function(peaks, genomicRegion, type, anno, sameStrand=FALSE) {
  hits <- getGenomicAnnotation.internal(peaks, genomicRegion, type, sameStrand=sameStrand)
  if (length(hits) > 1) {
    hitIndex <- hits$queryIndex
    anno[["annotation"]][hitIndex] <- hits$annotation
    anno[["detailGenomicAnnotation"]][hitIndex, type] <- TRUE
  }
  return(anno)
}

##' @import BiocGenerics S4Vectors IRanges
getGenomicAnnotation.internal <- function(peaks, genomicRegion, type, sameStrand=FALSE){
  GRegion <- unlist(genomicRegion)
  GRegionLen <- elementNROWS(genomicRegion)
  
  names(GRegionLen) <- names(genomicRegion)
  GRegion$gene_id <- rep(names(genomicRegion), times=GRegionLen)
  
  
  if (type == "Intron") {
    gr2 <- GRegion[!duplicated(GRegion$gene_id)]
    strd <- as.character(strand(gr2))
    len <- GRegionLen[GRegionLen != 0]
    
    GRegion$intron_rank <- unlist(lapply(seq_along(strd), function(i) {
      rank <- seq(1, len[i])
      if (strd[i] == '-')
        rank <- rev(rank)
      return(rank)
    }))
  }
  
  if (type == "Intron" || type =="Exon") {
    nn <- TXID2EG(names(genomicRegion))
    names(GRegionLen) <- nn
    GRegion$gene_id <- rep(nn, times=GRegionLen)
  }
  
  ## find overlap
  if (sameStrand) {
    GRegionHit <- findOverlaps(peaks, GRegion)
  } else {
    GRegionHit <- findOverlaps(peaks, unstrand(GRegion))
  }
  
  if (length(GRegionHit) == 0) {
    return(NA)
  }
  qh <- queryHits(GRegionHit)
  hit.idx <- getFirstHitIndex(qh)
  GRegionHit <- GRegionHit[hit.idx]
  queryIndex <- queryHits(GRegionHit)
  subjectIndex <- subjectHits(GRegionHit)
  
  hits <- GRegion[subjectIndex]
  geneID <- hits$gene_id
  
  if (type == "Intron") {
    anno <- paste(type, " (", geneID, ", intron ", hits$intron_rank,
                  " of ", GRegionLen[geneID], ")", sep="")
  } else if (type == "Exon") {
    anno <- paste(type, " (", geneID, ", exon ", hits$exon_rank,
                  " of ", GRegionLen[geneID], ")", sep="")
  } else if (type == "fiveUTR") {
    anno <- "5' UTR"
  } else if (type == "threeUTR") {
    anno <- "3' UTR"
  } else {
    anno <- type
  }
  res <- list(queryIndex=queryIndex, annotation=anno, gene=geneID)
  return(res)
}


TXID2EG <- function(txid, geneIdOnly=FALSE) {
  txid <- as.character(txid)
  if (geneIdOnly == TRUE) {
    res <- TXID2EGID(txid)
  } else {
    res <- TXID2TXEG(txid)
  }
  return(res)
}

##' @importFrom GenomicFeatures transcripts
TXID2TXEG <- function(txid) {
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if (exists("txid2geneid", envir=ChIPseekerEnv, inherits=FALSE)) {
    txid2geneid <- get("txid2geneid", envir=ChIPseekerEnv)
  } else {
    txdb <- get("TXDB", envir=ChIPseekerEnv)
    txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    idx <- which(sapply(txidinfo$gene_id, length) == 0)
    txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
    txid2geneid <- paste(mcols(txidinfo)[["tx_name"]],
                         mcols(txidinfo)[["gene_id"]],
                         sep="/")
    txid2geneid <- sub("/NA", "", txid2geneid)
    
    names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
    assign("txid2geneid", txid2geneid, envir=ChIPseekerEnv)
  }
  return(as.character(txid2geneid[txid]))
}

TXID2EGID <- function(txid) {
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  
  if (exists("txid2eg", envir=ChIPseekerEnv, inherits=FALSE)) {
    txid2geneid <- get("txid2eg", envir=ChIPseekerEnv)
  } else {
    txdb <- get("TXDB", envir=ChIPseekerEnv)
    txidinfo <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    idx <- which(sapply(txidinfo$gene_id, length) == 0)
    txidinfo[idx,]$gene_id <- txidinfo[idx,]$tx_name
    txid2geneid <- as.character(mcols(txidinfo)[["gene_id"]])
    
    names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
    assign("txid2eg", txid2geneid, envir=ChIPseekerEnv)
  }
  return(as.character(txid2geneid[txid]))
}

## according to: https://support.bioconductor.org/p/70432/#70545
## contributed by Hervé Pagès
getFirstHitIndex <- function(x) {
  ## sapply(unique(x), function(i) which(x == i)[1])
  which(!duplicated(x))
}

##' @importFrom GenomicFeatures exonsBy
get_exonList <- function(ChIPseekerEnv) {
  TxDb <- get("TXDB", envir=ChIPseekerEnv)
  if ( exists("exonList", envir=ChIPseekerEnv, inherits=FALSE) ) {
    exonList <- get("exonList", envir=ChIPseekerEnv)
  } else {
    exonList <- exonsBy(TxDb)
    assign("exonList", exonList, envir=ChIPseekerEnv)
  }
  return(exonList)
}

getGenomicAnnoStat <- function(peakAnno) {
  if ( class(peakAnno) == "GRanges" )
    peakAnno <- as.data.frame(peakAnno)
  anno <- peakAnno$annotation
  ## anno <- sub(" \\(.+", "", anno)
  
  e1 <- getOption("ChIPseeker.ignore_1st_exon")
  i1 <- getOption("ChIPseeker.ignore_1st_intron")
  ids <- getOption("ChIPseeker.ignore_downstream")
  
  if (is.null(e1) || !e1) {
    e1lab <- "1st Exon"
    anno[grep("exon 1 of", anno)] <- e1lab
    exonlab <- "Other Exon"
  } else {
    e1lab <- NULL
    exonlab <- "Exon"
  }
  
  if (is.null(i1) || !i1) {
    i1lab <- "1st Intron"
    anno[grep("intron 1 of", anno)] <- i1lab
    intronlab <- "Other Intron"
  } else {
    i1lab <- NULL
    intronlab <- "Intron"
  }
  
  anno[grep("Exon \\(", anno)] <- exonlab
  anno[grep("Intron \\(", anno)] <- intronlab
  
  if (is.null(ids) || !ids) {
    dsd <- getOption("ChIPseeker.downstreamDistance")
    if (is.null(dsd))
      dsd <- 3000 ## downstream 3k by default
    if (dsd > 1000) {
      dsd <- round(dsd/1000, 1)
      dsd <- paste0(dsd, "kb")
    }
    dslab <- paste0("Downstream (<=", dsd, ")")
    
    anno[grep("Downstream", anno)] <- dslab
    iglab <- "Distal Intergenic"
  } else {
    dslab <- NULL
    iglab <- "Intergenic"
    anno[grep("Downstream", anno)] <- iglab
  }
  anno[grep("^Distal", anno)] <- iglab
  
  lvs <- c(
    "5' UTR",
    "3' UTR",
    e1lab,
    exonlab,
    i1lab,
    intronlab,
    dslab,
    iglab
  )
  
  promoter <- unique(anno[grep("Promoter", anno)])
  ip <- getOption("ChIPseeker.ignore_promoter_subcategory")
  if ((is.null(ip) || !ip) && (length(promoter) > 0)) {
    plab <- sort(as.character(promoter))
  } else {
    plab <- "Promoter"
    anno[grep("^Promoter", anno)] <- plab
  }
  lvs <- c(plab, lvs)
  
  ## count frequency
  anno.table <- table(anno)
  
  ## calculate ratio
  anno.ratio <- anno.table/ sum(anno.table) * 100
  anno.df <- as.data.frame(anno.ratio)
  colnames(anno.df) <- c("Feature", "Frequency")
  
  anno.df$Feature <- factor(anno.df$Feature, levels=lvs[lvs %in% anno.df$Feature])
  anno.df <- anno.df[order(anno.df$Feature),]
  return(anno.df)
}

##' @importFrom AnnotationDbi get
.ChIPseekerEnv <- function(TxDb) {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists("ChIPseekerEnv", envir=.GlobalEnv)) {
    assign("ChIPseekerEnv", new.env(), envir = envir)
  }
  
  ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
  if (!exists("TXDB", envir=ChIPseekerEnv, inherits=FALSE)) {
    ## first run
    assign("TXDB", TxDb, envir=ChIPseekerEnv)
  } else {
    TXDB <- get("TXDB", envir=ChIPseekerEnv)
    m1 <- tryCatch(unlist(metadata(TXDB)), error=function(e) NULL)
    
    m2 <- unlist(metadata(TxDb))
    
    if (!is.null(m1)) {
      m1 <- m1[!is.na(m1)]
    }
    m2 <- m2[!is.na(m2)]
    
    if ( is.null(m1) || length(m1) != length(m2) || any(m1 != m2) ) {
      rm(ChIPseekerEnv)
      assign("ChIPseekerEnv", new.env(), envir = envir)
      ChIPseekerEnv <- get("ChIPseekerEnv", envir=.GlobalEnv)
      assign("TXDB", TxDb, envir=ChIPseekerEnv)
    }
  }
  
}



##' Class "csAnno"
##' This class represents the output of ChIPseeker Annotation
##'
##'
##' @name csAnno-class
##' @aliases csAnno-class
##' show,csAnno-method vennpie,csAnno-method
##' plotDistToTSS,csAnno-method plotAnnoBar,csAnno-method
##' plotAnnoPie,csAnno-method upsetplot,csAnno-method
##'
##' @docType class
##' @slot anno annotation
##' @slot tssRegion TSS region
##' @slot level transcript or gene
##' @slot hasGenomicAnnotation logical
##' @slot detailGenomicAnnotation Genomic Annotation in detail
##' @slot annoStat annotation statistics
##' @slot peakNum number of peaks
##' @exportClass csAnno
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{annotatePeak}}
##' @keywords classes
setClass("csAnno",
         representation=representation(
           anno = "GRanges",
           tssRegion = "numeric",
           level = "character",
           hasGenomicAnnotation = "logical",
           detailGenomicAnnotation="data.frame",
           annoStat="data.frame",
           peakNum="numeric"
         ))


##' convert csAnno object to GRanges
##'
##'
##' @title as.GRanges
##' @param x csAnno object
##' @return GRanges object
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @export
as.GRanges <- function(x) {
  if (!is(x, "csAnno"))
    stop("not supported...")
  return(x@anno)
}


getAnnoStat <- function(x) {
  if (!is(x, "csAnno"))
    stop("not supported...")
  return(x@annoStat)
}

##' convert csAnno object to data.frame
##'
##'
##' @title as.data.frame.csAnno
##' @param x csAnno object
##' @param row.names row names
##' @param optional should be omitted.
##' @param ... additional parameters
##' @return data.frame
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @method as.data.frame csAnno
##' @export
as.data.frame.csAnno <- function(x, row.names=NULL, optional=FALSE, ...) {
  y <- as.GRanges(x)
  if (!(is.null(row.names) || is.character(row.names)))
    stop("'row.names' must be NULL or a character vector")
  df <- as.data.frame(y)
  rownames(df) <- row.names
  return(df)
}
