#' Functions related to allele model


#' Set allele count matrices, creates in-silico pooled single cells as bulk reference if none provided
#' 
#' @name setAlleleMats
#' @param r.init SNP site alternate allele count matrix for single cells
#' @param n.sc.init SNP site coverage count matrix for single cells
#' @param l.init SNP site alternate allele counts for bulk reference. If NULL, in silico bulk will be created from single cells.
#' @param n.bulk.init SNP site coverage counts for bulk reference. If NULL, in silico bulk will be created from single cells.
#' @param het.deviance.threshold Deviation from expected 0.5 heterozygous fraction 
#' @param min.cell Minimum number of cells a SNP must have coverage observed in
#' @param n.cores Number of cores
#' @param verbose Verbosity
#' 
#' @examples 
#' data(r)
#' data(cov.sc)
#' allele.mats <- setAlleleMats(r, cov.sc)
#'
setAlleleMats=function(r.init, n.sc.init, l.init=NULL, n.bulk.init=NULL, filter=TRUE, het.deviance.threshold=0.05, min.cell=3, n.cores=1, verbose=TRUE) {
        if(verbose) {  
            cat("Initializing allele matrices ... \n")
        }
        
        if(is.null(l.init) | is.null(n.bulk.init)) {
            if(verbose) { 
                cat("Creating in-silico bulk ... \n")
                cat(paste0("using ", ncol(r.init), " cells ... \n"))
            }
            l <<- rowSums(r.init>0)
            n.bulk <<- rowSums(n.sc.init>0)
        } else {
            l <- l.init
            n.bulk <- n.bulk.init
        }
        
        if(filter) {
            if(verbose) { 
                cat("Filtering for putative heterozygous snps ... \n")
                cat(paste0("allowing for a ", het.deviance.threshold, " deviation from the expected 0.5 heterozygous allele fraction ... \n"))
            }
            E <- l/n.bulk
            vi <- names(which(E > het.deviance.threshold & E < 1-het.deviance.threshold))
            ##cat(paste0(length(vi), " heterozygous SNPs identified \n"))
            if(length(vi) < 0.01*length(l)) {
                cat("WARNING! CLONAL DELETION OR LOH POSSIBLE! \n")
            }
            r <- r.init[vi,]
            n.sc <- n.sc.init[vi,]
            l <- l[vi]
            n.bulk <- n.bulk[vi]
            
            ## must have coverage in at least 5 cells
            if(verbose) { 
                cat(paste0("must have coverage in at least ", min.cell, " cells ... \n"))
            }
            vi <- rowSums(n.sc > 0) >= min.cell
            cat(paste0(length(vi), " heterozygous SNPs identified \n"))
            r <- r[vi,]
            n.sc <- n.sc[vi,]
            l <- l[vi]
            n.bulk <- n.bulk[vi]
            
        } else {
            r <- r.init
            n.sc <- n.sc.init
        }
        if(!is.null(r.init) | !is.null(l.init)) {
            if(verbose) { 
                cat("Setting composite lesser allele count ... \n")
            }
            E <- l/n.bulk
            n <- nrow(r)
            m <- ncol(r)
            mat <- do.call(rbind, parallel::mclapply(1:n, function(i) {
                do.call(cbind, lapply(1:m, function(j) {
                    ri <- r[i,j]
                    n.sci <- n.sc[i,j]
                    Ei <- E[i]
                    if(is.na(Ei)) {
                        mut.frac <- 0
                    }
                    else if(Ei <= 0.5) {
                        mut.frac <- ri
                    }
                    else if(Ei > 0.5) {
                        mut.frac <- n.sci-ri
                    }
                    else {
                        mut.frac <- 0
                    }
                    
                    ## f will be high if inconsistent
                    ## f will be low if consistent
                    ## f will be NaN if no coverage
                    ## use colorRamp from green to red
                    f <- mut.frac
                    return(f)
                }))
            }, mc.cores=n.cores))
            rownames(mat) <- rownames(r)
            colnames(mat) <- colnames(r)
            r.maf <- mat
            
            mat <- sapply(1:n, function(i) {
                li <- l[i]
                n.bulki <- n.bulk[i]
                Ei <- E[i]
                if(is.na(Ei)) {
                    mut.frac <- 0
                }
                else if(Ei <= 0.5) {
                    mut.frac <- li
                }
                else if(Ei > 0.5) {
                    mut.frac <- n.bulki-li
                }
                else {
                    mut.frac <- 0
                }
                
                ## f will be high if inconsistent
                ## f will be low if consistent
                ## f will be NaN if no coverage
                ## use colorRamp from green to red
                f <- mut.frac
                return(f)
            })
            l.maf <<- mat
        }
        
        snps.df <- rownames(r)
        snps.df <- data.frame(do.call(rbind,strsplit(snps.df,":|-| ")), stringsAsFactors=F)
        if(ncol(snps.df)==2) {
            snps.df <- cbind(snps.df, snps.df[,2])
        }
        colnames(snps.df) <- c('chr','start','end');
        if(!grepl('chr', snps.df[1,1])) {
            snps.df[,1] <- paste0('chr', snps.df[,1])
        }
        snps <- with(snps.df, GenomicRanges::GRanges(chr, IRanges::IRanges(as.numeric(as.character(start)), as.numeric(as.character(end)))))
        names(snps) <- rownames(r) <- rownames(r.maf) <- rownames(n.sc) <- names(l) <- names(l.maf) <- names(n.bulk) <- paste0(snps.df[,1], ':', snps.df[,2], '-', snps.df[,3])
        
        if(verbose) { 
            cat("Done setting initial allele matrices! \n")
        }
        
        return(list(
            snps=snps,
            r=r,
            r.maf=r.maf,
            n.sc=n.sc,
            l=l,
            l.maf=l.maf,
            n.bulk=n.bulk
        ))
    }


#' Maps snps to genes
#'
#' @name setGeneFactors
#' @param txdb TxDb object (ex. TxDb.Hsapiens.UCSC.hg19.knownGene). 
#' @param fill SNPs mapping to genes not annotated in txdb will be given unique IDs
#' @param verbose Verbosity
#'
#' @examples 
#' data(r)
#' data(cov.sc)
#' allele.mats <- setAlleleMats(r, cov.sc)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' geneFactor <- setGeneFactors(allele.mats$snps, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
setGeneFactors=function(snps, txdb, fill=TRUE, verbose=TRUE) {
        if(verbose) {
            cat("Mapping snps to genes ... \n")
        }
        gf <- ChIPseeker::annotatePeak(peak=snps, TxDb=txdb, verbose=verbose)
        gf.df <- data.frame(gf)$geneId
        names(gf.df) <- names(snps)
        if(!fill) {
            gf.df[is.na(data.frame(gf)$annotation)] <- NA
            gf.df <- na.omit(gf.df)
        }
        if(verbose) {
            cat("Done mapping snps to genes! \n")
        }
        
        return(gf.df)
}


#' Plot allele profile
#'
#' @name HoneyBADGER_plotAlleleProfile
#' @param r.sub SNP lesser allele count matrix for single cells. If NULL, object's r.maf will be used
#' @param n.sc.sub SNP coverage count matrix for single cells. If NULL, object's n.sc will be used 
#' @param l.sub SNP lesser allele count matrix for bulk refernece. If NULL, object's l.maf will be used
#' @param n.bulk.sub SNP coverage count matrix for bulk refernece. If NULL, object's n.bulk will be used
#' @param region Limit plotting to particular GenomicRanges regions
#' @param chrs Limit plotting to select chromosomes. Default autosomes only. (default: paste0('chr', c(1:22))) 
#' @param widths Widths of chromosomes in plot. If 'set' will depend on number of SNPs in region. Else will be equal.
#' @param cellOrder Order of cells. If 'set' will be automatically ordered by clustering. Else will be same order as input.
#' @param filter Remove sites with no coverage
#' @param max.ps Maximum point size for plot. 
#' 
#' @examples 
#' data(r)
#' data(cov.sc)
#' allele.mats <- setAlleleMats(r, cov.sc)
#' plotAlleleProfile(allele.mats$r.maf, allele.mats$n.sc, allele.mats$l.maf, allele.mats$n.bulk, widths=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 51304566, 48129895)/1e7) 
#' 
plotAlleleProfile=function(r.maf, n.sc, l.maf, n.bulk, region=NULL, chrs=paste0('chr', c(1:22)), widths=NULL, cellOrder=NULL, filter=FALSE, max.ps=3) {
        if(!is.null(region)) {
            overlap <- IRanges::findOverlaps(region, snps)
            ## which of the ranges did the position hit
            hit <- rep(FALSE, length(snps))
            hit[S4Vectors::subjectHits(overlap)] <- TRUE
            if(sum(hit) < 10) {
                cat(paste0("WARNING! ONLY ", sum(hit), " SNPS IN REGION! \n"))
            }
            vi <- hit
            r.maf <- r.maf[vi,]
            n.sc <- n.sc[vi,]
            l.maf <- l.maf[vi]
            n.bulk <- n.bulk[vi]
            
            chrs <- region@seqnames@values 
        }
        if(filter) {
            ## filter out snps without coverage
            vi <- rowSums(n.sc) > 0
            r.maf <- r.maf[vi,]
            n.sc <- n.sc[vi,]
            l.maf <- l.maf[vi]
            n.bulk <- n.bulk[vi]
        }
        
        mat <- r.maf/n.sc
        ## organize into chromosomes
        gos <- as.data.frame(snps)
        tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
            mat[rownames(gos)[ii[order(gos[ii,]$start,decreasing=F)]],,drop=FALSE]
        })
        ## only care about these chromosomes
        tl <- tl[chrs]
        if(!is.null(region)) {
            ## remove empty; need more than 1 gene
            vi <- unlist(lapply(tl, function(x) {
                if(is.null(x)) { return(FALSE) }
                else { return(nrow(x)>1) }
            }))
            tl <- tl[vi]
        }
        
        if(is.null(cellOrder)) {
            cellOrder <- colnames(r.maf)
        } else if (cellOrder[1]=='set') {
            avgd <- do.call(rbind, lapply(names(tl),function(nam) {
                d <- tl[[nam]]
                d <- colMeans(d, na.rm=TRUE)
                d
            }))
            hc <- hclust(dist(t(avgd)))
            cellOrder <- hc$order
        }
        
        r.maf <- r.maf[,cellOrder]
        n.sc <- n.sc[,cellOrder]
        
        r.tot <- cbind(r.maf/n.sc, 'Bulk'=l.maf/n.bulk)
        n.tot <- cbind(n.sc, 'Bulk'=n.bulk)
        
        require(ggplot2)
        require(reshape2)
        
        if(is.null(widths)) {
            widths <- rep(1, length(chrs))
        } else if (widths[1]=='set'){
            widths <- sapply(chrs, function(chr) {
                sum(grepl(paste0('^',chr,':'), rownames(r)))
            }); widths <- widths/max(widths)*100
        }                
        
        plist <- lapply(chrs, function(chr) {
            vi <- grepl(paste0('^',chr,':'), rownames(r.tot))
            m <- melt(t(r.tot[vi,]))
            colnames(m) <- c('cell', 'snp', 'alt.frac')
            rownames(m) <- paste(m$cell, m$snp)
            m$alt.frac[is.nan(m$alt.frac)] <- NA
            n <- melt(t(n.tot[vi,]))
            colnames(n) <- c('cell', 'snp', 'coverage')
            rownames(n) <- paste(n$cell, n$snp)
            n$coverage[n$coverage>30] <- 30  # max for visualization purposes
            ##n$coverage <- log10(n$coverage+1)
            n$coverage <- n$coverage^(1/3) # cube root for visualization purposes only
            n$coverage[n$coverage==0] <- NA # if no coverage, just don't show
            dat <- cbind(m, coverage=n$coverage)
            
            p <- ggplot(dat, aes(snp, cell)) +
                ## geom_tile(alpha=0) +
                geom_point(aes(colour = alt.frac, size = coverage), na.rm=TRUE) +
                scale_size_continuous(range = c(0, max.ps)) +
                ## scale_colour_gradientn(colours = rainbow(10)) +
                scale_colour_gradient2(mid="yellow", low = "turquoise", high = "red", midpoint=0.5) +
                theme(
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    legend.position="none",
                    plot.margin=unit(c(0,0,0,0), "cm"),
                    panel.border = element_rect(fill = NA, linetype = "solid", colour = "black"),
                    plot.title = element_text(hjust = 0.5)
                ) + labs(title = chr)
            ## theme(
            ##     ## axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(0.5),lineheight=1),
            ##     ## axis.text.y=element_blank(),
            ##     axis.title.y=element_blank(),
            ##     axis.ticks.y=element_blank(),
            ##     ##axis.text.y=element_text(size=rel(0.5))
            ##     legend.position="bottom"
            ##     ##panel.margin=unit(0 , "lines")
            ## )
            return(p)
        })
        
        require(gridExtra)
        do.call("grid.arrange", c(plist, list(ncol=length(plist), widths=widths)))

    }


#' Calculate posterior probability of CNVs using allele data
#'
#' @name HoneyBADGER_calcAlleleCnvProb
#' @param r.sub Optional matrix of alt allele count in single cells. If not provided, internal r.sc matrix is used.  
#' @param n.sub Optional matrix of site coverage count in single cells. If not provided, internal n.sc matrix is used.  
#' @param l.sub Optional vector of alt allele count in pooled single cells or bulk. If not provided, internal l vector is used.  
#' @param n.bulk.sub Optional vector of site coverage count in pooled single cells or bulk. If not provided, internal n.bulk vector is used.  
#' @param region GenomicRanges region of interest such as expected CNV boundaries. 
#' @param filter Boolean for whether to filter out SNP sites with no coverage. (default: TRUE)
#' @param pe Effective error rate to capture error from sequencing, etc. (default: 0.01)
#' @param mono Rate of mono-allelic expression. (default: 0.7)
#' @param quiet Boolean of whether to suppress progress bar. (default: TRUE)
#' @param verbose Verbosity(default: FALSE)
#' 
#' @examples 
#' data(r)
#' data(cov.sc)
#' allele.mats <- setAlleleMats(r, cov.sc)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' geneFactor <- setGeneFactors(allele.mats$snps, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' results <- calcAlleleCnvProb(allele.mats$r.maf, allele.mats$n.sc, allele.mats$l.maf, allele.mats$n.bulk, allele.mats$snps, geneFactor, region=GenomicRanges::GRanges('chr10', IRanges::IRanges(0,1e9)), verbose=TRUE)
#' 
calcAlleleCnvProb=function(r.maf, n.sc, l.maf, n.bulk, snps, geneFactor, region=NULL, filter=FALSE, pe=0.1, mono=0.7, verbose=FALSE) {
    quiet = !verbose
    
        if(!is.null(region)) {
            overlap <- IRanges::findOverlaps(region, snps)
            ## which of the ranges did the position hit
            hit <- rep(FALSE, length(snps))
            hit[S4Vectors::subjectHits(overlap)] <- TRUE
            if(sum(hit) <= 1) {
                cat(paste0("ERROR! ONLY ", sum(hit), " SNPS IN REGION! \n"))
                return();
            }
            if(sum(hit) < 10) {
                cat(paste0("WARNING! ONLY ", sum(hit), " SNPS IN REGION! \n"))
            }
            vi <- hit
            r.maf <- r.maf[vi,]
            n.sc <- n.sc[vi,]
            l.maf <- l.maf[vi]
            n.bulk <- n.bulk[vi]
            geneFactor <- geneFactor[vi]
            snps <- snps[vi,]
        }
        if(filter) {
            ## filter out snps without coverage
            vi <- rowSums(n.sc) > 0
            r.maf <- r.maf[vi,]
            n.sc <- n.sc[vi,]
            l.maf <- l.maf[vi]
            n.bulk <- n.bulk[vi]
            geneFactor <- geneFactor[vi]
            snps <- snps[vi,]
        }
        if(verbose) {
            cat('Assessing posterior probability of CNV in region ... \n')
            cat(paste0('with ', length(n.bulk), ' snps ... '))
        }
        genes.of.interest <- unique(geneFactor)
        if(verbose) {
            cat(paste0('within ', length(genes.of.interest), ' genes ... \n'))
        }
        
        ## associate each gene factor with a set of snps
        genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
            names(geneFactor)[which(geneFactor %in% genes.of.interest[i])]
        })
        names(genes2snps.dict) <- genes.of.interest
        
        ## Model
        if(verbose) {
            cat('converting to multi-dimensional arrays ... ')
        }
        
        ## Convert to multi-dimensions based on j
        I.j <- unlist(lapply(genes2snps.dict, length))
        numGenes <- length(genes2snps.dict)
        numSnpsPerGene <- max(I.j)
        numCells <- ncol(r.maf)
        ## j, i, k
        r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
        for(i in seq_len(numGenes)) {
            snpst <- genes2snps.dict[[i]]
            for(s in seq_along(snpst)) {
                r.array[i,s,] <- r.maf[snpst[s],]
            }
        }
        n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
        for(i in seq_len(numGenes)) {
            snpst <- genes2snps.dict[[i]]
            for(s in seq_along(snpst)) {
                n.sc.array[i,s,] <- n.sc[snpst[s],]
            }
        }
        l.array <- array(0, c(numGenes, numSnpsPerGene))
        for(i in seq_len(numGenes)) {
            snpst <- genes2snps.dict[[i]]
            for(s in seq_along(snpst)) {
                l.array[i,s] <- l.maf[snpst[s]]
            }
        }
        n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
        for(i in seq_len(numGenes)) {
            snpst <- genes2snps.dict[[i]]
            for(s in seq_along(snpst)) {
                n.bulk.array[i,s] <- n.bulk[snpst[s]]
            }
        }
        if(verbose) {
            cat('Aggregating data to list ... \n')
        }
        data <- list(
            'l' = l.array,
            'r' = r.array,
            'n.bulk' = n.bulk.array,
            'n.sc' = n.sc.array,
            'J' = length(I.j),  # how many genes
            'K' = ncol(r.maf),  # how many cells
            'I.j' = I.j,
            'pseudo' = pe,
            'mono' = mono)
        
        modelFile <- system.file("bug", "snpModel.bug", package = "HoneyBADGER")
        
        if(verbose) {
            cat('Running model ... \n')
        }
        require(rjags)
        # 4 random chains
        model <- rjags::jags.model(modelFile, data=data, n.chains=4, n.adapt=100, quiet=quiet)
        update(model, 100, progress.bar=ifelse(quiet,"none","text"))
        if(verbose) {
            cat('Done modeling!')
        }
        
        parameters <- 'S'
        samples <- coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine samples across chains
        pm <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
            c(mean(samples[,paste("S[",ci,"]",sep="")]))
        }))
        colnames(pm) <- colnames(r.maf)
        return(pm)
    }


