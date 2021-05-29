#' Functions related to gene expression model


#' Set gene expression matrices, normalizes, and maps genes to genomic coordinates
#'
#' @param gexp.sc.init Single cell gene expression matrix
#' @param gexp.ref.init Reference gene expression matrix such as from GTEX or a match normal
#' @param mart.obj Biomart object used for mapping genes to genomic positions
#' @param filter Boolean of whether or not to filter genes (default: TRUE)
#' @param minMeanBoth Minimum mean gene expression in both the single cell and reference matrices (default: 4.5)
#' @param minMeanTest Minimum mean gene expression for the single cell expression matrix (default: 6)
#' @param minMeanRef Minimum mean gene expression for the reference expression matrix (default: 8)
#' @param scale Boolean of whether or not to scale by library size (default: TRUE)
#' @param id biomaRt ID for genes c(default: 'hgnc_symbol')
#' @param verbose Verbosity (default: TRUE)
#'
#' @examples 
#' data(gexp)
#' data(ref)
#' require(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#'    dataset = 'hsapiens_gene_ensembl', 
#'    host = "jul2015.archive.ensembl.org")
#' gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#' 
#' @export
#' 
setGexpMats=function(gexp.sc.init, gexp.ref.init, mart.obj, filter=TRUE, minMeanBoth=0, minMeanTest=mean(gexp.sc.init[gexp.sc.init!=0]), minMeanRef=mean(gexp.ref.init[gexp.ref.init!=0]), scale=TRUE, id="hgnc_symbol", verbose=TRUE) {
        if(verbose) {
            cat("Initializing expression matrices ... \n")
        }
        
        if(class(gexp.ref.init)!='Matrix') {
            gexp.ref.init <- as.matrix(gexp.ref.init)
        }
        
        vi <- intersect(rownames(gexp.sc.init), rownames(gexp.ref.init))
        if(length(vi) < 10) {
            cat('WARNING! GENE NAMES IN EXPRESSION MATRICES DO NOT SEEM TO MATCH! \n')
        }
        gexp.sc <- gexp.sc.init[vi,]
        gexp.ref <- gexp.ref.init[vi,,drop=FALSE]
        
        if(filter) {
            vi <- (rowMeans(gexp.sc) > minMeanBoth & rowMeans(gexp.ref) > minMeanBoth) | rowMeans(gexp.sc) > minMeanTest | rowMeans(gexp.ref) > minMeanRef
            if(verbose) {
                cat(paste0(sum(vi), " genes passed filtering ... \n"))
            }
            gexp.sc <- gexp.sc[vi,]
            gexp.ref <- gexp.ref[vi,,drop=FALSE]
        }
        if(scale) {
            if(verbose) {
                cat("Scaling coverage ... \n")
            }
            ## library size
            gexp.sc <- scale(gexp.sc)
            gexp.ref <- scale(gexp.ref)
        }
        
        if(verbose) {
            cat(paste0("Normalizing gene expression for ", nrow(gexp.sc), " genes and ", ncol(gexp.sc), " cells ... \n"))
        }
        refmean <- rowMeans(gexp.ref)
        gexp.norm <- gexp.sc - refmean
        
        gos <- getBM(values=rownames(gexp.norm),attributes=c(id, "chromosome_name","start_position","end_position"),filters=c(id),mart=mart.obj)
        ##gos$pos <- (gos$start_position + gos$end_position)/2
        rownames(gos) <- make.unique(gos[[id]])
        gos <- gos[rownames(gexp.norm),]
        
        if(nrow(gos)==0) {
            cat('ERROR! WRONG BIOMART GENE IDENTIFIER. Use ensembl_gene_id instead of hgnc_symbol? \n')
        }        
        
        if(length(grep('chr', gos$chromosome_name))==0) {
            gos$chromosome_name <- paste0('chr', gos$chromosome_name)
        }
        gos <- na.omit(gos)
        gs <- with(gos, GenomicRanges::GRanges(as.character(chromosome_name), IRanges::IRanges(as.numeric(as.character(start_position)), as.numeric(as.character(end_position)))))
        names(gs) <- rownames(gos)
        gvi <- intersect(rownames(gexp.norm), names(gs))
        genes <- gs[gvi]
        
        ## remove genes with no position information
        gexp.norm <- gexp.norm[gvi,]
        
        if(verbose) {
            cat("Done setting initial expression matrices! \n")
        }
        
        return(
            list(
                gexp.sc = gexp.sc,
                gexp.ref = gexp.ref,
                gexp.norm = gexp.norm,
                genes = genes
            )
        )
}


#' Plot gene expression profile
#'
#' @param gexp.norm Normalized gene expression matrix
#' @param genes GRanges annotation of gene names and coordinates
#' @param chrs Chromosomes to be plotted (default: paste0('chr', c(1:22, 'X')))
#' @param region Optional GenomicRanges region of interest such as expected CNV boundaries. (default: NULL)
#' @param window.size Window size for sliding window mean. Must be odd number. (default: 101)
#' @param zlim Limit for plotting heatmap (default: c(-2,2))
#' @param widths Widths of chromosomes in plot. If 'set' will depend on number of genes in region. Else will be equal.
#' @param cellOrder Order of cells. If 'set' will be automatically ordered by clustering. Else will be same order as input.
#' 
#' @examples 
#' data(gexp)
#' data(ref)
#' require(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#'     dataset = 'hsapiens_gene_ensembl', 
#'     host = "jul2015.archive.ensembl.org")
#' gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#' ##Set by known chromosome size widths:
#' ##https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
#' gexp.plot <- plotGexpProfile(gexp.mats$gexp.norm, gexp.mats$genes, 
#'     widths=c(249250621, 243199373, 198022430, 191154276, 180915260, 
#'     171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 
#'     133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 
#'     78077248, 59128983, 63025520, 51304566, 48129895)/1e7) 
#' 
#' @export
#' 
#' @import grDevices graphics stats
#' 
plotGexpProfile=function(gexp.norm, genes, chrs=paste0('chr', c(1:22)), region=NULL, window.size=101, zlim=c(-2,2), cellOrder=NULL, widths=NULL) {
        genes <- genes[rownames(gexp.norm)]

        if(!is.null(region)) {
            overlap <- IRanges::findOverlaps(region, genes)
            ## which of the ranges did the position hit
            hit <- rep(FALSE, length(genes))
            hit[S4Vectors::subjectHits(overlap)] <- TRUE
            if(sum(hit) < 10) {
                cat(paste0("WARNING! ONLY ", sum(hit), " GENES IN REGION! \n"))
            }
            vi <- hit
            gexp.norm <- gexp.norm[vi,]
            genes <- genes[rownames(gexp.norm)]
        }
        
        gos <- as.data.frame(genes)
        rownames(gos) <- names(genes)
        mat <- gexp.norm
        
        ## organize into chromosomes
        tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
            na.omit(mat[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
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
        
        if(is.null(widths)) {
            widths <- rep(1, length(tl))
        } else if(widths[1]=='set') {
            widths <- sapply(tl, nrow); widths <- widths/max(widths)*100
        } 
        l <- layout(matrix(seq(1,length(tl)),1,length(tl),byrow=TRUE), widths=widths)
        
        if(is.null(cellOrder)) {
            cellOrder <- colnames(gexp.norm)
        } else if(cellOrder[1]=='set') {
            avgd <- do.call(rbind, lapply(names(tl),function(nam) {
                d <- tl[[nam]]
                d <- apply(d,2,caTools::runmean,k=window.size, align="center")
                d
            }))
            hc <- hclust(dist(t(avgd)))
            cellOrder <- hc$order
        } 
        
        tlsub <- tl
        pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256)
        ## plot chromosomes
        tlsmooth <- lapply(names(tlsub),function(nam) {
            d <- tlsub[[nam]]
            d <- apply(d,2,caTools::runmean,k=window.size, align="center")
            d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
            d <- d[, cellOrder]
            par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
            image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=FALSE, main=nam)
            box()
            return(d)
        })
        
        return(
            list(
                tl=tl,
                tlsmooth=tlsmooth,
                cellOrder=cellOrder
            )
        )
    }


#' Model expected gene expression variance as a function of number of genes
#'
#' @param gexp.norm Normalized gene expression matrix
#' @param num.genes Number of random genes sampled (default: seq(5, 100, by=5))
#' @param rep Number of repeats/resampling (default: 50)
#' @param plot Whether to plot (default: FALSE)
#' @param verbose Verbosity (default: TRUE)
#' 
#' @examples 
#' data(gexp)
#' data(ref)
#' require(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#'     dataset = 'hsapiens_gene_ensembl', 
#'     host = "jul2015.archive.ensembl.org")
#' gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#' mvFit <- setMvFit(gexp.mats$gexp.norm)
#' 
#' @export
#' 
setMvFit=function(gexp.norm, num.genes = seq(5, 100, by=5), rep = 50, plot=FALSE, verbose=TRUE) {
        if(verbose) {
            cat('Modeling expected variance ... ')
        }
        mean.var.comp <- lapply(num.genes, function(ng) {
            set.seed(0)
            m <- do.call(rbind, lapply(1:rep, function(i) {
                nrmchr.sub <- gexp.norm[sample(1:nrow(gexp.norm), ng),]
                nm <- apply(nrmchr.sub, 2, mean)
                nm
            }))
            return(m)
        })
        names(mean.var.comp) <- num.genes
        
        fits <- lapply(1:ncol(gexp.norm), function(k) {
            ## for one cell
            mean.comp <- do.call(cbind, lapply(mean.var.comp, function(x) x[,k]))
            
            if(plot) {
                par(mfrow=c(1,3), mar=rep(2,4))
                perf.test <- function(mat) {
                    m <- reshape2::melt(mat)
                    p <- ggplot2::ggplot(m) + ggplot2::geom_boxplot(ggplot2::aes(x = factor(Var2), y = value))
                    return(p)
                }
                perf.test(mean.comp)
            }
            
            if(plot) {
                plot(log10(num.genes),log10(apply(mean.comp, 2, var)), type="l")
            }
            
            df <- data.frame('x'=num.genes, 'y'=apply(mean.comp, 2, var))
            fit <- lm(log10(y)~log10(x), data=df)
            
            if(plot) {
                x2 <- log10(num.genes)
                y2 <- predict(fit, x=x2, interval="predict")[, 'fit']
                plot(x2, y2)
                plot(log10(df$x), log10(df$y), type="l")
                points(x2, y2, type="l", col="red")
                plot(df$x, df$y, type="l")
                points(10^x2, 10^y2, type="l", col="red")
            }
            
            return(fit)
        })
        names(fits) <- colnames(gexp.norm)

        if(verbose) {
            cat('Done!')
        }
        
        return(fits)
    }


#' Set needed absolute gene expression deviance to be able to distinguish neutral from amplified or deletion regions
#' 
#' @param gexp.norm Normalized gene expression matrix
#' @param alpha Alpha level (default: 0.05)
#' @param n Number of repeats for estimating parameter (default: 100)
#' @param seed Random seed
#' @param plot Plotting
#' @param verbose Verbosity
#' 
#' @examples 
#' data(gexp)
#' data(ref)
#' require(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#'     dataset = 'hsapiens_gene_ensembl', 
#'     host = "jul2015.archive.ensembl.org")
#' gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#' dev <- setGexpDev(gexp.mats$gexp.norm)
#' 
#' @export
#' 
setGexpDev=function(gexp.norm, alpha=0.25, n=100, seed=0, plot=FALSE, verbose=FALSE) {
        k = 101
        set.seed(seed)
        gexp.sd <- sd(gexp.norm)
        devs <- seq_len(10)/10
        pvs <- unlist(lapply(devs, function(dev) {
            mean(unlist(lapply(seq_len(n), function(i) {
                pv <- ks.test(rnorm(k, 0, gexp.sd), rnorm(k, dev, gexp.sd))
                pv$p.value
            })))
        }))
        if(plot) {
            plot(pvs, devs, xlab="p-value", ylab="deviation", xlim=c(0,1))
        }
        fit <- lm(devs ~ pvs)
        optim.dev <- predict(fit, newdata=data.frame(pvs=alpha))
        if(verbose) {
            cat('Optimal deviance: ')
            cat(optim.dev)
        }
        return(optim.dev)
    }


#' Calculate posterior probability of CNVs using normalized expression data
#'
#' @param gexp.norm Normalized gene expression matrix.
#' @param genes GRanges annotation of gene names and coordinates
#' @param mvFit Mean variance fit
#' @param m Expression deviation due to copy number change (default: 0.15)
#' @param region Optional GenomicRanges region of interest such as expected CNV boundaries. (default: NULL)
#' @param verbose Verbosity (default: FALSE)
#'
#' @examples
#' data(gexp)
#' data(ref)
#' require(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#'     dataset = 'hsapiens_gene_ensembl', 
#'     host = "jul2015.archive.ensembl.org")
#' gexp.mats <- setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#' mvFit <- setMvFit(gexp.mats$gexp.norm)
#' results <- calcGexpCnvProb(gexp.mats$gexp.norm, gexp.mats$genes, 
#'     mvFit, region=GenomicRanges::GRanges('chr10', IRanges::IRanges(0,1e9)), verbose=TRUE)
#' 
#' @export
#' 
#' @import rjags
#' 
calcGexpCnvProb=function(gexp.norm, genes, mvFit, m=0.15, region=NULL, verbose=FALSE) {
        gexp <- gexp.norm
        gos <- genes[rownames(gexp.norm)]
        fits <- mvFit[colnames(gexp.norm)]
        quiet <- !verbose
        
        if(!is.null(region)) {
            overlap <- IRanges::findOverlaps(region, gos)
            ## which of the ranges did the position hit
            hit <- rep(FALSE, length(gos))
            names(hit) <- names(gos)
            hit[S4Vectors::subjectHits(overlap)] <- TRUE
            if(sum(hit) <= 1) {
                cat(paste0("ERROR! ONLY ", sum(hit), " GENES IN REGION! \n"))
                return();
            }
            if(sum(hit) < 3) {
                cat(paste0("WARNING! ONLY ", sum(hit), " GENES IN REGION! \n"))
            }
            vi <- hit
            if(verbose) {
                cat(paste0("restricting to ", sum(vi), " genes in region \n"))
            }
            if(sum(vi) <= 1) {
                pm <- rep(NA, ncol(gexp))
                names(pm) <- colnames(gexp)
                return(list(pm, pm, pm))
            }
            gexp <- gexp[vi,]
        }
        
        ## smooth
        mu0 <- apply(gexp, 2, mean)
        ng <- nrow(gexp)
        sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))
        
        ## Model
        if(verbose) {
            cat('Aggregating data to list ... \n')
        }
        data <- list(
            'K' = length(mu0),
            'JJ' = nrow(gexp),
            'gexp' = gexp,
            'sigma0' = sigma0,
            'mag0' = m
        )
        modelFile <-  system.file("bug", "expressionModel.bug", package = "HoneyBADGER")
        
        if(verbose) {
            cat('Initializing model ... \n')
        }
        ##model <- jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
        ##update(model, 1000, progress.bar=ifelse(quiet,"none","text"))
        inits <- list(
            list(S = rep(0, ncol(gexp)), dd = 0),
            list(S = rep(1, ncol(gexp)), dd = 0),
            list(S = rep(0, ncol(gexp)), dd = 1),
            list(S = rep(1, ncol(gexp)), dd = 1)
        )
        model <- rjags::jags.model(modelFile, data=data, inits=inits, n.chains=4, n.adapt=100, quiet=quiet)
        update(model, 100, progress.bar=ifelse(quiet,"none","text"))
        
        parameters <- c('S', 'dd')
        samples <- coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine chains
        
        if(verbose) {
            cat('...Done!')
        }
        
        snpLike <- samples
        v <- colnames(snpLike)
        S <- snpLike[,grepl('S', v)]
        dd <- snpLike[,grepl('dd', v)]
        delcall <- apply(S*(1-dd), 2, mean)
        ampcall <- apply(S*dd, 2, mean)
        names(ampcall) <- names(delcall) <- colnames(gexp)
        
        return(list('posterior probability of amplification'=ampcall,
                    'posterior probability of deletion'=delcall))
    }



#' HMM to identify CNV boundaries using normalized gene expression data
#' 
#' @param gexp.norm Normalized gene expression matrix
#' @param genes GRanges annotation of gene names and coordinates
#' @param m Expression magnitude deviation needed to distinguish CNV from neutral
#' @param chrs List of chromosome names. Genes not mapping to these chromosomes will be excluded. Default autosomes only: paste0('chr', c(1:22))
#' @param min.traverse Depth traversal to look for subclonal CNVs. Higher depth, potentially smaller subclones detectable. (default: 2)
#' @param min.num.genes Minimum number of genes within a CNV. (default: 3)
#' @param trim Trim boundary SNPs
#' @param t HMM transition parameter. Higher number, more transitions. (default: 1e-6)
#' @param verbose Verbosity (default: FALSE)
#' 
#' @export
#' 
#' @import stats
#' 
calcGexpCnvBoundaries=function(gexp.norm, genes, m=0.15, chrs=paste0('chr', c(1:22)), min.traverse=3, t=1e-6, min.num.genes=3, trim=0.1, verbose=FALSE) {

        genes <- genes[rownames(gexp.norm)]

        ## order
        gos <- as.data.frame(genes)
        rownames(gos) <- names(genes)
        gos <- gos[rownames(gexp.norm),]
        tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
            na.omit(gexp.norm[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],])
        })
        tl <- tl[chrs]
        gexp.norm <- do.call(rbind, lapply(tl, function(x) x))
        
        ## smooth
        k = 101
        mat.smooth <- apply(gexp.norm, 2, caTools::runmean, k)
        d <- dist(t(mat.smooth))
        d[is.na(d)] <- 0
        d[is.nan(d)] <- 0
        d[is.infinite(d)] <- 0
        hc <- hclust(d, method="ward.D2")
        
        ## iterative HMM
        heights <- seq_len(min(min.traverse, ncol(gexp.norm)))
        ## cut tree at various heights to establish groups
        boundgenes.pred <- lapply(heights, function(h) {
            ct <- cutree(hc, k = h)
            cuts <- unique(ct)
            ## look at each group, if deletion present
            boundgenes.pred <- lapply(cuts, function(group) {
                if(sum(ct==group)>1) {
                    mat.smooth <- apply(gexp.norm[, ct==group], 1, mean)
                    
                    ## change point
                    delta <- c(0, 1, 0)
                    t <- t
                    pd <- -m
                    pn <- 0
                    pa <- m
                    sd <- sd(mat.smooth)
                    z <- HiddenMarkov::dthmm(mat.smooth, matrix(c(1-2*t, t, t, t, 1-2*t, t, t, t, 1-2*t), byrow=TRUE, nrow=3), delta, "norm", list(mean=c(pd, pn, pa), sd=c(sd,sd,sd)))
                    results <- HiddenMarkov::Viterbi(z)
                    
                    ampgenes <- names(mat.smooth)[which(results==3)]
                    delgenes <- names(mat.smooth)[which(results==1)]
                    boundgenes <- list('amp'=ampgenes, 'del'=delgenes)
                    
                    return(boundgenes)
                }
            })
            
        })
        boundgenes.pred <- unlist(boundgenes.pred, recursive=FALSE)
        
        getTbv <- function(boundgenes.pred) {
          foo <- rep(0, nrow(gexp.norm)); names(foo) <- rownames(gexp.norm)
          foo[unique(unlist(boundgenes.pred))] <- 1
          ## vote
          vote <- rep(0, nrow(gexp.norm))
          names(vote) <- rownames(gexp.norm)
          lapply(boundgenes.pred, function(b) {
            vote[b] <<- vote[b] + 1
          })

          if(verbose) {
            cat(paste0('max vote:', max(vote), '\n'))
          }
          if(max(vote)==0) {
            if(verbose) {
              cat('Exiting; no genes found.\n')
            }
            return() ## exit iteration, no more bound genes found
          }
          
          vote[vote > 0] <- 1
          mv <- 1 ## at least 1 vote
          cs <- 1
          bound.genes.cont <- rep(0, length(vote))
          names(bound.genes.cont) <- names(vote)
          for(i in 2:length(vote)) {
            ##if(vote[i] == mv & vote[i] == vote[i-1]) {
            if(vote[i] >= mv & vote[i] == vote[i-1]) {
              bound.genes.cont[i] <- cs
            } else {
              cs <- cs + 1
            }
          }
          tb <- table(bound.genes.cont)
          tbv <- as.vector(tb); names(tbv) <- names(tb)
          tbv <- tbv[-1] # get rid of 0
          
          ## all detected deletions have fewer than 5 genes...reached the end
          tbv[tbv < min.num.genes] <- NA
          tbv <- na.omit(tbv)
          if(length(tbv)==0) {
            if(verbose) {
              cat(paste0('Exiting; fewer than ', min.num.genes, ' new bound genes found.\n'))
            }
            return()
          }
          
          boundgenes.info <- lapply(names(tbv), function(ti) {
            bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == ti]
            
            ## trim
            bound.genes.new <- bound.genes.new[1:round(length(bound.genes.new)-length(bound.genes.new)*trim)]
          })
          boundgenes.region <- do.call("c", lapply(boundgenes.info, function(bs) range(genes[bs])))
          
          return(list(info=boundgenes.info, regions=boundgenes.region))
        }
        
        amp.info <- getTbv(lapply(boundgenes.pred, function(x) x[['amp']]))
        del.info <- getTbv(lapply(boundgenes.pred, function(x) x[['del']]))
        
        if(verbose) {
          print(paste0("Identified ", length( amp.info$regions ), " potential amplifications"))
          print(paste0("Identified ", length( del.info$regions ), " potential deletions"))
        }
        return(list(amp=amp.info, del=del.info))
        
    }
