#' A Reference Class to represent single-cell RNA-seq data for HoneyBADGER analysis
#'
#' @field r single cell alternative allele count matrix
#' @field n.sc single cell snp coverage matrix
#' @field l bulk alternative allele count vector
#' @field n.bulk bulk snp coverage vector
#' @field r.maf single cell lesser allele count matrix
#' @field l.maf bulk lesser allele count vector
#' @field gexp.sc gene expression matrix
#' @field gexp.ref reference gene expression matrix
#' @field gexp.norm normalized gene expression matrix
#' @field mvFit estimated expression magnitude variance as a function of number of genes
#' @field dev expected absolute gene expression deviance due to CNV
#' @field snps GenomicRanges representation of snps in r
#' @field genes GenomicRanges representation of gene positions for genes in gexp.sc
#' @field geneFactor mapping of snps to genes
#' @field pred.snps.r matrix of cells and snps affected by cnv
#' @field pred.genes.r matrix of cells and snps affected by cnv
#' @field bound.snps.old temporary list of snps within cnv region used during recursion
#' @field bound.snps.final list of snps within cnv region
#' @field bound.genes.old temporary list of snps within cnv region used during recursion
#' @field bound.genes.final list of snps within cnv region
#' @field results retested posterior probabilities
#' @field summary summary of posterior probabilities
#' @field cnvs GenomicRanges of identified CNVs
#'
#' @export HoneyBADGER
#' @exportClass HoneyBADGER
#'
HoneyBADGER <- setRefClass(

    "HoneyBADGER",

    fields=c(
        'name', ## project name
        'r', ## r single cell alternative allele count matrix
        'n.sc', ## n.sc single cell snp coverage matrix
        'l', ## l bulk alternative allele count vector
        'n.bulk', ## n.bulk bulk snp coverage vector
        'r.maf', ## r.maf single cell lesser allele count matrix
        'l.maf', ## l.maf bulk lesser allele count vector
        'gexp.sc', ## gexp.sc gene expression matrix
        'gexp.ref', ## gexp.ref reference gene expression matrix
        'gexp.norm', ## gexp.norm normalized gene expression matrix
        'mvFit', ## mvFit estimated expression magnitude variance as a function of number of genes
        'dev', ## dev expected absolute gene expression deviance due to CNV
        'snps', ## snps GenomicRanges representation of snps in r
        'genes', ## genes GenomicRanges representation of gene positions for genes in gexp.sc
        'geneFactor', ## geneFactor mapping of snps to genes
        'pred.snps.r', ## pred.r matrix of cells and snps affected by cnv
        'pred.genes.r', ## pred.r matrix of cells and snps affected by cnv
        'bound.snps.old', ## bound.snps.old temporary list of snps within cnv region used during recursion
        'bound.snps.final', ## bound.snps.final list of snps within cnv region
        'bound.genes.old', ## bound.snps.old temporary list of snps within cnv region used during recursion
        'bound.genes.final', ## bound.snps.final list of snps within cnv region
        'results', ## results retested posterior probabilities
        'summary', ## summary of posterior probabilities
        'cnvs' ## GenomicRanges of identified CNVs
    ),

    methods = list(

        initialize=function(x=NULL, name='project') {

            name <<- name;
            r <<- NULL;
            n.sc <<- NULL;
            l <<- NULL;
            n.bulk <<- NULL;
            r.maf <<- NULL;
            l.maf <<- NULL;
            gexp.sc <<- NULL;
            gexp.ref <<- NULL;
            gexp.norm <<- NULL;
            mvFit <<- NULL;
            dev <<- NULL;
            snps <<- NULL;
            genes <<- NULL;
            geneFactor <<- NULL;

            pred.snps.r <<- NULL
            pred.genes.r <<- NULL
            bound.snps.old <<- c()
            bound.snps.final <<- list()
            bound.genes.old <<- c()
            bound.genes.final <<- list()

            results <<- list()
            summary <<- list()
            cnvs <<- list()
        }

    )
)


#' Set gene expression matrices, normalizes, and maps genes to genomic coordinates
#'
#' @name HoneyBADGER_setGexpMats
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
#'     dataset = 'hsapiens_gene_ensembl',
#'     host = "jul2015.archive.ensembl.org")
#' hb <- HoneyBADGER$new()
#' hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#'
HoneyBADGER$methods(
    setGexpMats=function(gexp.sc.init, gexp.ref.init, mart.obj, filter=TRUE, minMeanBoth=4.5, minMeanTest=6, minMeanRef=8, scale=TRUE, id="hgnc_symbol", verbose=TRUE) {
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
        gexp.sc <<- gexp.sc.init[vi,]
        gexp.ref <<- gexp.ref.init[vi,,drop=FALSE]

        if(filter) {
            vi <- (rowMeans(gexp.sc) > minMeanBoth & rowMeans(gexp.ref) > minMeanBoth) | rowMeans(gexp.sc) > minMeanTest | rowMeans(gexp.ref) > minMeanRef
            if(verbose) {
                cat(paste0(sum(vi), " genes passed filtering ... \n"))
            }
            gexp.sc <<- gexp.sc[vi,]
            gexp.ref <<- gexp.ref[vi,,drop=FALSE]
        }
        if(scale) {
            if(verbose) {
                cat("Scaling coverage ... \n")
            }
            ## library size
            gexp.sc <<- scale(gexp.sc)
            gexp.ref <<- scale(gexp.ref)
        }

        if(verbose) {
            cat(paste0("Normalizing gene expression for ", nrow(gexp.sc), " genes and ", ncol(gexp.sc), " cells ... \n"))
        }
        refmean <- rowMeans(gexp.ref)
        gexp.norm <<- gexp.sc - refmean

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
        genes <<- gs[gvi]

        ## remove genes with no position information
        gexp.norm <<- gexp.norm[gvi,]

        if(verbose) {
            cat("Done setting initial expression matrices! \n")
        }
    }
)


#' Plot gene expression profile
#'
#' @name HoneyBADGER_plotGexpProfile
#' @param gexp.norm.sub Optional normalized gene expression matrix. If not provided, internal normalized gene expression matrix is used.
#' @param chrs Chromosomes to be plotted (default: paste0('chr', c(1:22, 'X')))
#' @param window.size Window size for sliding window mean. Must be odd number. (default: 101)
#' @param zlim Limit for plotting heatmap (default: c(-2,2))
#' @param setOrder Boolean for whether or not to order cells (default: FALSE)
#' @param setWidths Boolean for whether or not to adjust widths of chomosomes based on number of genes. Otherwise uniform size. (default: FALSE)
#' @param order Order of cells (default: NULL)
#' @param defailt Boolean for whether to return detailed smoothed profiles (default: FALSE)
#'
#' @examples
#' data(gexp)
#' data(ref)
#' require(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#'     dataset = 'hsapiens_gene_ensembl',
#'     host = "jul2015.archive.ensembl.org")
#' hb <- HoneyBADGER$new()
#' hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE)
#' hb$plotGexpProfile()
#'
HoneyBADGER$methods(
    plotGexpProfile=function(gexp.norm.sub=NULL, chrs=paste0('chr', c(1:22)), region=NULL, window.size=101, zlim=c(-2,2), setOrder=FALSE, setWidths=FALSE, cellOrder=NULL, returnPlot=FALSE) {
        if(!is.null(gexp.norm.sub)) {
            gexp.norm <- gexp.norm.sub
            genes <- genes[rownames(gexp.norm.sub)]
        }
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

        if(!is.null(gexp.norm.sub) | !is.null(region)) {
            ## remove empty; need more than 1 gene
            vi <- unlist(lapply(tl, function(x) {
                if(is.null(x)) { return(FALSE) }
                else { return(nrow(x)>1) }
            }))
            tl <- tl[vi]
        }

        if(setWidths) {
            widths <- sapply(tl, nrow); widths <- widths/max(widths)*100
            ##Can also set by known chromosome size widths:
            ##https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
            ##widths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 51304566, 48129895)/1e7
        } else {
            widths <- rep(1, length(tl))
        }
        l <- layout(matrix(seq(1,length(tl)),1,length(tl),byrow=TRUE), widths=widths)

        if(setOrder) {
            avgd <- do.call(rbind, lapply(names(tl),function(nam) {
                d <- tl[[nam]]
                d <- colMeans(d)
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
            if(!is.null(cellOrder)) {
                d <- d[, cellOrder]
            }
            par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
            image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=FALSE, main=nam)
            box()
            return(d)
        })

        if(returnPlot) {
            return(tlsmooth)
        }
    }
)


#' Model expected gene expression variance as a function of number of genes
#'
#' @name HoneyBADGER_setMvFit
#' @param num.genes Number of random genes sampled (default: seq(5, 100, by=5))
#' @param rep Number of repeats/resampling (default: 50)
#' @param plot Whether to plot (default: FALSE)
#' @param verbose Verbosity (default: TRUE)
#'
HoneyBADGER$methods(
    setMvFit=function(num.genes = seq(5, 100, by=5), rep = 50, plot=FALSE, verbose=TRUE) {
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
                    require(ggplot2)
                    require(reshape2)
                    m <- melt(mat)
                    p <- ggplot(m) + geom_boxplot(aes(x = factor(Var2), y = value))
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
        mvFit <<- fits

        if(verbose) {
            cat('Done!')
        }
    }
)


#' Calculate posterior probability of CNVs using normalized expression data
#'
#' @name HoneyBADGER_calcGexpCnvProb
#' @param gexp.norm.sub Optional normalized gene expression matrix. If not provided, internal normalized gene expression matrix is used.
#' @param m Expected mean deviation due to copy number change (default: 0.15)
#' @param region Optional GenomicRanges region of interest such as expected CNV boundaries. (default: NULL)
#' @param quiet Boolean for whether to suppress progress display (default: TRUE)
#' @param verbose Verbosity (default: FALSE)
#'
HoneyBADGER$methods(
    calcGexpCnvProb=function(gexp.norm.sub=NULL, m=0.15, region=NULL, quiet=TRUE, verbose=FALSE) {
        if(!is.null(gexp.norm.sub)) {
            gexp.norm <- gexp.norm.sub
            genes <- genes[rownames(gexp.norm.sub)]
            mvFit <- mvFit[colnames(gexp.norm.sub)]
        }

        gexp <- gexp.norm
        gos <- genes
        fits <- mvFit

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
        ## mat <- apply(gexp, 2, runmean, k=window.size)
        mu0 <- apply(gexp, 2, mean)
        ng <- nrow(gexp)
        sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))
        ##gexp <- rbind(apply(gexp, 2, mean), apply(gexp, 2, median))

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
        require(rjags)
        model <- jags.model(modelFile, data=data, inits=inits, n.chains=4, n.adapt=100, quiet=quiet)
        update(model, 100, progress.bar=ifelse(quiet,"none","text"))

        parameters <- c('S', 'dd', 'mu')
        samples <- coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine chains

        if(verbose) {
            cat('...Done!')
        }

        snpLike <- samples
        v <- colnames(snpLike)
        S <- snpLike[,grepl('S', v)]
        dd <- snpLike[,grepl('dd', v)]
        mu <- snpLike[,grepl('mu', v)]
        ##plot(mu0, colMeans(mu))
        delcall <- apply(S*(1-dd), 2, mean)
        ampcall <- apply(S*dd, 2, mean)
        ##plot(mu0, delcall)
        ##plot(mu0, ampcall)
        names(ampcall) <- names(delcall) <- colnames(gexp)

        return(list('posterior probability of amplification'=ampcall,
                    'posterior probability of deletion'=delcall,
                    'estimated mean normalized expression deviation'=mu0))
    }
)


#' Set allele count matrices, creates in-silico pooled single cells as bulk reference if none provided
#'
#' @name HoneyBADGER_setAlleleMats
#' @param r.init SNP site alternate allele count matrix for single cells
#' @param n.sc.init SNP site coverage count matrix for single cells
#' @param l.init SNP site alternate allele counts for bulk reference. If NULL, in silico bulk will be created from single cells.
#' @param n.bulk.init SNP site coverage counts for bulk reference. If NULL, in silico bulk will be created from single cells.
#' @param het.deviance.threshold Deviation from expected 0.5 heterozygous fraction
#' @param min.cell Minimum number of cells a SNP must have coverage observed in
#' @param n.cores Number of cores
#' @param verbose Verbosity
#'
HoneyBADGER$methods(
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
            l <<- l.init
            n.bulk <<- n.bulk.init
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
            r <<- r.init[vi,]
            n.sc <<- n.sc.init[vi,]
            l <<- l[vi]
            n.bulk <<- n.bulk[vi]

            ## must have coverage in at least 5 cells
            if(verbose) {
                cat(paste0("must have coverage in at least ", min.cell, " cells ... \n"))
            }
            vi <- rowSums(n.sc > 0) >= min.cell
            cat(paste0(length(vi), " heterozygous SNPs identified \n"))
            r <<- r[vi,]
            n.sc <<- n.sc[vi,]
            l <<- l[vi]
            n.bulk <<- n.bulk[vi]

        } else {
            r <<- r.init
            n.sc <<- n.sc.init
        }
        if(is.null(r.maf) | is.null(l.maf) | !is.null(r.init) | !is.null(l.init)) {
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
            r.maf <<- mat

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
        snps <<- with(snps.df, GenomicRanges::GRanges(chr, IRanges::IRanges(as.numeric(as.character(start)), as.numeric(as.character(end)))))
        names(snps) <<- rownames(r) <<- rownames(r.maf) <<- rownames(n.sc) <<- names(l) <<- names(l.maf) <<- names(n.bulk) <<- apply(snps.df, 1, paste0, collapse=":")

        if(verbose) {
            cat("Done setting initial allele matrices! \n")
        }
    }
)


#' Maps snps to genes
#'
#' @name HoneyBADGER_setGeneFactors
#' @param txdb TxDb object ex. TxDb.Hsapiens.UCSC.hg19.knownGene
#' @param fill SNPs mapping to genes not annotated in txdb will be given unique IDs
#' @param verbose Verbosity
#'
HoneyBADGER$methods(
    setGeneFactors=function(txdb, fill=TRUE, verbose=TRUE) {
        if(verbose) {
            cat("Mapping snps to genes ... \n")
        }
        gf <- annotatePeak(peak=snps, TxDb=txdb)
        gf.df <- as.data.frame.csAnno(gf)$geneId
        names(gf.df) <- names(snps)
        if(!fill) {
            gf.df[is.na(as.data.frame.csAnno(gf)$annotation)] <- NA
            gf.df <- na.omit(gf.df)
        }
        geneFactor <<- gf.df
        if(verbose) {
            cat("Done mapping snps to genes! \n")
        }
    }
)


#' Plot allele profile
#'
#' @name HoneyBADGER_plotAlleleProfile
#' @param r.sub SNP lesser allele count matrix for single cells. If NULL, object's r.maf will be used
#' @param n.sc.sub SNP coverage count matrix for single cells. If NULL, object's n.sc will be used
#' @param l.sub SNP lesser allele count matrix for bulk refernece. If NULL, object's l.maf will be used
#' @param n.bulk.sub SNP coverage count matrix for bulk refernece. If NULL, object's n.bulk will be used
#' @param region Limit plotting to particular GenomicRanges regions
#' @param chrs Limit plotting to select chromosomes. Default autosomes only. (default: paste0('chr', c(1:22)))
#' @param setWidths Set widths on chromosomes based on SNP density. Otherwise will be uniform.
#' @param cellOrder Order of cells. Otherwise will be same order as in input matrix
#' @param filter Remove sites with no coverage
#' @param returnPlot Whether to return ggplot object
#' @param max.ps Maximum point size for plot
#'
#'
HoneyBADGER$methods(
    plotAlleleProfile=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, region=NULL, chrs=paste0('chr', c(1:22)), setWidths=FALSE, cellOrder=NULL, filter=FALSE, returnPlot=FALSE, max.ps=3) {
        if(!is.null(r.sub)) {
            r.maf <- r.sub
        }
        if(!is.null(n.sc.sub)) {
            n.sc <- n.sc.sub
        }
        if(!is.null(l.sub)) {
            l.maf <- l.sub
        }
        if(!is.null(n.bulk.sub)) {
            n.bulk <- n.bulk.sub
        }
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
        if(!is.null(cellOrder)) {
            r.maf <- r.maf[,cellOrder]
            n.sc <- n.sc[,cellOrder]
        }
        if(filter) {
            ## filter out snps without coverage
            vi <- rowSums(n.sc) > 0
            r.maf <- r.maf[vi,]
            n.sc <- n.sc[vi,]
            l.maf <- l.maf[vi]
            n.bulk <- n.bulk[vi]
        }

        r.tot <- cbind(r.maf/n.sc, 'Bulk'=l.maf/n.bulk)
        n.tot <- cbind(n.sc, 'Bulk'=n.bulk)

        require(ggplot2)
        require(reshape2)

        if(setWidths) {
            widths <- sapply(chrs, function(chr) {
                sum(grepl(paste0('^',chr,':'), rownames(r)))
            }); widths <- widths/max(widths)*100
        } else {
            widths <- rep(1, length(chrs))
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

        if(returnPlot) {
            return(plist)
        } else {
            require(gridExtra)
            do.call("grid.arrange", c(plist, ncol=length(plist)))
            ##print(p)
        }
    }
)


#' Plot smoothed allele profile
#'
#' @name HoneyBADGER_plotSmoothedAlleleProfile
#' @param inheritParams HoneyBADGER_plotAlleleProfile
#' @param window.size Size of sliding window for smoothing
#'
HoneyBADGER$methods(
  plotSmoothedAlleleProfile=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, region=NULL, chrs=paste0('chr', c(1:22)), setWidths=FALSE, setOrder=FALSE, cellOrder=NULL, filter=FALSE, returnPlot=FALSE, window.size=51) {
    if(!is.null(r.sub)) {
      r.maf <- r.sub
    }
    if(!is.null(n.sc.sub)) {
      n.sc <- n.sc.sub
    }
    if(!is.null(l.sub)) {
      l.maf <- l.sub
    }
    if(!is.null(n.bulk.sub)) {
      n.bulk <- n.bulk.sub
    }
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
    if(!is.null(cellOrder)) {
      r.maf <- r.maf[,cellOrder]
      n.sc <- n.sc[,cellOrder]
    }
    if(filter) {
      ## filter out snps without coverage
      vi <- rowSums(n.sc) > 0
      r.maf <- r.maf[vi,]
      n.sc <- n.sc[vi,]
      l.maf <- l.maf[vi]
      n.bulk <- n.bulk[vi]
    }

    r.tot <- cbind(r.maf/n.sc, 'Bulk'=l.maf/n.bulk)
    n.tot <- cbind(n.sc, 'Bulk'=n.bulk)

    mat <- r.maf/n.sc

    ## organize into chromosomes
    gos <- as.data.frame(snps)
    tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
      mat[rownames(gos)[ii[order(gos[ii,]$start,decreasing=F)]],,drop=FALSE]
    })
    ## only care about these chromosomes
    tl <- tl[chrs]

    if(!is.null(r.sub) | !is.null(region)) {
      ## remove empty; need more than 1 gene
      vi <- unlist(lapply(tl, function(x) {
        if(is.null(x)) { return(FALSE) }
        else { return(nrow(x)>1) }
      }))
      tl <- tl[vi]
    }

    ## https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
    ##chr.sizes <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 51304566, 48129895)
    ##l <- layout(matrix(seq(1, length(tl)),1,length(tl),byrow=T), widths=chr.sizes/1e7)
    if(setWidths) {
      widths <- sapply(tl, nrow); widths <- widths/max(widths)*100
    } else {
      widths <- rep(1, length(tl))
    }
    l <- layout(matrix(seq(1,length(tl)),1,length(tl),byrow=TRUE), widths=widths)

    if(setOrder) {
      avgd <- do.call(rbind, lapply(names(tl),function(nam) {
        d <- tl[[nam]]
        d <- colMeans(d, na.rm=TRUE)
        d
      }))
      hc <- hclust(dist(t(avgd)))
      cellOrder <- hc$order
    }

    tlsub <- tl

    pcol <- colorRampPalette(c('turquoise', 'yellow', 'red'))(256)

    ## plot chromosomes
    tlsmooth <- lapply(names(tlsub),function(nam) {
      d <- tlsub[[nam]]
      d <- apply(d,2,caTools::runmean,k=window.size, align="center")
      #d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
      if(!is.null(cellOrder)) {
        d <- d[, cellOrder]
      }
      par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
      image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=c(0,1), xlab="", ylab="", axes=FALSE, main=nam)
      box()
    })

    if(returnPlot) {
      return(tlsmooth)
    }

  }
)


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
#' @param n.iter Number of iterations in MCMC. (default: 1000)
#' @param quiet Boolean of whether to suppress progress bar. (default: TRUE)
#' @param verbose Verbosity(default: FALSE)
#'
HoneyBADGER$methods(
    calcAlleleCnvProb=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, region=NULL, filter=FALSE, pe=0.1, mono=0.7, n.iter=1000, quiet=TRUE, verbose=FALSE) {
        if(!is.null(r.sub)) {
            r.maf <- r.sub
            geneFactor <- geneFactor[rownames(r.sub)]
            snps <- snps[rownames(r.sub),]
        }
        if(!is.null(n.sc.sub)) {
            n.sc <- n.sc.sub
        }
        if(!is.null(l.sub)) {
            l.maf <- l.sub
        }
        if(!is.null(n.bulk.sub)) {
            n.bulk <- n.bulk.sub
        }
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
        model <- rjags::jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
        update(model, 300, progress.bar=ifelse(quiet,"none","text"))
        if(verbose) {
            cat('Done modeling!')
        }

        parameters <- 'S'
        samples <- coda.samples(model, parameters, n.iter=n.iter, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine samples across chains
        pm <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
            c(mean(samples[,paste("S[",ci,"]",sep="")]))
        }))
        colnames(pm) <- colnames(r.maf)
        return(pm)
    }
)


#' Set needed absolute gene expression deviance to be able to distinguish neutral from amplified or deletion regions
#'
#' @name HoneyBADGER_setGexpDev
#' @param alpha Alpha level (default: 0.05)
#' @param n Number of repeats for estimating parameter (default: 100)
#' @param seed Random seed
#' @param plot Plotting
#' @param verbose Verbosity
#'
HoneyBADGER$methods(
  setGexpDev=function(alpha=0.05, n=100, seed=0, plot=FALSE, verbose=FALSE) {
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
      dev <<- optim.dev
  }
)

#' Recursive HMM to identify CNV boundaries using normalized gene expression data
#'
#' @name HoneyBADGER_calcGexpCnvBoundaries
#' @param gexp.norm.sub Optional normalized gene expression matrix. Useful for manual testing of restricted regions. If NULL, gexp.norm within the HoneyBADGER object will be used.
#' @param chrs List of chromosome names. Genes not mapping to these chromosomes will be excluded. Default autosomes only: paste0('chr', c(1:22))
#' @param min.traverse Depth traversal to look for subclonal CNVs. Higher depth, potentially smaller subclones detectable. (default: 3)
#' @param min.num.genes Minimum number of genes within a CNV. (default: 5)
#' @param t HMM transition parameter. Higher number, more transitions. (default: 1e-6)
#' @param init Initialize recursion (default: FALSE)
#' @param verbose Verbosity (default: FALSE)
#' @param ... Additional parameters for \code{\link{calcGexpCnvProb}}
#'
HoneyBADGER$methods(
    calcGexpCnvBoundaries=function(gexp.norm.sub=NULL, chrs=paste0('chr', c(1:22)), min.traverse=5, min.num.genes=10, t=1e-6, init=FALSE, verbose=FALSE, ...) {
        if(!is.null(gexp.norm.sub)) {
            gexp.norm <- gexp.norm.sub
            genes <- genes[rownames(gexp.norm)]
        }
        if(init) {
            pred.genes.r <<- matrix(0, nrow(gexp.norm), ncol(gexp.norm))
            rownames(pred.genes.r) <<- rownames(gexp.norm)
            colnames(pred.genes.r) <<- colnames(gexp.norm)
            bound.genes.old <<- c()
            bound.genes.final <<- list()
        }
        if(is.null(bound.genes.final)) {
            cat('ERROR! USE init=TRUE! ')
            return()
        }

        ## remove old bound genes
        vi <- !(rownames(gexp.norm) %in% bound.genes.old)
        gexp.norm <- gexp.norm[vi,]

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
                    pd <- -dev
                    pn <- 0
                    pa <- dev
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
            vote[bound.genes.old] <- 0 ## do not want to rediscover old bounds

            if(verbose) {
                cat(paste0('max vote:', max(vote), '\n'))
            }
            if(max(vote)==0) {
                if(verbose) {
                  cat('Exiting; no new bound genes found.\n')
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

            ## test each of these highly confident deletions
            prob.info <- lapply(names(tbv), function(ti) {
                bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == ti]
                if(verbose) {
                    cat('GENES POTENTIALLY AFFECTED BY CNV: ')
                    cat(bound.genes.new)
                }
                ## now that we have boundaries, run on all cells
                prob <- calcGexpCnvProb(gexp.norm.sub=gexp.norm[bound.genes.new, ], m=dev, verbose=verbose, ...)

                if(verbose) {
                    cat("AMPLIFICATION PROBABILITY: ")
                    cat(prob[[1]])
                    cat("\n")

                    cat("DELETION PROBABILITY: ")
                    cat(prob[[2]])
                    cat("\n")
                }

                return(list('ap'=prob[[1]], 'dp'=prob[[2]], 'bs'=bound.genes.new))
            })
            ##cat(prob.info)
            return(prob.info)
        }
        amp.prob.info <- getTbv(lapply(boundgenes.pred, function(x) x[['amp']]))
        del.prob.info <- getTbv(lapply(boundgenes.pred, function(x) x[['del']]))

        del.prob <- do.call(rbind, lapply(seq_along(del.prob.info), function(i) del.prob.info[[i]]$dp))
        amp.prob <- do.call(rbind, lapply(seq_along(amp.prob.info), function(i) amp.prob.info[[i]]$ap))
        del.bound.genes.list <- lapply(seq_along(del.prob.info), function(i) del.prob.info[[i]]$bs)
        amp.bound.genes.list <- lapply(seq_along(amp.prob.info), function(i) amp.prob.info[[i]]$bs)

        prob.bin <- prob <- rbind(del.prob, amp.prob)
        bound.genes.list <- c(del.bound.genes.list, amp.bound.genes.list)

        if(sum(prob < 0.25) > 0) {
            prob.bin[prob < 0.25] <- 0
        }
        if(sum(prob > 0.75) > 0) {
            prob.bin[prob > 0.75] <- 1
        }
        if(sum(prob <= 0.75 & prob >= 0.25) > 0) {
            prob.bin[prob <= 0.75 & prob >= 0.25] <- NA
        }

        if(is.null(prob.bin)) {
            return()
        }
        if(ncol(prob.bin)>1) {
            dps <- rowSums(prob.bin, na.rm=TRUE)
            dpsi <- which(dps == max(dps))[1] ## pick one of the most clonal
            prob.fin <- t(as.matrix(prob[dpsi,]))
            bound.genes.new <- bound.genes.list[[dpsi]]
        } else {
            prob.fin <- prob
            bound.genes.new <- bound.genes.list
        }

        if(verbose) {
            print("CNV SNPS:")
            print(bound.genes.new)
            print(prob.fin)
        }

        ## need better threshold
        g1 <- colnames(prob.fin)[prob.fin > 0.75]
        if(verbose) {
          print("GROUP 1:")
          print(g1)
        }
        g2 <- colnames(prob.fin)[prob.fin <= 0.25]
        if(verbose) {
          print("GROUP 2:")
          print(g2)
        }
        ##clafProfile(r[, g1], n.sc[, g1], l, n.bulk)
        ##clafProfile(r[, g2], n.sc[, g2], l, n.bulk)

        ## record
        if(length(g1) > 0) {
            pred.genes.r[bound.genes.new, g1] <<- 1
            bound.genes.final[[length(bound.genes.final)+1]] <<- list(bound.genes.new)
        }

        ## all discovered
        bound.genes.old <<- unique(c(bound.genes.new, bound.genes.old))

        ## Recursion
        ##print('Recursion for Group1')
        if(length(g1)>=3) {
            tryCatch({
                calcGexpCnvBoundaries(gexp.norm.sub=gexp.norm[, g1])
            }, error = function(e) { cat(paste0("ERROR: ", e)) })
        }
        ##print('Recursion for Group2')
        if(length(g2)>=3) {
            tryCatch({
                calcGexpCnvBoundaries(gexp.norm.sub=gexp.norm[, g2])
            }, error = function(e) { cat(paste0("ERROR: ", e)) })
        }


    }
)


#' Recursive HMM to identify CNV boundaries using allele data
#'
#' @name HoneyBADGER_calcAlleleCnvBoundaries
#' @param r.sub Optional matrix of alt allele count in single cells. If not provided, internal r.sc matrix is used.
#' @param n.sub Optional matrix of site coverage count in single cells. If not provided, internal n.sc matrix is used.
#' @param l.sub Optional vector of alt allele count in pooled single cells or bulk. If not provided, internal l vector is used.
#' @param n.bulk.sub Optional vector of site coverage count in pooled single cells or bulk. If not provided, internal n.bulk vector is used.
#' @param min.traverse Depth traversal to look for subclonal CNVs. Higher depth, potentially smaller subclones detectable. (default: 3)
#' @param t HMM transition parameter. Higher number, more transitions. (default: 1e-6)
#' @param pd Probability of lesser allele detection in deleted region (ie. due to error)
#' @param pn Probability of lesser allele detection in neutral region (ie. 0.5 - error rate)
#' @param min.num.snps Minimum number of snps in candidate CNV
#' @param trim Trim boundary SNPs
#' @param init Boolean whether to initialize
#' @param verbose Verbosity(default: FALSE)
#' @param ... Additional parameters to pass to calcAlleleCnvProb
#'
HoneyBADGER$methods(
    calcAlleleCnvBoundaries=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, min.traverse=3, t=1e-6, pd=0.1, pn=0.45, min.num.snps=5, trim=0.1, init=FALSE, verbose=FALSE, ...) {

        if(!is.null(r.sub)) {
            r.maf <- r.sub
            snps <- snps[rownames(r.maf)]
            geneFactor <- geneFactor[rownames(r.maf)]
            snps <- rownames(r.maf)
        }
        if(!is.null(n.sc.sub)) {
            n.sc <- n.sc.sub
        }
        if(!is.null(l.sub)) {
            l.maf <- l.sub
        }
        if(!is.null(n.bulk.sub)) {
            n.bulk <- n.bulk.sub
        }

        if(init) {
            pred.snps.r <<- matrix(0, nrow(r.maf), ncol(r.maf))
            rownames(pred.snps.r) <<- rownames(r.maf)
            colnames(pred.snps.r) <<- colnames(r.maf)
            bound.snps.old <<- c()
            bound.snps.final <<- list()
        }

        if(is.null(bound.snps.final)) {
            cat('ERROR! USE init=TRUE! ')
            return()
        }

        if(verbose) {
            cat('ignore previously identified CNVs ... ')
        }
        ## remove old bound snps
        vi <- !(rownames(r.maf) %in% bound.snps.old)
        r.maf <- r.maf[vi,]
        n.sc <- n.sc[vi,]
        l.maf <- l.maf[vi]
        n.bulk <- n.bulk[vi]

        ## lesser allele fraction
        mat.tot <- r.maf/n.sc

        mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
        d <- dist(t(mat.smooth))
        d[is.na(d)] <- 0
        d[is.nan(d)] <- 0
        d[is.infinite(d)] <- 0
        hc <- hclust(d, method="ward.D2")

        if(verbose) {
            cat('iterative HMM ... ')
        }

        ## iterative HMM
        heights <- 1:min(min.traverse, ncol(r.maf))
        ## cut tree at various heights to establish groups
        boundsnps.pred <- lapply(heights, function(h) {

            ct <- cutree(hc, k = h)

            cuts <- unique(ct)

            ## look at each group, if deletion present
            boundsnps.pred <- lapply(cuts, function(group) {
                if(sum(ct==group)>1) {
                    mafl <- rowSums(r.maf[, ct==group]>0)
                    sizel <- rowSums(n.sc[, ct==group]>0)

                    ## change point
                    delta <- c(0, 1)
                    z <- HiddenMarkov::dthmm(mafl, matrix(c(1-t, t, t, 1-t), byrow=TRUE, nrow=2), delta, "binom", list(prob=c(pd, pn)), list(size=sizel), discrete=TRUE)
                    results <- HiddenMarkov::Viterbi(z)

                    ## Get boundaries from states
                    boundsnps <- rownames(r.maf)[results == 1]
                    return(boundsnps)
                }
            })
        })

        foo <- rep(0, nrow(r.maf)); names(foo) <- rownames(r.maf)
        foo[unique(unlist(boundsnps.pred))] <- 1
        ## vote
        vote <- rep(0, nrow(r.maf))
        names(vote) <- rownames(r.maf)
        lapply(boundsnps.pred, function(b) {
            vote[b[[1]]] <<- vote[b[[1]]] + 1
        })
        vote[bound.snps.old] <- 0 ## do not want to rediscover old bounds

        if(verbose) {
            cat(paste0('max vote:', max(vote), '\n'))
        }

        if(max(vote)==0) {
            if(verbose) {
              cat('Exiting; no new bound SNPs found.\n')
            }
            return() ## exit iteration, no more bound SNPs found
        }

        vote[vote > 0] <- 1
        ##bound.snps.new <- names(vote)[vote==max(vote)] ## newly discovered
        ## continous only
        ##mv <- max(vote)
        ##mv <- quantile(vote[vote!=0], 0.75)
        mv <- 1 ## at least 1 vote
        cs <- 1
        bound.snps.cont <- rep(0, length(vote))
        names(bound.snps.cont) <- names(vote)
        for(i in 2:length(vote)) {
            ##if(vote[i] == mv & vote[i] == vote[i-1]) {
            if(vote[i] >= mv & vote[i] == vote[i-1]) {
                bound.snps.cont[i] <- cs
            } else {
                cs <- cs + 1
            }
        }
        tb <- table(bound.snps.cont)
        tbv <- as.vector(tb); names(tbv) <- names(tb)
        tbv <- tbv[-1] # get rid of 0

        ## all detected deletions have fewer than 5 SNPs...reached the end
        tbv[tbv < min.num.snps] <- NA
        tbv <- na.omit(tbv)
        if(length(tbv)==0) {
            if(verbose) {
              cat(paste0('Exiting; less than ', min.num.snps, ' new bound SNPs found.\n'))
            }
            return()
        }

        ## test each of these highly confident deletions
        del.prob.info <- lapply(names(tbv), function(ti) {
            bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]

            ## trim
            bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]

            if(verbose) {
                cat('SNPS AFFECTED BY DELETION/LOH: ')
                cat(bound.snps.new)
            }

            ##clafProfile(r[bound.snps.new, hc$labels[hc$order]], n.sc[bound.snps.new, hc$labels[hc$order]], l[bound.snps.new], n.bulk[bound.snps.new])

            ## now that we have boundaries, run on all cells
            del.prob <- calcAlleleCnvProb(r.maf[bound.snps.new, ], n.sc[bound.snps.new, ], l.maf[bound.snps.new], n.bulk[bound.snps.new], region=NULL, n.iter=100, filter=FALSE, pe=pd, verbose=verbose, ...)

            if(verbose) {
                cat("DELETION/LOH PROBABILITY:")
                cat(del.prob)
            }

            return(list('dp'=del.prob, 'bs'=bound.snps.new))
        })
        ##print(del.prob.info)

        del.prob <- do.call(rbind, lapply(seq_along(del.prob.info), function(i) del.prob.info[[i]]$dp))
        bound.snps.list <- lapply(seq_along(del.prob.info), function(i) del.prob.info[[i]]$bs)

        del.prob.bin <- del.prob
        if(sum(del.prob < 0.25) > 0) {
            del.prob.bin[del.prob < 0.25] <- 0
        }
        if(sum(del.prob > 0.75) > 0) {
            del.prob.bin[del.prob > 0.75] <- 1
        }
        if(sum(del.prob <= 0.75 & del.prob >= 0.25) > 0) {
            del.prob.bin[del.prob <= 0.75 & del.prob >= 0.25] <- NA
        }

        if(length(tbv) > 1) {
            dps <- rowSums(del.prob.bin, na.rm=TRUE)
            dpsi <- which(dps == max(dps))[1] ## pick one of the most clonal
            del.prob.fin <- t(as.matrix(del.prob[dpsi,]))
            bound.snps.new <- bound.snps.list[[dpsi]]
        } else {
            del.prob.fin <- del.prob
            bound.snps.new <- unlist(bound.snps.list)
        }
        ##results.order <- order(del.prob)
        ##heatmap(cbind(del.prob[results.order], del.prob[results.order]), Rowv=NA, Colv=NA, scale="none", col=colorRampPalette(c("black", "red"))(100))
        ##clafProfile(r[, results.order], n.sc[, results.order], l, n.bulk)

        if(verbose) {
            cat("DELETION SNPS:")
            cat(bound.snps.new)
            cat(del.prob.fin)
        }

        ## need better threshold
        g1 <- colnames(del.prob.fin)[del.prob.fin > 0.75]
        if(verbose) {
            cat("GROUP1:")
            cat(g1)
        }
        g2 <- colnames(del.prob.fin)[del.prob.fin <= 0.25]
        if(verbose) {
            cat("GROUP 2:")
            cat(g2)
        }
        ##clafProfile(r[, g1], n.sc[, g1], l, n.bulk)
        ##clafProfile(r[, g2], n.sc[, g2], l, n.bulk)

        ## record
        if(length(g1) > 0) {
            pred.snps.r[bound.snps.new, g1] <<- 1
            bound.snps.final[[length(bound.snps.final)+1]] <<- list(bound.snps.new)
        }

        ## all discovered
        bound.snps.old <<- unique(c(bound.snps.new, bound.snps.old))
        ##sbs <- as.numeric(unlist(lapply(bound.snps, function(x) { strsplit(x, ':')[[1]][2] }))) ## sort
        ##names(sbs) <- bound.snps
        ##bound.snps <- names(sort(sbs))

        ## Recursion
        ##cat('Recursion for Group1')
        if(length(g1)>=3) {
            tryCatch({
                calcAlleleCnvBoundaries(r.sub=r.maf[, g1],
                                        n.sc.sub=n.sc[, g1],
                                        l.sub=rowSums(r.maf[, g1]>0),
                                        n.bulk.sub=rowSums(n.sc[, g1]>0)
                                        )
            }, error = function(e) { cat(paste0("ERROR: ", e)) })
        }
        ##cat('Recursion for Group2')
        if(length(g2)>=3) {
            tryCatch({
                calcAlleleCnvBoundaries(r.sub=r.maf[, g2],
                                        n.sc.sub=n.sc[, g2],
                                        l.sub=rowSums(r.maf[, g2]>0),
                                        n.bulk.sub=rowSums(n.sc[, g2]>0)
                                        )
            }, error = function(e) { cat(paste0("ERROR: ", e)) })
        }
    }
)


#' Calculate posterior probability of CNVs using normalized expression data and allele data
#'
#' @name HoneyBADGER_calcCombCnvProb
#' @inheritParams HoneyBADGER_calcGexpCnvProb
#' @inheritParams HoneyBADGER_calcAlleleCnvProb
#'
HoneyBADGER$methods(
    calcCombCnvProb=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, gexp.norm.sub=NULL, m=0.15, region=NULL, filter=FALSE, pe=0.1, mono=0.7, n.iter=1000, quiet=FALSE, verbose=FALSE) {
        if(!is.null(r.sub)) {
            r.maf <- r.sub
            geneFactor <- geneFactor[rownames(r.sub)]
            snps <- snps[rownames(r.sub),]
        }
        if(!is.null(n.sc.sub)) {
            n.sc <- n.sc.sub
        }
        if(!is.null(l.sub)) {
            l.maf <- l.sub
        }
        if(!is.null(n.bulk.sub)) {
            n.bulk <- n.bulk.sub
        }
        if(!is.null(gexp.norm.sub)) {
            gexp.norm <- gexp.norm.sub
            genes <- genes[rownames(gexp.norm.sub)]
            mvFit <- mvFit[colnames(gexp.norm.sub)]
        }
        gexp <- gexp.norm
        gos <- genes
        fits <- mvFit

        if(!is.null(region)) {
            ## limit gexp
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

            ## limit snp mats
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

        gexp <- gexp[, colnames(r.maf)]

        if(verbose) {
          cat('Assessing posterior probability of CNV in region ... \n')
          cat(paste0('with ', nrow(gexp), ' genes ... '))
          cat(paste0('and ', length(n.bulk), ' snps ... '))
        }

        genes.of.interest <- unique(geneFactor)
        if(verbose) {
            cat(paste0('... within ', length(genes.of.interest), ' genes ... \n'))
        }

        ## associate each gene factor with a set of snps
        genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
            names(geneFactor)[which(geneFactor %in% genes.of.interest[i])]
        })
        names(genes2snps.dict) <- genes.of.interest

        ## smooth
        ## mat <- apply(gexp, 2, runmean, k=window.size)
        mu0 <- apply(gexp, 2, mean)
        ng <- nrow(gexp)
        sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))
        ##gexp <- rbind(apply(gexp, 2, mean), apply(gexp, 2, median))

        if(verbose) {
          cat('Converting to multi-dimensional arrays...')
        }

        ## Convert to multi-dimensions based on j
        I.j <- unlist(lapply(genes2snps.dict, length))
        numGenes <- length(genes2snps.dict)
        numSnpsPerGene <- max(I.j)
        numCells <- ncol(r)
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
          cat('Aggregating data to list...')
        }
        data <- list(
            'l' = l.array,
            'r' = r.array,
            'n.bulk' = n.bulk.array,
            'n.sc' = n.sc.array,
            'J' = length(I.j),  # how many genes
            'K' = ncol(r),  # how many cells
            'I.j' = I.j,
            'pseudo' = pe,
            'mono' = mono,
            'gexp' = gexp,
            'JJ' = nrow(gexp),
            'sigma0' = sigma0,
            'mag0' = m
        )

        modelFile <- system.file("bug", "combinedModel.bug", package = "HoneyBADGER")

        if(verbose) {
            cat('Initializing model...')
        }

        model <- rjags::jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
        update(model, 300, progress.bar=ifelse(quiet,"none","text"))

        parameters <- c('S', 'dd')
        samples <- coda.samples(model, parameters, n.iter=300, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine chains

        snpLike <- samples
        v <- colnames(snpLike)
        S <- snpLike[,grepl('S', v)]
        dd <- snpLike[,grepl('dd', v)]
        ##plot(mu0, colMeans(mu))
        delcall <- apply(S*(1-dd), 2, mean)
        ##delcall
        ampcall <- apply(S*dd, 2, mean)
        ##ampcall
        ##plot(mu0, delcall)
        ##plot(mu0, ampcall)
        names(ampcall) <- names(delcall) <- colnames(gexp)

        return(list('posterior probability of amplification'=ampcall,
                    'posterior probability of deletion'=delcall)
               )
    }
)


#' Calculate posterior probability of CNVs for regions identified by the recursive HMM approach
#'
#' @name HoneyBADGER_retestIdentifiedCnvs
#' @param retestBoundGenes Boolean of whether to retest using expression model
#' @param retestBoundSnps Boolean of whether to retest using allele model
#' @param intersect If true will intersect regions identified from allele and expression HMMs. Otherwise, will union
#' @param verbose Verbosity
#' @param ... Additional parameters to pass to calcGexpCnvProb or calcAlleleCnvProb
#'
HoneyBADGER$methods(
    retestIdentifiedCnvs=function(retestBoundGenes=TRUE, retestBoundSnps=FALSE, intersect=FALSE, verbose=FALSE, ...) {
        if(retestBoundGenes) {
            if(verbose) {
                cat('Retesting bound genes ... ')
            }
            if(length(bound.genes.final)==0) {
                cat('ERROR NO GENES AFFECTED BY CNVS IDENTIFIED! Run calcGexpCnvBoundaries()? ')
            } else {
                ## retest <- lapply(seq_along(bound.genes.final), function(i) {
                ##     bgs <- bound.genes.final[[i]][[1]]
                ##     bgs <- bgs[trimGenes:(length(bgs)-trimGenes)]
                ##     x <- calcGexpCnvProb(gexp.norm[bgs,])
                ##     list(x[[1]], x[[2]])
                ## })
                rgs <- range(genes[unlist(bound.genes.final),])
                retest <- lapply(seq_len(length(rgs)), function(i) {
                    x <- calcGexpCnvProb(region=rgs[i], m=dev, verbose=verbose, ...)
                    list(x[[1]], x[[2]])
                })
                cnvs[['gene-based']] <<- list()
                cnvs[['gene-based']][['all']] <<- rgs
                results[['gene-based']] <<- retest
            }
        }

        if(retestBoundSnps) {
            if(verbose) {
                cat('Retesting bound snps ... ')
            }
            if(length(bound.snps.final)==0) {
                cat('ERROR NO SNPS AFFECTED BY CNVS IDENTIFIED! Run calcAlleleCnvBoundaries()? ')
            } else {
                ## retest <- lapply(seq_along(bound.snps.final), function(i) {
                ##     bgs <- bound.snps.final[[i]][[1]]
                ##     bgs <- bgs[trimSnps:(length(bgs)-trimSnps)]
                ##     x <- calcAlleleCnvProb(r.sub=r[bgs,], n.sc.sub=n.sc[bgs,], l.sub=l[bgs], n.bulk.sub=n.bulk[bgs])
                ##     x
                ## })
                rgs <- range(snps[unlist(bound.snps.final),])
                retest <- lapply(seq_len(length(rgs)), function(i) {
                    x <- calcAlleleCnvProb(region=rgs[i], verbose=verbose, ...)
                    x
                })
                cnvs[['allele-based']] <<- list()
                cnvs[['allele-based']][['all']] <<- rgs
                results[['allele-based']] <<- retest
            }
        }

        if(retestBoundSnps & retestBoundGenes) {
            if(verbose) {
                cat('Retesting bound snps and genes using joint model')
            }
            if(length(bound.genes.final)==0) {
                cat('ERROR NO GENES AFFECTED BY CNVS IDENTIFIED! Run calcGexpCnvBoundaries()? ')
                return();
            }
            if(length(bound.snps.final)==0) {
                cat('ERROR NO SNPS AFFECTED BY CNVS IDENTIFIED! Run calcAlleleCnvBoundaries()? ')
                return();
            }

            gr <- range(genes[unlist(bound.genes.final),])
            sr <- range(snps[unlist(bound.snps.final),])
            if(intersect) {
                rgs <- intersect(gr, sr)
            } else {
                rgs <- union(gr, sr)
            }

            retest <- lapply(seq_len(length(rgs)), function(i) {
                x <- calcCombCnvProb(region=rgs[i], m=dev, verbose=verbose, ...)
                list(x[[1]], x[[2]])
            })
            cnvs[['combine-based']] <<- list()
            cnvs[['combine-based']][['all']] <<- rgs
            results[['combine-based']] <<- retest
        }
    }
)


#' Summarize results
#'
#' @name HoneyBADGER_summarizeResults
#' @param geneBased Boolean of whether to summarize gene-based results
#' @param alleleBased Boolean of whether to summarize allele-based results
#' @param min.num.cells CNV must be present in this minimum number of cells
#' @param t Minimum probability threshold for call
#'
HoneyBADGER$methods(
    summarizeResults=function(geneBased=TRUE, alleleBased=FALSE, min.num.cells=2, t=0.75) {
        if(geneBased & !alleleBased) {
            rgs <- cnvs[['gene-based']][['all']]
            retest <- results[['gene-based']]
            amp.gexp.prob <- do.call(rbind, lapply(retest, function(x) x[[1]]))
            del.gexp.prob <- do.call(rbind, lapply(retest, function(x) x[[2]]))
            df <- cbind(as.data.frame(rgs),
                        avg.amp.gexp=rowMeans(amp.gexp.prob),
                        avg.del.gexp=rowMeans(del.gexp.prob),
                        amp.gexp.prob,
                        del.gexp.prob)

            ## filter to regions with at least some highly confident cells
            vi1 <- rowSums(amp.gexp.prob > t) > min.num.cells
            amp.gexp.prob <- amp.gexp.prob[vi1,] ## amplifications
            vi2 <- rowSums(del.gexp.prob > t) > min.num.cells
            del.gexp.prob <- del.gexp.prob[vi2,] ## amplifications

            names <- apply(as.data.frame(rgs), 1, paste0, collapse=":")
            rownames(amp.gexp.prob) <- paste0('amp', names[vi1])
            rownames(del.gexp.prob) <- paste0('del', names[vi2])
            ret <- rbind(del.gexp.prob, amp.gexp.prob)
            cnvs[['gene-based']][['amp']] <<- rgs[vi1]
            cnvs[['gene-based']][['del']] <<- rgs[vi2]
            summary[['gene-based']] <<- ret

            colnames(amp.gexp.prob) <- paste0('amp.gexp.', colnames(amp.gexp.prob))
            colnames(del.gexp.prob) <- paste0('del.gexp.', colnames(del.gexp.prob))

        }
        if(alleleBased & !geneBased) {
            rgs <- cnvs[['allele-based']][['all']]
            retest <- results[['allele-based']]
            del.loh.allele.prob <- do.call(rbind, lapply(retest, function(x) x))
            df <- cbind(as.data.frame(rgs), avg.del.loh.allele=rowMeans(del.loh.allele.prob), del.loh.allele.prob)

            ## filter to regions with at least some highly confident cells
            vi1 <- rowSums(del.loh.allele.prob > t) > min.num.cells
            del.loh.allele.prob <- del.loh.allele.prob[vi1,] ## amplifications

            names <- apply(as.data.frame(rgs), 1, paste0, collapse=":")
            rownames(del.loh.allele.prob) <- paste0('del.loh.', names[vi1])
            cnvs[['allele-based']][['del.loh']] <<- rgs[vi1]
            summary[['allele-based']] <<- del.loh.allele.prob

            colnames(del.loh.allele.prob) <- paste0('del.loh.allele.', colnames(del.loh.allele.prob))
        }
        if(alleleBased & geneBased) {
            rgs <- cnvs[['combine-based']][['all']]

            retest <- results[['combine-based']]
            amp.comb.prob <- do.call(rbind, lapply(retest, function(x) x[[1]]))
            del.comb.prob <- do.call(rbind, lapply(retest, function(x) x[[2]]))

            df <- cbind(as.data.frame(rgs), avg.amp.comb=rowMeans(amp.comb.prob), avg.del.comb=rowMeans(del.comb.prob), amp.comb.prob, del.comb.prob)

            ## filter to regions with at least some highly confident cells
            vi1 <- rowSums(amp.comb.prob > t) > min.num.cells
            amp.comb.prob <- amp.comb.prob[vi1,] ## amplifications
            vi2 <- rowSums(del.comb.prob > t) > min.num.cells
            del.comb.prob <- del.comb.prob[vi2,] ## amplifications

            names <- apply(as.data.frame(rgs), 1, paste0, collapse=":")
            rownames(amp.comb.prob) <- paste0('amp', names[vi1])
            rownames(del.comb.prob) <- paste0('del', names[vi2])
            ret <- rbind(del.comb.prob, amp.comb.prob)

            cnvs[['combine-based']][['amp']] <<- rgs[vi1]
            cnvs[['combine-based']][['del']] <<- rgs[vi2]
            summary[['combine-based']] <<- ret

            colnames(amp.comb.prob) <- paste0('amp.comb.', colnames(amp.comb.prob))
            colnames(del.comb.prob) <- paste0('del.comb.', colnames(del.comb.prob))
        }
        return(df)
    }
)

#' Plot posterior probability as heatmap
#'
#' @name HoneyBADGER_visualizeResults
#' @param geneBased Boolean whether to use gene-based results summary
#' @param alleleBased Boolean whether to use allele-based results summary
#' @param hc hclust object for cells
#' @param vc hclust object for CNVs
#' @param power Parameter to tweak clustering by weighing posterior probabilities
#' @param details Boolean whether to return details
#' @param ... Additional parameters to pass to heatmap
#'
HoneyBADGER$methods(
  visualizeResults=function(geneBased=TRUE, alleleBased=FALSE, hc=NULL, vc=NULL, power=1, details=FALSE, ...) {
    if(geneBased & !alleleBased) {
      df <- summary[['gene-based']]
    }
    if(alleleBased & !geneBased) {
      df <- summary[['allele-based']]
    }
    if(alleleBased & geneBased) {
      df <- summary[['combine-based']]
    }

    ## visualize as heatmap
    if(is.null(hc)) {
      hc <- hclust(dist(t(df^(power))), method='ward.D')
    }
    if(is.null(vc)) {
      vc <- hclust(dist(df^(power)), method='ward.D')
    }
    heatmap(t(df), Colv=as.dendrogram(vc), Rowv=as.dendrogram(hc), scale="none", col=colorRampPalette(c('beige', 'grey', 'black'))(100), ...)
    if(details) {
      return(list(hc=hc, vc=vc))
    }
  }
)
