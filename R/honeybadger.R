honeybadger <- setRefClass(

    "honeybadger",

    fields=c(
        'r', ## r single cell alternative allele count matrix
        'n.sc', ## n.sc single cell snp coverage matrix
        'l', ## l bulk alternative allele count vector
        'n.bulk', ## n.bulk bulk snp coverage vector
        'r.maf', ## r.maf single cell minor allele count matrix
        'l.maf', ## l.maf bulk minor allele count vector
        'gexp.sc', ## gexp.sc gene expression matrix
        'gexp.ref', ## gexp.ref reference gene expression matrix
        'gexp.norm', ## gexp.norm normalized gene expression matrix
        'mvFit', ## mvFit estimated expression magnitude variance as a function of number of genes
        'snps', ## snps GenomicRanges representation of snps in r
        'genes', ## genes GenomicRanges representation of gene positions for genes in gexp.sc
        'geneFactor', ## geneFactor mapping of snps to genes
        'pred.snps.r', ## pred.r matrix of cells and snps affected by cnv
        'pred.genes.r', ## pred.r matrix of cells and snps affected by cnv
        'bound.snps.old', ## bound.snps.old temporary list of snps within cnv region used during recursion
        'bound.snps.final', ## bound.snps.final list of snps within cnv region
        'bound.genes.old', ## bound.snps.old temporary list of snps within cnv region used during recursion
        'bound.genes.final' ## bound.snps.final list of snps within cnv region
    ),

    methods = list(

        initialize=function() {
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
            snps <<- NULL;
            genes <<- NULL;
            geneFactor <<- NULL;

            pred.snps.r <<- NULL
            pred.genes.r <<- NULL
            bound.snps.old <<- c()
            bound.snps.final <<- list()
            bound.genes.old <<- c()
            bound.genes.final <<- list()
        },
        

        setGexpMats=function(gexp.sc.init, gexp.ref.init, mart.obj, filter=TRUE, minMinBoth=4.5, minMeanTest=6, minMeanRef=8, scale=TRUE) {
            cat("Initializing expression matrices ... \n")

            genes <- intersect(rownames(gexp.sc.init), rownames(gexp.ref.init))
            if(length(genes) < 10) {
                cat('WARNING! GENE NAMES IN EXPRESSION MATRICES DO NOT SEEM TO MATCH! \n')
            }
            gexp.sc <<- gexp.sc.init[genes,]
            gexp.ref <<- gexp.ref.init[genes,]

            if(filter) {
                vi <- (rowMeans(gexp.sc) > minMeanBoth & rowMeans(gexp.ref) > minMeanBoth) | rowMeans(gexp.sc) > minMeanTest | rowMeans(gexp.ref) > minMeanRef
                cat(paste0(sum(vi), " genes passed filtering ... \n"))
                gexp.sc <- gexp.sc[vi,]
                gexp.ref <- gexp.ref[vi,]
            }
            if(scale) {
                cat("scaling coverage ... \n")
                ## library size
                mat <- scale(mat)
                mat.ref <- scale(mat.ref)
            }

            cat(paste0("normalizing gene expression for ", length(genes), " genes and ", ncol(gexp.sc), " cells ... \n"))
            refmean <- rowMeans(gexp.ref)
            gexp.norm <<- gexp.sc - refmean

            gos <- getBM(values=rownames(gexp.norm),attributes=c("hgnc_symbol", "chromosome_name","start_position","end_position"),filters=c("hgnc_symbol"),mart=mart.obj)
            ##gos$pos <- (gos$start_position + gos$end_position)/2
            rownames(gos) <- make.unique(gos$hgnc_symbol)
            gos <- gos[rownames(gexp.norm),]
            
            require(GenomicRanges)
            if(!grepl('chr', gos$chromosome_name)) {
                gos$chromosome_name <- paste0('chr', gos$chromosome_name)
            }
            gos <- na.omit(gos)
            genes <- with(gos, GRanges(chromosome_name, IRanges(as.numeric(start_position), as.numeric(end_position)), strand=NULL))
            names(genes) <- rownames(gos)
            genes <<- genes

            ## remove genes with no position information
            gexp.norm <<- gexp.norm[names(genes),]

            cat("Done setting initial expression matrices! \n")
        },


        plotGexpProfile=function(chrs=paste0('chr', c(1:22, 'X')), window.size=101, zlim=c(-2,2), setOrder=FALSE, order=NULL, details=FALSE) {
            gos <- as.data.frame(genes)
            pos <- (gos$start+gos$end)/2
            mat <- gexp.norm
            
            ## organize into chromosomes
            tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
                na.omit(mat[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],])
            })
            ## only care about these chromosomes
            tl <- tl[chrs]

            ## https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
            #chr.sizes <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 51304566, 48129895)
            #l <- layout(matrix(seq(1, length(tl)),1,length(tl),byrow=T), widths=chr.sizes/1e7)
            l <- layout(matrix(seq(1,length(tlsub)),1,length(tlsub),byrow=T))

            if(setOrder) {
                avgd <- do.call(rbind, lapply(names(tl),function(nam) {
                    d <- tl[[nam]]
                    d <- colMeans(d)
                    d
                }))
                hc <- hclust(dist(t(avgd)))
                order <- hc$order
            }

            tlsub <- tl
            require(RColorBrewer)
            pcol <- rev(brewer.pal(11, 'RdBu'))
            ## plot chromosomes
            tlsmooth <- lapply(names(tlsub),function(nam) {
                d <- tlsub[[nam]]
                require(caTools)
                d <- apply(d,2,caTools::runmean,k=window.size, align="center")
                d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
                if(!is.null(order)) {
                    d <- d[, order]
                }
                par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
                image(1:nrow(d),1:ncol(d),d,col=pcol,zlim=zlim,xlab="",ylab="",axes=F,main=nam); box()
            })

            if(details) {
                return(tlsmooth)
            }
        },

        
        setMvFit=function(num.genes = seq(5, 100, by=5), rep = 50, plot=FALSE) {
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
        },

        calcGexpCnvProb=function(gexp.norm.sub=NULL, m=0.15, region=NULL, quiet=TRUE) {
            if(!is.null(gexp.norm.sub)) {
                gexp.norm <- gexp.norm.sub
                genes <- genes[rownames(gexp.norm.sub)]
                mvFit <- mvFit[colnames(gexp.norm.sub)]
            }
            
            gexp <- gexp.norm
            gos <- genes
            fits <- mvFit

            if(!is.null(region)) {
                require(GenomicRanges)
                overlap <- GenomicRanges::findOverlaps(region, gos)
                # which of the ranges did the position hit
                hit <- rep(FALSE, length(gos))
                names(hit) <- names(gos)
                hit[GenomicRanges::subjectHits(overlap)] <- TRUE
                if(sum(hit) < 3) {
                    cat(paste0("WARNING! ONLY ", sum(hit), " GENES IN REGION! \n"))
                }
                vi <- hit                
                cat(paste0("restricting to ", sum(vi), " genes in region \n"))
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
            cat('aggregating data to list ... \n')
            data <- list(
                'K' = length(mu0),
                'JJ' = nrow(gexp),
                'gexp' = gexp,
                'sigma0' = sigma0,
                'mag0' = m
            )
            modelFile <-  system.file("bug", "expressionModel.bug", package = "badger")

            cat('Initializing model ... \n')
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

            print('...Done!')

            snpLike <- samples
            v <- colnames(snpLike)
            S <- snpLike[,grepl('S', v)]
            dd <- snpLike[,grepl('dd', v)]
            mu <- snpLike[,grepl('mu', v)]
            ##plot(mu0, colMeans(mu))
            delcall <- apply(S*(1-dd), 2, mean)
            delcall
            ampcall <- apply(S*dd, 2, mean)
            ampcall
            ##plot(mu0, delcall)
            ##plot(mu0, ampcall)
            names(ampcall) <- names(delcall) <- colnames(gexp)

            return(list('posterior probability of amplification'=ampcall,
                        'posterior probability of deletion'=delcall,
                        'estimated mean normalized expression deviation'=mu0))
        },

        setAlleleMats=function(r.init, n.sc.init, l.init=NULL, n.bulk.init=NULL, filter=TRUE, het.deviance.threshold=0.05, min.cell=3, n.cores=1) {
            cat("Initializing allele matrices ... \n")

            if(is.null(l.init) | is.null(n.bulk.init)) {
                cat("creating in-silico bulk ... \n")
                cat(paste0("using ", ncol(r.init), " cells ... \n"))
                l <<- rowSums(r.init>0)
                n.bulk <<- rowSums(n.sc.init>0)
            } else {
                l <<- l.init
                n.bulk <<- n.bulk.init
            }

            if(filter) {
                cat("filtering for putative heterozygous snps ... \n")
                cat(paste0("allowing for a ", het.deviance.threshold, " deviation from the expected 0.5 heterozygous allele fraction ... \n"))
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
                cat(paste0("must have coverage in at least ", min.cell, " cells ... \n"))
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
            if(is.null(r.maf) | is.null(l.maf)) {
                cat("setting composite minor allele count ... \n")
                E <- l/n.bulk
                n <- nrow(r)
                m <- ncol(r)
                require(parallel)
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
                l.maf <<- rowSums(r.maf>0)
            }

            snps.df <- rownames(r)
            snps.df <- data.frame(do.call(rbind,strsplit(snps.df,":|-")), stringsAsFactors=F)
            if(ncol(snps.df)==2) {
                snps.df <- cbind(snps.df, snps.df[,2])
            }
            colnames(snps.df) <- c('chr','start','end');
            require(GenomicRanges)
            if(!grepl('chr', snps.df[1,1])) {
                snps.df[,1] <- paste0('chr', snps.df[,1])
            }
            snps <<- with(snps.df, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand=NULL))
            rownames(r) <<- rownames(r.maf) <<- rownames(n.sc) <<- names(l) <<- names(l.maf) <<- names(n.bulk) <<- paste0(snps)

            cat("Done setting initial allele matrices! \n")
        },

        setGeneFactors=function(txdb, fill=TRUE, gene=TRUE) {
            cat("Mapping snps to genes ... \n")
            require(ChIPseeker)
            gf <- ChIPseeker::annotatePeak(peak=snps, TxDb=txdb, verbose = FALSE, genomicAnnotationPriority="Exon")
            gf.df <- data.frame(gf)$geneId
            names(gf.df) <- paste0(snps)
            if(!fill) {
                gf.df[is.na(data.frame(gf)$annotation)] <- NA
                gf.df <- na.omit(gf.df)
            }
            geneFactor <<- gf.df
            cat("Done mapping snps to genes! \n")
        },

        plotClafProfile=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, region=NULL, order=NULL, filter=FALSE, return.plot=FALSE) {
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
                require(GenomicRanges)
                overlap <- GenomicRanges::findOverlaps(gr, snps)
                # which of the ranges did the position hit
                hit <- rep(FALSE, length(snps))
                hit[GenomicRanges::subjectHits(overlap)] <- TRUE
                if(sum(hit) < 10) {
                    cat(paste0("WARNING! ONLY ", sum(hit), " SNPS IN REGION! \n"))
                }
                vi <- hit
                r.maf <- r.maf[vi,]
                n.sc <- n.sc[vi,]
                l.maf <- l.maf[vi]
                n.bulk <- n.bulk[vi]
            }
            if(!is.null(order)) {
                r.maf <- r.maf[,order]
                n.sc <- n.sc[,order]
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

            m <- melt(t(r.tot))
            colnames(m) <- c('cell', 'snp', 'alt.frac')
            rownames(m) <- paste(m$cell, m$snp)
            m$alt.frac[is.nan(m$alt.frac)] <- NA
            n <- melt(t(n.tot))
            colnames(n) <- c('cell', 'snp', 'coverage')
            rownames(n) <- paste(n$cell, n$snp)
            n$coverage[n$coverage>30] <- 30  # max for visualization purposes
            ##n$coverage <- log10(n$coverage+1)
            n$coverage <- n$coverage^(1/3) # cube root for visualization purposes only
            dat <- cbind(m, coverage=n$coverage)

            p <- ggplot(dat, aes(snp, cell)) +
                ## geom_tile(alpha=0) +
                geom_point(aes(colour = alt.frac, size = coverage)) +
                scale_size_continuous(range = c(0,3)) +
                ## scale_colour_gradientn(colours = rainbow(10)) +
                scale_colour_gradient2(mid="yellow", low = "turquoise", high = "red", midpoint=0.5) +
                theme(
                    ## axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(0.5),lineheight=1),
                    ## axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    ##axis.text.y=element_text(size=rel(0.5))
                    legend.position="bottom"
                    ##panel.margin=unit(0 , "lines")
                )
            if(return.plot) {
                return(p)
            } else {
                print(p)
            }
        },

        calcAlleleCnvProb=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, region=NULL, filter=FALSE, pe=0.1, mono=0.7, n.iter=1000, quiet=FALSE) {
            if(!is.null(r.sub)) {
                r.maf <- r.sub
                geneFactor <- geneFactor[rownames(r.sub)]
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
                require(GenomicRanges)
                overlap <- GenomicRanges::findOverlaps(region, snps)
                # which of the ranges did the position hit
                hit <- rep(FALSE, length(snps))
                hit[GenomicRanges::subjectHits(overlap)] <- TRUE
                if(sum(hit) < 10) {
                    cat(paste0("WARNING! ONLY ", sum(hit), " SNPS IN REGION! \n"))
                }
                vi <- hit
                r.maf <- r.maf[vi,]
                n.sc <- n.sc[vi,]
                l.maf <- l.maf[vi]
                n.bulk <- n.bulk[vi]
                geneFactor <- geneFactor[vi]
            }
            if(filter) {
                ## filter out snps without coverage
                vi <- rowSums(n.sc) > 0
                r.maf <- r.maf[vi,]
                n.sc <- n.sc[vi,]
                l.maf <- l.maf[vi]
                n.bulk <- n.bulk[vi]
                geneFactor <- geneFactor[vi]
            }
            cat('Assessing posterior probability of CNV in region ... \n')
            cat(paste0('with ', length(n.bulk), ' snps ... '))

            genes.of.interest <- unique(geneFactor)
            cat(paste0('within ', length(genes.of.interest), ' genes ... \n'))

            ## associate each gene factor with a set of snps
            genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
                names(geneFactor)[which(geneFactor %in% genes.of.interest[i])]
            })
            names(genes2snps.dict) <- genes.of.interest

            ## Model
            cat('converting to multi-dimensional arrays ... ')

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
            cat('aggregating data to list ... \n')
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

            ##modelFile <- system.file("bug", "snpModel.bug", package = "badger")
            modelFile <- "/home/jfan/Projects/badger/inst/bug/snpModel.bug"

            cat('Running model ... \n')
            require('rjags')
            model <- rjags::jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
            update(model, 300, progress.bar=ifelse(quiet,"none","text"))
            cat('Done modeling!')

            parameters <- 'S'
            samples <- coda.samples(model, parameters, n.iter=n.iter, progress.bar=ifelse(quiet,"none","text"))
            samples <- do.call(rbind, samples) # combine samples across chains
            pm <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
                c(mean(samples[,paste("S[",ci,"]",sep="")]))
            }))
            colnames(pm) <- colnames(r.maf)
            return(pm)
        },

        
        calcGexpCnvBoundaries=function(gexp.norm.sub=NULL, min.traverse=3, min.num.genes=5, init=FALSE) {
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
            chrs=paste0('chr', c(1:22, 'X'))
            tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
                na.omit(mat[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],])
            })
            tl <- tl[chrs]
            gexp.norm <- do.call(rbind, lapply(tl, function(x) x))

            ## smooth
            mat.smooth <- apply(gexp.norm, 2, caTools::runmean, k=101)
            d <- dist(t(mat.smooth))
            d[is.na(d)] <- 0
            d[is.nan(d)] <- 0
            d[is.infinite(d)] <- 0
            hc <- hclust(d, method="ward.D2")

            ## iterative HMM
            heights <- 1:min(min.traverse, ncol(gexp.norm))
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
                        require(HiddenMarkov)
                        t <- 1e-6
                        pd <- -1
                        pn <- 0
                        pa <- 1
                        sd <- sd(mat.smooth)
                        z <- HiddenMarkov::dthmm(mat.smooth, matrix(c(1-t, t, t, t, 1-t, t, t, t, 1-t), byrow=TRUE, nrow=3), delta, "norm", list(mean=c(pd, pn, pa), sd=c(sd,sd,sd)))
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
                
                print(paste0('max vote:', max(vote)))
                if(max(vote)==0) {
                    return() ## exit iteration, no more bound SNPs found
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
                    return()
                }
                
                ## test each of these highly confident deletions
                prob.info <- lapply(names(tbv), function(ti) {
                    bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == ti]
                    
                    cat('GENES AFFECTED BY CNV: ')
                    cat(bound.genes.new)
                    
                    ## now that we have boundaries, run on all cells
                    prob <- calcGexpCnvProb(gexp.norm.sub=gexp.norm[bound.genes.new, ])
                    
                    cat("AMPLIFICATION PROBABILITY: ")
                    print(prob[[1]])
                    
                    cat("DELETION PROBABILITY: ")
                    print(prob[[2]])
                    
                    return(list('ap'=prob[[1]], 'dp'=prob[[2]], 'bs'=bound.genes.new))
                })
                print(prob.info)
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

            if(ncol(prob.bin>1)) {
                dps <- rowSums(prob.bin, na.rm=TRUE)
                dpsi <- which(dps == max(dps))[1] ## pick one of the most clonal
                prob.fin <- t(as.matrix(prob[dpsi,]))
                bound.genes.new <- bound.genes.list[[dpsi]]
            } else {
                prob.fin <- prob
                bound.genes.new <- bound.genes.list
            }

            print("CNV SNPS:")
            print(bound.genes.new)
            print(prob.fin)

            ## need better threshold
            g1 <- colnames(prob.fin)[prob.fin > 0.75]
            print("GROUP1:")
            print(g1)
            g2 <- colnames(prob.fin)[prob.fin <= 0.25]
            print("GROUP2:")
            print(g2)
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
            print('Recursion for Group1')
            if(length(g1)>=3) {
                tryCatch({
                    calcGexpCnvBoundaries(gexp.norm.sub=gexp.norm[, g1])
                }, error = function(e) { print("error detected"); print(e) })
            }
            print('Recursion for Group2')
            if(length(g2)>=3) {
                tryCatch({
                    calcGexpCnvBoundaries(gexp.norm.sub=gexp.norm[, g2])
                }, error = function(e) { print("error detected"); print(e) })
            }            

            
        },

        calcAlleleCnvBoundaries=function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, min.traverse=3, t=1e-5, pd=0.1, pn=0.45, min.num.snps=5, init=FALSE) {
            if(!is.null(r.sub)) {
                r.maf <- r.sub
                snps <- snps[rownames(r.maf)]
                geneFactor <- geneFactor[rownames(r.maf)]
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

            ## remove old bound snps
            vi <- !(rownames(r.maf) %in% bound.snps.old)
            r.maf <- r.maf[vi,]
            n.sc <- n.sc[vi,]
            l.maf <- l.maf[vi]
            n.bulk <- n.bulk[vi]

            ## Minor allele fraction
            mat.tot <- r.maf/n.sc

            mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
            d <- dist(t(mat.smooth))
            d[is.na(d)] <- 0
            d[is.nan(d)] <- 0
            d[is.infinite(d)] <- 0
            hc <- hclust(d, method="ward.D2")

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
                        require(HiddenMarkov)
                        z <- HiddenMarkov::dthmm(mafl, matrix(c(1-t, t, t, 1-t), byrow=TRUE, nrow=2), delta, "binom", list(prob=c(pd, pn)), list(size=sizel), discrete=TRUE)
                        results <- HiddenMarkov::Viterbi(z)

                        ## Get boundaries from states
                        boundsnps <- rownames(r.maf)[results == 1]
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

            print(paste0('max vote:', max(vote)))
            if(max(vote)==0) {
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

            ## all detected deletions have fewer than 5 genes...reached the end
            tbv[tbv < min.num.snps] <- NA
            tbv <- na.omit(tbv)
            if(length(tbv)==0) {
                return()
            }
            cat(tbv)

            ## test each of these highly confident deletions
            del.prob.info <- lapply(names(tbv), function(ti) {
                bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]

                print('SNPS AFFECTED BY DELETION')
                print(bound.snps.new)

                ##clafProfile(r[bound.snps.new, hc$labels[hc$order]], n.sc[bound.snps.new, hc$labels[hc$order]], l[bound.snps.new], n.bulk[bound.snps.new])

                ## now that we have boundaries, run on all cells
                del.prob <- calcAlleleCnvProb(r.maf[bound.snps.new, ], n.sc[bound.snps.new, ], l.maf[bound.snps.new], n.bulk[bound.snps.new], region=NULL, n.iter=100, filter=FALSE, pe=pd, quiet=FALSE)

                print("DELETION PROBABILITY")
                print(del.prob)

                return(list('dp'=del.prob, 'bs'=bound.snps.new))
            })
            print(del.prob.info)

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

            print("DELETION SNPS:")
            print(bound.snps.new)
            print(del.prob.fin)

            ## need better threshold
            g1 <- colnames(del.prob.fin)[del.prob.fin > 0.75]
            print("GROUP1:")
            print(g1)
            g2 <- colnames(del.prob.fin)[del.prob.fin <= 0.25]
            print("GROUP2:")
            print(g2)
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
            print('Recursion for Group1')
            if(length(g1)>=3) {
                tryCatch({
                    calcAlleleCnvBoundaries(r.sub=r.maf[, g1],
                                 n.sc.sub=n.sc[, g1],
                                 l.sub=rowSums(r.maf[, g1]>0),
                                 n.bulk.sub=rowSums(n.sc[, g1]>0)
                                 )
                }, error = function(e) { print("error detected"); print(e) })
            }
            print('Recursion for Group2')
            if(length(g2)>=3) {
                tryCatch({
                    calcAlleleCnvBoundaries(r.sub=r.maf[, g2],
                                 n.sc.sub=n.sc[, g2],
                                 l.sub=rowSums(r.maf[, g2]>0),
                                 n.bulk.sub=rowSums(n.sc[, g2]>0)
                                 )
                }, error = function(e) { print("error detected"); print(e) })
            }
        }

    ) ## end methods
)
