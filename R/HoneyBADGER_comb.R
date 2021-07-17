#' Calculate posterior probability of CNVs using normalized expression data and allele data 
#'
#' @inheritParams calcGexpCnvProb 
#' @inheritParams calcAlleleCnvProb 
#' 
#' @export
#' 
calcCombCnvProb=function(gexp.norm, genes, mvFit, m=0.15, r.maf, n.sc, l.maf, n.bulk, snps, geneFactor, region=NULL, filter=FALSE, pe=0.1, mono=0.7, verbose=FALSE) {
  
  quiet <- !verbose
  
  geneFactor <- geneFactor[rownames(r.maf)]
  snps <- snps[rownames(r.maf),]
  
  genes <- genes[rownames(gexp.norm)]
  mvFit <- mvFit[colnames(gexp.norm)]
             
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
  numCells <- ncol(r.maf)## original name error --Rongting
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
    'K' = ncol(r.maf),  # how many cells # name error Rongting
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
