# March 8, 2023


#' Gene set analysis using pre-computed test statistic
#'
#' Perform gene set analysis on the result of a pre-computed test statistic.  Test whether statistics in a gene set are larger/smaller than statistics not in the set.
#'
#' @param statistics pre-computed test statistics
#' @param ids name of gene for each entry in \code{statistics}
#' @param geneSets \code{GeneSetCollection} 
#' @param use.ranks do a rank-based test \code{TRUE} or a parametric test \code{FALSE}? default: FALSE
#' @param n_genes_min minumum number of genes in a geneset
#' @param progressbar if TRUE, show progress bar
#' @param inter.gene.cor correlation of test statistics with in gene set
#'
#' @details
#' This is the same as \code{zenith_gsa()}, but uses pre-computed test statistics.  Note that \code{zenithPR_gsa()} may give slightly different results for small samples sizes, if \code{zenithPR_gsa()} is fed t-statistics instead of z-statistics.
#'
#' @return
#' \itemize{
#'   \item \code{NGenes}: number of genes in this set
#'   \item \code{Correlation}: mean correlation between expression of genes in this set
#'   \item \code{delta}: difference in mean t-statistic for genes in this set compared to genes not in this set
#'   \item \code{se}: standard error of \code{delta}
#'   \item \code{p.less}: p-value for hypothesis test of \code{H0: delta < 0}
#'   \item \code{p.greater}: p-value for hypothesis test of \code{H0: delta > 0}
#'   \item \code{PValue}:  p-value for hypothesis test \code{H0: delta != 0}
#'   \item \code{Direction}: direction of effect based on sign(delta)
#'   \item \code{FDR}: false discovery rate based on Benjamini-Hochberg method in \code{p.adjust}
#'   \item \code{coef.name}: name for pre-computed test statistics. Default: \code{zenithPR}
#' }
#'
#' @seealso \code{zenith_gsa()}
#' @importFrom Rdpack reprompt 
#' @importFrom stats runif
#'
#' @export
zenithPR_gsa = function(statistics, ids, geneSets, use.ranks = FALSE, n_genes_min = 10, progressbar=TRUE, inter.gene.cor = 0.01, coef.name = "zenithPR"){

	if( length(statistics) != length(ids) ){
		stop("statsitics and ids must be the same length")
	}

	if( any(duplicated(ids)) ){
		stop("All entries in ids must be unique")
	}

	if( any(is.null(ids) | is.na(ids)) ){
		stop("All entries in ids characters and not NULL or NA")
	}

	# Global statistics
	meanStat <- mean(statistics)
	varStat <- var(statistics)
	G = length( statistics )
	df.camera = G - 2

	allow.neg.cor = FALSE

	# convert GeneSetCollection to list
	geneSets.lst = geneIds( geneSets )

	# Map from genes to gene sets
	index = ids2indices( geneSets.lst, ids)
	   
	# filter by size of gene set
	index = index[vapply(index, length, FUN.VALUE=numeric(1)) >= n_genes_min]
	nsets = length(index)

	if( nsets == 0 ){
		stop("No sets contain genes matching entries in ids")
	}

	# set up progress bar
	pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta", total = nsets, width= 60, clear=FALSE)

	# Analysis for each gene set
	tab = lapply(seq_len(nsets), function(i){   

		if( progressbar & runif(1) < .01 ){
			pb$update(ratio = i / nsets )
		} 

	    iset <- index[[i]]
	    if(is.character(iset)) iset <- which(ids %in% iset)

	    StatInSet <- statistics[iset]

	    m <- length(StatInSet)
	    m2 <- G - m

	    # convert correlation to VIF
	    correlation = inter.gene.cor
	    vif = 1 + inter.gene.cor * (m - 1) 

	    if(use.ranks) {

	      corr.use = correlation
	      if( ! allow.neg.cor ) corr.use <- max(0,corr.use)

	      res = .rankSumTestWithCorrelation(iset, statistics=statistics, correlation=corr.use, df=df.camera)

	      df = data.frame(  NGenes      = m,
	                        Correlation = correlation,
	                        delta       = res$effect,
	                        se          = res$se,
	                        p.less      = res$less,
	                        p.greater   = res$greater )
	    }else{ 

	      if( ! allow.neg.cor ) vif <- max(1,vif)

	      meanStatInSet <- mean(StatInSet)
	      delta <- G/m2*(meanStatInSet-meanStat)
	      varStatPooled <- ( (G-1)*varStat - delta^2*m*m2/G ) / (G-2)
	      delta.se = sqrt( varStatPooled * (vif/m + 1/m2) )
	      two.sample.t = delta / delta.se

	      df = data.frame(NGenes      = m,
	                      Correlation = correlation,
	                      delta       = delta,
	                      se          = delta.se,
	                      p.less      = pt(two.sample.t,df=df.camera),
	                      p.greater   = pt(two.sample.t,df=df.camera,lower.tail=FALSE) )
	    }
	    df 
	  })
	tab = do.call(rbind, tab)
	rownames(tab) <- names(index)

	if( progressbar & ! pb$finished ) pb$update(ratio = 1)
	pb$terminate()

  	# Post-process results
	tab$PValue = 2*pmin(tab$p.less, tab$p.greater)
	tab$Direction = ifelse(tab$p.less < tab$p.greater, "Down", "Up")
	tab$FDR = p.adjust(tab$PValue, "BH")

	if( progressbar & ! pb$finished ) pb$update( 1.0 )
	pb$terminate() 

	# Sort by p-value
	o <- order(tab$PValue)
	tab <- tab[o,]

	# Make results compatible with plotZenithResults
	tab$Geneset <- rownames(tab)
	tab$coef <- coef.name

	tab
}
