

#' Gene set analysis following differential expression with dream
#'
#' Perform gene set analysis on the result of differential expression using linear (mixed) modeling with \code{variancePartition::dream} by considering the correlation between gene expression traits.  This package is a slight modification of \code{limma::camera} to 1) be compatible with dream, and 2) allow identification of gene sets with log fold changes with mixed sign.
#'
#' @param fit result of differential expression with dream
#' @param coef coefficient to test using topTable(fit, coef)
#' @param index an index vector or a list of index vectors.  Can be any vector such that 'fit[index,]'  selects the rows corresponding to the test set.  The list can be made using 'ids2indices'.
#' @param  use.ranks do a rank-based test ('TRUE') or a parametric test ('FALSE')?
#' @param allow.neg.cor should reduced variance inflation factors be allowed for negative correlations?
#' @param squaredStats Test squared test statstics to identify gene sets with log fold change of mixed sign.
#'
#' @details
#' \code{zenith} gives the same results as `\code{camera(..., inter.gene.cor=NA)} which estimates the correlation with each gene set.
#'
#' For differential expression with dream using linear (mixed) models see Hoffman and Roussos (2020).  For the original camera gene set test see Wu and Smyth (2012).
#' 
#' @references{
#'   \insertRef{hoffman2020dream}{zenith}
#' 
#'   \insertRef{wu2012camera}{zenith}
#' }
#' @importFrom Rdpack reprompt
#' @import limma stats utils methods
#'
#' @export
zenith <- function( fit, coef, index, use.ranks=FALSE, allow.neg.cor=FALSE, squaredStats=FALSE ){

  if( ! is(fit, 'MArrayLM') ){
    stop("fit must be of class MArrayLM from variancePartition::dream")
  }

  if( is.null(fit$residuals) ){
    stop("fit must be result of dream(..., computeResiduals=TRUE)")
  }

  if( ! (coef %in% colnames(coef(fit))) ){
    stop("coef must be in colnames(coef(fit))")
  }

  # Check index
  if(!is.list(index)) index <- list(set1=index)
  nsets <- length(index)
  if(nsets==0L) stop("index is empty")
  
  # extract test statistics
  Stat = topTable(fit, coef, number=Inf, sort.by="none")$t

  if( ! use.ranks ){
    Stat <- zscoreT( Stat, df= fit$df.total[1], approx=TRUE, method="hill")
  }

  # use squared test statistics
  if( squaredStats ){
    Stat = Stat^2
  }

  # get number of statistics
  G = length( Stat )
  ID = rownames(fit)
  
  df.camera <- min(fit$df.residual[1], G - 2L)

  # Global statistics
  meanStat <- mean(Stat)
  varStat <- var(Stat)

  tab <- matrix(0,nsets,5)
  rownames(tab) <- names(index)
  colnames(tab) <- c("NGenes","Correlation","Down","Up","TwoSided")
  for (i in 1:nsets) {

    iset <- index[[i]]
    if(is.character(iset)) iset <- which(ID %in% iset)

    StatInSet <- Stat[iset]

    m <- length(StatInSet)
    m2 <- G-m

    # Compute correlation within geneset
    res = corInGeneSet( fit, iset, squaredStats)
    correlation = res$correlation
    vif = res$vif

    tab[i,1] <- m
    tab[i,2] <- correlation

    if(use.ranks) {
      if(!allow.neg.cor) correlation <- max(0,correlation)
      tab[i,3:4] <- rankSumTestWithCorrelation(iset,statistics=Stat,correlation=correlation,df=df.camera)
    }else{  
      if(!allow.neg.cor) vif <- max(1,vif)
      meanStatInSet <- mean(StatInSet)
      delta <- G/m2*(meanStatInSet-meanStat)
      varStatPooled <- ( (G-1)*varStat - delta^2*m*m2/G ) / (G-2)
      two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
      tab[i,3] <- pt(two.sample.t,df=df.camera)
      tab[i,4] <- pt(two.sample.t,df=df.camera,lower.tail=FALSE)
    }
  }
  tab[,5] <- 2*pmin(tab[,3],tab[,4])

  # New column names (Jan 2013)
  tab <- data.frame(tab,stringsAsFactors=FALSE)
  Direction <- rep_len("Up",length.out=nsets)
  Direction[tab$Down < tab$Up] <- "Down"
  tab$Direction <- Direction
  tab$PValue <- tab$TwoSided
  tab$Down <- tab$Up <- tab$TwoSided <- NULL

  # Add FDR
  if(nsets>1) tab$FDR <- p.adjust(tab$PValue,method="BH")

  # Sort by p-value
  o <- order(tab$PValue)
  tab <- tab[o,]

  tab
}



