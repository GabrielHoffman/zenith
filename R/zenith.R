

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
#' @param progressbar if TRUE, show progress bar
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
#' @import limma stats utils methods progress
#'
#' @export
zenith <- function( fit, coef, index, use.ranks=FALSE, allow.neg.cor=FALSE, squaredStats=FALSE, progressbar=TRUE ){

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

  # Only keep residuals for genes present in the main part of fit
  # Currently fit[1:10,] subsets objects in a standard MArrayLM
  #   but residuals is not standard, so it is not subsetted
  fit$residuals = fit$residuals[rownames(fit),,drop=FALSE]
  
  if( fit$method == "ls" ){

    # extract test statistics
    Stat = topTable(fit, coef, number=Inf, sort.by="none")$t

    if( ! use.ranks ){
      df = fit$df.total[1]
      Stat <- zscoreT( Stat, df=df, approx=TRUE, method="hill")
    }
  }else if( fit$method == "lmer"){
    # extract test statistics
    Stat = topTable(fit, coef, number=Inf, sort.by="none")$z.std
  }else{
    stop("Model method must be either 'ls' or 'lmer'")
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

  # setup progressbar
  if( progressbar ){
    # since time is quadratic in the size of the gene set
    total_work = sum(sapply(index, length)^2)

    pb <- progress_bar$new(
      format = " [:bar] :percent eta: :eta",
      clear = FALSE,
      total = total_work, width = 60)
  }
  cumulative_work = cumsum(sapply(index, length)^2)

  # tab <- matrix(0,nsets,5)
  # rownames(tab) <- names(index)
  # colnames(tab) <- c("NGenes","Correlation","Down","Up","TwoSided")
  tab = lapply(seq_len(nsets), function(i){    

    iset <- index[[i]]
    if(is.character(iset)) iset <- which(ID %in% iset)

    StatInSet <- Stat[iset]

    m <- length(StatInSet)
    m2 <- G-m

    # cumulative_work <<- cumulative_work + m^2

    # Compute correlation within geneset
    res = corInGeneSet( fit, iset, squaredStats)
    correlation = res$correlation
    vif = res$vif

    # tab[i,1] <- m
    # tab[i,2] <- correlation

    if(use.ranks) {

      corr.use = correlation
      if( ! allow.neg.cor ) corr.use <- max(0,corr.use)

      res = .rankSumTestWithCorrelation(iset, statistics=Stat, correlation=corr.use, df=df.camera)

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

      # two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
      # tab[i,3] <- pt(two.sample.t,df=df.camera)
      # tab[i,4] <- pt(two.sample.t,df=df.camera,lower.tail=FALSE)

      df = data.frame(NGenes      = m,
                      Correlation = correlation,
                      delta       = delta,
                      se          = delta.se,
                      p.less      = pt(two.sample.t,df=df.camera),
                      p.greater   = pt(two.sample.t,df=df.camera,lower.tail=FALSE) )
    }

    if( progressbar & (i %% 100 == 0) ) pb$update( cumulative_work[i] / total_work )

    df 
  })
  tab = do.call(rbind, tab)
  rownames(tab) <- names(index)

  # p-value for two sided test
  if( squaredStats ){
    tab$PValue = tab$p.greater
    tab$Direction = "Up"
  }else{    
    tab$PValue = 2*pmin(tab$p.less, tab$p.greater)
    tab$Direction = ifelse(tab$p.less < tab$p.greater, "Down", "Up")
  }

  tab$FDR = p.adjust(tab$PValue, "BH")

  if( progressbar ){
    pb$update( 1.0 )
    pb$terminate() 
  }

  # tab[,5] <- 2*pmin(tab[,3],tab[,4])

  # # New column names (Jan 2013)
  # tab <- data.frame(tab,stringsAsFactors=FALSE)
  # Direction <- rep_len("Up",length.out=nsets)
  # Direction[tab$Down < tab$Up] <- "Down"
  # tab$Direction <- Direction
  # tab$PValue <- tab$TwoSided
  # tab$Down <- tab$Up <- tab$TwoSided <- NULL

  # # Add FDR
  # if(nsets>1) tab$FDR <- p.adjust(tab$PValue,method="BH")

  # Sort by p-value
  o <- order(tab$PValue)
  tab <- tab[o,]

  tab
}




#' Two Sample Wilcoxon-Mann-Whitney Rank Sum Test Allowing For Correlation
#' 
#' Same as \code{limma::.rankSumTestWithCorrelation}, but returns effect size.
#'
#' @param index any index vector such that 'statistics[index]' contains the values of the statistic for the test group.
#' @param statistics numeric vector giving values of the test statistic.
#' @param correlation numeric scalar, average correlation between cases in the test group.  Cases in the second group are assumed independent of each other and other the first group.
#' @param df degrees of freedom which the correlation has been estimated.
#'
#' @details See \code{limma::.rankSumTestWithCorrelation}
#'
.rankSumTestWithCorrelation = function (index, statistics, correlation = 0, df = Inf){
    n <- length(statistics)
    r <- rank(statistics)
    r1 <- r[index]
    n1 <- length(r1)
    n2 <- n - n1
    U <- n1 * n2 + n1 * (n1 + 1)/2 - sum(r1)
    mu <- n1 * n2/2
    if (correlation == 0 || n1 == 1) {
        sigma2 <- n1 * n2 * (n + 1)/12
    }
    else {
        sigma2 <- asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 - 
            1) + asin(correlation/2) * n1 * (n1 - 1) * n2 * (n2 - 
            1) + asin((correlation + 1)/2) * n1 * (n1 - 1) * 
            n2
        sigma2 <- sigma2/2/pi
    }
    TIES <- (length(r) != length(unique(r)))
    if (TIES) {
        NTIES <- table(r)
        adjustment <- sum(NTIES * (NTIES + 1) * (NTIES - 1))/(n * 
            (n + 1) * (n - 1))
        sigma2 <- sigma2 * (1 - adjustment)
    }
    zlowertail <- (U + 0.5 - mu)/sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu)/sqrt(sigma2)
    
    data.frame( effect  = -1*(U - mu), 
                se      = sqrt(sigma2),
                less    = pt(zuppertail, df = df, lower.tail = FALSE), 
                greater = pt(zlowertail, df = df))
}



















