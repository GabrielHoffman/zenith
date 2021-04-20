

#' Evaluate mean correlation between residuals in gene set
#'
#' Evaluate mean correlation between residuals in gene set based on results from dream
#'
#' @param fit result of differential expression with dream
#' @param idx indeces or rownames to extract
#' @param squareCorr compute the mean squared correlation instead
#'
#' @importFrom Rfast cora
corInGeneSet <- function( fit, idx, squareCorr = FALSE ){

  if( is.null(fit$residuals) ){
    stop("fit must be result of dream(..., computeResiduals=TRUE)")
  }

  # Extract residuals 
  resid = fit$residuals[idx,,drop=FALSE]

  m = nrow(resid)

  # evaluate correlation between residuals
  # This is the correlation between test statistics (almost, consider eBayes??)
  # C = cor(t(resid))
  C = cora(t(resid)) # this is a faster version of correlation in Rfast

  if( squareCorr ){
    # Correlation between squared test statistics
    # note that the covariance is 2*C^2
    C = C^2
  }

  # evaluate mean correlation
  correlation = mean(C[lower.tri(C)])

  vif = 1 + correlation*(m-1) 

  # if gene set has only 1 gene, correlation will be NAN
  # and VIF should be set to 1
  if( is.nan(correlation) ){
    correlation = NA
    vif = 1
  }

  list(vif=vif, correlation=correlation)
}




