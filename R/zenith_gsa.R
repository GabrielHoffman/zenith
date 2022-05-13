# Gabriel Hoffman
# May 13, 2022
# 
# Apply zenith to dream result


#' Perform gene set analysis using zenith
#' 
#' Perform a competitive gene set analysis accounting for correlation between genes.
#'
#' @param fit results from \code{dream()}
#' @param geneSets \code{GeneSetCollection} 
#' @param coefs coefficients to test using \code{topTable(fit, coef=coefs[i])}
#' @param use.ranks do a rank-based test \code{TRUE} or a parametric test \code{FALSE}? default: FALSE
#' @param n_genes_min minumum number of genes in a geneset
#' @param inter.gene.cor if NA, estimate correlation from data.  Otherwise, use specified value
#' @param progressbar if TRUE, show progress bar
#' @param ... other arguments
#'  
#' @return \code{data.frame} of results for each gene set and cell type 
#'
#' @details This code adapts the widely used \code{camera()} analysis \insertCite{wu2012camera}{zenith} in the \code{limma} package \insertCite{ritchie2015limma}{zenith} to the case of linear (mixed) models used by \code{variancePartition::dream()}.
#'
#' @references{
#' \insertAllCited{}
#' }
#' @rdname zenith_gsa-methods
#' @export
#' @seealso \code{limma::camera}
setGeneric('zenith_gsa', function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){
	standardGeneric("zenith_gsa")
	})





#' @examples
#' # Load packages
#' library(edgeR)
#' library(variancePartition)
#' library(tweeDEseqCountData)
#' 
#' # Load RNA-seq data from LCL's
#' data(pickrell)
#' 
#' # Filter genes
#' # Note this is low coverage data, so just use as code example
#' dsgn = model.matrix(~ gender, pData(pickrell.eset))
#' keep = filterByExpr(exprs(pickrell.eset), dsgn, min.count=5)
#' 
#' # Compute library size normalization
#' dge = DGEList(counts = exprs(pickrell.eset)[keep,])
#' dge = calcNormFactors(dge)
#' 
#' # Estimate precision weights using voom
#' vobj = voomWithDreamWeights(dge, ~ gender, pData(pickrell.eset))
#' 
#' # Apply dream analysis
#' fit = dream(vobj, ~ gender, pData(pickrell.eset))
#' fit = eBayes(fit)
#' 
#' # Load Gene Ontology database 
#' # use gene 'SYMBOL', or 'ENSEMBL' id
#' # use get_MSigDB() to load MSigDB
#' go.gs = get_GeneOntology("BP", to="ENSEMBL")
#'    
#' # Run zenith analysis
#' res.gsa = zenith_gsa(fit, go.gs, 'gendermale', progressbar=FALSE )
#' 
#' # Show top gene sets
#' head(res.gsa, 2)
#' 
#' # for each cell type select 3 genesets with largest t-statistic
#' # and 1 geneset with the lowest
#' # Grey boxes indicate the gene set could not be evaluted because
#' #    to few genes were represented
#' plotZenithResults(res.gsa)
#' 
#' @importFrom limma ids2indices
#' @importFrom stats p.adjust
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,MArrayLM,GeneSetCollection,ANY-method
#' @export
setMethod("zenith_gsa", signature(fit="MArrayLM", geneSets = "GeneSetCollection", coefs="ANY"),
	function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){

	# convert GeneSetCollection to list
	geneSets.lst = recodeToList( geneSets )

	# Map from genes to gene sets
	index = ids2indices( geneSets.lst, rownames(fit))
	   
	# filter by size of gene set
	index = index[vapply(index, length, FUN.VALUE=numeric(1)) >= n_genes_min]

	# for each coefficient selected
	df_zenith = lapply( coefs, function(coef){
		# run zenith on dream fits
		df_res = zenith(fit, coef, index, use.ranks=use.ranks, inter.gene.cor=inter.gene.cor, progressbar=progressbar)
		
		data.frame(coef = coef, Geneset = rownames(df_res), df_res)
	})
	df_zenith = do.call(rbind, df_zenith)
	
	df_zenith
})





