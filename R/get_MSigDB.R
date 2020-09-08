
#' Load MSigDB genesets
#'
#' Load MSigDB genesets
#'
#' @param cat array of categories to load
#' @param to convert gene names to this type using EnrichmentBrowser::idMap().  See EnrichmentBrowser::idTypes(org="hsa") for valid types
#'
#' @details
#' This function loads the MSigDB gene sets using the packages EnrichmentBrowser and msigdbr.  It can take a mintute to load because converting gene name type is slow.   
#'
#' @return Gene sets stored as GeneSetCollection
#' @examples
#' # load Hallmark gene sets
#' gs = get_MSigDB('H')
#' 
#' # load all gene sets
#' # gs = get_MSigDB()
#' 
#' @import EnrichmentBrowser GSEABase msigdbr
#' @export
get_MSigDB = function(cat = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7"), to = "ENSEMBL"){

	gs.list = lapply( cat, function(x){
		gs <- getGenesets(org="hsa", db="msigdb", return.type='GeneSetCollection', cat=x, subcat=NA)

		# Convert gene identifiers from the default Entrez to Ensembl
		idMap(gs, org = "hsa", from = "ENTREZID", to = to) 
	})

	# combine gene sets and convert to GeneSetCollection
	GeneSetCollection( do.call(c, gs.list ) )
}





