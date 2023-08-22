
#' Load MSigDB genesets
#'
#' Load MSigDB genesets
#'
#' @param cat array of categories to load.  Defaults to array of all MSigDB categories
#' @param to convert gene names to this type using \code{EnrichmentBrowser::idMap()}.  See \code{EnrichmentBrowser::idTypes(org="hsa")} for valid types
#' @param org organism.  human ('hsa'), mouse ('mmu')
#'
#' @details
#' This function loads the MSigDB gene sets using the packages \code{EnrichmentBrowser} and \code{msigdbr}.  It can take a mintute to load because converting gene name type is slow.   
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
get_MSigDB = function(cat = unique(msigdbr_collections()$gs_cat), to = "ENSEMBL", org="hsa"){

	gs.list = lapply( cat, function(x){
		# get gene sets and convert to gene.id.type "to"
		getGenesets(org=org, db="msigdb", gene.id.type = to, return.type='GeneSetCollection', cat=x, subcat=NA)
	})

	# combine gene sets and convert to GeneSetCollection
	GeneSetCollection( do.call(c, gs.list ) )
}






