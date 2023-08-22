

#' Load Gene Ontology genesets
#'
#' Load Gene Ontology genesets
#'
#' @param onto array of categories to load
#' @param to convert gene names to this type using \code{EnrichmentBrowser::idMap()}.  See \code{EnrichmentBrowser::idTypes(org="hsa")} for valid types
#' @param includeOffspring if TRUE, follow the GO hierarchy down and include all genes in offspring sets for a given gene set
#' @param org organism.  human (\code{'hsa'}), mouse (\code{'mmu'}), etc
#'
#' @details
#' This function loads the GO gene sets using the packages \code{EnrichmentBrowser} and \code{GO.db}  It can take a mintute to load because converting gene name type is slow.   
#'
#' @return Gene sets stored as GeneSetCollection
#' @examples
#' # load GO Biological Process
#' # gs = get_GeneOntology('BP')
#' 
#' # load all gene sets
#' # gs = get_GeneOntology()
#' 
#' @import EnrichmentBrowser GSEABase
#' @export
get_GeneOntology = function( onto = c("BP", "MF", "CC"), to = 'ENSEMBL', includeOffspring=TRUE, org="hsa" ){

	gs.list = lapply( onto, function(x){
		getGenesets(org=org, db="go", gene.id.type = to, return.type='GeneSetCollection', onto=x, hierarchical=includeOffspring)
		})

	# combine gene sets and convert to GeneSetCollection
	gs.go = GeneSetCollection( do.call(c,gs.list) )

	# now superceded by hierarchical=TRUE
	# if( includeOffspring ){
	# 	# follow the GO hierarchy down and include all genes in offspring sets for a given gene set
	# 	gs.go = aggregate_GO_offspring( gs.go )
	# }

	# summary(sapply(geneIds(gs.go), length))

	# modify gene set names to include description
	gs.go_rename = lapply( gs.go, function(gs){
	   setName(gs) = paste0( setName(gs), ': ', description(gs))
	   gs
	  } )
	GeneSetCollection( gs.go_rename )
}



