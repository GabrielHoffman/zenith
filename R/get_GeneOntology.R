

#' Load Gene Ontology genesets
#'
#' Load Gene Ontology genesets
#'
#' @param onto array of categories to load
#' @param to convert gene names to this type using \code{EnrichmentBrowser::idMap()}.  See \code{EnrichmentBrowser::idTypes(org="hsa")} for valid types
#' @param includeOffspring if TRUE, follow the GO hierarchy down and include all genes in offspring sets for a given gene set
#'
#' @details
#' This function loads the GO gene sets using the packages \code{EnrichmentBrowser} and \code{GO.db}  It can take a mintute to load because converting gene name type is slow.   
#'
#' @return Gene sets stored as GeneSetCollection
#' @examples
#' # load GO Biological Process
#' gs = get_GeneOntology('BP')
#' 
#' # load all gene sets
#' # gs = get_GeneOntology()
#' 
#' @import EnrichmentBrowser GSEABase
#' @export
get_GeneOntology = function( onto = c("BP", "MF", "CC"), to = 'ENSEMBL', includeOffspring=TRUE ){

	gs.list = lapply( onto, function(x){
		getGenesets(org="hsa", db="go", gene.id.type = to, return.type='GeneSetCollection', onto=x)
		})

	# combine gene sets and convert to GeneSetCollection
	gs.go = GeneSetCollection( do.call(c,gs.list) )

	if( includeOffspring ){
		# follow the GO hierarchy down and include all genes in offspring sets for a given gene set
		gs.go = aggregate_GO_offspring( gs.go )
	}

	# summary(sapply(geneIds(gs.go), length))

	# modify gene set names to include description
	gs.go_rename = lapply( gs.go, function(gs){
	   setName(gs) = paste0( setName(gs), ': ', description(gs))
	   gs
	  } )
	GeneSetCollection( gs.go_rename )
}


#' Aggregate GO offspring genesets and genes
#'
#' Follow the GO hierarchy down and include all genes in offspring sets for a given gene set
#'
#' @param gs.GO GeneSetCollection from getGenesets()
#'
#' @return GeneSetCollection with aggregated gene sets
#'
#' @import EnrichmentBrowser GSEABase GO.db
#' @importFrom data.table data.table setkey
#'
aggregate_GO_offspring = function( gs.GO ){

	if( ! is(gs.GO, 'GeneSetCollection') ){
		stop("gsGO must be a GeneSetCollection")
	}

	# create a data.table of gene set names vs index for fast searching
	dt = data.table(GOID = gsub("GO", "GO:", names(gs.GO)), idx = seq_len(length(gs.GO)))
	setkey(dt, 'GOID')

	# summary(sapply(geneIds(gs.GO), length))

	# get combined of offspring for each gene set
	df_all = c(as.list(GOMFOFFSPRING), as.list(GOBPOFFSPRING), as.list(GOCCOFFSPRING))

	# for each gene set in GO, extract genes in this set and all offspring sets
	gs.GO.all_offspring = lapply(gs.GO, function(x){
		
		# get GO ID for this gene set
		goid = ids(collectionType(x))

		# get GO ID's or this and offspring gene sets
		query = c(goid, df_all[[goid]])
		query = query[!is.na(query)]

		# get the index through a fast look up in dt
		idx = dt[query,]$idx
		idx = idx[!is.na(idx)]

		gs.offspring = gs.GO[idx]

		# get gene array to current and offspring genes
		geneIds(x) = unique(unlist(geneIds(gs.offspring)))

		x
		})
	GeneSetCollection(gs.GO.all_offspring)
}




