
#' Load MSigDB genesets
#'
#' Load MSigDB genesets
#'
#' @param cat array of categories to load.  
#' @param to return genes names as \code{'ENSEMBL'} or \code{'SYMBOL'}
#' @param organism organism: human (\code{'HS'}) or mouse (\code{'MS'})
#'
#' @details
#' This function loads the MSigDB gene sets using the packages  and \code{msigdbr}.  It can take a mintute to load because converting gene name type is slow.   
#'
#' @return Gene sets stored as GeneSetCollection
#' @examples
#' # load Hallmark gene sets
#' gs = get_MSigDB('H')
#' 
#' # load all gene sets
#' # gs = get_MSigDB()
#' 
#' @importFrom GSEABase GeneSet BroadCollection ENSEMBLIdentifier SymbolIdentifier EntrezIdentifier
#' @importFrom dplyr `%>%` bind_rows
#' @importFrom msigdbr msigdbr
#' @export
get_MSigDB = function(cat, to = c("ENSEMBL", "SYMBOL", "ENTREZ"), organism = c("HS", "MM")){

	to <- match.arg(to)
	organism <- match.arg(organism)

	species <- switch(organism, 
						HS = "Homo sapiens",
						MM = "Mus musculus")

	cats <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "H")
	if( missing(cat) ){
		cat <- cats
	}

	if( ! all(cat %in% cats) ){
		stop("cat is not valid")
	}

	switch(to, 
		"ENSEMBL" = {	
			geneKey <- "db_ensembl_gene"
			geneIdType <- ENSEMBLIdentifier()
			},
		"SYMBOL" = {			
			geneKey <- "db_gene_symbol"
			geneIdType <- SymbolIdentifier()
			},
		"ENTREZ" = {
			geneKey <- "db_ncbi_gene"
			geneIdType <- EntrezIdentifier()
			})		

	# pass R CMD check
	gs_name <- gs_collection <- NULL

	if( ! requireNamespace("msigdbdf", quietly = TRUE) ){
		txt <- "! Please run the following command to install the 'msigdbdf' package:\ninstall.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')"
		warning(txt)
	}

	# filter by category
	df1 <- lapply(cat, function(x){
		msigdbr(species = species, 
				db_species = organism,
				collection = x)		
		})
	df1 <- bind_rows(df1)

	# loop thru gene sets
	gs.list <- lapply( unique(df1$gs_name), function(ID){

		df2 <- df1 %>%
			dplyr::filter(gs_name == ID)

		# create gene set
		GeneSet(geneIds = unique(df2[[geneKey]]),
			setName = ID,
			shortDescription = df2$gs_description[1],
			organism = organism,
			geneIdType = geneIdType, 
			collectionType = BroadCollection(tolower(unique(df2$gs_collection)), df2$gs_subcollection[1]))
		})
	# combine gene sets into GeneSetCollection
	GeneSetCollection( do.call(c, gs.list ) )
}











