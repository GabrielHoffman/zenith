
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
#' @importFrom dplyr `%>%`
#' @export
get_MSigDB = function(cat, to = c("ENSEMBL", "SYMBOL", "ENTREZ"), organism = c("HS", "MM")){

	to <- match.arg(to)
	organism <- match.arg(organism)

	cats = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "H")
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

	requireNamespace("msigdbdf", quietly=TRUE)
	if( ! isNamespaceLoaded("msigdbdf") ){
		txt = "Please run the following command to install the 'msigdbdf' package:\ninstall.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')"
		stop(txt)
	}

	# pass R CMD check
	gs_name <- gs_collection <- NULL

	# filter by category
	df1 <- msigdbdf::msigdbdf(organism) %>%
			dplyr::filter(gs_collection %in% cat) 

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
			collectionType = BroadCollection(tolower(cat), df2$gs_subcollection[1]))
		})
	# combine gene sets into GeneSetCollection
	GeneSetCollection( do.call(c, gs.list ) )
}











