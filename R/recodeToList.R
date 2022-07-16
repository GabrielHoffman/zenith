# #' Recode GeneSetCollection to list used by limma
# #'
# #' Recode \code{GeneSetCollection} to \code{list} used by \code{limma}
# #'
# #' @param gsc GeneSetCollection object
# #'
# #' @return list storing data from GeneSetCollection
# #' @examples
# #' 
# #' # load Hallmark gene sets
# #' gs = get_MSigDB('H')
# #' 
# #' # recode GeneSetCollection as a list
# #' gs.list = recodeToList(gs)
# #'
# #' @rdname recodeToList
# #' @export
# setGeneric("recodeToList", function(gsc)
# 	standardGeneric("recodeToList"))

# #' Recode GeneSetCollection to list used by limma
# #'
# #' Recode GeneSetCollection to list used by limma
# #'
# #' @param gsc GeneSetCollection object
# #'
# #' @return list storing data from GeneSetCollection
# #' @rdname recodeToList-GeneSetCollection
# #' @importFrom GSEABase geneIds 
# #' @export
# setMethod("recodeToList", c("GeneSetCollection"),
#   function( gsc ){

#   	# convert to list  
# 	gsList = lapply( gsc, geneIds) 
# 	names(gsList) = names(gsc)

# 	gsList
# })


