library(RUnit)

test_corInGeneSet = function(){

	library(variancePartition)

	set.seed(1)

	y <- matrix(rnorm(1000*6),1000,6)
	colnames(y) = paste0("sample_", 1:6)
	rownames(y) = paste0("gene_", 1:1000)

	design <- data.frame(Intercept=1,Group=c(0,0,0,1,1,1))
	rownames(design) = paste0("sample_", 1:6)

	# First set of 20 genes are genuinely differentially expressed
	index1 <- 1:20
	y[index1,4:6] <- y[index1,4:6]+1

	# Second set of 20 genes are not DE
	index2 <- 21:40

	# gene indeces
	lst = list(set1=index1,set2=index2)

	# differential expression with dream	
	fit = dream(y, ~ Group, design, computeResiduals=TRUE)
	fit = eBayes(fit)

	# Calculate correlation between residuals
	res1 = limma::interGeneCorrelation(y[index2,], design)

	res2 = zenith:::corInGeneSet(fit, rownames(fit)[index2])

	checkEquals(res1, res2)
}


test_zenith = function(){

	library(variancePartition)

	set.seed(1)

	y <- matrix(rnorm(1000*6),1000,6)
	colnames(y) = paste0("sample_", 1:6)
	rownames(y) = paste0("gene_", 1:1000)

	design <- data.frame(Intercept=1,Group=c(0,0,0,1,1,1))
	rownames(design) = paste0("sample_", 1:6)

	# First set of 20 genes are genuinely differentially expressed
	index1 <- 1:20
	y[index1,4:6] <- y[index1,4:6]+1

	# Second set of 20 genes are not DE
	index2 <- 21:40

	# gene indeces
	lst = list(set1=index1,set2=index2)

	# differential expression with dream	
	fit = dream(y, ~ Group, design, computeResiduals=TRUE)
	fit = eBayes(fit)

	# Parametric test
	res1 = camera(y, lst, design, inter.gene.cor=NA)
	res2 = zenith( fit, "Group", lst)

	# only check columns in common
	cols = colnames(res1)[1:5]

	test1 = checkEquals(res1[,cols], res2[,cols])

	# Non-parametric test
	res1 = camera(y, lst, design, inter.gene.cor=NA, use.ranks=TRUE)
	res2 = zenith( fit, "Group", lst, use.ranks=TRUE)

	test2 = checkEquals(res1[,cols], res2[,cols])

	test1 & test2
}


