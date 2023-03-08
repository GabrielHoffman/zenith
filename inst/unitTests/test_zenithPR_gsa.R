
library(RUnit)

# test that zenith_gsa() and zenithPR_gsa() give
# equivalent results
test_zenithPR_gsa = function(){

	library(edgeR)
	library(variancePartition)
	library(tweeDEseqCountData)
	library(zenith)

	# Load RNA-seq data from LCL's
	data(pickrell)
	geneCounts = exprs(pickrell.eset)
	df_metadata = pData(pickrell.eset)

	# Filter genes
	# Note this is low coverage data, so just use as code example
	dsgn = model.matrix(~ gender, df_metadata)
	keep = filterByExpr(geneCounts, dsgn, min.count=5)

	# Compute library size normalization
	dge = DGEList(counts = geneCounts[keep,])
	dge = calcNormFactors(dge)

	# Estimate precision weights using voom
	vobj = voomWithDreamWeights(dge, ~ gender, df_metadata)

	# Apply dream analysis
	fit = dream(vobj, ~ gender, df_metadata)
	fit = eBayes(fit)

	gs = get_MSigDB("H", to="ENSEMBL")

	res.gsa = zenith_gsa(fit, gs, 'gendermale', progressbar=FALSE )

	tab = topTable(fit, coef="gendermale", number=Inf )

	stat = zscoreT( tab$t, df=fit$df.residual[1], approx=TRUE, method="hill")
	df_gsa = zenithPR_gsa( stat, rownames(tab), gs, progressbar=FALSE)

	i = match(rownames(df_gsa), rownames(res.gsa))

	checkEqualsNumeric(res.gsa$delta[i], df_gsa$delta, tol=1e-2)
}






