# Zenith

Perform gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package implements the [camera](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/camera) method from the [limma](https://bioconductor.org/packages/limma/) package proposed by [Wu and Smyth (2012)](https://doi.org/10.1093/nar/gks461).  `zenith()` is a simple extension of `camera()` to be compatible with linear mixed models implemented in `dream()`.


### Note
A future version will allow identification of gene sets with log fold changes with mixed sign. But this is not currently supported.

### Usage
For compatability with `limma::camera`, use 
```r
zenith(fit, coef, index, use.ranks=FALSE, inter.gene.cor=0.01)
```

## Install
```r
# Install package and dependencies
devtools::install_github("GabrielHoffman/zenith")
```

### Dependencies
In case code aboves doesn't install these automatically
```r
# install EnrichmentBrowser >= v2.22.0
devtools::install_github("lgeistlinger/EnrichmentBrowser")

# install variancePartition >= v1.19.16
devtools::install_github("GabrielHoffman/variancePartition")
```

### Issues
If you encounter the error:

```
Error in BroadCollection(category = tolower(cat), subCategory = tolower(subcat)) :
	invalid BroadCollection category: ‘c8’
```

This is due to a recent addition of [cell type signature gene sets](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C8) in [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) v7.2.  In order to be compatible with this new addition, you must update the R package `GSEABase`: 

```r
remotes::install_github("Bioconductor/GSEABase")
```
In the future, this will be installed by default


<sub><sub>Nov 1, 2021</sub></sub>

