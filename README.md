# Zenith

Perform gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package implements the [camera](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/camera) method from the [limma](https://bioconductor.org/packages/limma/) package proposed by [Wu and Smyth (2012)](https://doi.org/10.1093/nar/gks461).  `zenith()` is a simple extension of `camera()` to be compatible with linear mixed models implemented in `dream()`.


### Note
A future version will allow identification of gene sets with log fold changes with mixed sign. But this is not currently supported.

## Install

```r
library(devtools)

# install EnrichmentBrowser >= v2.22.0
install_github("lgeistlinger/EnrichmentBrowser")

# install variancePartition >= v1.19.16
install_github("GabrielHoffman/variancePartition")

install_github("GabrielHoffman/zenith")
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


<sub><sub>Updated January 11, 2021</sub></sub>

