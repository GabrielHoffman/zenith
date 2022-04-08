# Zenith

Perform gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package implements the [camera](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/camera) method from the [limma](https://bioconductor.org/packages/limma/) package proposed by [Wu and Smyth (2012)](https://doi.org/10.1093/nar/gks461).  `zenith()` is a simple extension of `limma::camera()` to be compatible with linear mixed models implemented in `variancePartition::dream()`.


## Install
```r
# Install package and dependencies
devtools::install_github("DiseaseNeuroGenomics/zenith")
```

## Supported gene set databases
[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser/) provides access to a range of gene set databases.  Here are some shortcuts to load common databases:  

```r
library(zenith)

# MSigDB
gs.msigdb = get_MSigDB()

# Gene Ontology
gs.go = get_GeneOntology()
```

For more detailed control including alternate gene identifiers (i.e. ENSEMBL, ENTREZ) or species (i.e. hsa, mmu)

```r
library(EnrichmentBrowser)

# KEGG
gs.kegg = getGenesets(  org = "hsa", 
	   					db = "kegg", 
	   					gene.id.type = "ENSEMBL", 
	   					return.type = "GeneSetCollection")

## ENRICHR resource
df = showAvailableCollections( org = "hsa", db = "enrichr")

# Allen_Brain_Atlas_10x_scRNA_2021
gs.allen = getGenesets(  org = "hsa", 
	   					db = "enrichr", 
	   					lib = "Allen_Brain_Atlas_10x_scRNA_2021",
	   					gene.id.type = "ENSEMBL", 
	   					return.type = "GeneSetCollection")
```



<sub><sub>April 1, 2022</sub></sub>

