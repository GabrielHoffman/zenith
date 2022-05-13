
<br>

## Gene set testing after dream analysis

<div style="text-align: justify">
`zenith()` performs gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package implements the [camera](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/camera) method from the [limma](https://bioconductor.org/packages/limma/) package proposed by [Wu and Smyth (2012)](https://doi.org/10.1093/nar/gks461).  `zenith()` is a simple extension of `limma::camera()` to be compatible with linear mixed models implemented in `variancePartition::dream()`.
</div>

## Install
```r
# Install package and dependencies
devtools::install_github("DiseaseNeuroGenomics/zenith")
```

<sub><sub>May 13, 2022</sub></sub>

