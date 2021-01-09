# Zenith

Perform gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package implements the [camera](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/camera) method from the [limma](https://bioconductor.org/packages/limma/) package proposed by [Wu and Smyth (2012)](https://doi.org/10.1093/nar/gks461).  `zenith()` is a simple extension of `camera()` to be compatible with linear mixed models implemented in `dream()`.

See extensive [documentation](https://hoffmg01.u.hpc.mssm.edu/software/zenith/zenith-manual.pdf) and [vignette](https://hoffmg01.u.hpc.mssm.edu/software/zenith/geuvadis.html).


### Note
A future version will allow identification of gene sets with log fold changes with mixed sign. But this is not currently supported.

## Install

```r
library(devtools)

# install EnrichmentBrowser >= v2.19.13
install_github("lgeistlinger/EnrichmentBrowser")

# install variancePartition >= v1.19.16
install_github("GabrielHoffman/variancePartition")

install_github("GabrielHoffman/zenith")
```


<sub><sub>Updated January 8, 2021</sub></sub>

