

# Zenith
## [ACTIVE DEVELOPMENT, USE WITH CAUTION]


Perform gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package is a slight modification of [limma::camera](https://doi.org/10.1093/nar/gks461) to 1) be compatible with dream, and 2) allow identification of gene sets with log fold changes with mixed sign.

## See [vignette](https://hoffmg01.u.hpc.mssm.edu/software/zenith/geuvadis.html)

## Install

```r
library(devtools)

# install EnrichmentBrowser >= v2.19.13
install_github("lgeistlinger/EnrichmentBrowser")

# install variancePartition >= v1.19.16
install_github("GabrielHoffman/variancePartition")

install_github("GabrielHoffman/zenith")
```
