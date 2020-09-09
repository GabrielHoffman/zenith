
 # Zenith

 Perform gene set analysis on the result of differential expression using linear (mixed) modeling with [dream](https://doi.org/10.1093/bioinformatics/btaa687) by considering the correlation between gene expression traits.  This package is a slight modification of [limma::camera](https://doi.org/10.1093/nar/gks461) to 1) be compatible with dream, and 2) allow identification of gene sets with log fold changes with mixed sign.

## Install

````
git clone https://github.com/GabrielHoffman/zenith.git
R CMD INSTALL zenith
```