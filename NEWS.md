# zenith 1.5.3
 - March 8, 2024
 - add `coef` to `zenithPR_gsa()`

 # zenith 1.5.2
 - Feb 2, 2024
 - resolve issue with renameing gene sets when multiple coefs are used

# zenith 1.5.1
 - throw error when zenith is run on multiple coefs

# zenith 1.3.1
 - user can specify organism

# zenith 1.1.2
 - April 20, 2023
 - fix warnings for limma

# zenith 1.1.1
 - fix warnings

# zenith 1.0.2
 - add `zenithPR_gsa()`

# zenith 1.0.1
 - fixes to Bioconductor

# zenith 0.99.13
 - fix bug in degrees of freedom for `zenith()` null

# zenith 0.99.12
 - update dependencies
    
# zenith 0.99.10
 - replace `recodeToList()` since this functionality already exists
    
# zenith 0.99.9
 - Update to address issues from https://github.com/Bioconductor/Contributions/issues/2630#issuecomment-1171369001

# zenith 0.99.8
 - `plotZenithResults()` as default

# zenith 0.99.7
 - update `plotZenithResults_gg()`

# zenith 0.99.6
 - add `plotZenithResults_gg()`

# zenith 0.99.2
 - fix documentation
 - add `zenith_gsa()`

# zenith 0.99.0
 - submit to BioC

# zenith 1.0.7
 - update docs

# zenith 1.0.6
 - in `zenith()` set `inter.gene.cor=0.01` to be default to be consistent with `limma::camera`

# zenith 1.0.5
 - bug fix in `zenith()` when `progressbar=FALSE`

# zenith 1.0.3
 - fix issue with `corInGeneSet()` when some residuals are NA

# zenith 1.0.2
 - flag to disable correlation in zenith()

# zenith 1.0.1
 - improve documentation
 - get_GeneOntology() uses getGenesets(...,hierarchical=TRUE)

# zenith 1.0.0
 - Improve documetation
 - fix bug with residuals

# zenith 0.99.10
 - Upgrade t- MSigDB 7.2

# zenith 0.99.10
 - Update documentation

# zenith 0.99.9
 - zenith report error when squaredStats is TRUE

# zenith 0.99.8
 - Much faster loading of cached gene sets
 - enforce: 
  - variancePartition (>= 1.19.16)
  - EnrichmentBrowser (>= 2.19.13)

# zenith 0.99.7
 - Enforce subsetting of residuals from MArrayLM

# zenith 0.99.6
 - zenith() now returns effect size and standard error
 - if squaredStats==TRUE, return PValue = p.greater

# zenith 0.99.5
 - Fix small error caused by fit$method

# zenith 0.99.4
 - zenith uses t statistical for fixed effect models and z.std for mixed effects

# zenith 0.99.3
 - zenith now works when dream uses random effects
 - more accurate time in progress bar