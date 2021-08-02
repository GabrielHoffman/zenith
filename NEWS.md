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