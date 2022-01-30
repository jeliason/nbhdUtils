* Spatial statistics
	* Moran's I, Geary's C, G (local or not), ELSA (elsa)
	* Joint count (spdep)
	* Biodiversity analyses (functional summaries), e.g. Shannon index, etc. (spatialsegregation)
	* Ripley's K,L,pcf for determining interaction ranges
	* Spatial point process models (Gibbs, area-interaction, Geyer saturation) for interaction effect size and direction (joint or conditional on location)
	* mobsim - spatial simulation of communities
	* SAR,CAR models for tumor cell # prediction in neighborhoods


Neighborhood/niche questions:
* Can we predict the number of tumor cells in a given neighborhood, based on the counts of other cells (and maybe their cell markers?)?
* Or maybe just predict whether tumor cells will be present?
* Or whether we can predict PCs or some other low-dim manifold of tumor cell markers per neighborhood, based on loadings in that neighborhood?
* 

* Predict density of tumor cells in neighborhoods, based on mean of several different markers and the interaction with different immune cell types