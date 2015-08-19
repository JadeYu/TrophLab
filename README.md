# TrophLab
Simulate food webs to study macroecological patterns (between biophysical constraints and trophic structure).

Each new species chooses a position in the food web that correponds to the highest competitive advantage (in terms of total resource acquired).
The competition model is based on maximizing resource allocation microstates (Zhang and Harte 2015;http://dx.doi.org/10.1016/j.tpb.2015.07.003).
Please refer to this paper for the derivation and definitions of the key parameters, e.g. relative individual distinguishability.

In these simulations relative individual distinguishability is assumed to be positively related to the species resource requirement (assumed to 
be proportional to species mean metabolic rate) to account for the effect of lower reproduction success due to higher within-species variation. 
The linearity of this relationship is regulated by the dispersion parameter: when dispersion -> 0, relative individual distinguishability is linearly 
related to resource requirement; the bigger the dispersion, the more extreme (very small or very big) relative individual distinguishability is for
any given resource requirement.

Generalization cost is assumed to be negatively related to relative individual distinguishability to account for the interpretation that species 
with higher within-species variation are better generalists. A similar dispersion parameter is used to regulate the linearity of this relationship.

For more detailed explanations of the underlying process and optimization algorithm, please refer to the pdf write-up:FWA_switch_cost.
