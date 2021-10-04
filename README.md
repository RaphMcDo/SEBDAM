## SEBDAM (Functional but still under construction)
- SEBDAM is an R package for fitting a Spatially-Explicit Biomass Dynamics Assessment Model in a state-space framework to fishery-independent data. Built for the assessment of scallop fisheries, it can easily be applied to any marine animals that exhibits spatially heterogeneous area utilization inside of a given stock.
- The SEBDAM package contains 2 models: SEBDAM, previously described, and TLM (Tow-Level Model), a strictly temporal biomass dynamics models that can be directly compared to SEBDAM.

## Installation
- Installation requires the installation of the *devtools* package in R and uses the function *devtools::install_github("RaphMcDo/SEBDAM")*. Other packages that need to be installed are *sf*, *dplyr*, *stringr*, *optimx*, *parallel* (for non-Windows users) and *INLA*. The first five packages are on CRAN, while INLA has to be installed following the instructions at https://www.r-inla.org/download-install.

## Background
- SEBDAM was designed to incorporate spatial variability and autocorrelations in fishery-independent observations directly into traditional biomass dynamics model and model the underyling processes (biomass, recruitment and natural mortality) as explicitly spatial processes.
- This approach builds on different spatial and temporal techniques, incorporating advances related to the use of Gaussian Markov Random Fields through the Stochastic Partial Differential Equations approach (Lindgren and Rue, 2011), the use of "knots" to create a predictive framework (Thorson et al., 2015), delta models (Martin et al., 2005), and of a delay-difference stock assessment model (Deriso, 1980, Schnute, 1985, Smith and Lundy 2002).
- The use of Gaussian Markov Random Fields assumes that spatial correlations decrease according to distance. This decorrelation based on distance can further account for geometric anisotropy (Lindgren et al., 2012) or the presence of land through a barrier model (Bakka et al., 2019).

# Resources

- *User Manual*: Please see **Under Construction** for the model equations in relations to the input and output as seen when fitting the model, along with all possible options that can be picked by the user.
- *Examples*: Please see the **Examples** folder for examples of both models, both for real-life cases or simulation experiments.
- *Papers using TLM*: 
  -  McDonald, R.R., Keith, D.M., Sameoto, J.A., Hutchings, J.A., Flemming, J.M. 2021. Incorporating intra-annual variability in fisheries abundance data to better capture population dynamics. *Fisheries Research*. In press.
- *Papers using SEBDAM*:
  - McDonald, R.R., Keith, D.M., Sameoto, J.A., Hutchings, J.A., Flemming, J.M. 2021. Explicit incorporation of spatial variability in a biomass dynamics assessment model. *ICES Journal of Marine Science*. In press.

# References

Bakka, H., Vanhatalo, J., Illian, J.B., Simpson, D., Rue, H. Non-stationary Gaussian models with physical barriers. Spatial Statistics. 29: 268-288. 2019.

Deriso, R.B. Harvesting Strategies and Parameter Estimation for an Age-Structured Model. Canadian Journal of Fisheries and Aquatic Sciences, 37:268–282, 1980.

Lindgren, F., Rue, H. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B, 73:423–498, 2011

Martin, T.G., Wintle, B.A., Rhodes, J.R., Kuhnert, P.M., Field, S.A., Low-Choy, S.J., Tyre, A.J, Possingham, H.P. Zero tolerance ecology: Improving ecological inference by modelling the source of zero observations. Ecology Letters, 8 (11): 1235-1246. 2005.

Schnute, J. A General Theory for Analysis of Catch and Effort Data. Canadian Journal of Fisheries and Aquatic Sciences, 42:414–429, 1985.

Smith, S.J., Lundy, M.J. Scallop Production Area 4 in the Bay of Fundy: Stock status and
672 forecast. Technical Report 2002/018, Fisheries and Oceans Canada, 2002.

Thorson, J.T., Shelton, A.O., Ward, E.J., Skaug, H.J. Geostatistical delta-generalized linear mixed models improve precision for estimated abundance indices for West Coast groundfishes. ICES Journal of Marine Science, 72(5):1297–1310, 2015.

