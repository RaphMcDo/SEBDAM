## SEBDAM
- SEBDAM is an R package for fitting a Spatially-Explicit Biomass Dynamics Assessment Model in a state-space framework to fishery-independent data. Built for the assessment of scallop fisheries, it can easily be applied to any marine animals that exhibits spatially heterogeneous area utilization inside of a given stock.
- The SEBDAM package contains 2 models: SEBDAM, previously described, and TLM (Tow-Level Model), a strictly temporal biomass dynamics models that can be directly compared to SEBDAM.

## Background
- SEBDAM was designed to incorporate spatial variability and autocorrelations in fishery-independent observations directly into traditional biomass dynamics model and model the underyling processes (biomass, recruitment and natural mortality) as explicitly spatial processes.
- This approach builds on different spatial and temporal techniques, incorporating advances related to the use of Gaussian Markov Random Fields through the Stochastic Partial Differential Equations approach (Lindgren and Rue, 2012), the use of "knots" to create a predictive framework (Thorson et al., 2015), delta models (Martin et al., 2005), and of a delay-difference stock assessment model (Deriso, 1980, Schnute, 1985, Smith and Lundy 2002).
- The use of Gaussian Markov Random Fields assumes that spatial correlations decrease according to distance. This decorrelation based on distance can further account for geometric anisotropy (Lindgren et al., 2012) or the presence of land through a barrier model (Bakka et al., 2019).

# Resources

- *User Manual*: Please see **Under Construction** for the model equations in relations to the input and output as seen when fitting the model, along with all possible options that can be picked by the user.
- *Examples*: Please see the **Under Construction** folder for many simple to complex examples, both for real-life cases or simulation experiments.
