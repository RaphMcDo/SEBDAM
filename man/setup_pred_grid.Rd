% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_pred_grid.R
\name{setup_pred_grid}
\alias{setup_pred_grid}
\title{Setting up prediction grid for spatial model}
\usage{
setup_pred_grid(knots = NULL, model_bound = NULL)
}
\arguments{
\item{knots}{kmeans: knot object obtained from setup_mesh function}

\item{model_bound}{sf object: provided by user, modelling area and bounds}
}
\value{
List containing gridded bound for modelling purposes and area associated with each knot
}
\description{
Setting up prediction grid for spatial model
}
