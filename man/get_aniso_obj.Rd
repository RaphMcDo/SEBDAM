% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_aniso_obj}
\alias{get_aniso_obj}
\title{Obtain the necessary matrices and parameters to include geometric anisotropy in SEBDAM}
\usage{
get_aniso_obj(mesh)
}
\arguments{
\item{mesh}{INLA mesh object, preferably obtained from setup_mesh function}
}
\value{
list of matrices and parameters for SPDE approach with geometric anisotropy
}
\description{
Obtain the necessary matrices and parameters to include geometric anisotropy in SEBDAM
}
