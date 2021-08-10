
#' Find the EPGS code for UTM zone code of given WGS 84 longitude and latitude
#'
#' @param lonlat Vector of WGS 84 with entries of longitude and latitude
#'
#' @return EPSG code UTM zone number
#' @export
#'
#' @examples
find_utm_code <- function(lonlat) {
  ## Special Cases for Norway & Svalbard
  if (lonlat[2]> 55 & lonlat[2]< 64 & lonlat[1]> 2 & lonlat[1]< 6){
    utm <- 32
  } else if (lonlat[2]> 71 & lonlat[1]>= 6 & lonlat[1]< 9){
    utm <- 31
  } else if (lonlat[2]> 71 & lonlat[1]>= 9 & lonlat[1]< 12){
    utm <- 33
  } else if (lonlat[2]> 71 & lonlat[1]>= 18 & lonlat[1]< 21){
    utm <- 33
  } else if (lonlat[2]> 71 & lonlat[1]>= 21 & lonlat[1]< 24){
    utm <- 35
  } else if (lonlat[2]> 71 & lonlat[1]>= 30 & lonlat[1]< 33){
    utm <- 35
  } else utm<-(floor((lonlat[1]+ 180)/6) %% 60) + 1

  if(lonlat[2] > 0) {
    return(utm + 32600)
  } else {
    return(utm + 32700)
  }

}

#' Logit
#'
#' @param x Value of which the logit is desired
#'
#' @return Logit of x
#' @export
#'
#' @examples
logit<-function(x) {log(x/(1-x))}

#' Obtain the necessary matrices and parameters to include geometric anisotropy in SEBDAM
#'
#' @param mesh INLA mesh object, preferably obtained from setup_mesh function
#'
#' @return list of matrices and parameters for SPDE approach with geometric anisotropy
#' @export
#'
#' @examples
get_aniso_obj<-function(mesh){
  spde = INLA::inla.spde2.matern(mesh)
  Dset = 1:2
  TV = mesh$graph$tv           # Triangle to vertex indexing
  V0 = mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
  V1 = mesh$loc[TV[,2],Dset]
  V2 = mesh$loc[TV[,3],Dset]
  E0 = V2 - V1                      # E = edge for each triangle
  E1 = V0 - V2
  E2 = V1 - V0
  # Calculate Areas
  TmpFn = function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
  Tri_Area = rep(NA, nrow(E0))
  for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
  spde_aniso <- list(
    "n_s"      = spde$n.spde,
    "n_tri"    = nrow(TV),
    "Tri_Area" = Tri_Area,
    "E0"       = E0,
    "E1"       = E1,
    "E2"       = E2,
    "TV"       = TV - 1,
    "G0"       = spde$param.inla$M0,
    "G0_inv"   = as(diag(1/diag(as.matrix(spde$param.inla$M0))), "dgTMatrix"))

  return(spde_aniso)
}

#' Finding INLA mesh vertices located on land for barrier model
#'
#' @param mesh INLA mesh, preferably obtained from setup_mesh function
#' @param bound sf objects: modelling bounds including islands/land for barrier model
#'
#' @return object for barrier model indicating triangles/vertices on land
#' @export
#'
#' @examples
get_triangles<-function(mesh,bound) {
  #FInd which triangles are on land
  tl = length(mesh$graph$tv[,1])
  # - the number of triangles in the mesh
  posTri = matrix(0, tl, 2)
  for (t in 1:tl){
    temp = mesh$loc[mesh$graph$tv[t, ], ]
    posTri[t,] = colMeans(temp)[c(1,2)]
  }
  posTri = sp::SpatialPoints(posTri)
  # - the positions of the triangle centres

  barrier = sp::over(sf::as_Spatial(bound), sp::SpatialPoints(posTri), returnList=T)
  # - checking which mesh triangles are inside the barrier area
  barrier = unlist(barrier)
  return(barrier)
}

#' Randomly create vector of integers that sum up to specified total
#'
#' @param min minimum number desired
#' @param max max number desired
#' @param n_grp number of integers to be simulated
#' @param total sum of integers desired
#'
#' @return vector of integers
#' @export
#'
#' @examples
rand_int_vect<-function(min,max,n_grp,total) {

  if (max*n_grp < total | min*n_grp > total) stop("Values provided cannot work together")

  x<-sample(min:max,n_grp,replace=T)

  while (sum(x) < total) {
    loc<-sample(1:n_grp,1)
    if (x[loc] < max) x[loc]<-x[loc]+1
  }

  while (sum(x) > total) {
    loc<-sample(1:n_grp,1)
    if (x[loc]>0 & x[loc]>min) x[loc]<-x[loc]-1
  }

  return(x)

}

