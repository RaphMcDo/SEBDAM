#' Create an area, random locations, knots and corresponding mesh for the purpose of simulations
#'
#' @param n_obs Number of simulated observations
#' @param n_knots Number of desired knots
#' @param seed Seed for reproducibility of knots
#' @param ... Given to sf::st_sample for different types of random samplings, and to INLA::inla.mesh.2d for additional parameters, especially if not using defaults
#' @param x_coord x coordinates for creation of simulation area
#' @param y_coord y coordinates for creation of simulation area
#'
#' @return list containing sf object for simulation area, kmeans object for knots, inla mesh, and sfc object with location of every simulated points
#' @export
#'
#' @examples
simulate_area<-function(n_obs=NULL, n_knots=NULL,
                     seed=NULL,
                     ...,
                     x_coord=c(175,225,225,175),
                     y_coord=c(175,175,225,225)) {

  if (is.null(n_obs) | is.null(n_knots)) stop("Provide desired number of observation and desired number of knots")
  if (!is.numeric(n_obs) | !is.numeric(n_knots)) stop("Number of observations and knots need to be numeric values")

  n_obs<-as.integer(n_obs)
  n_knots<-as.integer(n_knots)

  #Create modelling area, default is a 100X100km square
  xy<-cbind(x_coord,y_coord)
  points_xy<-sf::st_multipoint(xy)
  poly_xy<-sf::st_sf(sf::st_sfc(sf::st_cast(points_xy,"POLYGON")))
  colnames(poly_xy)<-"geometry"
  sf::st_geometry(poly_xy)<-"geometry"

  #Sample locations, default is random distribution across space
  rand_loc<-sf::st_sample(poly_xy,size=n_obs,...)

  #Create knots from these random locations
  if (is.numeric(seed)) set.seed(as.integer(seed)) else warning("No seed set, knot creation not reproducible")
  knots<-tryCatch(stats::kmeans(sf::st_coordinates(rand_loc),n_knots,nstart=25),
                  warning=function(w) stop("kmeans ",w,", try different seed"))

  knots.loc<-as.data.frame(knots[[2]])
  knots.loc<-dplyr::rename(knots.loc,lon=X,lat=Y)
  knots.loc<-sf::st_as_sf(knots.loc,coords=c("lon","lat"))

  #Create mesh based on the knots
  if (suppressWarnings(all(x_coord==c(175,225,225,175)) & all(y_coord==c(175,175,225,225)))){
    mesh<-INLA::inla.mesh.2d(sf::st_coordinates(knots.loc),
                           max.edge=c(8,30),cutoff=2,
                           boundary=INLA::inla.sp2segment(sf::as_Spatial(poly_xy)),
                           ...)
    cellsize<-1
  } else {
    mesh<-tryCatch(INLA::inla.mesh.2d(sf::st_coordinates(knots.loc),
                             boundary=INLA::inla.sp2segment(sf::as_Spatial(poly_xy)),
                             ...),
    error=function(e) stop(e,"for inla.mesh.2d when not using defaults"))
    cellsize<-10^(nchar(as.integer(mean(max(x_coord)-min(x_coord))))-3)
  }


  #Setup a 1x1 grid based on the bounds provided
  grid<-sf::st_make_grid(poly_xy,cellsize = cellsize)
  #Cut out the parts that we don't want modelled
  gridded_bound<-sf::st_intersection(poly_xy,grid)
  #Pick out centers of each 1x1 grid cell to attribute them to knots
  centroids<-sf::st_centroid(gridded_bound)

  #Get distance between centroids of grid and knots
  distmat<-sf::st_distance(centroids, knots.loc)
  polyknotID<-c()
  for (i in 1:length(distmat[,1])) {
    polyknotID<-c(polyknotID,which(distmat[i,]==min(distmat[i,])))
  }
  gridded_bound$knotID<-polyknotID
  gridded_bound$area<-sf::st_area(gridded_bound)

  #Area of each knots
  stratarea<-stats::aggregate(area~knotID,FUN=sum,data=sf::st_drop_geometry(gridded_bound))

  listy<-list(area=poly_xy,knots=knots,mesh=mesh,
              sim_obs=rand_loc,stratarea=stratarea,knots_sf=knots.loc,
              grid=gridded_bound)
  return(listy)

}
