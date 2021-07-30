
#' Setting up the predictive knots and mesh based on the knots
#'
#' @param data sf object containing data and locations of all tows
#' @param seed specify seed number for reproducibility
#' @param nknot integer: number of knots wanted for number of knots, default set at 25
#' @param model_bound sf object for boundaries of modelling area
#' @param crs_input Specify crs that want to model in, but not tested for anything that is not UTM at the moment
#' @param max.edge inla.mesh.2d control
#' @param cutoff inla.mesh.2d control
#' @param ... other parameters fed to inla.mesh.2d
#'
#' @return
#' @export
#'
#' @examples
setup_mesh<-function(data=NULL,seed=NULL,nknot=25,model_bound=NULL,
                     crs_input=NULL, max.edge=c(5,30),cutoff=2, ...) {

  #Verify if data is an sf object and has a crs
  if (is.null(data) | !("sf" %in% attr(data,"class")) | !("sf" %in% attr(model_bound,"class"))) {
    stop("Requires sf object for data and/or modelling boundaries")
  }
  if (is.na(sf::st_crs(data))) stop("Data needs a coordinate reference system")

  #Transform to WGS 84 to standardize no matter what the crs is, then into UTM
  wgs_data<-sf::st_transform(data,crs=4326)
  wgs_bound<-sf::st_transform(model_bound,crs=4326)
  if (!is.null(crs_input)) {
    crs_change<-crs_input
    warning("Non UTM CRS not tested, weird results likely")
  } else {
    crs_change<-find_utm_code(colMeans(sf::st_coordinates(wgs_data)))
  }
  utm_data<-sf::st_transform(wgs_data,crs=crs_change)

  #From the lat long, create the knots
  if (is.numeric(seed)) set.seed(as.integer(seed)) else warning("No seed set, knot creation not reproducible")
  knots<-tryCatch(stats::kmeans(sf::st_coordinates(utm_data),nknot,nstart=25),
                  warning=function(w) stop("kmeans ",w,", try different seed"))

  knots.loc<-as.data.frame(knots[[2]])
  knots.loc<-dplyr::rename(knots.loc,lon=X,lat=Y)
  knots.loc<-sf::st_as_sf(knots.loc,coords=c("lon","lat"),crs=crs_change)

  if (is.null(crs_input)) knots.loc<-sf::st_sf(knots.loc/1000,crs=crs_change)
  if (is.null(crs_input)) knots[[2]]<-knots[[2]]/1000

  utm_bound<-sf::st_transform(wgs_bound,crs=crs_change)
  utm_bound$geometry<-utm_bound$geometry/1000

  # Create mesh
  mesh<-INLA::inla.mesh.2d(sf::st_coordinates(knots.loc),
                     max.edge=max.edge, cutoff=cutoff,
                     boundary=INLA::inla.sp2segment(sf::as_Spatial(utm_bound)),
                     ...)

  end_list<-list(knots=knots,mesh=mesh,crs=crs_change,utm_bound=utm_bound)
  return(end_list)
}
