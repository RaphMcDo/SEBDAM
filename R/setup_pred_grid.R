#' Setting up prediction grid for spatial model
#'
#' @param knots kmeans: knot object obtained from setup_mesh function
#' @param model_bound sf object: provided by user, modelling area and bounds
#'
#' @return List containing gridded bound for modelling purposes and area associated with each knot
#' @export
#'
#' @examples
setup_pred_grid<-function(knots=NULL,model_bound=NULL) {

  if(is.null(knots) | is.null(model_bound) | !("sf" %in% attr(model_bound,"class")) | !("kmeans" %in% attr(knots,"class"))) {
    stop("Either knot object or modelling bound object don't exist or are wrong type of objects")
  }

  #Setup a 1x1 grid based on the bounds provided
  grid<-sf::st_make_grid(model_bound,cellsize = 1)
  #Cut out the parts that we don't want modelled
  gridded_bound<-sf::st_intersection(model_bound,grid)
  #Pick out centers of each 1x1 grid cell to attribute them to knots
  centroids<-sf::st_centroid(gridded_bound)

  #Make the points for the knots
  temp_knots<-as.data.frame(knots[[2]])
  colnames(temp_knots)<-c("x","y")
  knots_sf<- sf::st_as_sf(temp_knots,coords=c("x","y"),crs=sf::st_crs(model_bound))

  #Get distance between centroids of grid and knots
  distmat<-sf::st_distance(centroids, knots_sf)
  polyknotID<-c()
  for (i in 1:length(distmat[,1])) {
    polyknotID<-c(polyknotID,which(distmat[i,]==min(distmat[i,])))
  }
  gridded_bound$knotID<-polyknotID
  gridded_bound$area<-sf::st_area(gridded_bound)

  #Area of each knots
  stratarea<-stats::aggregate(area~knotID,sum,data=gridded_bound)

  listy<-list(grid=gridded_bound,area=stratarea)
  return(listy)

}
