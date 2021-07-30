#' Spread raw catch data to appropriate knots
#'
#' @param catch sf object with columns Catch, Year, geometry
#' @param knots kmeans object obtained from setup_mesh function
#' @param area area covered by each knot as obtained from setup_pred_grid function
#' @param crs_input Desired crs, currently only UTM is accepted
#'
#' @return List containing properly formatted catch, either sums per at each knot or density (catch/km^2) at each knot
#' @export
#'
#' @examples
catch_spread<-function(catch=NULL,knots=NULL,area=NULL,crs_input=NULL) {

  if (!("sf" %in% attr(catch,"class")) | !("kmeans" %in% attr(knots,"class")) | colnames(catch)!=c("Catch","Year","geometry")) {
    stop("Either catch or knots are wrong formats")
  }
  if (is.na(sf::st_crs(catch))) stop("Catch needs a coordinate reference system")

  #Transform to WGS 84 to standardize no matter what the crs is, then into UTM
  wgs_catch<-sf::st_transform(catch,crs=4326)
  if (!is.null(crs_input)) {
    crs_change<-crs_input
    warning("Non UTM CRS not tested, weird results likely")
  } else {
    crs_change<-find_utm_code(colMeans(sf::st_coordinates(wgs_catch)))
  }
  utm_catch<-sf::st_transform(wgs_catch,crs=crs_change)

  temp_knots<-as.data.frame(knots[[2]])
  colnames(temp_knots)<-c("x","y")
  knots_sf<- sf::st_as_sf(temp_knots,coords=c("x","y"),crs=sf::st_crs(utm_catch))


  #Attribute each fishing trip to closest knot
  distmat<-sf::st_distance(catch,knots_sf)
  polyknotID<-c()
  for (i in 1:length(distmat[,1])) {
    polyknotID<-c(polyknotID,which(distmat[i,]==min(distmat[i,])))
  }
  utm_catch$knotID<-polyknotID

  #By year and knot
  agcatch<-aggregate(Catch~Year+knotID,data=utm_catch,FUN=sum)

  catch_frame<-data.frame("knotID"=c(1:length(unique(polyknotID))))
  for (i in min(utm_catch$Year):max(utm_catch$Year)){
    catchtemp <-subset(agcatch,Year==i)
    catchtemp <-catchtemp[,-1]
    catch_frame<-left_join(catch_frame,catchtemp,by="knotID")
  }
  colnames(catch_frame)<-c("knotID",paste(1998:2019))
  #Change empty entries to 0, cause 0 catches there that year
  catch_frame[is.na(catch_frame)]<-0

  sum_catches<-catch_frame
  density_catches<-catch_frame/area$area

  listy<-list(sum_catches=sum_catches,density_catches=density_catches)
  return(listy)

}
