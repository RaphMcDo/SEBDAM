#' Fix specific parameters to user-specified values
#'
#' @param obj TMB object ready to be fit, obtained from data_setup; same object that would be given to fit_model
#' @param pars list:

#' @return
#' @export
#'
#' @examples


fix_param<-function(obj=NULL,pars=NULL){

  if (!(is.list(pars) & !(is.data.frame(pars)) & any(names(pars) %in% names(obj$par)))){
    stop("Either no or wrong parameters provided, or wrong format")
  }

  if (!(all(names(pars) %in% names(obj$par)))) {
    warning("Some of the parameters indicated do not match with the model parameters")
  }

  #Subsetting only the parameters that are relevant and ignoring the others
  good_names<-names(pars)[which(names(pars) %in% names(obj$par))]

  for (params in good_names){
    obj$map[params]<-pars[params]
  }

  return(obj)

}
