#' Fix specific parameters to user-specified values
#'
#' @param obj TMB object ready to be fit, obtained from data_setup; same object that would be given to fit_model
#' @param pars list:

#' @return
#' @export
#'
#' @examples


fix_param<-function(tmb_obj=NULL,pars=NULL,par_vals=NULL){

  if (!(is.list(pars) & !(is.data.frame(pars)) & any(names(pars) %in% names(obj$par)))){
    stop("Either no or wrong parameters provided, or wrong format")
  }

  if (!(all(names(pars) %in% names(obj$par)))) {
    warning("Atleast some of the parameters indicated do not match with the model parameters")
  }

  if (!(length(pars)==length(par_vals))){
    stop("Number of parameters does not match up with number of values provided")
  }

  #Subsetting only the parameters that are relevant and ignoring the others
  good_names<-names(pars)[which(names(pars) %in% names(obj$par))]
  good_vals<-par_vals[which(names(pars) %in% names(obj$par))]

  for (params in good_names){
    obj$map[params]<-as.factor(NA)
    obj$par[params]<-good_vals[which(names(good_names)==params)]
  }

  return(obj)

}
