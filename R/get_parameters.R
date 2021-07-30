#' Obtain data frames of parameter and their standard errors
#'
#' @param return_obj Fitted modelling output for fit_model function
#'
#' @return dataframe: parameter estimates and standard errors
#' @export
#'
#' @examples
get_parameters<-function(return_obj) {

  par_names<-names(return_obj$obj$par)
  par_names<-stringr::str_remove(par_names,"log_")
  par_names<-stringr::str_remove(par_names,"logit_")

  pars<-return_obj$sdrep$value[which(names(return_obj$sdrep$value) %in% par_names)]
  pars_se<-return_obj$sdrep$sd[which(names(return_obj$sdrep$value) %in% par_names)]

  par_frame<-as.data.frame(cbind(pars,pars_se))
  colnames(par_frames)<-c("Estimate","SE")

  return(par_frame)
}
