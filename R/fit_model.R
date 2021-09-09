#' Fit TLM or SEBDAM to fisheries data
#'
#' @param tmb_obj List: properly formatted data and parameter starting values from data_setup function
#' @param optim Character: choice of optimizers to use, choices are nlminb, optimr, or parallel (use multiple cores)
#' @param control list: control optimizer parameters
#' @param optim_method if using optimr, choice of optimization method
#' @param cores detects cores for parallel optimization unless otherwise specified
#' @param bias.correct for TMB::sdreport, if want to apply bias correction
#' @param silent for TMB::MakeADFun and
#'
#' @return Object containing modelling output
#' @export
#'
#' @examples
fit_model<-function(tmb_obj,optim="optimr",
                    control=NULL,optim_method="nlminb",
                    cores=parallel::detectCores(),
                    bias.correct=F,
                    silent=T) {
  if (length(tmb_obj$data$options_vec)>2){
    if (tmb_obj$data$options_vec[3]!=1) warning("Standard Errors for random fields and densities (B, R, and m) will not be calculated, if desired rerun data_setup() with all_se=T")
  }

  #Create object
  obj<-TMB::MakeADFun(data=tmb_obj$data,parameters=tmb_obj$par,random=tmb_obj$random,map=tmb_obj$map,DLL="SEBDAM",silent=silent)

  if (!(optim %in% c("nlminb","optimr","parallel"))) stop("Incorrect optimizer specification, options: nlminb, optimr, parallel (last does not work on Windows currently)")

  if (optim == "nlminb") {
    Opt<-try(stats::nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                control=control),T)
  } else if (optim == "optimr") {
    Opt<-optimx::optimr(obj$par,obj$fn,obj$gr,control=control,method=optim_method)
  } else if(optim == "parallel"){
    cl <- parallel::makeCluster(cores,type="FORK")
    parallel::setDefaultCluster(cl=cl)
    Opt = optimParallel::optimParallel(obj$par,obj$fn,obj$gr,control=control)
    parallel::setDefaultCluster(cl=NULL)
    parallel::stopCluster(cl)
  }

  rep<-TMB::sdreport(obj,bias.correct=bias.correct)

  return_obj<-list()
  return_obj$obj<-obj
  return_obj$opt<-Opt
  return_obj$report<-obj$report()
  return_obj$sdrep<-rep
  return_obj$convergence<-Opt$convergence
  return_obj$message<-Opt$message
  return_obj$log_likelihood<--(sum(return_obj$report$nll_comp))

  print(return_obj$convergence)
  print(return_obj$message)
  return(return_obj)

}
