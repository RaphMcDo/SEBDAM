#' Obtain data frame of underlying processes and standard errors
#'
#' @param return_obj Fitted modelling output from fit_model function
#'
#' @return dataframe: underlying processes and standard errors
#' @export
#'
#' @examples
get_processes<-function(return_obj) {

  if (return_obj$obj$env$data$model == "TLM"){

    B<-return_obj$report$B
    se_B<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="B")])

    R<-return_obj$report$R
    se_R<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="R")])

    m<-return_obj$report$m
    if (return_obj$obj$env$data$options_vec[2]==1) se_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="m")])
    else se_m<-rep(NA,length(m))

    proc_frame<-data.frame(B=B,se_B=se_B,R=R,se_R=se_R,m=m,se_m=se_m)

    listy<-list(processes=proc_frame)

  } else if (return_obj$obj$env$data$model == "SEBDAM"){

    if (return_obj$obj$env$data$options_vec[3]!=1) warning("Standard errors for processes were not calculated as part of fitting process, if desired please rerun with options_vec[3]==1")

    #Obtaining density processes
    B<-return_obj$report$B
    if (return_obj$obj$env$data$options_vec[3]==1) se_B<-matrix(unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="B")]),ncol=return_obj$obj$env$data$n_t+1)
    else se_B<-matrix(rep(NA,length(B)),ncol=return_obj$obj$env$data$n_t+1)

    R<-return_obj$report$R
    R<-c(R,rep(NA,return_obj$obj$env$data$n_s))
    if (return_obj$obj$env$data$options_vec[3]==1) {
      se_R<-matrix(unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="R")]),ncol=return_obj$obj$env$data$n_t)
      se_R[,(return_obj$obj$env$data$n_t)]<-rep(NA,return_obj$obj$env$data$n_s)
    } else se_R<-matrix(rep(NA,length(R)),ncol=return_obj$obj$env$data$n_t+1)

    m<-return_obj$report$m
    if (return_obj$obj$env$data$options_vec[3]==1) se_m<-matrix(unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="m")]),ncol=return_obj$obj$env$data$n_t+1)
    else se_m<-matrix(rep(NA,length(m)),ncol=return_obj$obj$env$data$n_t+1)

    dens_list<-list(B=B,se_B=se_B,R=R,se_R=se_R,m=m,se_m=se_m)

    #Obtaining predictions for whole area
    totB<-return_obj$report$totB
    se_totB<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="totB")])

    totR<-return_obj$report$totR
    totR<-c(totR,NA)
    se_totR<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="totR")])
    se_totR<-c(se_totR,NA)

    mean_m<-return_obj$report$mean_m
    se_mean_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="mean_m")])

    tot_frame<-data.frame(totB=totB,se_totB=se_totB,totR=totR,se_totR=se_totR,mean_m=mean_m,se_mean_m=se_mean_m)

    listy<-list(densities=dens_list,totals=tot_frame)
  }

  return(listy)

}
