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

    exp_rates<-return_obj$report$exp_rates
    exp_rates<-c(exp_rates,NA,NA)

    proc_frame<-data.frame(B=B,se_B=se_B,R=R,se_R=se_R,m=m,se_m=se_m,exploitation_rates=exp_rates)

    log_B<-return_obj$report$log_B
    se_log_B<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="log_B")])

    log_R<-return_obj$report$log_R
    se_log_R<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="log_R")])

    log_m<-unname(return_obj$sdrep$value[which(names(return_obj$sdrep$value)=="log_m")])
    if (return_obj$obj$env$data$options_vec[2]==1) se_log_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="log_m")])
    else se_log_m<-rep(NA,length(log_m))

    log_proc_frame<-data.frame(log_B=log_B,se_log_B=se_log_B,
                               log_R=log_R,se_log_R=se_log_R,
                               log_m=log_m,se_log_m=se_log_m)

    listy<-list(processes=proc_frame,log_processes=log_proc_frame)

  } else if (return_obj$obj$env$data$model == "SEBDAM"){

    if (return_obj$obj$env$data$options_vec[3]!=1) warning("Standard errors for density processes were not calculated as part of fitting process, if desired please rerun data_setup() with all_se=T")

    #Obtaining density processes
    B<-return_obj$report$B
    if (return_obj$obj$env$data$options_vec[3]==1) se_B<-matrix(unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="B")]),ncol=return_obj$obj$env$data$n_t+1)
    else se_B<-matrix(rep(NA,length(B)),ncol=return_obj$obj$env$data$n_t+1)

    R<-return_obj$report$R
    R<-c(R,rep(NA,return_obj$obj$env$data$n_s))
    R<-matrix(R,ncol=return_obj$obj$env$data$n_t+1,nrow=return_obj$obj$env$data$n_s)
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

    exp_rates<-return_obj$report$exp_rates
    exp_rates<-c(exp_rates,NA,NA)

    tot_frame<-data.frame(totB=totB,se_totB=se_totB,totR=totR,se_totR=se_totR,mean_m=mean_m,se_mean_m=se_mean_m,exploitation_rates=exp_rates)

    log_totB<-return_obj$report$log_totB
    se_log_totB<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="log_totB")])

    log_totR<-return_obj$report$log_totR
    log_totR<-c(log_totR,NA)
    se_log_totR<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="log_totR")])
    se_log_totR<-c(se_log_totR,NA)

    log_mean_m<-unname(return_obj$sdrep$value[which(names(return_obj$sdrep$value)=="log_mean_m")])
    se_log_mean_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$value)=="log_mean_m")])

    log_tot_frame<-data.frame(log_totB=log_totB,se_log_totB=se_log_totB,
                              log_totR=log_totR,se_log_totR=se_log_totR,
                              log_mean_m=log_mean_m,se_log_mean_m=se_log_mean_m)

    listy<-list(densities=dens_list,totals=tot_frame,log_tot_frame=log_tot_frame)
  }

  return(listy)

}
