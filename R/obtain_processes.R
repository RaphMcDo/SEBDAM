#' Obtain data frame of underlying processes and standard errors
#'
#' @param return_obj Fitted modelling output from fit_model function
#'
#' @return dataframe: underlying processes and standard errors
#' @export
#'
#' @examples
obtain_processes<-function(return_obj) {

  if (return_obj$obf$env$data$model == "TLM"){

    B<-return_obj$report$B
    se_B<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="B")])

    R<-return_obj$report$R
    R<-c(R,rep(NA,return_obj$obj$env$data$n_s))
    se_R<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="R")])
    se_R<-c(se_R,rep(NA,report_obj$obj$env$data$n_s))

    m<-return_obj$report$m
    se_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="m")])

    proc_frame<-data.frame(B=B,se_B=se_B,R=R,se_R=se_R,m=m,se_m=se_m)

    listy<-list(processes=proc_frame)

  } else if (return_obj$obj$env$data$model == "SEBDAM"){

    if (return_obj$obj$env$data$options_vec[3]!=1) warning("Standard errors for processes were not calculated as part of fitting process, if desired please rerun with options_vec[3]==1")

    #Obtaining density processes
    B<-return_obj$report$B
    if (return_obj$obj$env$data$options_vec[3]==1) se_B<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="B")])
    else se_B<-rep(NA,length(B))

    R<-return_obj$report$R
    R<-c(R,rep(NA,return_obj$obj$env$data$n_s))
    if (return_obj$obj$env$data$options_vec[3]==1) {
      se_R<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="R")])
      se_R<-c(se_R,rep(NA,report_obj$obj$env$data$n_s))
    } else se_R<-rep(NA,length(R))

    m<-return_obj$report$m
    if (return_obj$obj$env$data$options_vec[3]==1) se_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="m")])
    else se_m<-rep(NA,length(m))

    dens_frame<-data.frame(B=B,se_B=se_B,R=R,se_R=se_R,m=m,se_m=se_m)

    #Obtaining predictions for whole area
    totB<-return_obj$report$totB
    se_totB<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="totB")])

    totR<-return_obj$report$totR
    totR<-c(totR,NA)
    se_totR<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="totR")])
    se_totR<-c(se_totR,NA)

    mean_m<-return_obj$report$mean_m
    se_mean_m<-unname(return_obj$sdrep$sd[which(names(return_obj$sdrep$sd)=="mean_m")])

    tot_frame<-data.frame(totB=totB,se_totB=se_totB,totR=totR,se_totR=se_totR,mean_m=mean_m,se_mean_m=se_mean_m)

    listy<-list(densities=dens_frame,totals=tot_frame)
  }

  return(listy)

}
