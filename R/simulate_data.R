#' Simulate data and underlying processes using either TLM or SEBDAM
#'
#' @param simul_obj List: object for simulations obtained from simul_TLM_params or simul_SEBDAM_params function
#' @param seed integer: to set the seed for data simulation
#' @param format character: choose format of output, "raw" simply returns the list directly created by obj$simulate, "formatted" returns an organized list separating data, growths, catches, and simulated processes
#' @param sim_obs_loc sfc objecT: tow locations obtained from simul_area function, only needed if using SEBDAM
#'
#' @return List containing simulated data and processes
#' @export
#'
#' @examples
simulate_data<-function(simul_obj=NULL,seed=NULL,format="formatted",
                        sim_obs_loc=NULL) {

  #Setup simulation object
  obj <- MakeADFun( data=simul_obj$simul_data, parameters=simul_obj$simul_par,
                    random=simul_obj$random, map = simul_obj$map )
  if (!is.null(seed) & is.numeric(seed)) set.seed(seed)
  simdata <- obj $ simulate(complete=T)

  if (format == "raw") {
    return_obj<-simdata
  } else if (format == 'formatted') {
    temp_data<-data.frame(I=exp(simdata$logI),IR=exp(simdata$logIR),
                          Year=simdata$t_i,Knot=simdata$s_i)
    temp_data$I[which(is.na(temp_data$I))]<-0
    temp_data$IR[which(is.na(temp_data$IR))]<-0
    if (!is.null(simdata$L)){
      temp_data$L<-simdata$L
      temp_data$N<-simdata$n_bin-sim_data$L
    }
    temp_growths<-data.frame(g=simdata$g,gR=simdata$gR)
    temp_catch<-simdata$C

    check_frame<-data.frame(n_tows=simdata$n_tows,pos_tows_I=simdata$pos_tows_I,
                            pos_tows_IR=simdata$pos_tows_IR)

    temp_proc<-list()

    if (simdata$model=="SEBDAM"){
      temp_data$geometry<-sim_obs_loc
      temp_data<-st_sf(temp_data)

      temp_proc$density_B<-simdata$B
      temp_proc$density_R<-simdata$R
      temp_proc$totals<-data.frame(B=simdata$totB,R=c(simdata$totR,NA))
      temp_proc$omega_B<-simdata$omega_B
      temp_proc$omega_R<-simdata$omega_R

      if (simdata$options_vec[4]==1){
        temp_proc$m<-simdata$m
        temp_proc$totals$m<-simdata$mean_m
        temp_proc$omega_m<-simdata$omega_m
      }

    } else if (simdata$model=="TLM") {
      temp_proc<-data.frame(B=simdata$B,R=simdata$R,m=simdata$m)
    }

    return_obj<-list(data=temp_data,growths=test_growths,
                     catch=test_catch,processes=temp_proc)

  }

  return(return_obj)

}
