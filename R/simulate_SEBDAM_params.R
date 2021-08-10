#' Setup parameters for simulations using SEBDAM
#'
#' @param simul_data_obj List: object obtained from simulate_obj function
#' @param sigma_epsilon double: biomass observation variance parameter
#' @param sigma_upsilon double: recruitment observation variance parameter
#' @param R0 double: recruitment density mean value in initial year
#' @param B0 double: biomass density mean value in initial year
#' @param m0 double: if using clapper approach, natural mortality mean value in initial year, if not then value at which natural mortality will be fixed
#' @param p_I double: biomass probability of capture
#' @param p_IR double: recruit probability of capture
#' @param qR double: recruitment catchability
#' @param range_B double: range at which two points become decorrelated for biomass
#' @param range_R double: range at which two points become decorrelated for recruitment
#' @param sigma_B double: marginal variance of biomass random fields
#' @param sigma_R double: marginal variance of recruit random fields
#' @param H_input_B double: anisotropy parameters for biomass random fields
#' @param H_input_R double: anisotropy parameters for recruitment random fields
#' @param qI double: biomass catchability
#' @param range_m double: range at which two points become decorrelated for natural mortality
#' @param sigma_m double: marginal variance of natural mortality random fields
#' @param H_input_m double: anisotropy parameters for natural mortality random fields
#' @param S double: clapper catchability
#'
#' @return Formatted list for simulation purposes using SEBDAM
#' @export
#'
#' @examples
simulate_SEBDAM_params<-function(simul_data_obj=NULL,
                                 sigma_epsilon=0.1,sigma_upsilon=0.1,
                                 R0=NULL,B0=NULL,m0=0.1, p_I=0.9,
                                 p_IR=0.4,qR=0.2,
                                 range_B=40,range_R=30,
                                 sigma_B=0.1,sigma_R=0.1,
                                 H_input_B=c(-1.2,-0.5),
                                 H_input_R=c(-1,-0.8),
                                 qI=0.45,
                                 range_m=20,sigma_m=0.1,
                                 H_input_m=c(-0.5,0.4),
                                 S=0.4) {

  if (is.null(simul_data_obj)) stop("Missing simulated data object, run simulate_obj() and use its output for this function")

  if (simul_data_obj$model=="SEBDAM") {

    if (is.null(R0)) R0 = mean(exp(simul_data_obj$logIR)*10,na.rm=T)
    if (is.null(B0)) B0 = mean(exp(simul_data_obj$logI)*10,na.rm=T)

    if(!all(c(sigma_epsilon>0,sigma_upsilon>0,R0>0,B0>0,
              m0>0,p_I>0,p_IR>0,qR>0,range_B>0,range_R>0,
              sigma_B>0,sigma_R>0,qI>0,range_m>0,sigma_m>0,
              S>0))) stop("All parameters (except anisotropy) must be positive")

    temp_par_list<-list()
    temp_par_list$log_sigma_epsilon<-log(sigma_epsilon)
    temp_par_list$log_sigma_upsilon<-log(sigma_upsilon)
    temp_par_list$log_R0<-log(R0)
    temp_par_list$log_B0<-log(B0)
    temp_par_list$log_m0<-log(m0)
    temp_par_list$logit_p_I<-logit(p_I)
    temp_par_list$logit_p_IR<-logit(p_IR)
    temp_par_list$log_qR<-log(qR)

    temp_par_list$omega_B<-matrix(rep(0,(simul_data_obj$n_t+1)*simul_data_obj$n_m),ncol=(simul_data_obj$n_t+1))
    temp_par_list$omega_R<-matrix(rep(0,simul_data_obj$n_t*simul_data_obj$n_m),ncol=simul_data_obj$n_t)

    if (simul_data_obj$options_vec[6] %in% c(0,1)){
      temp_par_list$log_kappa_B<-log(sqrt(8)/range_B)
      temp_par_list$log_tau_B<-log(1/((sigma_B)*sqrt(4*pi)*exp(temp_par_list$log_kappa_B)))
      temp_par_list$log_kappa_R<-log(sqrt(8)/range_R)
      temp_par_list$log_tau_R<-log(1/((sigma_R)*sqrt(4*pi)*exp(temp_par_list$log_kappa_R)))
      if (simul_data_obj$options_vec[6] == 1) {
        temp_par_list$log_H_input_B<-H_input_B
        temp_par_list$log_H_input_R<-H_input_R
      }
    } else if (simul_data_obj$options_vec[6] == 2) {
      temp_par_list$log_range_B<-c(log(range_B),log(range_B*0.1))
      temp_par_list$log_sigma_B<-log(sigma_B)
      temp_par_list$log_range_R<-c(log(range_R),log(range_R*0.1))
      temp_par_list$log_sigma_R<-log(sigma_R)
    }

    temp_random<-c("omega_B","omega_R")
    temp_map<-list()

    if (simul_data_obj$options_vec[2] == 1) {
      if(!(length(qI)>1)) stop("Length of qI needs to be the same as number of knots, needs to be specified manually if more than 1")
      temp_par_list$log_qI<-log(qI)
    }
    else if (simul_data_obj$options_vec[2] == 0) temp_par_list$log_qI<-log(qI)

    if (simul_data_obj$options_vec[4] == 1){
      temp_par_list$omega_m<-matrix(rep(0,(simul_data_obj$n_t+1)*simul_data_obj$n_m),ncol=(simul_data_obj$n_t+1))
      if (simul_data_obj$options_vec[6] %in% c(0,1)){
        temp_par_list$log_kappa_m<-log(sqrt(8)/range_m)
        temp_par_list$log_tau_m<-log(1/((sigma_m)*sqrt(4*pi)*exp(temp_par_list$log_kappa_m)))
        if (simul_data_obj$options_vec[6] == 1) temp_par_list$log_H_input_m<-H_input_m
      } else if (simul_data_obj$options_vec[6] == 2){
        temp_par_list$log_range_m<-c(log(range_m),log(range_m*0.1))
        temp_par_list$log_sigma_m<-log(sigma_m)
      }
      temp_random[3]<-"omega_m"
      temp_par_list$log_S<-log(S)
    } else if (obs_mort == FALSE) temp_map<-list(log_m0=as.factor(NA))

    simul_obj<-list(simul_data=simul_data_obj,simul_par=temp_par_list,random=temp_random,map=temp_map)

  } else stop("Function is for SEBDAM, use simul_TLM_params() for TLM")

  return(simul_obj)

}
