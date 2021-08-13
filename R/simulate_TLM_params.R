#' Setup parameters for simulations using TLM
#'
#' @param simul_data_obj List: object obtained from simulate_obj function
#' @param sigma_tau double: biomass process variance parameter
#' @param sigma_phi double: recruitment process variance parameter
#' @param sigma_epsilon double: biomass observation variance parameter
#' @param sigma_upsilon double: recruitment observation variance parameter
#' @param q_R double: recruitment catchability
#' @param q_I double: biomass catchability
#' @param p_I double: probability of biomass capture
#' @param p_IR double: probability of recruit capture
#' @param sigma_m double: natural mortality process variance parameter
#' @param S double: clapper catchability
#'
#' @return Formatted list for simulation purposes using TLM
#' @export
#'
#' @examples
simulate_TLM_params<-function(simul_data_obj=NULL, sigma_tau=0.1,
                           sigma_phi=0.1, sigma_epsilon=0.1,
                           sigma_upsilon=0.1,q_R=0.2,
                           q_I=0.45,p_I=0.9,p_IR=0.4,
                           sigma_m=0.1,S=0.4) {

  if (is.null(simul_data_obj)) stop("Missing simulated data object, run simulate_data_obj() and use its output for this function")

  if (!(all(c(sigma_tau>0,sigma_phi>0,sigma_epsilon>0,
              sigma_upsilon>0,q_R>0,q_I>0,p_I>0,p_IR>0,
              sigma_m>0,S>0)))) stop("All parameters must be positive")

  if (simul_data_obj$model=="TLM") {

    temp_par_list<-list()
    temp_par_list$log_sigma_tau<-log(sigma_tau)
    temp_par_list$log_sigma_phi<-log(sigma_phi)
    temp_par_list$log_sigma_epsilon<-log(sigma_epsilon)
    temp_par_list$log_sigma_upsilon<-log(sigma_upsilon)
    temp_par_list$log_q_R<-log(q_R)
    temp_par_list$log_q_I<-log(q_I)
    temp_par_list$logit_p_I<-logit(p_I)
    temp_par_list$logit_p_IR<-logit(p_IR)

    temp_par_list$log_B<-rep(1,length(unique(simul_data_obj$t_i))+1)
    temp_par_list$log_R<-rep(1,length(unique(simul_data_obj$t_i))+1)

    temp_random<-c("log_B","log_R")

    if (simul_data_obj$options_vec[2]==1){
      temp_par_list$log_sigma_m<-log(sigma_m)
      temp_par_list$log_S<-log(S)
      temp_par_list$log_input_m<-rep(log(0.1),length(unique(simul_data_obj$t_i))+1)

      temp_random[3]<-"log_input_m"
    }

    simul_obj<-list(simul_data=simul_data_obj,simul_par=temp_par_list,random=temp_random, map=list())

  } else stop("Function is for TLM, use simul_SEBDAM_params() for SEBDAM")

  return(simul_obj)

}
