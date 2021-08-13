#' Create object used for simulations using SEBDAM or TLM
#'
#' @param model Character: Choose which model to use, options are TLM or SEBDAM
#' @param n_years Integer: Number of years desired
#' @param obs_mort Boolean: True to simulate clapper observations and processes, false to fix natural mortality
#' @param n_obs Integer: required for TLM and ignored for SEBDAM (as number is chosen in simulate_area function), choose total number of observations
#' @param even_spread Boolean: choice to spread observations evenly across years, or randomly split them up
#' @param fix_m Double: only matters for TLM, if obs_mort is FALSE, choice of value to fix the natural mortality (default of 0.1)
#' @param prior_q boolean: TRUE to use a prior distribution on q_I in later fitting processes
#' @param catch_spread Character: ignored if using TLM, choice of how to simulate catches, choices are prop, aggregated or aggregated_extra
#' @param simul_area_obj List: ignored if using TLM, object obtained from running simulate_area function for SEBDAM
#' @param mult_qI boolean: ignored if using TLM, choice of using a single q_I (FALSE) or one per knot (TRUE)
#' @param spat_approach character: ignored if using TLM, choice of spatial approach, options are spde, spde_aniso and barrier
#' @param separate_R_aniso boolean: ignored if using TLM or using spde or barrier spatial approaches, choice of having separate anisotropy parameters for R or the same as for B
#' @param bound sf object: boundaries of modelling only when using barrier model, must includes islands
#'
#' @return list: object used for simulations in simulate_data function
#' @export
#'
#' @examples
simulate_obj<-function(model=NULL, n_years=1L, obs_mort=FALSE,
                        n_obs=NULL,even_spread=TRUE, fix_m=0.1,
                        prior_q=FALSE, catch_spread=NULL,
                        simul_area_obj=NULL,
                        mult_qI=FALSE,spat_approach="spde",
                        separate_R_aniso=TRUE,bound=NULL) {

  if (!(model %in% c("TLM","SEBDAM"))) stop("Incorrect Model Choice")

  if (model=="TLM") {

    temp_data_list<-list()

    temp_data_list$model<-"TLM"

    temp_data_list$options_vec<-c(0,0)
    if (prior_q == TRUE) temp_data_list$options_vec[1]<-1
    if (obs_mort == TRUE) temp_data_list$options_vec[2]<-1

    temp_data_list$logI<-rep(1,n_obs)
    temp_data_list$logIR<-rep(1,n_obs)

    temp_data_list$C<-rep(1,n_years)

    temp_data_list$g<-rep(1,n_years)
    temp_data_list$gR<-rep(1,n_years)

    temp_div<-length(temp_data_list$logI)/n_years
    rand_vec<-rand_int_vect(min=round(temp_div/2),max=round(temp_div+0.5*temp_div),n_grp=n_years,total=length(temp_data_list$logI))
    if (even_spread==T){
      if (length(temp_data_list$logI)/n_years %% 2 %in% c(0,1)) {
        temp_data_list$n_tows<-rep(length(temp_data_list$logI)/n_years,n_years)
      } else {
        temp<-rep(length(temp_data_list$logI)/n_years,n_years-1)
        temp_data_list$n_tows<-c(temp,length(temp_data_list$logI)-sum(temp))
      }
    } else if (even_spread==F) {
      temp_data_list$n_tows<-rand_vec
    }

    temp_data_list$pos_tows_I<-rep(1,n_years)
    temp_data_list$pos_tows_IR<-rep(1,n_years)

    if (even_spread==T) {
      if (length(temp_data_list$logI)/n_years %% 2 %in% c(0,1)){
        temp_data_list$t_i<-rep(0:(n_years-1),each=length(temp_data_list$logI)/n_years)
      } else {
        temp<-rep(0:(n_years-1),each=floor(length(temp_data_list$logI)/n_years))
        temp_data_list$t_i<-c(temp,0:(length(temp_data_list$logI)-length(temp_data_list$t_i)))
      }
    } else if (even_spread==F) {
      temp_data_list$t_i<-rep(0:(n_years-1),times=rand_vec)
    }

    if (obs_mort == TRUE) {
      temp_data_list$L<-rep(1,length(temp_data_list$logI))
      temp_data_list$n_bin<-rep(1,length(temp_data_list$logI))
    } else if (obs_mort == FALSE) temp_data_list$set_m<-fix_m

    sim_data<-temp_data_list

  } else if (model=="SEBDAM") {

    if (is.null(simul_area_obj)) stop("SEBDAM requires running simul_area first and using its output")

    temp_data_list<-list()

    temp_data_list$model<-"SEBDAM"

    if (is.null(catch_spread)) warning("No appropriate specification of how to simulate catch, therefore none will be simulated")
    else if (!(catch_spread %in% c("prop","aggregated","aggregated_extra"))) warning("No appropriate specification of how to simulate catch, therefore none will be simulated")

    if (!(spat_approach %in% c("spde","spde_aniso","barrier"))) warning("Incorrect specification of spatial approach, spde used")

    temp_data_list$options_vec<-c(0,0,0,0,0,0,0)
    if (!is.null(catch_spread)){
      if (catch_spread == "prop") temp_data_list$options_vec[1]<-1
      else if (catch_spread == "aggregated") temp_data_list$options_vec[1]<-2
      else if (catch_spread == "aggregated_extra") temp_data_list$options_vec[1]<-3
    }
    if (prior_q == TRUE) temp_data_list$options_vec[2]<-1
    if (obs_mort == TRUE) temp_data_list$options_vec[4]<-1
    if (mult_qI == TRUE) temp_data_list$options_vec[5]<-1
    if (spat_approach == "spde_aniso") {
      temp_data_list$options_vec[6]<-1
      if (separate_R_aniso==T) temp_data_list$options_vec[7]<-1
    }
    else if (spat_approach == "barrier") temp_data_list$options_vec[6]<-2

    temp_data_list$logI<-rep(1,length(simul_area_obj$sim_obs))
    temp_data_list$logIR<-rep(1,length(simul_area_obj$sim_obs))
    temp_data_list$area<-simul_area_obj$stratarea$area

    temp_data_list$C<-matrix(rep(1,(n_years+1)*length(unique(simul_area_obj$knots$cluster))),ncol=(n_years+1))

    temp_div<-length(temp_data_list$logI)/n_years
    rand_vec<-rand_int_vect(min=round(temp_div/2),max=round(temp_div+0.5*temp_div),n_grp=n_years,total=length(temp_data_list$logI))
    if (even_spread==T){
      if (length(temp_data_list$logI)/n_years %% 2 %in% c(0,1)) {
        temp_data_list$n_tows<-rep(length(temp_data_list$logI)/n_years,n_years)
      } else {
        temp<-rep(length(temp_data_list$logI)/n_years,n_years-1)
        temp_data_list$n_tows<-c(temp,length(temp_data_list$logI)-sum(temp))
      }
    } else if (even_spread==F) {
      temp_data_list$n_tows<-rand_vec
    }

    temp_data_list$pos_tows_I<-rep(1,n_years)
    temp_data_list$pos_tows_IR<-rep(1,n_years)

    temp_data_list$n_i<-length(temp_data_list$logI)
    temp_data_list$n_t<-n_years
    temp_data_list$n_s<-length(unique(simul_area_obj$stratarea$knotID))
    temp_data_list$n_m<-simul_area_obj$mesh$n

    temp_data_list$s_i<-unname(simul_area_obj$knots$cluster-min(simul_area_obj$knots$cluster))

    #
    if (even_spread==T) {
      if (length(temp_data_list$logI)/n_years %% 2 %in% c(0,1)){
        temp_data_list$t_i<-rep(0:(n_years-1),each=length(temp_data_list$logI)/n_years)
      } else {
        temp<-rep(0:(n_years-1),each=floor(length(temp_data_list$logI)/n_years))
        temp_data_list$t_i<-c(temp,0:(length(temp_data_list$logI)-length(temp_data_list$t_i)))
      }
    } else if (even_spread==F) {
      temp_data_list$t_i<-rep(0:(n_years-1),times=rand_vec)
    }

    temp_data_list$v_i<-simul_area_obj$mesh$idx$loc-1

    temp_data_list$gI<-rep(1.1,n_years)
    temp_data_list$gR<-rep(1.1,n_years)

    if (spat_approach %in% c("spde","spde_aniso")){
      if (spat_approach == "spde") temp_data_list$mesh_obj<-(INLA::inla.spde2.matern(simul_area_obj$mesh))$param.inla[c("M0","M1","M2")]
      else if (spat_approach == "spde_aniso") {
        temp_data_list$mesh_obj<-get_aniso_obj(simul_area_obj$mesh)
      }
    } else if (spat_approach == "barrier") {
      warning("INLA function for barrier model can give warning message coming from their use of Matrix::sparseMatrix function depending on Matrix version. If it sets repr=T, can be safely ignored")
      temp_data_list$mesh_obj<-INLA::inla.barrier.fem(simul_area_obj$mesh,barrier.triangles = get_triangles(simul_area_obj$mesh,bound))
    }

    if (obs_mort == TRUE){
      temp_data_list$L<-rep(1,length(temp_data_list$logI))
      temp_data_list$n_bin<-rep(1,length(temp_data_list$logI))
    }

    sim_data<-temp_data_list

  }

  return(sim_data)

}
