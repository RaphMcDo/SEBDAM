simulate_data<-function(model=NULL, n_years=1, even_spread=TRUE,
                        simul_area_obj=NULL, obs_mort=FALSE,
                        mult_qI=FALSE,spat_approach="spde",
                        separate_R_aniso=TRUE) {

  if (model=="TLM") {



  }

  if (model=="SEBDAM") {

    if (is.null(simul_area_obj)) stop("SEBDAM requires running simul_area first and using its output")

    temp_data_list<-list()

    temp_data_list$model<-"SEBDAM"

    temp_data_list$options_vec<-c(0,0,0,0,0,0,0)
    if (obs_mort == TRUE) temp_data_list$options_vec[4]<-1
    if (mult_qI == TRUE) temp_data_list$options_vec[5]<-1
    if (spat_approach == "spde_aniso") {
      temp_data_list$options_vec[6]<-1
      if (separate_R_aniso==F) temp_data_list$options_vec[7]<-1
    }
    else if (spat_approach == "barrier") temp_data_list$options_vec[6]<-2

    temp_data_list$logI<-rep(1,length(simul_area_obj$sim_obs))
    temp_data_list$logIR<-rep(1,length(simul_area_obj$sim_obs))
    temp_data_list$area<-simul_area_obj$area$area

    temp_data_list$C<-matrix(rep(1,n_year*length(unique(simul_area_obj$knots$cluster))))

    temp_div<-temp_data_list$logI/n_years
    rand_vec<-rand_int_vect(min=round(temp_div)/2,max=round(temp_div)+0.5*round(temp_div),n_grp=n_years,total=length(temp_data_list$logI))
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
    temp_data_list$n_s<-length(unique(simul_area_obj$area$knotID))
    temp_data_list$n_m<-simul_area_obj$mesh$n

    temp_data_list$s_i<-unname(simul_area_obj$knots$cluster-min(simul_area_obj$knots$cluster))

    #
    if (even_spread==T) {
      if (length(temp_data_list$logI)/n_years %% 2 %in% c(0,1)){
        data$t_i<-rep(0:(n_years-1),each=temp_data_list$logI/n_years)
      } else {
        temp<-rep(0:(n_years-1),each=floor(temp_data_list$logI/n_years))
        data$t_I<-c(temp,0:(length(temp_data_list$logI)-length(data$t_i)))
      }
    } else if (even_spread==F) {
      data$t_I<-rep(0:(n_years-1),times=rand_vec)
    }

    temp_data_list$v_i<-simul_area_obj$mesh$idx$loc-1

    temp_data_list$gI<-rep(1.1,n_years)
    temp_data_list$gR<-rep(1.1,n_years)

    if (spat_approach %in% c("spde","spde_aniso")){
      if (spat_approach == "spde") temp_data_list$mesh_obj<-(INLA::inla.spde2.matern(mesh))$param.inla[c("M0","M1","M2")]
      else if (spat_approach == "spde_aniso") {
        temp_data_list$mesh_obj<-get_aniso_obj(mesh)
      }
    } else if (spat_approach == "barrier") {
      warning("INLA function for barrier model can give warning message coming from their use of Matrix::sparseMatrix function depending on Matrix version. If it sets repr=T, can be safely ignored")
      temp_data_list$mesh_obj<-INLA::inla.barrier.fem(mesh,barrier.triangles = get_triangles(mesh,bound))
    }

    if (obs_mort == TRUE){
      temp_data_list$L<-rep(1,length(temp_data_list$logI))
      temp_data_list$n_bin<-rep(1,length(temp_data_list$logI))
    }

    sim_data<-temp_data_list

  }

  return(sim_data)

}
