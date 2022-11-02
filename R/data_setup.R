#' Organize a dataframe containing the individual observations to work with TMB model
#'
#' @param data list: containing individual observations: I for commercial size, IR for recruit size, Year for year, if spatial: geometry to show object is an sf object, optional: L for clappers, N for total numbers.
#' @param growths dataframe: one column for commercial size growth rate g, one for recruit growth rate gR
#' @param catch if non-spatial: vector of yearly catches, if spatial: data frame of catches with knots as rows and years as columns
#' @param model character: choice of model, options: SEBDAM or TLM
#' @param mesh INLA mesh: used for spatial approach, if non-spatial leave null, obtained from setup_mesh function
#' @param bound sf object: used to get triangles for barrier model, if not barrier model leave null
#' @param obs_mort boolean: true if using observations of mortality (e.g. clappers), false if fixing mortality
#' @param prior boolean: true if using a prior distribution to inform q_I, false to freely estimate
#' @param prior_pars numeric: parameters for beta distribution used as a prior, default set to 10 and 12
#' @param fix_m double: value used for fixing mortality when not using any observations, ignore if using observations
#' @param mult_qI boolean: if not using spatial approach, ignore, if using spatial: true when using a separate catchability at each knot, false for using a single overall catchability
#' @param spat_approach character: if using TLM, ignore, if using SEBDAM, choice of spatial model: "spde" for normal SPDE, "spde_aniso" for SPDE approach with geometric anisotropy, and "barrier" for barrier model
#' @param knot_obj kmeans: objects that set knot for spatial approach, can be ignored if using TLM, obtained from setup_mesh function
#' @param knot_area dataframe: obtained from setup_pred_grid function, used to know area used for each knot in spatial approach, ignore if using TLM
#' @param separate_R_aniso boolean: if using SPDE with geometric anisotropy, set to F is want biomass and recruitment to have same anisotropy
#' @param all_se boolean: if using SEBDAM, FALSE to speed up fitting process but not get all standard errors for random effects, TRUE to get all standard errors for random effects but slow down fitting process
#'
#' @return Properly formatted lists for model fitting
#' @export
#'
#' @examples
data_setup<-function(data=NULL, growths=NULL, catch=NULL, model=NULL, mesh=NULL, bound=NULL,
                     obs_mort=FALSE, prior=FALSE, prior_pars=c(10,12), fix_m=0.1,
                     mult_qI=FALSE, s_a=NULL, spat_approach=NULL, knot_obj=NULL,
                     knot_area=NULL,separate_R_aniso=T,all_se=F) {

  #Data tables do not always work the same way as data.frames, so ensure they are data.frames
  if (is.null(data) | is.null(growths)) {
    stop("Needs data")
  } else if(attr(data,"class")[1] != "data.frame") {
    data <- as.data.frame(data)
  }

  #Check that model is a string and appropriate model choice
  if (!(model %in% c("SEBDAM","TLM"))) {
    stop("Incorrect model choice, options are currently SEBDAM or TLM")
  }

  if (model == "TLM") {

    #Check
    if (obs_mort == TRUE) {
      req_names<-c("I","IR","L","N","Year")
    } else if (obs_mort == FALSE) {
      req_names<-c("I","IR","Year")
    }

    #Verify that data contains the required columns
    if (!all(req_names %in% names(data))) {
      miss_names<-req_names[which(!(req_names %in% names(data)))]
      stop(paste(c("Missing data columns:",miss_names),collapse=" "))
    }

    #Check if data elements are right lengths and numbers
    if (!is.numeric(data$I) | !is.numeric(data$IR) | !is.numeric(data$Year) | !is.numeric(growths$g) | !is.numeric(growths$gR)) {
      stop("I, IR, g, and gR and Year need to be numeric values")
    }

    if (!is.numeric(catch) & !(length(catch) %in% c(length(unique(data$Year)),length(unique(data$Year))+1))) {
      stop("Catch is either not a vector or has a different length than data")
    }

    if (ncol(growths)!=2 | !(any(colnames(growths) %in% c("g"))) | !(any(colnames(growths) %in% c("gR")))){
      stop("Issue with growths data frame, needs to have 2 columns with names g and gR")
    }

    temp_data_list<-list()

    temp_data_list$model<-"TLM"

    temp_data_list$options_vec<-c(0,0)
    if (prior == TRUE) {
      temp_data_list$options_vec[1]<-1
      temp_data_list$prior_pars<-prior_pars
    }
    if (obs_mort == TRUE) temp_data_list$options_vec[2]<-1

    #For adult/commercial size animals, can be directly used
    temp_data_list$logI<-data$I
    temp_data_list$logI[temp_data_list$logI==0]<-NA
    temp_data_list$logI<-log(temp_data_list$logI)

    #For recruits, can be directly used
    temp_data_list$logIR<-data$IR
    temp_data_list$logIR[temp_data_list$logIR==0]<-NA
    temp_data_list$logIR<-log(temp_data_list$logIR)

    #Catch
    temp_data_list$C<-as.matrix(catch)

    #Growth rates
    temp_data_list$g<-growths$g
    temp_data_list$gR<-growths$gR

    temp_tows<-c()
    for (i in unique(data$Year)){
      temp<-subset(data,Year==i)
      temp_tows<-c(temp_tows,length(temp$Year))
    }
    temp_data_list$n_tows<-temp_tows

    #Positive Tows:
    non_zeroes_I<-c()
    non_zeroes_IR<-c()
    for(i in unique(data$Year)) {
      non_zeroes_I[i] <- length(which(subset(data,Year==(i))$I!=0))
      non_zeroes_IR[i] <- length(which(subset(data,Year==(i))$IR!=0))
    }
    temp_data_list$pos_tows_I<-non_zeroes_I
    temp_data_list$pos_tows_IR<-non_zeroes_IR

    temp_data_list$t_i<-(data$Year-min(data$Year))

    temp_data_list$plug_exploit<-0.1

    if (obs_mort == FALSE) temp_data_list$set_m<-fix_m

    temp_par_list<-list()
    temp_par_list$log_sigma_tau<--1
    temp_par_list$log_sigma_phi<--1
    temp_par_list$log_sigma_epsilon<--1
    temp_par_list$log_sigma_upsilon<--1
    temp_par_list$log_q_R<--1
    temp_par_list$log_q_I<--1
    temp_par_list$logit_p_I<-logit(0.5)
    temp_par_list$logit_p_IR<-logit(0.5)

    temp_par_list$log_B<-rep(log(max(data$I,na.rm=T)*10),length(unique(data$Year))+1)
    temp_par_list$log_R<-rep(log(max(data$IR,na.rm=T)*10),length(unique(data$Year))+1)

    temp_random<-c("log_B","log_R")

    if (obs_mort == TRUE) {
      if (!is.numeric(data$L) | !is.numeric(data$N)) stop("L and N need to be numeric values")
      #Round behaves weirdly, so have to check if even, just in case the values given are standardized numbers and not integers
      temp_bin<-data$L+data$N
      temp_L<-data$L
      for (i in 1:length(temp_bin)) {
        if (temp_bin[i] %% 2 >= 1) temp_bin[i]<-round(temp_bin[i])
        ##If number odd, then add 0.01 to avoid rounding down at exactly .5
        else temp_bin[i]<-round(temp_bin[i]+0.01)

        if (temp_L[i] %% 2 >= 1) temp_L[i]<-round(temp_L[i])
        else temp_L[i]<-round(temp_L[i]+0.01)

      }
      temp_bin[which(temp_bin==0)]<-NA
      temp_L[which(is.na(temp_bin))]<-NA
      temp_data_list$L<-temp_L
      temp_data_list$n_bin<-temp_bin
      temp_data_list$n_bin<-temp_bin
      temp_data_list$L<-temp_L

      temp_par_list$log_sigma_m<--1
      temp_par_list$log_S<--1
      temp_par_list$log_input_m<-rep(log(0.3),length(unique(data$Year))+1)

      temp_random[3]<-"log_input_m"
    }

    tmb_obj<-list(data=temp_data_list,par=temp_par_list,random=temp_random,map=list())


    #Checking that the output of the data will work (e.g. sizes of everything works)
    if (obs_mort==T) {
      if (!all.equal(length(temp_data_list$logI),length(temp_data_list$logIR),length(temp_data_list$n_bin),length(temp_data_list$L),length(temp_data_list$t_i),length(temp_data_list$s_i))) {
        stop("One of I, IR, L, N, or indicators s_i or t_i is not the correct length, check data")
      }
    } else {
      if (!all.equal(length(temp_data_list$logI),length(temp_data_list$logIR),length(temp_data_list$t_i),length(temp_data_list$s_i))) {
        stop("One of I, IR, or indicators s_i or t_i is not the correct length, check data")
      }
    }

    if (!all.equal(length(temp_data_list$gI),length(temp_data_list$gR),temp_data_list$n_t,length(temp_data_list$pos_tows_I),length(temp_data_list$pos_tow_IR),length(temp_data_list$n_tows))) {
      stop("Growth rates do not match up with number of years, check data")
    }

  } else if (model == "SEBDAM") { #Choosing SEBDAM for spatial approach

    #Check that the mesh is an appropriate INLA mesh
    if (class(mesh)!="inla.mesh") {
      stop("Mesh is not an INLA mesh")
    }

    if (!(spat_approach %in% c("spde","spde_aniso","barrier"))) {
      stop("Incorrect spatial approach choice, options are currently spde, spde_aniso or barrier")
    }

    if (obs_mort == TRUE) {
      req_names<-c("I","IR","L","N","Year","geometry")
    } else if (obs_mort == FALSE) {
      req_names<-c("I","IR","Year","geometry")
    }

    #Verify that data contains the required columns
    if (!all(req_names %in% names(data))) {
      miss_names<-req_names[which(!(req_names %in% names(data)))]
      stop(paste(c("Missing data columns:",miss_names),collapse=" "))
    }

    #Check if data elements are right lengths and numbers
    if (!is.numeric(data$I) | !is.numeric(data$IR) | !is.numeric(data$Year) | !is.numeric(growths$g) | !is.numeric(growths$gR)) {
      stop("I, IR, g, and gR and Year need to be numeric values")
    }

    if (attr(catch,"class")[1]!="data.frame" | dim(catch)[1]!=length(mesh$idx$loc) | !(dim(catch)[2] %in% c(length(unique(data$Year)),length(unique(data$Year))+1))) {
      stop("Catch is either not a data frame or has incompatible dimensions")
    }

    #Setting up the data as required by SEBDAM
    temp_data_list<-list()

    temp_data_list$model<-"SEBDAM"

    temp_data_list$options_vec<-c(0,0,0,0,0,0,0)
    if (prior == TRUE) {
      temp_data_list$options_vec[2]<-1
      temp_data_list$prior_pars<-prior_pars
    }
    if (all_se == T) temp_data_list$options_vec[3]<-1
    if (obs_mort == TRUE) temp_data_list$options_vec[4]<-1
    if (mult_qI == TRUE) temp_data_list$options_vec[5]<-1
    if (spat_approach == "spde_aniso") {
      temp_data_list$options_vec[6]<-1
      if (separate_R_aniso==T) temp_data_list$options_vec[7]<-1
    }
    else if (spat_approach == "barrier") temp_data_list$options_vec[6]<-2

    #For adult/commercial size animals, can be directly used
    temp_data_list$logI<-data$I
    temp_data_list$logI[temp_data_list$logI==0]<-NA
    temp_data_list$logI<-log(temp_data_list$logI)

    #For recruits, can be directly used
    temp_data_list$logIR<-data$IR
    temp_data_list$logIR[temp_data_list$logIR==0]<-NA
    temp_data_list$logIR<-log(temp_data_list$logIR)

    temp_data_list$area<-knot_area$area

    temp_data_list$C<-as.matrix(catch/knot_area$area)

    temp_tows<-c()
    for (i in unique(data$Year)){
      temp<-subset(data,Year==i)
      temp_tows<-c(temp_tows,length(temp$Year))
    }
    temp_data_list$n_tows<-temp_tows

    #Positive Tows:
    non_zeroes_I<-c()
    non_zeroes_IR<-c()
    for(i in unique(data$Year)) {
      non_zeroes_I[i] <- length(which(subset(data,Year==(i))$I!=0))
      non_zeroes_IR[i] <- length(which(subset(data,Year==(i))$IR!=0))
    }
    temp_data_list$pos_tows_I<-non_zeroes_I
    temp_data_list$pos_tows_IR<-non_zeroes_IR

    temp_data_list$n_i<-length(temp_data_list$logI)
    temp_data_list$n_t<-length(unique(data$Year))
    temp_data_list$n_s<-length(unique(knot_area$knotID))
    temp_data_list$n_m<-mesh$n

    temp_data_list$s_i<-unname(knot_obj$cluster-min(knot_obj$cluster))
    temp_data_list$t_i<-(data$Year-min(data$Year))
    temp_data_list$v_i<-mesh$idx$loc-1

    if (is.data.frame(growths)){
      temp_data_list$gI<-t(as.matrix(growths$g,ncol=1))
      temp_data_list$gR<-t(as.matrix(growths$gR,ncol=1))
      temp_data_list$s_a<-rep(0,temp_data_list$n_s)
    } else if (is.list(growths) & (nrow(growths$gI)==nrow(growths$gR)) & (ncol(growths$gI)>=temp_data_list$n_t) & (ncol(growths$gR)>=temp_data_list$n_t)){
      if (length(s_a)!=temp_data_list$n_s){
        stop("Issue with s_a")
      }
      temp_data_list$gI<-as.matrix(growths$gI)
      temp_data_list$gR<-as.matrix(growths$gR)
      temp_data_list$s_a<-s_a
    } else (stop("Issue with growth data"))

    temp_data_list$plug_exploit<-0.1

    temp_par_list<-list()
    temp_par_list$log_sigma_epsilon<--1
    temp_par_list$log_sigma_upsilon<--1
    temp_par_list$log_R0<-mean(log(exp(temp_data_list$logIR)*10),na.rm=T)
    temp_par_list$log_B0<-mean(log(exp(temp_data_list$logI)*10),na.rm=T)
    temp_par_list$log_m0<-log(fix_m)
    temp_par_list$logit_p_I<-logit(0.5)
    temp_par_list$logit_p_IR<-logit(0.5)
    temp_par_list$log_qR<--1

    temp_par_list$omega_B<-matrix(rep(0,(temp_data_list$n_t+1)*temp_data_list$n_m),ncol=(temp_data_list$n_t+1))
    temp_par_list$omega_R<-matrix(rep(0,temp_data_list$n_t*temp_data_list$n_m),ncol=temp_data_list$n_t)

    if (spat_approach %in% c("spde","spde_aniso")){
      temp_par_list$log_kappa_B<--1
      temp_par_list$log_tau_B<-1
      temp_par_list$log_kappa_R<--1
      temp_par_list$log_tau_R<-1
      if (spat_approach == "spde") temp_data_list$mesh_obj<-(INLA::inla.spde2.matern(mesh))$param.inla[c("M0","M1","M2")]
      else if (spat_approach == "spde_aniso") {
        temp_data_list$mesh_obj<-get_aniso_obj(mesh)
        temp_par_list$log_H_input_B<-c(0,0)
        temp_par_list$log_H_input_R<-c(0,0)
      }
    } else if (spat_approach == "barrier") {
        warning("INLA function for barrier model can give warning message coming from their use of Matrix::sparseMatrix function depending on Matrix version. If it sets repr=T, can be safely ignored")
        temp_data_list$mesh_obj<-INLA::inla.barrier.fem(mesh,barrier.triangles = get_triangles(mesh,bound))
        temp_par_list$log_range_B<-c(log(40),log(2))
        temp_par_list$log_sigma_B<--1
        temp_par_list$log_range_R<-c(log(40),log(2))
        temp_par_list$log_sigma_R<--1
    }

    temp_random<-c("omega_B","omega_R")
    temp_map<-list()

    if (mult_qI == TRUE) temp_par_list$log_qI<-rep(-1,temp_data_list$n_s)
    else if (mult_qI == FALSE) temp_par_list$log_qI<--1

    if (obs_mort == TRUE){
      temp_par_list$omega_m<-matrix(rep(0,(temp_data_list$n_t+1)*temp_data_list$n_m),ncol=(temp_data_list$n_t+1))
      if (spat_approach %in% c("spde","spde_aniso")){
        temp_par_list$log_kappa_m<-0
        temp_par_list$log_tau_m<-0
        if (spat_approach == "spde_aniso") temp_par_list$log_H_input_m<-c(0,0)
      } else if (spat_approach == "barrier"){
        temp_par_list$log_range_m<-c(log(40),log(2))
        temp_par_list$log_sigma_m<--1
      }
      temp_random[3]<-"omega_m"
      if (!is.numeric(data$L) | !is.numeric(data$N)) stop("L and N need to be numeric values")
      #Round behaves weirdly, so have to check if even, just in case the values given are standardized numbers and not integers
      temp_bin<-data$L+data$N
      temp_L<-data$L
      for (i in 1:length(temp_bin)) {

        if (temp_L[i] %% 2 >= 1) temp_L[i]<-round(temp_L[i])
        else temp_L[i]<-round(temp_L[i]+0.01)

        temp_data_list$L<-temp_L
        if (temp_bin[i] %% 2 >= 1) temp_bin[i]<-round(temp_bin[i])
        ##If number odd, then add 0.01 to avoid rounding down at exactly .5
        else temp_bin[i]<-round(temp_bin[i]+0.01)
      }
      temp_bin[which(temp_bin==0)]<-NA
      temp_L[which(is.na(temp_bin))]<-NA
      temp_data_list$L<-temp_L
      temp_data_list$n_bin<-temp_bin
      temp_par_list$log_S<--1
    } else if (obs_mort == FALSE) temp_map$log_m0<-as.factor(NA)
    if (spat_approach == "spde_aniso" & separate_R_aniso==FALSE) temp_map$log_H_input_R<-as.factor(c(NA,NA))

    tmb_obj<-list(data=temp_data_list,par=temp_par_list,random=temp_random,map=temp_map)

    #Checking that the output of the data will work (e.g. sizes of everything works)
    if (obs_mort==T) {
      if (!all.equal(length(temp_data_list$logI),length(temp_data_list$logIR),length(temp_data_list$n_bin),length(temp_data_list$L),length(temp_data_list$t_i),length(temp_data_list$s_i))) {
        stop("One of I, IR, L, N, or indicators s_i or t_i is not the correct length, check data")
      }
    } else {
      if (!all.equal(length(temp_data_list$logI),length(temp_data_list$logIR),length(temp_data_list$t_i),length(temp_data_list$s_i))) {
        stop("One of I, IR, or indicators s_i or t_i is not the correct length, check data")
      }
  }

    if (!all.equal(length(temp_data_list$gI),length(temp_data_list$gR),temp_data_list$n_t,length(temp_data_list$pos_tows_I),length(temp_data_list$pos_tow_IR),length(temp_data_list$n_tows))) {
      stop("Growth rates do not match up with number of years, check data")
    }

    if (!all.equal(length(temp_data_list$area),length(temp_data_list$C[,1]),temp_data_list$n_s,temp_data_list$v_i)) {
      stop("Knot area, catches and number of knots do not all match, check data")
    }

  }

  return(tmb_obj)

}

