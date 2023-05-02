#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sehbam(objective_function <Type>* obj) {
  // Read in data, parameters, random effects
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  using namespace barrier_q;

  DATA_IVECTOR(options_vec); //Vector of multiple choices
  //Slot 0: Choice of simulated landings
  //  0 (default)= No landing s
  //  1 = Proportionally distributed at every knot
  //  2 = Proportionally distributed at knots with above average biomass
  //  3 = Proportionally distributed at knots with above average biomass with extra landings at other knots
  //Slot 1: Choice of using prior on catchability
  //  0 (default) = No prior
  //  1 = Use prior on catchabilities
  //Slot 2: Choice of obtaining computationally expensive standard errors
  // 0 (default) = Obtain standard errors only for parameters and overall processes
  // 1 = Obtain standard errors for knot-specific densities
  //Slot 3: Choice of mortality approach, given that most fisheries don't have any observations for mortality
  // 0 (default) = Fix mortality to user-provided value
  // 1 = use spatial approach
  // 2 = Temporal random walk to natural mortality
  //Slot 4: choice of single or multiple commercial size catchabilities
  // 0 (default) = Single commercial size catchability
  // 1 = different catchability per knot
  // 2 = different catchability per habitat type
  //Slot 5: Choice of spatial approach:
  //0 (default): Standard SPDE approach
  //1: SPDE approach with geometric anisotropy
  //2: Barrier model
  //Slot 6: Separate anisotropy for recruitment or not
  //0 (default): same as for biomass
  //1: separate anisotropy parameters

  //Data
  DATA_VECTOR(logI); //commercial biomass survey index by tow (dim n_i)
  DATA_VECTOR(logIR); //recruit survey index by tow (dim n_i)
  DATA_VECTOR(area); //area covered by each knot (dim n_s)
  DATA_MATRIX(C); //commercial catch (dim n_s,n_t)
  DATA_VECTOR(n_I); //number of commercial size scallops captured in tow (dim n_i)
  DATA_VECTOR(n_IR); //number of recruit size scallops captured in tow (dim n_i)
  DATA_MATRIX(X); //covariate matrix for estimating probability of extra zeroes

  //For residuals
  DATA_VECTOR_INDICATOR(keep_I, logI);
  DATA_VECTOR_INDICATOR(keep_IR, logIR);

  //Sizes
  DATA_INTEGER(n_i); //Number of observations
  DATA_INTEGER(n_t); //Number of years
  DATA_INTEGER(n_s); //Number of knots
  DATA_INTEGER(n_m); //number of vertex in mesh
  DATA_INTEGER(n_h); //number of unique habitat types

  //Indices to associate individual observations and proper mesh locations
  DATA_FACTOR(s_i); //Indexing for knot (dim n_i)
  DATA_FACTOR(t_i); //Indexing for year (dim n_i)
  DATA_FACTOR(v_i); //Indexing for specific mesh location to match knots with right vertex (dim n_m)
  DATA_FACTOR(h_i); //Indexing for habitat type (dim n_i)
  DATA_FACTOR(s_a); //Indexing knot for area for growth rates

  //Fixed covariates
  DATA_MATRIX(gI); //commercial biomass growth
  DATA_MATRIX(gR); //recruitment growth

  DATA_SCALAR(plug_exploit); //Desired exploitation rate for simulations, not used when fitting data

  //Parameters
  PARAMETER(log_sigma_epsilon); //obs sd survey index commercial biomass
  PARAMETER(log_sigma_upsilon); //obs sd survey index recruits
  PARAMETER(log_R0);//initial recruit mean value
  PARAMETER(log_B0);//initial commercial biomass mean value
  PARAMETER(log_m0);//initial mortality mean value
  PARAMETER(mu_I);
  PARAMETER_VECTOR(mu_IR);
  PARAMETER(log_var_I);
  PARAMETER_VECTOR(log_var_IR);
  PARAMETER_VECTOR(log_qR); //recruit catchability
  PARAMETER_VECTOR(beta_I); //intercept and slopes for regression to obtain excess zero probabilities for commercial size

  //Transform parameters
  Type pi = M_PI;
  Type sigma_epsilon = exp(log_sigma_epsilon);
  Type sigma_upsilon = exp(log_sigma_upsilon);
  vector <Type> qR = exp(log_qR);
  Type R0 = exp(log_R0);
  Type B0 = exp(log_B0);
  Type m0 = exp(log_m0);
  vector <Type> logit_zip_I(n_i);
  vector <Type> logit_zip_IR(n_i);
  vector <Type> zip_I(n_i);
  Type var_I = exp(log_var_I);
  vector <Type> var_IR = exp(log_var_IR);

  //Spatial Random Effects
  PARAMETER_ARRAY(omega_B);
  PARAMETER_ARRAY(omega_R);

  //Subsetting to only report the predicted field values at the knots
  matrix <Type> sub_omega_B(n_m,n_t);
  matrix <Type> sub_omega_R(n_m,n_t);

  // Set up matrices for processes of interest
  matrix <Type> log_B(n_s,(n_t+1));
  matrix <Type> B(n_s,(n_t+1));
  matrix <Type> areaB(n_s,(n_t+1));
  matrix <Type> areaR(n_s,(n_t));
  matrix <Type> log_R(n_s,(n_t));
  matrix <Type> R(n_s,(n_t));
  matrix <Type> log_m(n_s,(n_t+1));
  matrix <Type> m(n_s,(n_t+1));
  //

  //Setup for simulations and derived values
  vector <Type> totB(n_t+1); totB.setZero();
  vector <Type> totR(n_t); totR.setZero();
  vector <Type> mean_m(n_t+1); mean_m.setZero();
  matrix <Type> mean_pro_B(n_s,n_t+1);

  // ----------------------------------------------
  // nll
  vector <Type> nll_comp(14); nll_comp.setZero();


  //Setup spatial approach
  if (options_vec[5] == 0 || options_vec[5] == 1){
    //SPDE spatial parameter
    PARAMETER(log_kappa_B);//commercial size range parameter
    PARAMETER(log_tau_B); //commercial size spatio-temporal variability parameter
    PARAMETER(log_kappa_R); //recruit range parameter
    PARAMETER(log_tau_R); //recruit spatio-temporal variability parameter

    Type kappa_B = exp(log_kappa_B);
    Type tau_B = exp(log_tau_B);
    Type kappa_R = exp(log_kappa_R);
    Type tau_R = exp(log_tau_R);

    // Calculate marginal field variances and ranges based on relation described in Lindgren et al., 2012
    Type SigmaO_B = 1 / (sqrt(4*pi)*exp(log_tau_B)*exp(log_kappa_B));
    Type SigmaO_R = 1 / (sqrt(4*pi)*exp(log_tau_R)*exp(log_kappa_R));
    Type Range_B = sqrt(8)/exp(log_kappa_B);
    Type Range_R = sqrt(8)/exp(log_kappa_R);

    if (options_vec[5] == 0) {
      //Mesh object
      DATA_STRUCT(mesh_obj,spde_t);//SPDE

      //Set up GMRF for commercial size biomass
      SparseMatrix <Type> Q_B = Q_spde(mesh_obj, kappa_B);
      for (int t = 0; t < (n_t+1); t++){
        nll_comp(0) += SCALE(GMRF(Q_B),1/tau_B)(omega_B.col(t));
      }

      //Set up GMRF for recruits
      SparseMatrix <Type> Q_R = Q_spde(mesh_obj, kappa_R);
      for (int t = 0; t < (n_t); t++){
        nll_comp(1) += SCALE(GMRF(Q_R),1/tau_R)(omega_R.col(t));
      }

      SIMULATE {
        SparseMatrix <Type> Q_B = Q_spde(mesh_obj, kappa_B);
        density::GMRF_t<Type> gmrf1(Q_B);
        for (int t = 0; t < (n_t+1); t++){
          vector <Type> temp_omega(n_m);
          SCALE(gmrf1,1/tau_B).simulate(temp_omega);
          omega_B.col(t)=temp_omega;
        }
        REPORT(Q_B);
        REPORT(omega_B);
      }

      SIMULATE {
        SparseMatrix <Type> Q_R = Q_spde(mesh_obj, kappa_R);
        density::GMRF_t<Type> gmrf2(Q_R);
        for (int t = 0; t < n_t; t++){
          vector <Type> temp_omega(n_m);
          SCALE(gmrf2,1/tau_R).simulate(temp_omega);
          omega_R.col(t)=temp_omega;
        }
        REPORT(Q_R);
        REPORT(omega_R);
      }

    } else if (options_vec[5] == 1) {
      //Mesh object
      DATA_STRUCT(mesh_obj,spde_aniso_t);//SPDE with anisotropy

      PARAMETER_VECTOR(log_H_input_B);
      PARAMETER_VECTOR(log_H_input_R);

      //Create anisotropy matrices
      matrix <Type> H_B(2,2); H_B = create_H(log_H_input_B);
      matrix <Type> H_R(2,2);
      if (options_vec[6] == 1){
        H_R = create_H(log_H_input_R);
      } else if (options_vec[6]==0) {H_R = H_B;}

      //Set up GMRF for commercial size biomass
      SparseMatrix <Type> Q_B = Q_spde(mesh_obj, kappa_B, H_B);
      for (int t = 0; t < (n_t+1); t++){
        nll_comp(0) += SCALE(GMRF(Q_B),1/tau_B)(omega_B.col(t));
      }

      //Set up GMRF for recruits
      SparseMatrix <Type> Q_R = Q_spde(mesh_obj, kappa_R, H_R);
      for (int t = 0; t < (n_t); t++){
        nll_comp(1) += SCALE(GMRF(Q_R),1/tau_R)(omega_R.col(t));
      }
      //Report anisotropy parameters
      REPORT(log_H_input_B);
      REPORT(log_H_input_R);
      REPORT(H_B);
      REPORT(H_R);

      SIMULATE {
        SparseMatrix <Type> Q_B = Q_spde(mesh_obj, kappa_B, H_B);
        density::GMRF_t<Type> gmrf1(Q_B);
        for (int t = 0; t < (n_t+1); t++){
          vector <Type> temp_omega(n_m);
          SCALE(gmrf1,1/tau_B).simulate(temp_omega);
          omega_B.col(t)=temp_omega;
        }
        REPORT(Q_B);
        REPORT(H_B);
        REPORT(omega_B);
      }

      SIMULATE {
        SparseMatrix <Type> Q_R = Q_spde(mesh_obj, kappa_R, H_R);
        density::GMRF_t<Type> gmrf2(Q_R);
        for (int t = 0; t < n_t; t++){
          vector <Type> temp_omega(n_m);
          SCALE(gmrf2,1/tau_R).simulate(temp_omega);
          omega_R.col(t)=temp_omega;
        }
        REPORT(Q_R);
        REPORT(H_R);
        REPORT(omega_R);
      }
    }
    //Report SPDE parameters
    REPORT(kappa_B);
    REPORT(tau_B);
    REPORT(kappa_R);
    REPORT(tau_R);
    //REPORT derived parameters for SPDE approach
    REPORT(SigmaO_B);
    REPORT(SigmaO_R);
    REPORT(Range_B);
    REPORT(Range_R);
    ADREPORT(kappa_B);
    ADREPORT(tau_B);
    ADREPORT(kappa_R);
    ADREPORT(tau_R);
    //ADREPORT derived parameters for SPDE approach to get SEs
    ADREPORT(SigmaO_B);
    ADREPORT(SigmaO_R);
    ADREPORT(Range_B);
    ADREPORT(Range_R);

  } else if (options_vec[5] == 2){
    //Barrier object
    DATA_STRUCT(mesh_obj, fem_barrier_t); //INLA FEM object

    PARAMETER_VECTOR(log_range_B);//commercial size range parameter
    PARAMETER(log_sigma_B); //commercial size spatio-temporal variability parameter
    PARAMETER_VECTOR(log_range_R); //recruit range parameter
    PARAMETER(log_sigma_R); //recruit spatio-temporal variability parameter

    vector <Type> range_B; range_B = exp(log_range_B);
    vector <Type> range_R; range_R = exp(log_range_R);
    Type sigma_B = exp(log_sigma_B);
    Type sigma_R = exp(log_sigma_R);

    //Set up GMRF for commercial size biomass
    SparseMatrix <Type> Q_B = Q_barrier(mesh_obj, range_B, sigma_B);
    for (int t = 0; t < (n_t+1); t++){
      nll_comp(0) += GMRF(Q_B)(omega_B.col(t));
    }

    //Set up GMRF for recruits
    SparseMatrix <Type> Q_R = Q_barrier(mesh_obj, range_R, sigma_R);
    for (int t = 0; t < (n_t); t++){
      nll_comp(1) += GMRF(Q_R)(omega_R.col(t));
    }

    REPORT(range_B);
    REPORT(range_R);
    REPORT(sigma_B);
    REPORT(sigma_R);
    ADREPORT(range_B);
    ADREPORT(range_R);
    ADREPORT(sigma_B);
    ADREPORT(sigma_R);

    SIMULATE {
      SparseMatrix <Type> Q_B = Q_barrier(mesh_obj, range_B, sigma_B);
      density::GMRF_t<Type> gmrf1(Q_B);
      for (int t = 0; t < (n_t+1); t++){
        vector <Type> temp_omega(n_m);
        gmrf1.simulate(temp_omega);

        omega_B.col(t)=temp_omega;
      }
      REPORT(Q_B);
      REPORT(omega_B);
    }

    SIMULATE {
      SparseMatrix <Type> Q_R = Q_barrier(mesh_obj, range_R, sigma_R);
      density::GMRF_t<Type> gmrf2(Q_R);
      for (int t = 0; t < (n_t); t++){
        vector <Type> temp_omega(n_m);
        gmrf2.simulate(temp_omega);
        omega_R.col(t)=temp_omega;
      }
      REPORT(Q_R);
      REPORT(omega_R);
    }

  }

  //Subsetting the fields:
  for (int s=0; s<n_s;s++){
    for (int t=0;s<n_t;t++){
      sub_omega_B(s,t) = omega_B(v_i(s),t);
      sub_omega_R(s,t) = omega_R(v_i(s),t);
    }
  }

  if (options_vec[3] == 0){
    //Derived values
    //For mortality
    for (int s = 0; s < n_s; s++){
      log_m(s,0) = log(m0);
      m(s,0) = exp(log_m(s,0));
      for (int t = 1; t < (n_t); t++){
        log_m(s,t) = log(m0);
        m(s,t) = exp(log_m(s,t));
      }
      //Project 1 year ahead
      log_m(s,n_t) = log(m0);
      m(s,n_t) = exp(log_m(s,n_t));
    }

  } else if (options_vec[3] == 2){

    PARAMETER(log_sigma_m);
    Type sigma_m = exp(log_sigma_m);
    PARAMETER_VECTOR(log_temporal_m);
    vector<Type>temporal_m = exp(log_temporal_m);

    nll_comp(2) -= dnorm(log_temporal_m(0),log_m0-(sqr(sigma_m)/Type(2.0)),sigma_m,true);
    for (int t = 1; t < n_t; t++){
      nll_comp(2) -= dnorm(log_temporal_m(t),log_temporal_m(t-1)-(sqr(sigma_m)/Type(2.0)),sigma_m,true);
    }

    SIMULATE{
      log_temporal_m(0) = log_m0 + rnorm(Type(0.0),sigma_m);
      temporal_m(0) = exp(log_temporal_m(0));
      for (int t = 1; t < n_t; t++){
        log_temporal_m(t) = log_temporal_m(t-1)+rnorm(Type(0.0),sigma_m);
        temporal_m(t) = exp(log_temporal_m(t));
      }
      REPORT(log_temporal_m);
      REPORT(temporal_m);
    }

    //So as to not change the other equations
    for (int s = 0; s < n_s; s++){
      for (int t = 0; t < (n_t); t++){
        log_m(s,t) = log_temporal_m(t);
        m(s,t) = exp(log_m(s,t));
      }
      //Project 1 year ahead
      log_m(s,n_t) = log_temporal_m(n_t-1);
      m(s,n_t) = exp(log_m(s,n_t));
    }

  } else if (options_vec[3] == 1){

    PARAMETER_ARRAY(omega_m);
    matrix <Type> sub_omega_m(n_m,n_t);

    if (options_vec[5] == 0 || options_vec[5] == 1){
      //Load in parameters
      PARAMETER(log_kappa_m); //mortality range parameter
      PARAMETER(log_tau_m); //mortality spatio-temporal variability parameter

      //Transform parameter
      Type kappa_m = exp(log_kappa_m);
      Type tau_m = exp(log_tau_m);
      Type SigmaO_m = 1 / (sqrt(4*pi)*exp(log_tau_m)*exp(log_kappa_m));
      Type Range_m = sqrt(8)/exp(log_kappa_m);

      if (options_vec[5] == 0) {
        //SPDE
        //Mesh object
        DATA_STRUCT(mesh_obj,spde_t);//SPDE

        //Set up GMRF for natural mortality
        SparseMatrix <Type> Q_m = Q_spde(mesh_obj, kappa_m);
        for (int t = 0; t < (n_t+1); t++){
          nll_comp(0) += SCALE(GMRF(Q_m),1/tau_m)(omega_m.col(t));
        }

        SIMULATE {
          SparseMatrix <Type> Q_m = Q_spde(mesh_obj, kappa_m);
          density::GMRF_t<Type> gmrf3(Q_m);
          for (int t = 0; t < (n_t+1); t++){
            vector <Type> temp_omega(n_m);
            SCALE(gmrf3,1/tau_m).simulate(temp_omega);
            omega_m.col(t)=temp_omega;
          }
          REPORT(Q_m);
          REPORT(omega_m);
        }

      } else if (options_vec[5] == 1) {
        //SPDE with Anisotropy
        //Mesh object
        DATA_STRUCT(mesh_obj,spde_aniso_t);//SPDE with anisotropy

        PARAMETER_VECTOR(log_H_input_m); //mortality anisotropy parameters

        //mortality anisotropy matrix
        matrix<Type> H_m(2,2); H_m = create_H(log_H_input_m);

        //Set up GMRF for mortality
        SparseMatrix <Type> Q_m = Q_spde(mesh_obj, kappa_m, H_m);
        for (int t = 0; t < (n_t+1); t++){
          nll_comp(2) += SCALE(GMRF(Q_m),1/tau_m)(omega_m.col(t));
        }

        SIMULATE {
          SparseMatrix <Type> Q_m = Q_spde(mesh_obj, kappa_m, H_m);
          density::GMRF_t<Type> gmrf3(Q_m);
          for (int t = 0; t < (n_t+1); t++){
            vector <Type> temp_omega(n_m);
            SCALE(gmrf3,1/tau_m).simulate(temp_omega);
            omega_m.col(t)=temp_omega;
          }
          REPORT(Q_m);
          REPORT(H_m);
          REPORT(omega_m);
        }

        REPORT(log_H_input_m);
        REPORT(H_m);
      }

      REPORT(kappa_m);
      REPORT(tau_m);
      REPORT(SigmaO_m);
      REPORT(Range_m);
      ADREPORT(kappa_m);
      ADREPORT(tau_m);
      ADREPORT(SigmaO_m);
      ADREPORT(Range_m);

    } else if (options_vec[5] == 2) {
      //Barrier object
      DATA_STRUCT(mesh_obj, fem_barrier_t); //INLA FEM object

      PARAMETER_VECTOR(log_range_m); //mortality range parameter
      PARAMETER(log_sigma_m); //mortality spatio-temporal variability parameter

      vector <Type> range_m; range_m = exp(log_range_m);
      Type sigma_m = exp(log_sigma_m);

      //Set up GMRF for mortality
      SparseMatrix <Type> Q_m = Q_barrier(mesh_obj, range_m, sigma_m);
      for (int t = 0; t < (n_t+1); t++){
        nll_comp(2) += GMRF(Q_m)(omega_m.col(t));
      }

      REPORT(range_m);
      REPORT(sigma_m);
      ADREPORT(range_m);
      ADREPORT(sigma_m);

      SIMULATE {
        SparseMatrix <Type> Q_m = Q_barrier(mesh_obj, range_m, sigma_m);
        density::GMRF_t<Type> gmrf3(Q_m);
        for (int t = 0; t < (n_t+1); t++){
          vector <Type> temp_omega(n_m);
          gmrf3.simulate(temp_omega);
          omega_m.col(t)=temp_omega;
        }
        REPORT(Q_m);
        REPORT(omega_m);
      }
    }

    //Derived values
    //For mortality
    for (int s = 0; s < n_s; s++){
      log_m(s,0) = log(m0 * exp(omega_m(v_i(s),0)));
      m(s,0) = exp(log_m(s,0));
      for (int t = 1; t < (n_t); t++){
        log_m(s,t) = log( m(s,t-1) * exp(omega_m(v_i(s),t)) );
        m(s,t) = exp(log_m(s,t));
      }
      //Project 1 year ahead
      log_m(s,n_t) = log( m(s,n_t-1)*exp(omega_m(v_i(s),n_t)));
      m(s,n_t) = exp(log_m(s,n_t));
    }

    REPORT(omega_m);
    for (int s=0; s<n_s;s++){
      for (int t=0;s<n_t;t++){
        sub_omega_m(s,t) = omega_m(v_i(s),t);
      }
    }

  }

  //Recruit derivation
  for (int s = 0; s < n_s; s++){
    log_R(s,0) = log(R0 * exp(omega_R(v_i(s),0)));
    R(s,0) = exp(log_R(s,0));
    for (int t = 1; t < (n_t); t++){
      log_R(s,t) = log(R(s,t-1)*exp(omega_R(v_i(s),t)));
      R(s,t) = exp(log_R(s,t));
    }
  }

  //Biomass derivation
  for (int s = 0; s < n_s; s++) {
    log_B(s,0) = log(B0*exp(omega_B(v_i(s),0)));
    B(s,0) = exp(log_B(s,0));
    for (int t = 1; t < (n_t); t++) {
      Type mean_pro = (exp(-m(s,t))*gI(s_a(s),t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(s_a(s),t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
      log_B(s,t) = log(mean_pro);
      B(s,t) = exp(log_B(s,t));
    }
    //Project 1 year ahead
    Type mean_pro = (exp(-m(s,n_t))*gI(s_a(s),n_t-1)*(B(s, (n_t - 1)) - C(s,(n_t-1))) + exp(-m(s,n_t))*gR(s_a(s),n_t-1)*R(s,(n_t-1)))*exp(omega_B(v_i(s),n_t));
    log_B(s,n_t) = log(mean_pro);
    B(s,n_t) = exp(log_B(s,n_t));
  }

  //Simulate biomass and commercial catch
  SIMULATE {
    //Simulate biomass based on having 0 landings
    if (options_vec[0]==0) {
      vector <Type> counter_B(n_t); counter_B.setZero();//Calculate biomass in each year to simulate landings
      for (int t = 0; t < n_t; t++){
        for (int s = 0; s < n_s; s++){
          if (t == 0){
            mean_pro_B(s,t) = B0 * exp(omega_B(v_i(s),t));
            // mean_pro_B(s,t) = init_B(s) * exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
          else {
            C(s,t-1) = Type(0.0);
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(s_a(s),t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(s_a(s),t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
        }
      }
      //Simulating 1-year projection
      for (int s = 0; s < n_s; s++){
        C(s,n_t-1) = Type(0.0);
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(s_a(s),n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(s_a(s),n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
        B(s,n_t) = mean_pro_B(s,n_t);
        log_B(s,n_t) = log(B(s,n_t));
        C(s,n_t) = 0;
      }
    }
    //Simulating biomass while having landings distributed at each knot proportionally to biomass
    if (options_vec[0]==1) {
      vector <Type> counter_B(n_t); counter_B.setZero();//Calculate biomass in each year to simulate landings
      for (int t = 0; t < n_t; t++){
        for (int s = 0; s < n_s; s++){
          if (t == 0){
            mean_pro_B(s,t) = B0 * exp(omega_B(v_i(s),t));
            // mean_pro_B(s,t) = init_B(s) * exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
          else {
            //Calculate proportion of total biomass at knot s
            Type prop_B; prop_B = (B(s,t-1)*area(s))/counter_B(t-1);
            //Simulate total landings from total biomass on log scale
            Type exploitation =  exp(log(((counter_B(t-1)/1000)*plug_exploit)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
            //Simulate landings per km^2 at knot s
            C(s,t-1) = exp(log((exploitation*prop_B)/area(s)) + rnorm(Type(0.0),Type(0.2)));
            //Simulate biomass density
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(s_a(s),t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(s_a(s),t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
        }
      }
      for (int s = 0; s < n_s; s++){
        //Calculate proportion of total biomass at knot s
        Type prop_B; prop_B = (B(s,n_t-1)*area(s))/counter_B(n_t-1);
        //Simulate total landings from total biomass on log scale
        Type exploitation = exp(log(((counter_B(n_t-1)/1000)*plug_exploit)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
        //Simulate landings per km^2 at knot s
        C(s,n_t-1) = exp(log((exploitation*prop_B)/area(s)) + rnorm(Type(0.0),Type(0.2)));
        //Simulate biomass density
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(s_a(s),n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(s_a(s),n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
        B(s,n_t) = mean_pro_B(s,n_t);
        log_B(s,n_t) = log(B(s,n_t));
        C(s,n_t) = 0;
      }
    }
    //Simulating biomass based on having landings aggregated in above average knots
    if (options_vec[0]==2) {
      vector <Type> counter_B(n_t); counter_B.setZero();
      for (int t = 0; t < n_t; t++){
        for (int s = 0; s < n_s; s++){
          if (t == 0){
            mean_pro_B(s,t) = B0 * exp(omega_B(v_i(s),t));
            // mean_pro_B(s,t) = init_B(s) * exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
          else {
            //Calculate total biomass only in above average knots
            Type sum_pop; sum_pop = 0;
            for (int i = 0; i < n_s; i++){
              if ((B(i,t-1)*area(i)) >  (counter_B(t-1)/n_s)) {
                sum_pop = sum_pop + (B(i,t-1)*area(i));
              }
            }
            //Proportion of above average biomass at knot
            Type prop_B; prop_B = (B(s,t-1)*area(s))/sum_pop;
            Type exploitation =  exp(log(((counter_B(t-1)/1000)*plug_exploit)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
            if ( (B(s,t-1)*area(s)) > (counter_B(t-1)/n_s) ) C(s,t-1) = (exp(log((exploitation/area(s))*prop_B))+ rnorm(Type(0.0),Type(0.2)));
            else C(s,t-1) = Type(0.0);
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(s_a(s),t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(s_a(s),t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
        }
      }
      for (int s = 0; s < n_s; s++){
        Type sum_pop; sum_pop = 0;
        for (int i = 0; i < n_s; i++){
          if ((B(i,n_t-1)*area(i)) >  (counter_B(n_t-1)/n_s)) {
            sum_pop = sum_pop + (B(i,n_t-1)*area(i));
          }
        }
        Type prop_B; prop_B = (B(s,n_t-1)*area(s))/sum_pop;
        Type exploitation = exp(log(((counter_B(n_t-1)/1000)*plug_exploit)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
        if ( (B(s,n_t-1)*area(s)) > (counter_B(n_t-1)/n_s) ) C(s,n_t-1) = exp(log((exploitation/area(s))*prop_B) + rnorm(Type(0.0),Type(0.2)));
        else C(s,n_t-1) = Type(0.0);
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(s_a(s),n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(s_a(s),n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
        B(s,n_t) = mean_pro_B(s,n_t);
        log_B(s,n_t) = log(B(s,n_t));
        C(s,n_t) = 0;
      }
    }
    //Simulating biomass based on having landings aggregated in above average knots
    //with extra 2% catches at all other knots
    if (options_vec[0]==3) {
      vector <Type> counter_B(n_t); counter_B.setZero();
      for (int t = 0; t < n_t; t++){
        for (int s = 0; s < n_s; s++){
          if (t == 0){
            mean_pro_B(s,t) = B0 * exp(omega_B(v_i(s),t));
            // mean_pro_B(s,t) = init_B(s) * exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
          else {
            Type sum_pop; sum_pop = 0;
            for (int i = 0; i < n_s; i++){
              if ((B(i,t-1)*area(i)) >  (counter_B(t-1)/n_s)) {
                sum_pop = sum_pop + (B(i,t-1)*area(i));
              }
            }

            Type prop_B; prop_B = (B(s,t-1)*area(s))/sum_pop;
            Type exploitation =  exp(log(((counter_B(t-1)/1000)*plug_exploit)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
            if ( (B(s,t-1)*area(s)) > (counter_B(t-1)/n_s) ) C(s,t-1) = (exp(log((exploitation/area(s))*prop_B))+ rnorm(Type(0.0),Type(0.2)));
            else C(s,t-1) = B(s,t-1)*0.02*exp(rnorm(Type(0.0),Type(0.2)));
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(s_a(s),t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(s_a(s),t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
        }
      }
      for (int s = 0; s < n_s; s++){
        Type sum_pop; sum_pop = 0;
        for (int i = 0; i < n_s; i++){
          if ((B(i,n_t-1)*area(i)) >  (counter_B(n_t-1)/n_s)) {
            sum_pop = sum_pop + (B(i,n_t-1)*area(i));
          }
        }
        Type prop_B; prop_B = (B(s,n_t-1)*area(s))/sum_pop;
        Type exploitation = exp(log(((counter_B(n_t-1)/1000)*plug_exploit)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
        if ( (B(s,n_t-1)*area(s)) > (counter_B(n_t-1)/n_s) ) C(s,n_t-1) = exp(log((exploitation/area(s))*prop_B) + rnorm(Type(0.0),Type(0.2)));
        else C(s,n_t-1)=B(s,n_t-1)*0.02*exp(rnorm(Type(0.0),Type(0.2)));
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(s_a(s),n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(s_a(s),n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
        B(s,n_t) = mean_pro_B(s,n_t);
        log_B(s,n_t) = log(B(s,n_t));
        C(s,n_t) = 0;
      }
    }
    REPORT(B);
    REPORT(log_B);
    REPORT(C);
  }

  //Calculating predicted biomass and recruitment over area covered by each knots
  //Commercial biomass
  matrix <Type> areaC(n_s,n_t);
  areaC.setZero();
  for (int s = 0; s < n_s; s++){
    for (int t = 0; t < (n_t+1); t++){
      areaB(s,t) = B(s,t) * area(s);
      areaC(s,t) = C(s,t) * area(s);
    }
  }

  //Recruits
  for (int s = 0; s < n_s; s++){
    for (int t = 0; t < (n_t); t++){
      areaR(s,t) = R(s,t) * area(s);
    }
  }

  //Calculate mean natural mortality, and total biomass and recruitment
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < (n_s); s++){
      mean_m(t) = mean_m(t) + m(s,t);
    }
    mean_m(t) = mean_m(t) / n_s;
  }
  vector <Type> log_mean_m = log(mean_m);

  //Divide by 1000 to represent metric tonnes
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < n_s; s++){
      totB(t) = totB(t) + (areaB(s,t)/1000);
    }
  }
  vector <Type> log_totB = log(totB);

  vector <Type> totC(n_t);
  totC.setZero();
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < n_s; s++){
      totC(t) = totC(t) + (areaC(s,t)/1000);
    }
  }

  vector <Type> exp_rates(n_t-1);
  exp_rates.setZero();
  for (int t = 0; t < (n_t-1); t++){
    exp_rates(t) = totC(t+1)/totB(t);
  }

  for (int t = 0; t < (n_t); t++){
    for (int s = 0; s < n_s; s++){
      totR(t) = totR(t)+(areaR(s,t)/1000);
    }
  }
  vector <Type> log_totR = log(totR);

  vector <Type> log_bio_weighted_mean_m(n_t+1);
  for (int t = 0; t < (n_t+1); t++){
    Type temp_m = 0;
    for (int s = 0; s < (n_s); s++){
      Type temp_prop = areaB(s,t)/(totB(t)*1000);
      log_bio_weighted_mean_m(t) = temp_m + (m(s,t)*temp_prop);
    }
    log_bio_weighted_mean_m(t) = exp(temp_m);
  }

  // Observation equations

  // Probability of capturing commercial biomass
  for (int i = 0; i < n_i; i++){
    logit_zip_I(i) = (vector<Type>(X.row((i)))*beta_I).sum();
    zip_I(i) = invlogit(logit_zip_I(i));
    nll_comp[5] -= dzinbinom2(n_I(i), mu_I, var_I, zip_I(i), true);
  }

  Type zero_prob_I = pow((mu_I/var_I),(sqr(mu_I)/(var_I-mu_I)));
  REPORT(zero_prob_I);

  vector <Type> report_prob_I(n_h); report_prob_I.setZero();
  for (int h = 0; h < n_h; h++){
    for (int i = 0; i < n_i; i++){
      if (h_i(i)==h){
        report_prob_I(h) = (1-zip_I(i) +(1-zip_I(i))*zero_prob_I);
        if (report_prob_I(h) > 0) break;
      }
    }
  }
  ADREPORT(report_prob_I);

  SIMULATE{
    for (int i = 0; i < n_i; i++) {
      logit_zip_I(i) = (vector<Type>(X.row((i)))*beta_I).sum();
      zip_I(i) = invlogit(logit_zip_I(i));
      Type berny_boi = rbinom(Type(1.0),zip_I(i));
      if (berny_boi == 1) n_I(i) = Type(0.0);
      else n_I(i) = rnbinom2(mu_I,var_I);
    }
    REPORT(n_I);
  }


  //Commercial size observations
  if (options_vec[4] == 0){
    PARAMETER(log_qI);//Commercial biomass catchability
    Type qI = exp(log_qI);

    if (options_vec[1] == 1){
      DATA_VECTOR(prior_pars);
      nll_comp[8] -= dbeta(qI,prior_pars[0],prior_pars[1],true);
    }

    // Commercial Index
    for (int i = 0; i < n_i; i++){
      if( !isNA(logI(i) )) {
        Type mean_B = qI*B(s_i(i),t_i(i))/(1-zip_I(h_i(i))+(1-zip_I(h_i(i)))*zero_prob_I);
        nll_comp(4) -= keep_I(i) * dnorm(logI(i), log (mean_B)- sqr(sigma_epsilon)/Type(2.0), sigma_epsilon, true);
      }
    }

    SIMULATE{
      for (int i = 0; i < n_i; i++){
        Type mean_I = qI*B(s_i(i),t_i(i))/(1-zip_I(h_i(i))+(1-zip_I(h_i(i)))*zero_prob_I);
        if (n_I(i)>0) logI(i) = log(mean_I) - (sqr(sigma_epsilon)/2.0) + rnorm(Type(0.0),sigma_epsilon);
        else logI(i) = NA_REAL;
      }
      REPORT(logI);
    }
    REPORT(log_qI);
    REPORT(qI);
    ADREPORT(qI);
  } else if (options_vec[4] == 1) {
    PARAMETER_VECTOR(log_qI); //commercial biomass catchabilities
    vector <Type> qI(n_s); qI = exp(log_qI);
    if (options_vec[1] == 1){
      DATA_VECTOR(prior_pars_1);
      DATA_VECTOR(prior_pars_2);
      for (int s = 0; s < n_s; s++){
        nll_comp[8] -= dbeta(qI(s),prior_pars_1[s],prior_pars_2[s],true);
      }
    }

    // Commercial Index
    for (int i = 0; i < n_i; i++){
      if( !isNA(logI(i) )) {
        Type mean_B = qI(s_i(i))*B(s_i(i),t_i(i))/(1-zip_I(h_i(i))+(1-zip_I(h_i(i)))*zero_prob_I);
        nll_comp(4) -= keep_I(i) * dnorm(logI(i), log (mean_B)- sqr(sigma_epsilon)/Type(2.0), sigma_epsilon, true);
      }
    }

    SIMULATE{

      for (int i = 0; i < n_i; i++){
        Type mean_I = qI(s_i(i))*B(s_i(i),t_i(i))/(1-zip_I(h_i(i))+(1-zip_I(h_i(i)))*zero_prob_I);
        if (n_I(i)>0) logI(i) = log(mean_I) - (sqr(sigma_epsilon)/2.0) + rnorm(Type(0.0),sigma_epsilon);
        else logI(i) = NA_REAL;
      }
      REPORT(logI);
    }
    REPORT(log_qI);
    REPORT(qI);
    ADREPORT(qI);
  } else if (options_vec[4] == 2) {
    PARAMETER_VECTOR(log_qI); //commercial biomass catchabilities
    vector <Type> qI = exp(log_qI);
    if (options_vec[1] == 1){
      DATA_VECTOR(prior_pars);
      for (int h = 0; h < n_h; h++){
        nll_comp[8] -= dbeta(qI(h),prior_pars[(h*2)],prior_pars[(h*2)+1],true);
      }
    }

    // Commercial Index
    for (int i = 0; i < n_i; i++){
      if( !isNA(logI(i) )) {
        Type mean_B = qI(h_i(i))*B(s_i(i),t_i(i))/(1-zip_I(h_i(i))+(1-zip_I(h_i(i)))*zero_prob_I);
        nll_comp(4) -= keep_I(i) * dnorm(logI(i), log (mean_B)- sqr(sigma_epsilon)/Type(2.0), sigma_epsilon, true);
      }
    }

    SIMULATE{
      for (int i = 0; i < n_i; i++){
        Type mean_I = qI(h_i(i))*B(s_i(i),t_i(i))/(1-zip_I(h_i(i))+(1-zip_I(h_i(i)))*zero_prob_I);
        if (n_I(i)>0) logI(i) = log(mean_I) - (sqr(sigma_epsilon)/2.0) + rnorm(Type(0.0),sigma_epsilon);
        else logI(i) = NA_REAL;
      }
      REPORT(logI);
    }
    REPORT(log_qI);
    REPORT(qI);
    ADREPORT(qI);
  }

  //Probability of capturing recruits
  for (int i = 0; i < n_i; i++){
    nll_comp[11] -= dnbinom2(n_IR(i), mu_IR(h_i(i)), var_IR(h_i(i)), true);
  }

  SIMULATE {
    for (int i = 0; i < n_i; i++) {
      n_IR(i) = rnbinom2(mu_IR(h_i(i)), var_IR(h_i(i)));
    }
    REPORT(n_IR);
  }

  vector <Type> zero_prob_IR(n_h);
  for (int h = 0; h < n_h; h++){
    zero_prob_IR(h) = pow((mu_IR(h)/var_IR(h)),(sqr(mu_IR(h))/(var_IR(h)-mu_IR(h))));
  }
  REPORT(zero_prob_IR);

  vector <Type> report_prob_IR(n_h);
  for (int h = 0; h < n_h; h++){
    report_prob_IR(h) = 1-zero_prob_IR(h);
  }
  ADREPORT(report_prob_IR);

  // Recruit Index
  for (int i = 0; i < n_i; i++){
    if ( !isNA(logIR(i) )) {
      Type mean_R = qR(h_i(i))*R(s_i(i),t_i(i))/(1-zero_prob_IR(h_i(i)));
      nll_comp(6) -= keep_IR(i) * dnorm(logIR(i), log(mean_R)-sqr(sigma_upsilon)/Type(2.0), sigma_upsilon, true);
    }
  }

  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_IR = qR(h_i(i))*R(s_i(i),t_i(i))/(1-zero_prob_IR(h_i(i)));
      if(n_IR(i)>0) logIR(i) = log(mean_IR) - (sqr(sigma_upsilon)/2.0) + rnorm(Type(0.0),sigma_upsilon);
      else logIR(i) = NA_REAL;
    }
    REPORT(logIR);
  }


  if (options_vec[3] == 1 || options_vec[3] == 2) {
    DATA_VECTOR(L); //number of clappers (dim n_i)
    DATA_VECTOR(n_bin); //number of shell caught in tow i (dim n_i)

    DATA_VECTOR_INDICATOR(keep_L, L);

    PARAMETER_VECTOR(log_S); //clapper catchability
    vector <Type> S = exp(log_S);

    // Observations to natural mortality
    for (int i = 0; i < n_i; i++){
      if (!isNA(L(i))) {
        nll_comp(7) -= keep_L(i) * dbinom_robust(L(i), n_bin(i),logit(m(s_i(i),t_i(i))*S(h_i(i))),true);
      }
    }

    SIMULATE {
      for (int i = 0; i < n_i; i++){
        n_bin(i) = rpois(Type(400));
        L(i) = rbinom(n_bin(i), m(s_i(i),t_i(i))*S(h_i(i)));
      }
      REPORT(L);
      REPORT(n_bin);
    }
    REPORT(S);
    ADREPORT(S);
  }

  //Reporting
  // ADREPORT all individual parameters to get their standard errors
  REPORT(sigma_epsilon);
  REPORT(sigma_upsilon);
  REPORT(R0);
  REPORT(B0);
  REPORT(m0);
  REPORT(qR);
  REPORT(mu_I);
  REPORT(mu_IR);
  REPORT(var_I);
  REPORT(var_IR);
  REPORT(beta_I);
  REPORT(zip_I);
  ADREPORT(zip_I);
  ADREPORT(sigma_epsilon);
  ADREPORT(sigma_upsilon);
  ADREPORT(R0);
  ADREPORT(B0);
  ADREPORT(m0);
  ADREPORT(qR);
  ADREPORT(mu_I);
  ADREPORT(mu_IR);
  ADREPORT(var_I);
  ADREPORT(var_IR);
  ADREPORT(beta_I);

  //Reporting processes
  REPORT(B);
  REPORT(R);
  REPORT(omega_B);
  REPORT(omega_R);
  REPORT(areaB);
  REPORT(areaR);

  //Reporting derived totals
  REPORT(totB);
  REPORT(totR);
  ADREPORT(totB);
  ADREPORT(totR);
  REPORT(exp_rates);

  REPORT(m);
  REPORT(mean_m);
  ADREPORT(mean_m);

  REPORT(log_totB);
  REPORT(log_totR);
  REPORT(log_mean_m);
  ADREPORT(log_totB);
  ADREPORT(log_totR);
  ADREPORT(log_mean_m);
  REPORT(log_bio_weighted_mean_m);
  ADREPORT(log_bio_weighted_mean_m);

  if (options_vec[2] == 1){
    ADREPORT(R);
    ADREPORT(B);
    ADREPORT(m);
  } else {}

  //Report individual negative log-likelihood contributions
  REPORT(nll_comp);

  Type nll = nll_comp.sum();
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
