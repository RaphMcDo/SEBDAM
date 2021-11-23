#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tlm(objective_function <Type>* obj) {
  using namespace density;
  using namespace Eigen;

  DATA_IVECTOR(options_vec);
  //Slot 0: use a prior for q or not
  //    0 (default): no prior
  //    1: use prior
  //Slot 1: Choice of mortality approach, given that most fisheries don't have any observations for mortality
  // 0 (default) = Fix mortality to user-provided value
  // 1 = use clapper approach

  //----------------------------------------------------------------------------
  // Data
  //----------------------------------------------------------------------------
  //(* Directly observable variable)
  DATA_VECTOR(logI); // observed survey commercial size biomass (dim n_i)
  DATA_VECTOR(logIR); // observed survey recruit biomass (dim n_i)
  DATA_VECTOR(C); // observed commercial landings (dim NY)
  DATA_VECTOR(g); // covariate: growth rate commercial size (dim NY)
  DATA_VECTOR(gR); // covariate: growth rate recruit size (dim NY)
  DATA_VECTOR(n_tows); //total number of tows in year (dim NY)
  DATA_VECTOR(pos_tows_I); //number of tows with positive catches (dim NY)
  DATA_VECTOR(pos_tows_IR); //number of tows with positive catches (dim NY)

  //For residuals
  DATA_VECTOR_INDICATOR(keep_I, logI);
  DATA_VECTOR_INDICATOR(keep_IR, logIR);

  //Indices
  DATA_FACTOR(t_i); //index to identify what year to attribute observation

  DATA_SCALAR(plug_exploit); // Exploitation rate desired for simulations, not used when fitting data

  //----------------------------------------------------------------------------
  // Parameters
  //----------------------------------------------------------------------------

  PARAMETER(log_sigma_tau); // proc sd biomass
  PARAMETER(log_sigma_phi); // proc sd recruits
  PARAMETER(log_sigma_epsilon); // obs sd survey tows commercial size
  PARAMETER(log_sigma_upsilon); // obs sd survey tows recruits
  PARAMETER(log_q_I); //catchability for commercial size
  PARAMETER(log_q_R); // catchability for recruits
  PARAMETER(logit_p_I); //probability of capturing commercial biomass
  PARAMETER(logit_p_IR); //probability of capturing recruits

  PARAMETER_VECTOR(log_B); // process: biomass (dim NY)
  PARAMETER_VECTOR(log_R); // process: recruits (dim NY)

  //----------------------------------------------------------------------------
  // Setup
  //----------------------------------------------------------------------------

  int NY = pos_tows_I.size(); // number of years
  int n_i = logI.size(); //number of observations

  // Transform input parameters to natural scale
  Type sigma_tau = exp(log_sigma_tau);
  Type sigma_phi = exp(log_sigma_phi);
  Type sigma_epsilon = exp(log_sigma_epsilon);
  Type sigma_upsilon = exp(log_sigma_upsilon);
  Type q_R = exp(log_q_R);
  Type q_I = exp(log_q_I);
  Type p_I = invlogit(logit_p_I);
  Type p_IR = invlogit(logit_p_IR);

  vector <Type> B = exp(log_B);
  vector <Type> R = exp(log_R);
  vector <Type> m(NY+1); m.setZero();
  vector <Type> log_m(NY+1); log_m.setZero();

  vector<Type> nll_comp(10); // nll components, break down contributions
  nll_comp.setZero(); // initialize

  // For simulating empty tows in data vectors I and IR
  vector <Type> bern_I(n_i);
  vector <Type> bern_IR(n_i);
  vector <Type> sum_tows(NY); sum_tows.setZero();

  //Simulate initial states and other datatypes
  SIMULATE{
    for (int i = 0; i < n_i; i++){
      sum_tows(t_i(i)) = sum_tows(t_i(i)) + 1;
    }

    for (int t = 0; t < NY; t++){
      n_tows(t) = sum_tows(t);
    }
    REPORT(n_tows);
  }

  //----------------------------------------------------------------------------
  // Proc eq
  //----------------------------------------------------------------------------

  // RW for recruits
  for (int t=1; t<(NY+1); t++){ // R[0] free, yet predicted
    Type mean_proc_R =  R[t-1];
    nll_comp[0] -= dnorm(log_R[t], log(mean_proc_R)-sqr(sigma_phi)/2.0, sigma_phi, true);
  }

  //Simulate recruitment
  SIMULATE {
    R(0) = Type(100);
    log_R(0) = log(R(0));
    for (int t=1; t<(NY+1); t++){
      Type mean_R =  R[t-1];
      R(t) = exp(log(mean_R)-(sqr(sigma_phi)/2) + rnorm(Type(0.0), sigma_phi));
      log_R(t) = log(R(t));
    }
    REPORT(log_R);
    REPORT(R);
  }

  if (options_vec[1] == 0){

    DATA_SCALAR(set_m);
    for (int t=0; t<(NY+1);t++){
      m(t) = set_m;
      log_m(t) = log(m(t));
    }
    REPORT(m);
    REPORT(log_m);

  } else if (options_vec[1] == 1){
    PARAMETER(log_sigma_m); // proc sd nat mort
    PARAMETER_VECTOR(log_input_m); // process: nat mort (dim NY)

    //Transform
    Type sigma_m = exp(log_sigma_m);
    vector<Type> input_m = exp(log_input_m);
    m(0) = input_m(0);
    // RW for nat mort
    for (int t=1; t<(NY+1); t++){ // m[0] free, yet predicted
      Type mean_proc_m = input_m[t-1];
      nll_comp[1] -= dnorm(log_input_m[t], log(mean_proc_m)-sqr(sigma_m)/2.0, sigma_m, true);
      m(t) = input_m(t);
    }

    //Simulate Natural Mortality
    SIMULATE {
      m(0) = Type(0.1);
      log_m(0) = log(m(0));
      for (int t=1; t<(NY+1); t++){
        m(t) = exp(log(m(t-1))-(sqr(sigma_m)/2) + rnorm(Type(0.0),sigma_m));
        log_m(t) = log(m(t));
      }
      input_m = m;
      log_input_m = log(input_m);
      REPORT(log_m);
      REPORT(m);
      REPORT(input_m);
      REPORT(log_input_m);
    }

    ADREPORT(sigma_m);
    ADREPORT(m);
    REPORT(m);
  }

  // Biomass Process
  for (int t=1; t<(NY+1); t++){ // B[0] free, yet predicted
    Type mean_proc_B = (exp(-m[t])*g[t-1]*(B[t-1]-C[t-1])+exp(-m[t])*gR[t-1]*R[t-1]);
    nll_comp[2] -= dnorm(log_B[t], log(mean_proc_B)-sqr(sigma_tau)/2.0, sigma_tau, true);
  }

  //Simulate Biomass and Commercial landings
  SIMULATE {
    B(0) = Type(1000);
    log_B(0) = log(B(0));
    C(0) = exp(log(B(0)*plug_exploit)+rnorm(Type(0.0),Type(0.2)));
    for (int t = 1; t < (NY+1); t++){
      Type mean_proc_B = (exp(-m[t])*g[t-1]*(B[t-1]-C[t-1])+exp(-m[t])*gR[t-1]*R[t-1]);
      B(t) = exp(log(mean_proc_B) - (sqr(sigma_tau)/2) + rnorm(Type(0.0),sigma_tau));
      log_B(t) = log(B(t));
      C(t) = exp(log(B(t)/Type(10.0))+rnorm(Type(0.0),Type(0.2))) ;
    }
    REPORT(log_B);
    REPORT(B);
    REPORT(C);
  }

  //----------------------------------------------------------------------------
  // Obs eq
  //----------------------------------------------------------------------------

  // Prob of capturing commercial biomass
  for (int t = 0; t < NY; t++){
    nll_comp[3] -= dbinom_robust(pos_tows_I(t), n_tows(t), logit(p_I),true);
  }

  //Simulate number of tows with positive commercial size catches
  SIMULATE{
    for (int t = 0; t < NY; t++){
      pos_tows_I(t) = rbinom(n_tows(t),p_I);
    }
    REPORT(pos_tows_I);
  }

  //If prior on q
  if (options_vec[0] == 1) {
    DATA_VECTOR(prior_pars);
    nll_comp[4] -= dbeta(q_I,prior_pars[0],prior_pars[1],true);
    }

  //Observation equations for observed survey commercial size biomass
  for (int i = 0; i < n_i; i++){
    if( !isNA(logI(i) )) {
      Type mean_B = q_I*B(t_i(i))/p_I;
      nll_comp[4] -= keep_I(i) * dnorm(logI(i), log(mean_B)-sqr(sigma_epsilon)/2.0, sigma_epsilon, true);
    }
  }

  //Simulate observed survey commercial size biomass and empty tows
  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_I = q_I*B(t_i(i))/p_I;
      bern_I(i) = rbinom(Type(1.0),p_I);
      if (bern_I(i)>0) logI(i) = log(mean_I) - (sqr(sigma_epsilon)/2) + rnorm(Type(0.0),sigma_epsilon);
      else logI(i) = NA_REAL;
    }
    REPORT(bern_I);
    REPORT(logI);
  }

  //Probability of capturing recruits
  for (int t = 0; t < NY; t++){
    nll_comp[5] -= dbinom_robust(pos_tows_IR(t), n_tows(t), logit(p_IR),true);
  }

  //Simulate number of tows with positive recruit catches
  SIMULATE {
    for (int t = 0; t < NY; t++){
      pos_tows_IR(t) = rbinom(n_tows(t),p_IR);
    }
    REPORT(pos_tows_IR);
  }

  //Observation Equation for Recruits
  for (int i = 0; i < n_i; i++){
    if( !isNA(logIR(i) )) {
      Type mean_R = q_R*R(t_i(i))/p_IR;
      nll_comp[6] -= keep_IR(i) * dnorm(logIR(i), log(mean_R)-sqr(sigma_upsilon)/2.0, sigma_upsilon, true);
    }
  }

  //Simulate positive and empty recruit survey tows
  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_IR = q_R*R(t_i(i))/p_IR;
      bern_IR(i) = rbinom(Type(1.0),p_IR);
      if(bern_IR(i)>0) logIR(i) = log(mean_IR) - (sqr(sigma_upsilon)/2) + rnorm(Type(0.0),sigma_upsilon);
      else logIR(i) = NA_REAL;
    }
    REPORT(bern_IR);
    REPORT(logIR);
  }

  if (options_vec[1] == 1) {
    DATA_VECTOR(n_bin); //total number of shells caught in tow (dim n_i)
    DATA_VECTOR(L); //number of clappers in tow (dim n_i)

    DATA_VECTOR_INDICATOR(keep_L, L);

    PARAMETER(log_S); // clapper catchability
    Type S = exp(log_S);

    // Clappers with replicates with binomial approach
    for (int i = 0; i < n_i; i++){
      if (!isNA(L(i))) {
        nll_comp[7] -= keep_L(i) * dbinom(L(i), n_bin(i),m(t_i(i))*S,true);
      }
    }

    //Simulate clapper count
    SIMULATE {

      for (int i = 0; i < n_i; i++){
        n_bin(i) = rpois(100);
      }
      for (int i = 0; i < n_i; i++){
        L(i) = rbinom(n_bin(i), m(t_i(i))*S);
      }
      REPORT(n_bin);
      REPORT(L);
    }
    ADREPORT(S);

  }


  //----------------------------------------------------------------------------
  // Outputs
  //----------------------------------------------------------------------------

  REPORT(nll_comp);

  ADREPORT(sigma_tau);
  ADREPORT(sigma_phi);
  ADREPORT(sigma_epsilon);
  ADREPORT(sigma_upsilon);
  ADREPORT(q_I);
  ADREPORT(q_R);
  ADREPORT(p_I);
  ADREPORT(p_IR);

  ADREPORT(B);
  REPORT(B);
  ADREPORT(R);
  REPORT(R);

  Type nll = nll_comp.sum();
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
