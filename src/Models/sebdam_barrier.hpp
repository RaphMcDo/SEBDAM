#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sebdam_barrier(objective_function <Type>* obj) {
  // Read in data, parameters, random effects
  using namespace barrier_q;
  using namespace density;
  using namespace Eigen;

  //Data
  DATA_VECTOR(I); //commercial biomass survey index by tow (dim n_i)
  DATA_VECTOR(IR); //recruit survey index by tow (dim n_i)
  DATA_VECTOR(area); //area covered by each knot (dim n_s)
  DATA_MATRIX(C); //commercial catch (dim n_s,n_t)
  DATA_VECTOR(L); //number of clappers (dim n_i)
  DATA_VECTOR(n_bin); //number of shell caught in tow (dim n_i)
  DATA_VECTOR(n_tows); //number of tows in year (dim n_t)
  DATA_VECTOR(pos_tows_I); //number of tows that captured commercial biomass (dim n_t)
  DATA_VECTOR(pos_tows_IR); //number of tows that captured recruits (dim n_t)

  //Sizes
  DATA_INTEGER(n_i); //Number of observations per year
  DATA_INTEGER(n_t); //Number of years
  DATA_INTEGER(n_s); //Number of knots
  DATA_INTEGER(n_m); //number of vertex in mesh

  //Indices
  DATA_FACTOR(s_i); //Indexing for knot (dim n_i)
  DATA_FACTOR(t_i); //Indexing for year (dim n_i)
  DATA_FACTOR(v_i); //Indexing for specific mesh location to match knots with right vertex (dim n_m)

  //SPDE objects
  DATA_STRUCT(fem, fem_barrier_t); //INLA FEM object

  //Fixed covariates
  DATA_VECTOR(gI); //commercial biomass growth
  DATA_VECTOR(gR); //recruitment growth

  DATA_IVECTOR(options_vec); //Vector of multiple choices
  //Slot 0: Choice of simulated landings
  //  0 (default)= No landings
  //  1 = Proportionally distributed at every knot
  //  2 = Proportionally distributed at knots with above average biomass
  //  3 = Proportionally distributed at knots with above average biomass with extra landings at other knots
  //Slot 1: Choice of using prior on catchability
  //  0 (default) = No prior
  //  1 = Use prior on catchabilities
  //Slot 2: Choice of obtaining computationally expensive standard errors
  // 0 (default) = Obtain standard errors only for parameters and overall processes
  // 1 = Obtain all standard errors

  //Parameters
  PARAMETER(log_sigma_epsilon); //obs sd survey index commercial biomass
  PARAMETER(log_sigma_upsilon); //obs sd survey index recruits
  PARAMETER(log_S); //clapper catchability
  PARAMETER_VECTOR(log_range_B);//commercial size range parameter
  PARAMETER(log_sigma_B); //commercial size spatio-temporal variability parameter
  PARAMETER_VECTOR(log_range_R); //recruit range parameter
  PARAMETER(log_sigma_R); //recruit spatio-temporal variability parameter
  PARAMETER_VECTOR(log_range_m); //mortality range parameter
  PARAMETER(log_sigma_m); //mortality spatio-temporal variability parameter
  PARAMETER(log_R0);//initial recruit mean value
  PARAMETER(log_B0);//initial commercial biomass mean value
  PARAMETER(log_m0);//initial mortality mean value
  PARAMETER(logit_p_I); //probability of capturing commercial biomass
  PARAMETER(logit_p_IR); //probability of capturing recruits
  PARAMETER_VECTOR(log_qI); //commerical biomass catchability
  PARAMETER(log_qR); //recruit catchability

  Type pi = M_PI;
  vector <Type> range_B; range_B = exp(log_range_B);
  vector <Type> range_R; range_R = exp(log_range_R);
  vector <Type> range_m; range_m = exp(log_range_m);
  Type sigma_B = exp(log_sigma_B);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_m = exp(log_sigma_m);
  Type sigma_epsilon = exp(log_sigma_epsilon);
  Type sigma_upsilon = exp(log_sigma_upsilon);
  Type S = exp(log_S);
  Type qR = exp(log_qR);
  vector <Type> qI(n_s); qI = exp(log_qI);
  Type R0 = exp(log_R0);
  Type B0 = exp(log_B0);
  Type m0 = exp(log_m0);
  Type p_I = invlogit(logit_p_I);
  Type p_IR = invlogit(logit_p_IR);

  //Random Effects
  PARAMETER_ARRAY(omega_B);
  PARAMETER_ARRAY(omega_R);
  PARAMETER_ARRAY(omega_m);

  // Set up matrices for processes of interest
  matrix <Type> log_B(n_s,(n_t+1));
  matrix <Type> B(n_s,(n_t+1));
  matrix <Type> areaB(n_s,(n_t+1));
  matrix <Type> areaR(n_s,(n_t));
  matrix <Type> log_R(n_s,(n_t));
  matrix <Type> R(n_s,(n_t));
  matrix <Type> log_m(n_s,(n_t+1));
  matrix <Type> m(n_s,(n_t+1));

  //Set up initial states and other elements
  SIMULATE{
    n_bin = rpois(Type(100));
    n_tows = Type(120);
    for (int t = 0; t < n_t; t++){
      gI(t) = 1.1;
      gR(t) = 1.5;
    }
    REPORT(gI);
    REPORT(gR);
    REPORT(n_bin);
    REPORT(n_tows);
  }

  //Setup for simulations and derived values
  vector <Type> bern_I(n_i);
  vector <Type> bern_IR(n_i);
  vector <Type> totB(n_t+1); totB.setZero();
  vector <Type> totR(n_t); totR.setZero();
  vector <Type> mean_m(n_t+1); mean_m.setZero();
  matrix <Type> mean_pro_B(n_s,n_t+1);

  // ----------------------------------------------
  // nll
  vector <Type> nll_comp(9); nll_comp.setZero();


  //Set up GMRF for commercial size biomass
  SparseMatrix <Type> Q_B = Q_barrier(fem, range_B, sigma_B);
  for (int t = 0; t < (n_t+1); t++){
    nll_comp(0) += GMRF(Q_B)(omega_B.col(t));
  }

  SIMULATE {
    SparseMatrix <Type> Q_B = Q_barrier(fem, range_B, sigma_B);
    density::GMRF_t<Type> gmrf1(Q_B);
    for (int t = 0; t < (n_t+1); t++){
      vector <Type> temp_omega(n_m);
      gmrf1.simulate(temp_omega);
      omega_B.col(t)=temp_omega;
    }
    REPORT(Q_B);
    REPORT(omega_B);
  }

  //Set up GMRF for recruits
  SparseMatrix <Type> Q_R = Q_barrier(fem, range_R, sigma_R);
  for (int t = 0; t < (n_t); t++){
    nll_comp(1) += GMRF(Q_R)(omega_R.col(t));
  }

  SIMULATE {
    SparseMatrix <Type> Q_R = Q_barrier(fem, range_R, sigma_R);
    density::GMRF_t<Type> gmrf2(Q_R);
    for (int t = 0; t < n_t; t++){
      vector <Type> temp_omega(n_m);
      gmrf2.simulate(temp_omega);
      omega_R.col(t)=temp_omega;
    }
    REPORT(Q_R);
    REPORT(omega_R);
  }

  //Set up GMRF for mortality
  SparseMatrix <Type> Q_m = Q_barrier(fem, range_m, sigma_m);
  for (int t = 0; t < (n_t+1); t++){
    nll_comp(2) += GMRF(Q_m)(omega_m.col(t));
  }

  SIMULATE {
    SparseMatrix <Type> Q_m = Q_barrier(fem, range_m, sigma_m);
    density::GMRF_t<Type> gmrf3(Q_m);
    for (int t = 0; t < (n_t+1); t++){
      vector <Type> temp_omega(n_m);
      gmrf3.simulate(temp_omega);
      omega_m.col(t)=temp_omega;
    }
    REPORT(Q_m);
    REPORT(omega_m);
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
    //Predict 1 year ahead, assume same spatial pattern for error but propagate error
    log_m(s,n_t) = log( m(s,n_t-1)*exp(omega_m(v_i(s),n_t)));
    m(s,n_t) = exp(log_m(s,n_t));
  }

  //Simulate Natural mortality
  SIMULATE {
    for (int s = 0; s < n_s; s++){
      m(s,0) = m0 * exp(omega_m(v_i(s),0));
      log_m(s,0) = log(m(s,0));
      for (int t = 1; t < (n_t); t++){
        m(s,t) = m(s,t-1) * exp(omega_m(v_i(s),t)) ;
        log_m(s,t) = log(m(s,t));
      }
      m(s,n_t) = m(s,n_t-1) * exp(omega_m(v_i(s),n_t));
      log_m(s,n_t) = log(m(s,n_t));
    }
    REPORT(m);
    REPORT(log_m);
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

  //Simulate recruitment
  SIMULATE{
    for (int s = 0; s < n_s; s++){
      R(s,0) = R0 * exp(omega_R(v_i(s),0));
      log_R(s,0) = log(R(s,0));
      for (int t = 1; t < (n_t); t++){
        R(s,t) = R(s,t-1)*exp(omega_R(v_i(s),t));
        log_R(s,t) = log(R(s,t));
      }
    }
    REPORT(R);
    REPORT(log_R);
  }

  //Biomass derivation
  for (int s = 0; s < n_s; s++) {
    log_B(s,0) = log(B0*exp(omega_B(v_i(s),0)));
    B(s,0) = exp(log_B(s,0));
    for (int t = 1; t < (n_t); t++) {
      Type mean_pro = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
      log_B(s,t) = log(mean_pro);
      B(s,t) = exp(log_B(s,t));
    }
    //Predict 1 year ahead, assume same spatial pattern as last year but propagate error
    Type mean_pro =(exp(-m(s,n_t))*gI(n_t-1)*(B(s, (n_t - 1)) - C(s,(n_t-1))) + exp(-m(s,n_t))*gR(n_t-1)*R(s,(n_t-1)))*exp(omega_B(v_i(s),n_t));
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
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
          else {
            C(s,t-1) = Type(0.0);
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
        }
      }
      //Simulating 1-year projection
      for (int s = 0; s < n_s; s++){
        C(s,n_t-1) = Type(0.0);
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
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
            B(s,t) = exp(log(mean_pro_B(s,t)));
            log_B(s,t) = log(B(s,t));
            counter_B(t) = counter_B(t) + (B(s,t)*area(s));
          }
          else {
            //Calculate proportion of total biomass at knot s
            Type prop_B; prop_B = (B(s,t-1)*area(s))/counter_B(t-1);
            //Simulate total landings from total biomass on log scale
            Type exploitation =  exp(log(((counter_B(t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
            //Simulate landings per km^2 at knot s
            C(s,t-1) = exp(log((exploitation*prop_B)/area(s)) + rnorm(Type(0.0),Type(0.2)));
            //Simulate biomass density
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
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
        Type exploitation = exp(log(((counter_B(n_t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
        //Simulate landings per km^2 at knot s
        C(s,n_t-1) = exp(log((exploitation*prop_B)/area(s)) + rnorm(Type(0.0),Type(0.2)));
        //Simulate biomass density
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
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
            Type exploitation =  exp(log(((counter_B(t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
            if ( (B(s,t-1)*area(s)) > (counter_B(t-1)/n_s) ) C(s,t-1) = (exp(log((exploitation/area(s))*prop_B))+ rnorm(Type(0.0),Type(0.2)));
            else C(s,t-1) = Type(0.0);
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
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
        Type exploitation = exp(log(((counter_B(n_t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
        if ( (B(s,n_t-1)*area(s)) > (counter_B(n_t-1)/n_s) ) C(s,n_t-1) = exp(log((exploitation/area(s))*prop_B) + rnorm(Type(0.0),Type(0.2)));
        else C(s,n_t-1) = Type(0.0);
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
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
            Type exploitation =  exp(log(((counter_B(t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
            if ( (B(s,t-1)*area(s)) > (counter_B(t-1)/n_s) ) C(s,t-1) = (exp(log((exploitation/area(s))*prop_B))+ rnorm(Type(0.0),Type(0.2)));
            else C(s,t-1) = B(s,t-1)*0.02*exp(rnorm(Type(0.0),Type(0.2)));
            mean_pro_B(s,t) = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
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
        Type exploitation = exp(log(((counter_B(n_t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
        if ( (B(s,n_t-1)*area(s)) > (counter_B(n_t-1)/n_s) ) C(s,n_t-1) = exp(log((exploitation/area(s))*prop_B) + rnorm(Type(0.0),Type(0.2)));
        else C(s,n_t-1)=B(s,n_t-1)*0.02*exp(rnorm(Type(0.0),Type(0.2)));
        mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
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
  for (int s = 0; s < n_s; s++){
    for (int t = 0; t < (n_t+1); t++){
      areaB(s,t) = B(s,t) * area(s);
    }
  }

  SIMULATE{
    for (int s = 0; s < n_s; s++){
      for (int t = 0; t < (n_t+1); t++){
        areaB(s,t) = B(s,t) * area(s);
      }
    }
    REPORT(areaB);
  }

  //Recruits
  for (int s = 0; s < n_s; s++){
    for (int t = 0; t < (n_t); t++){
      areaR(s,t) = R(s,t) * area(s);
    }
  }

  SIMULATE{
    for (int s = 0; s < n_s; s++){
      for (int t = 0; t < (n_t); t++){
        areaR(s,t) = R(s,t) * area(s);
      }
    }
    REPORT(areaR);
  }

  //Calculate mean natural mortality, and total biomass and recruitment
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < (n_s); s++){
      mean_m(t) = mean_m(t) + m(s,t);
    }
    mean_m(t) = mean_m(t) / n_s;
  }

  //Divide by 1000 to represent metric tonnes
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < n_s; s++){
      totB(t) = totB(t) + (areaB(s,t)/1000);
    }
  }

  for (int t = 0; t < (n_t); t++){
    for (int s = 0; s < n_s; s++){
      totR(t) = totR(t)+(areaR(s,t)/1000);
    }
  }

  // Observation equations

  //Probability of capturing commercial biomass
  for (int t = 0; t <n_t; t++){
    nll_comp[3] -= dbinom_robust(pos_tows_I(t), n_tows(t), logit(p_I),true);
  }

  SIMULATE{
    for (int t = 0; t < n_t; t++){
      pos_tows_I(t) = rbinom(n_tows(t),p_I);
    }
    REPORT(pos_tows_I);
  }

  // Set prior on q
  if (options_vec[1] == 1){
    for (int s = 0; s < n_s; s++){
      nll_comp[8] -= dbeta(qI(s),Type(10.0),Type(12.0),true);
    }
  }

  // Commercial Index
  for (int i = 0; i < n_i; i++){
    if( !isNA(I(i) )) {
      Type mean_B = qI(s_i(i))*B(s_i(i),t_i(i))/p_I;
      nll_comp(4) -= dnorm(log(I(i)), log (mean_B)- sqr(sigma_epsilon)/Type(2.0), sigma_epsilon, true);
    }
  }

  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_I = qI(s_i(i))*B(s_i(i),t_i(i))/p_I;
      bern_I(i) = rbinom(Type(1.0),p_I);
      if (bern_I(i)>0) I(i) = exp(log(mean_I) - (sqr(sigma_epsilon)/2.0) + rnorm(Type(0.0),sigma_epsilon));
      else I(i) = NA_REAL;
    }
    REPORT(bern_I);
    REPORT(I);
  }

  //Probability of capturing recruits
  for (int t = 0; t <n_t; t++){
    nll_comp[5] -= dbinom_robust(pos_tows_IR(t), n_tows(t), logit(p_IR),true);
  }

  SIMULATE {
    for (int t = 0; t < n_t; t++){
      pos_tows_IR(t) = rbinom(n_tows(t),p_IR);
    }
    REPORT(pos_tows_IR);
  }

  // Recruit Index
  for (int i = 0; i < n_i; i++){
    if ( !isNA(IR(i) )) {
      Type mean_R = qR*R(s_i(i),t_i(i))/p_IR;
      nll_comp(6) -= dnorm(log(IR(i)), log(mean_R)-sqr(sigma_upsilon)/Type(2.0), sigma_upsilon, true);
    }
  }

  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_IR = qR*R(s_i(i),t_i(i))/p_IR;
      bern_IR(i) = rbinom(Type(1.0),p_IR);
      if(bern_IR(i)>0) IR(i) = exp(log(mean_IR) - (sqr(sigma_upsilon)/2.0) + rnorm(Type(0.0),sigma_upsilon));
      else IR(i) = NA_REAL;
    }
    REPORT(bern_IR);
    REPORT(IR);
  }

  // Observations to natural mortality
  for (int i = 0; i < n_i; i++){
    if (!isNA(L(i))) {
      nll_comp(7) -= dbinom_robust(L(i), n_bin(i),logit(m(s_i(i),t_i(i))*S),true);
    }
  }

  SIMULATE {
    for (int i = 0; i < n_i; i++){
      L(i) = rbinom(n_bin(i), m(s_i(i),t_i(i))*S);
    }
    REPORT(L);
  }

  //Reporting
  // ADREPORT all individual parameters to get their standard errors
  REPORT(sigma_epsilon);
  REPORT(sigma_upsilon);
  REPORT(S);
  REPORT(R0);
  REPORT(B0);
  REPORT(m0);
  REPORT(qI);
  REPORT(qR);
  REPORT(p_I);
  REPORT(p_IR);
  REPORT(range_B);
  REPORT(range_R);
  REPORT(range_m);
  REPORT(sigma_B);
  REPORT(sigma_R);
  REPORT(sigma_m);

  ADREPORT(sigma_epsilon);
  ADREPORT(sigma_upsilon);
  ADREPORT(S);
  ADREPORT(R0);
  ADREPORT(B0);
  ADREPORT(m0);
  ADREPORT(qI);
  ADREPORT(qR);
  ADREPORT(p_I);
  ADREPORT(p_IR);
  ADREPORT(range_B);
  ADREPORT(range_R);
  ADREPORT(range_m);
  ADREPORT(sigma_B);
  ADREPORT(sigma_R);
  ADREPORT(sigma_m);

  //Reporting processes
  REPORT(B);
  REPORT(R);
  REPORT(m);
  REPORT(omega_B);
  REPORT(omega_R);
  REPORT(omega_m);
  REPORT(areaB);
  REPORT(areaR);

  //Reporting derived totals
  REPORT(totB);
  REPORT(totR);
  REPORT(mean_m);
  ADREPORT(totB);
  ADREPORT(totR);
  ADREPORT(mean_m);

  if (options_vec[2] == 1){
    ADREPORT(R);
    ADREPORT(m);
    ADREPORT(B);
    ADREPORT(omega_R);
    ADREPORT(omega_B);
    ADREPORT(omega_m);
    ADREPORT(areaB);
    ADREPORT(areaR);
  } else {}

  //Report individual negative log-likelihood contributions
  REPORT(nll_comp);

  Type nll = nll_comp.sum();
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
