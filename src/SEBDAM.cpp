#define TMB_LIB_INIT R_init_SEBDAM
#include <TMB.hpp>

#include "setups/helpers.hpp"
#include "setups/barrier_q.hpp"
#include "setups/create_H.hpp"

#include "Models/tlm.hpp"
#include "Models/sebdam.hpp"
#include "Models/sehbam.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);


  if (model == "SEBDAM") {
    return sebdam(this);
  }  else if (model =="TLM"){
      return tlm(this);
  } else if (model=="SEHBAM"){
    return sehbam(this);
  } else {
      error("Unknown model");
  }

    return 0;

}
