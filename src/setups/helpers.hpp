

template <class Type>
Type sqr(Type x) {
  return x * x;
}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}
