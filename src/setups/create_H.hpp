//Function to create the H anisotropy matrices utilizing two parameters
//Parameterize so that det(H)=1 to preserve volume as seen in Lindgren 2011

template <class Type>
matrix <Type> create_H(vector <Type> log_H_input) {
  matrix<Type> H(2,2);
  H(0,0) = exp(log_H_input(0));
  H(1,0) = log_H_input(1);
  H(0,1) = log_H_input(1);
  H(1,1) = (1+log_H_input(1)*log_H_input(1)) / exp(log_H_input(0));
  return H;
}
