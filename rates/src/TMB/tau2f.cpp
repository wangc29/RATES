#include <TMB.hpp>

//using namespace Eigen;
using namespace tmbutils;
/* List of matrices */
template<class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x){  /* x = List passed from R */
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

/* List of vectors */
template<class Type>
struct LOV_t : vector<vector<Type> > {
  LOV_t(SEXP x){  /* x = List passed from R */
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(sm);
    }
  }
};



template<class Type>
Type objective_function<Type>::operator() ()
{
  //Data from R
  // z2invsj list of matrices
  DATA_STRUCT(z2invsj, LOM_t);
  // dj list of vectors
  DATA_STRUCT(dj, LOV_t);
  DATA_VECTOR(gammaj);
  //Parameters from R
  PARAMETER_VECTOR(logTau2);
  vector<Type> tau2 = exp(logTau2);
  // make it diagonal
  int tau2_size = tau2.size();
  matrix<Type> tau2_diag(tau2_size,tau2_size);
  for(int i=0; i< tau2_size; i++){
    tau2_diag(i,i)=tau2(i);
  }
  matrix<Type> tau2inv = atomic::matinv(tau2_diag);
  //REPORT(tau2inv);
  Type nll=0; 
  // Get the first...
  //matrix<Type> m = z(1);
  //REPORT(m);
  //vector<Type> v = z2(1);
  //REPORT(v);
  for(int j=0; j< gammaj.size(); j++) {
    matrix<Type> mj = z2invsj(j) + tau2inv;
    matrix<Type> mj_inv = mj.inverse();
    vector<Type> qj = mj_inv * dj(j);
    nll -= gammaj(j) * ((dj(j)*qj).sum() - log(tau2.prod()) - atomic::logdet(mj)) ;
  }
  return nll;
}
