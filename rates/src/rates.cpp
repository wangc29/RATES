#include <Rcpp.h>
#include <RcppEigen.h>
#include "rates_types.h"

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double logsumexp(NumericVector x){
  return(max(x)+log(sum(exp(x-max(x)))));
}
 
// [[Rcpp::export]]
double bjR1(VectorXd betajk_j, VectorXd sdjk2_j, MatrixXd zkl,VectorXd tau2){
  static const double pi = 3.14159265358979323846;
  int k_j=betajk_j.size();
  MatrixXd sd_kk = sdjk2_j.asDiagonal();
  MatrixXd M_jp1 = zkl.adjoint() * sd_kk.inverse() * zkl;
  MatrixXd M_jp2 =tau2.cwiseInverse().asDiagonal();
  MatrixXd M_j = M_jp1 + M_jp2;
  VectorXd d_j = zkl.adjoint() * sd_kk.inverse() * betajk_j;
  double term1 = -0.5*k_j*log(2*pi) - 0.5*sdjk2_j.array().log().sum() - 0.5*tau2.array().log().sum() - 0.5*log(M_j.determinant());
  double term2 = -0.5* ((betajk_j.adjoint() * sd_kk.inverse() * betajk_j).sum() - (d_j.adjoint() * M_j.inverse() * d_j).sum());
  return term1+term2;
 }


// [[Rcpp::export]]
double bjR0(NumericVector betajk_j,NumericVector sdjk2_j,double alpha,double lambda){
  int k_j = betajk_j.size();
  NumericVector bjkR0 = NumericVector(k_j);
  for(int i=0;i< k_j; ++i){
    NumericVector llb0_inds(2L);
    llb0_inds[0]=log(1-lambda)+R::dnorm(betajk_j[i],0.0,sqrt(sdjk2_j[i]),true);
    llb0_inds[1]=log(lambda)+R::dnorm(betajk_j[i],0.0,sqrt(alpha*sdjk2_j[i]),true);
    bjkR0[i] = logsumexp(llb0_inds);
  }
  return sum(bjkR0);
}


// [[Rcpp::export]]
NumericVector deltis(NumericVector betajk_j, NumericVector sdjk2_j,double lambda, double alpha) {
  int k_j = betajk_j.size();
  NumericVector OjkR1(k_j);
  for(int i=0;i<k_j;++i){
    if(!Rcpp::NumericVector::is_na(betajk_j[i])){
      NumericVector llb0_inds(2L);
      llb0_inds[0]=log(1-lambda)+R::dnorm(betajk_j[i],0.0,sqrt(sdjk2_j[i]),true);
      llb0_inds[1]=log(lambda)+R::dnorm(betajk_j[i],0.0,sqrt(alpha*sdjk2_j[i]),true);
      OjkR1[i] = exp(llb0_inds[1])/exp(logsumexp(llb0_inds));
      if(Rcpp::NumericVector::is_na(OjkR1[i])){
        double m=which_max(llb0_inds) ; 
        OjkR1[i]= m;
        }
    } else{
      OjkR1[i]=NA_REAL;
      }
    }
  return OjkR1;
  }















