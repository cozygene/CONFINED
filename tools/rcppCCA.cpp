#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <tuple>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat centerCPP(arma::mat& x){
    arma::rowvec xcmeans = mean(x);
    arma::mat xhold(x.n_rows, x.n_cols);
    for(int i = 0 ; i < x.n_rows; i++){
      xhold.row(i) = xcmeans;
    }
    return(x-xhold);  
}

// [[Rcpp::export]]
arma::mat getUV(arma::mat& x, arma::mat& eigspc)
{
  arma::mat x_centered = centerCPP(x);
   return(x_centered * eigspc);
 }

//' @export
// [[Rcpp::export]]
arma::field<arma::mat> doRCCA(arma::mat& x, arma::mat& y){
	arma::mat Cxy = cov(x, y); // Amat
	arma::mat Cxx = cov(x); // Bmat
	arma::mat Cyy = cov(y); // Cmat
	int p = Cxx.n_rows;
	int q = Cyy.n_rows;
	Cxx = (Cxx + Cxx.t()) / 2.;
	Cyy = (Cyy + Cyy.t()) / 2.;
	arma::mat Cxx_f = chol(Cxx);
	arma::mat Cyy_f = chol(Cyy);
	arma::mat Cxx_fi = inv(Cxx_f);
	arma::mat Cyy_fi = inv(Cyy_f);
	arma::mat D = Cxx_fi.t() * Cxy * Cyy_fi;
	arma::mat left;
	arma::mat right;
	arma::vec eigs;
	arma::mat A;
	arma::mat B;
	if (p >= q){
		svd(left, eigs, right, D);
		A = Cxx_fi * left;
		B = Cyy_fi * right;
	} 
	else{
		svd(left, eigs, right, D.t());
		A = Cxx_fi * right;
		B = Cyy_fi * left;
	}
	//std::vector<arma::mat> ret_vec;
	//ret_vec.push_back(A);
	//ret_vec.push_back(B); 
	arma::field<arma::mat> ret_matrices(3);
	ret_matrices(0) = A;
	ret_matrices(1) = B;
	ret_matrices(2) = eigs;
	return ret_matrices;
}

// [[Rcpp::export]]
arma::mat getQCPP(arma::mat& x){
  arma::mat Q, R;
  qr_econ(Q, R, x);
  return Q;
}

// [[Rcpp::export]]
arma::field<arma::mat> svdCPP(arma::mat x){
 arma::mat left, right;
 arma::vec eigs;
 svd_econ(left, eigs, right, x);
 arma::field<arma::mat> ret_matrices(3);
 ret_matrices(0) = left;
 ret_matrices(1) = diagmat(eigs);
 ret_matrices(2) = right;
 return ret_matrices; 
}
      
