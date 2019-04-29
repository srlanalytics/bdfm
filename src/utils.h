#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
arma::mat QuickReg(arma::mat X, arma::mat Y);
List UVreg(arma::vec x, arma::vec y, arma::uword rm_outlier = 0);
arma::sp_mat MakeSparse(arma::mat A);
arma::sp_mat sp_rows(arma::sp_mat A, arma::uvec r);
arma::sp_mat sp_cols(arma::sp_mat A, arma::uvec r);
arma::sp_mat sprow(arma::sp_mat A, arma::mat a, arma::uword r);
arma::mat comp_form(arma::mat B);
arma::mat mvrnrm(int n, arma::vec mu, arma::mat Sigma);
arma::cube rinvwish(int n, int v, arma::mat S);
double invchisq(double nu, double scale);
arma:: mat stack_obs(arma::mat nn, arma::uword p, arma::uword r = 0);


#endif