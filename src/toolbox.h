#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <RcppArmadillo.h>
//#include "utils.h"
using namespace arma;
using namespace Rcpp;

arma::sp_mat J_MF(arma::uword days, arma::uword m, arma::uword ld, arma::uword sA);
List PrinComp(arma::mat Y, arma::uword m);
List BReg(arma::mat X, arma::mat Y, bool Int, arma::mat Bp, double lam, double nu,
          arma::uword reps = 1000, arma::uword burn = 1000);
List BReg_diag(arma::mat X,  arma::mat Y, bool Int, arma::mat Bp, double lam, arma::vec nu,
               arma::uword reps = 1000, arma::uword burn = 1000); 
List DSmooth(arma::mat B,  arma::sp_mat Jb, arma::mat q, arma::mat H, arma::mat R, arma::mat Y,     
             arma::uvec freq, arma::uvec LD);
arma::mat DSMF( arma::mat B,  arma::sp_mat Jb, arma::mat q,  arma::mat H,  arma::mat R,  arma::mat Y,
                arma::uvec freq, arma::uvec LD);
arma::field<arma::mat> FSimMF(arma::mat B, arma::sp_mat Jb,  arma::mat q,  arma::mat H,  
                              arma::mat R,  arma::mat Y,  arma::uvec freq, arma::uvec LD);
arma::field<arma::mat> Identify(arma::mat H, arma::mat q);




#endif