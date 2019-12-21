// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "utils.h"
#include "toolbox.h"
using namespace arma;
using namespace Rcpp;

// A note on dimensions: these programs are written so that the long axis of the data is in rows.
// That means for the complete data set time is in rows.
// When the program uses an individual observation (i.e. one period) it takes the transpose of the data so that companion matrices and the like keep their standard form.
// That is, for a single period series are also indexed by rows.

// v 2019.04.27

// //cumsum for mixed frequency data

// arma::vec mf_cumsum(arma::uword fq,
//                     arma::vec y){
//   for(uword j = fq; j<y.n_elem; j++){
//     y(j) = y(j) + y(j-fq);
//   }
//   return(y);
// }

// [[Rcpp::export]]
List EstDFM(      arma::mat B,     // transition matrix
                  arma::mat Bp,    // prior for B
                  arma::sp_mat Jb, // aggrigations for transition matrix
                  double lam_B,    // prior tightness on transition matrix
                  arma::mat q,     // covariance matrix of shocks to states
                  double nu_q,     // prior deg. of freedom for variance of shocks in trans. eq.
                  arma::mat H,     // measurement equation
                  arma::mat Hp,    //prior for H
                  double lam_H,    // prior tightness on obs. equation
                  arma::vec R,     // covariance matrix of shocks to observables; Y are observations
                  arma::vec nu_r,     //prior degrees of freedom for elements of R used to normalize
                  arma::mat Y,     // data
                  arma::uvec freq, // frequency denoted as number of high frequency periods in a low frequency period
                  arma::uvec LD,  // 0 for level data and 1 for first difference
                  bool store_Y, //Store distribution of Y?
                  arma::uvec store_idx, // index to store distribution of predicted values
                  arma::uword reps = 1000, //repetitions
                  arma::uword burn = 500,
                  bool verbose = true){ //burn in periods

  // preliminaries

  // preliminaries
  uword m  = B.n_rows;
  uword p  = Jb.n_cols/m;
  uword T  = Y.n_rows - p;
  uword k  = H.n_rows;
  uword sA    = m*p; // size A matrix
  uword sB    = B.n_cols; // columns of transition matrix B
  mat Lam_B   = lam_B*eye<mat>(sB,sB);
  mat Lam_H   = lam_H*eye<mat>(m,m);

  // ----- Priors -------
  //mat Bp(m,sA, fill::zeros);  //prior for B
  //mat Hp(k,m,fill::zeros); //prior for M and H (treated as the same parameters)

  // Initialize variables
  mat v_1, V_1, mu, Mu, Beta, scale, xx, yy, Zd, Zs, Zsim, Yd, Ys, Rmat;
  mat Ht(m,m,fill::zeros), aa;
  sp_mat Jh;
  vec Yt;
  uvec ind;
  uword count_reps;
  double scl;
  double ev = 2;
  cube Bstore(m,sB,reps); //store draws for B
  cube Hstore(k,m,reps);  //store draws for H
  cube Qstore(m,m,reps);  //store draws for Q
  mat  Rstore(k,reps);    //R is diagonal so only diagonals stored (hence matrix not cube)
  cube Ystore;
  mat  Ymedian;
  if(store_Y){
    Ystore = zeros<cube>(T+p, store_idx.n_elem, reps);
    Ymedian = zeros<mat>(T+p, store_idx.n_elem);
  }
  List Out;
  field<mat> FSim;

  vec eigval;
  mat eigvec;
  cx_vec eigval_cx;
  cx_mat eigvec_cx;

  mat Ytmp = Y;
  Ytmp.shed_rows(0,p-1); //shed initial values to match Z

  for(uword rep = 0; rep<burn; rep++){

    Rcpp::checkUserInterrupt();

    if(verbose){
      Rcpp::Rcout << "\rProgress: " << round(100*rep/(burn + reps)) << "% (burning)";
    }

    // --------- Sample Factors given Data and Parameters ---------

    // Sampling follows Durbin and Koopman 2002/2012

    Rmat  = diagmat(R); // Make a matrix out of R to plug in to DSimMF
    // Draw observations Y^star and Z^star
    FSim  = FSimMF(B, Jb, q, H, Rmat, Y, freq, LD);
    Zd    = FSim(0); //draw for Z
    Yd    = FSim(1); //draw for Y
    Ys    = Y-Yd;
    // Smooth using Ys (i.e. Y^star)
    Zs    = DSMF(B, Jb, q, H, Rmat, Ys, freq, LD);
    Zsim  = Zs + Zd; // Draw for factors

    Zsim.shed_rows(0,p-1); //shed initial values (not essential)

    // -------- Sample Parameters given Factors -------

    // For H, M, and R

    //For observations used to normalize
    for(uword j=0; j<m; j++){ //loop over variables
      Yt   = Ytmp.col(j);
      ind  = find_finite(Yt);   //find non-missing values
      yy   = Yt(ind);         //LHS variable
      Jh = J_MF(freq(j), m, LD(j), sA);
      xx = Zsim.rows(ind)*trans(Jh);
      V_1   = trans(xx)*xx+Lam_H;
      V_1   = (trans(V_1)+V_1)/2;
      V_1   = inv_sympd(V_1);
      mu    = V_1*(trans(xx)*yy+Lam_H*trans(Hp.row(j)));
      scl   = 1 + as_scalar(trans(yy-xx*mu)*(yy-xx*mu)+trans(mu-trans(Hp.row(j)))*Lam_H*(mu-trans(Hp.row(j)))); // prior scale is the 1, + as_scalar(junk) comes from the posterior
      R(j)  = invchisq(nu_r(j)+yy.n_elem,scl); //Draw for r
      Beta  = mvrnrm(1, mu, V_1*R(j));
      Ht.row(j) = trans(Beta.col(0));
      //H.row(j) = trans(Beta.col(0));
    }

    //Rotate and scale the factors to fit our normalization for H
    Zsim     = Zsim*kron(eye<mat>(p,p),trans(Ht));

    //For observations not used to normalize
    for(uword j=m; j<k; j++){ //loop over variables
      Yt   = Ytmp.col(j);
      ind  = find_finite(Yt);   //find non-missing values
      yy   = Yt(ind);         //LHS variable
      Jh = J_MF(freq(j), m, LD(j), sA);
      xx = Zsim.rows(ind)*trans(Jh);
      V_1   = trans(xx)*xx+Lam_H;
      V_1   = (trans(V_1)+V_1)/2;
      V_1   = inv_sympd(V_1);
      mu    = V_1*(trans(xx)*yy+Lam_H*trans(Hp.row(j)));
      scl   = 1 + as_scalar(trans(yy-xx*mu)*(yy-xx*mu)+trans(mu-trans(Hp.row(j)))*Lam_H*(mu-trans(Hp.row(j)))); // prior scale is the 1, + as_scalar(junk) comes from the posterior
      R(j)  = invchisq(nu_r(j)+yy.n_elem,scl); //Draw for r
      Beta  = mvrnrm(1, mu, V_1*R(j));
      H.row(j) = trans(Beta.col(0));
    }

    // For B and q

    yy    = Zsim.cols(0,m-1);
    yy.shed_row(0);
    xx    = Zsim*trans(Jb);
    xx.shed_row(T-1);
    v_1   = trans(xx)*xx+Lam_B;
    v_1   = (trans(v_1)+v_1)/2;
    v_1   = inv_sympd(v_1);
    Mu    = v_1*(trans(xx)*yy+Lam_B*trans(Bp));
    scale = eye(m,m)+trans(yy-xx*Mu)*(yy-xx*Mu)+trans(Mu-trans(Bp))*Lam_B*(Mu-trans(Bp)); // eye(k) is the prior scale parameter for the IW distribution and eye(k)+junk the posterior.
    scale = (scale+trans(scale))/2;
    q     = rinvwish(1,nu_q+T,scale); //Draw for q
    mu    = vectorise(trans(Mu));  //posterior mean for beta
    V_1   = kron(v_1,q);  //covariance of vectorized parameters beta
    ev    = 2;
    count_reps = 1;
    do{ //this loop ensures the fraw for B is stationary --- non stationary draws are rejected
      Beta  = mvrnrm(1,mu,V_1); //Draw for B
      B     = reshape(Beta,m,sB); //recovering original dimensions
      // Check wheter B is stationary and reject if not
      aa      = comp_form(B); //
      eig_gen(eigval_cx, eigvec_cx, aa);
      ev     = as_scalar(max(abs(eigval_cx)));
      Rcpp::checkUserInterrupt();
      if(count_reps == 10000){
        Rcpp::Rcout << "Draws Non-Stationary" << endl;
      }
      if(count_reps == 30000 || B.has_nan()){
        throw("Draws Non-Stationary"); //break program if still no stationary draws
      }
      count_reps = count_reps+1;
    } while(ev>1);
  }

  // ------------------ Sampling Loop ------------------------------------

  for(uword rep = 0; rep<reps; rep++){

    Rcpp::checkUserInterrupt();

    if(verbose){
      Rcpp::Rcout << "\rProgress: " << round(100*(burn + rep)/(burn + reps)) << "% (sampling)";
    }

    // --------- Sample Factors given Data and Parameters

    // Sampling follows Durbin and Koopman 2002/2012

    Rmat  = diagmat(R); // Make a matrix out of R to plug in to DSimMF
    // Draw observations Y^star and Z^star
    FSim  = FSimMF(B, Jb, q, H, Rmat, Y, freq, LD);
    Zd    = FSim(0); //draw for Z
    Yd    = FSim(1); //draw for Y
    Ys    = Y-Yd;
    // Smooth using Ys (i.e. Y^star)
    Zs    = DSMF(B, Jb, q, H, Rmat, Ys, freq, LD);
    Zsim  = Zs + Zd; // Draw for factors

    if(store_Y){
      for(uword j=0;j<store_idx.n_elem;j++){
        Jh        = J_MF(freq(store_idx[j]), m, LD(store_idx[j]), sA);
        Ymedian.col(store_idx[j]) = Zsim*trans(Jh)*trans(H.row(store_idx[j]));
      }
      Ystore.slice(rep) = Ymedian;
    }

    Zsim.shed_rows(0,p-1);

    // -------- Sample Parameters given Factors

    // For H, M, and R

    //For observations used to normalize
    for(uword j=0; j<m; j++){ //loop over variables
      Yt   = Ytmp.col(j);
      ind  = find_finite(Yt);   //find non-missing values
      yy   = Yt(ind);         //LHS variable
      Jh = J_MF(freq(j), m, LD(j), sA);
      xx = Zsim.rows(ind)*trans(Jh);
      V_1   = trans(xx)*xx+Lam_H;
      V_1   = (trans(V_1)+V_1)/2;
      V_1   = inv_sympd(V_1);
      mu    = V_1*(trans(xx)*yy+Lam_H*trans(Hp.row(j)));
      scl   = 1 + as_scalar(trans(yy-xx*mu)*(yy-xx*mu)+trans(mu-trans(Hp.row(j)))*Lam_H*(mu-trans(Hp.row(j)))); // prior scale is the 1, + as_scalar(junk) comes from the posterior
      R(j)  = invchisq(nu_r(j)+yy.n_elem,scl); //Draw for r
      Beta  = mvrnrm(1, mu, V_1*R(j));
      Ht.row(j) = trans(Beta.col(0));
      //H.row(j) = trans(Beta.col(0));
    }

    //Rotate and scale the factors to fit our normalization for H
    Zsim     = Zsim*kron(eye<mat>(p,p),trans(Ht));

    //For observations not used to normalize
    for(uword j=m; j<k; j++){ //loop over variables
      Yt   = Ytmp.col(j);
      ind  = find_finite(Yt);   //find non-missing values
      yy   = Yt(ind);         //LHS variable
      Jh = J_MF(freq(j), m, LD(j), sA);
      xx = Zsim.rows(ind)*trans(Jh);
      V_1   = trans(xx)*xx+Lam_H;
      V_1   = (trans(V_1)+V_1)/2;
      V_1   = inv_sympd(V_1);
      mu    = V_1*(trans(xx)*yy+Lam_H*trans(Hp.row(j)));
      scl   = 1 + as_scalar(trans(yy-xx*mu)*(yy-xx*mu)+trans(mu-trans(Hp.row(j)))*Lam_H*(mu-trans(Hp.row(j)))); // prior scale is the 1, + as_scalar(junk) comes from the posterior
      R(j)  = invchisq(nu_r(j)+yy.n_elem,scl); //Draw for r
      Beta  = mvrnrm(1, mu, V_1*R(j));
      H.row(j) = trans(Beta.col(0));
    }

    // For B and q

    yy    = Zsim.cols(0,m-1);
    yy.shed_row(0);
    xx    = Zsim*trans(Jb);
    xx.shed_row(T-1);
    v_1   = trans(xx)*xx+Lam_B;
    v_1   = (trans(v_1)+v_1)/2;
    v_1   = inv_sympd(v_1);
    Mu    = v_1*(trans(xx)*yy+Lam_B*trans(Bp));
    scale = eye(m,m)+trans(yy-xx*Mu)*(yy-xx*Mu)+trans(Mu-trans(Bp))*Lam_B*(Mu-trans(Bp)); // eye(k) is the prior scale parameter for the IW distribution and eye(k)+junk the posterior.
    scale = (scale+trans(scale))/2;
    q     = rinvwish(1,nu_q+T,scale); //Draw for q
    mu    = vectorise(trans(Mu));  //posterior mean for beta
    V_1   = kron(v_1,q);  //covariance of vectorized parameters beta
    ev    = 2;
    count_reps = 1;
    do{ //this loop ensures the fraw for B is stationary --- non stationary draws are rejected
      Beta  = mvrnrm(1,mu,V_1); //Draw for B
      B     = reshape(Beta,m,sB); //recovering original dimensions
      // Check wheter B is stationary and reject if not
      aa    = comp_form(B); //
      eig_gen(eigval_cx, eigvec_cx, aa);
      ev    = as_scalar(max(abs(eigval_cx)));
      Rcpp::checkUserInterrupt();
      if(count_reps == 10000){
        Rcpp::Rcout << "Draws Non-Stationary" << endl;
      }
      if(count_reps == 30000){
        throw("Draws Non-Stationary" || B.has_nan()); //break program if still no stationary draws
      }
      count_reps = count_reps+1;
    } while(ev>1);

    Bstore.slice(rep) = B;
    Qstore.slice(rep) = q;
    Hstore.slice(rep) = H;
    Rstore.col(rep)   = R;

  }

  Rcpp::Rcout << "\r                          \r";

  //Getting posterior medians

  //For B
  for(uword rw=0;rw<m;rw++){
    for(uword cl=0;cl<sB;cl++){
      B(rw,cl) = as_scalar(median(vectorise(Bstore.tube(rw,cl))));
    }
  }
  //For H
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<m;cl++){
      H(rw,cl) = as_scalar(median(vectorise(Hstore.tube(rw,cl))));
    }
  }
  //For q
  for(uword rw=0;rw<m;rw++){
    for(uword cl=0;cl<m;cl++){
      q(rw,cl) = as_scalar(median(vectorise(Qstore.tube(rw,cl))));
    }
  }
  //For R
  for(uword rw=0;rw<k;rw++){
    R(rw) = median(Rstore.row(rw));
  }
  //For Y
  if(store_Y){
    for(uword rw=0;rw<T;rw++){
      for(uword cl=0;cl<store_idx.n_elem;cl++)
        Ymedian(rw,cl) = as_scalar(median(vectorise(Ystore.tube(rw,cl))));
    }
  }else{
    Ystore = zeros<cube>(0,0,0);
    Ymedian = zeros<mat>(0,0);
  }

  Out["B"]  = B;
  Out["H"]  = H;
  Out["Q"]  = q;
  Out["R"]  = R;
  Out["Bstore"]  = Bstore;
  Out["Hstore"]  = Hstore;
  Out["Qstore"]  = Qstore;
  Out["Rstore"]  = Rstore;
  Out["Zsim"]  = Zsim;
  Out["Ystore"] = Ystore;
  Out["Y_median"] = Ymedian;

  return(Out);
}

//-------------------------------------------------
// --------- Maximum Likelihood Programs ----------
//-------------------------------------------------


// [[Rcpp::export]]
List Ksmoother(arma::sp_mat A,  // companion form of transition matrix
               arma::sp_mat Q,  // covariance matrix of shocks to states
               arma::sp_mat HJ, // measurement equation
               arma::mat R,     // covariance matrix of shocks to observables; Y are observations
               arma::mat Y){    //data
  // preliminaries
  uword T  = Y.n_rows;
  uword sA = A.n_rows;

  // specifying initial values
  mat P0, P1, S, C;
  P0 = 100000*eye<mat>(sA,sA);
  P1 = P0;

  //Declairing variables for the filter
  cube P0str(sA,sA,T);
  P0str.slice(0) = P0;
  cube P1str(sA,sA,T+1);
  P1str.slice(0) = P1;
  field<mat> Kstr(T,1);
  field<vec> PEstr(T);
  mat Z1, VarY, Z, Zp(sA, 1, fill::zeros), Lik, K, G, Rn;
  mat Yf = Y;
  sp_mat Hn;
  vec PE, Yt, Yn, Yp;
  uvec ind;
  Lik << 0;
  double tmp;
  double tmpp;
  Z.zeros(T,sA);
  Z1.zeros(T+1,sA);

  for(uword t=0; t<T; t++) {
    Rcpp::checkUserInterrupt();
    //Allowing for missing Y values
    Yt     = trans(Y.row(t));
    ind    = find_finite(Yt);
    Yn     = Yt(ind);
    // if nothing is observed
    if(Yn.is_empty()){
      Z.row(t) = trans(Zp);
      P0       = P1;
      P0str.slice(t) = P0;
    } else {
      Hn     = sp_rows(HJ,ind);
      Rn     = R.rows(ind);
      Rn     = Rn.cols(ind);
      Yp     = Hn*Zp; //prediction step for Y
      S      = Hn*P1*trans(Hn)+Rn; //variance of Yp
      S      = symmatu((S+trans(S))/2);
      C      = P1*trans(Hn); //covariance of Zp Yp
      K      = trans(solve(S,trans(C))); //Kalman Gain
      PE         = Yn-Yp; // prediction error
      PEstr(t)   = PE;
      Kstr(t,0)  = K;
      Z.row(t)   = trans(Zp+K*PE); //updating step for Z
      P0     = P1-C*solve(S,trans(C)); // variance Z(t+1)|Y(1:t+1)
      P0     = (P0+trans(P0))/2;
      P0str.slice(t) = P0;
      log_det(tmp,tmpp,S);
      Lik    = -.5*tmp-.5*trans(PE)*solve(S,PE)+Lik;
    }
      //next period variables
      Zp               = A*trans(Z.row(t)); //prediction for Z(t+1)
      Z1.row(t+1)      = trans(Zp);
      P1               = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
      P1               = (P1+trans(P1))/2;
      P1str.slice(t+1) = P1;
  }

  //Declairing additional variables
  mat Zs = Z;
  cube Ps(sA,sA,T);
  Ps = P0str;

  //Note indexing starts at 0 so the last period is T-1

  for(uword t=T-1; t>0; t = t-1) {
    G = P0str.slice(t-1)*trans(A)*inv(P1str.slice(t)); //inv_sympd
    Zs.row(t-1)   = Z.row(t-1) + (Zs.row(t)-Z1.row(t))*trans(G);
    P0            = P0str.slice(t-1)-G*(P1str.slice(t)-Ps.slice(t))*trans(G);
    Ps.slice(t-1) = (P0+trans(P0))/2;
  }

  mat Yhat   = Zs*trans(HJ);
  Yf.elem(find_nonfinite(Y)) = Yhat.elem(find_nonfinite(Y));

  List Out;
  Out["Lik"]  = Lik;
  Out["Yf"]   = Yf;
  Out["Ys"]   = Yhat;
  Out["Zz"]   = Z;
  Out["Z"]    = Zs;
  Out["Kstr"] = Kstr;
  Out["PEstr"]= PEstr;
  Out["Ps"]   = Ps;
  Out["P1str"] = P1str;
  return(Out);
}

// [[Rcpp::export]]
List KestExact(arma::sp_mat A,
               arma::sp_mat Q,
               arma::mat H,
               arma::mat R,
               arma::mat Y,
               arma::vec itc,
               arma::uword m,
               arma::uword p){


  uword T  = Y.n_rows;
  uword k  = Y.n_cols;

  // Helper matrix J

  sp_mat J(m,m*(p+1));
  J(span(0,m-1),span(0,m-1)) = speye<sp_mat>(m,m);

  sp_mat   HJ  = MakeSparse(H*J);

  // Removing intercept terms from Y

  mat Ytmp  = Y - kron(ones<mat>(T,1),trans(itc));

  List Smth = Ksmoother(A, Q, HJ, R, Ytmp);

  mat Z     = Smth["Z"];
  cube Ps   = Smth["Ps"];
  mat Lik   = Smth["Lik"];

  mat xx, Zx, axz, XZ, azz, ZZ, axx, tmp, B;

  //For B and qB

  xx  = Z.cols(span(0,m-1));
  Zx  = Z.cols(span(m,(p+1)*m-1));
  axz = sum(Ps(span(0,m-1),span(m,(p+1)*m-1),span::all),2);
  XZ  = trans(Zx)*xx + trans(axz);
  azz = sum(Ps(span(m,(p+1)*m-1),span(m,(p+1)*m-1),span::all),2);
  ZZ  = trans(Zx)*Zx + azz;
  B = trans(solve(ZZ,XZ));
  A(span(0,m-1),span(0,m*p-1)) = B;

  axx = sum(Ps(span(0,m-1),span(0,m-1),span::all),2);
  mat q = (trans(xx-Zx*trans(B))*(xx-Zx*trans(B)) + axx - axz*trans(B) - B*trans(axz) + B*azz*trans(B) )/T;

  Q(span(0,m-1),span(0,m-1))   = q;

  //For H, R

  uword n_elm;
  uvec ind;
  vec yy, y, h;
  mat x, XX;

  xx     =  join_horiz( ones<mat>(T,1), Z*trans(J) ); // ones are for the intercept term

  for(uword j=0; j<k; j++) {
    yy    =  Y.col(j);
    ind   = find_finite(yy);
    y     = yy(ind);
    x     = xx.rows(ind);
    n_elm = y.size();
    // This section is super clunky... should be able to clean it up somehow but .slices does not accept non-contiguous indexes
    axx   = zeros<mat>(m+1,m+1);
    for(uword i = 0; i<n_elm; i++){
      axx(span(1,m),span(1,m))   = J*Ps.slice(ind(i))*trans(J)+axx(span(1,m),span(1,m));
    }
    // End clunky bit
    XX        = trans(x)*x+axx;
    h         = solve(XX,trans(x)*y);
    itc(j)    = h(0);
    H.row(j)  = trans(h(span(1,m)));
    R(j,j)    = as_scalar( trans(y-x*h)*(y-x*h)  + trans(h(span(1,m)))*axx(span(1,m),span(1,m))*h(span(1,m)) )/n_elm;
  }

  //Normalization --- cholesky ordering
  mat Thet, ThetI;
  Thet  = chol(q, "lower");
  ThetI = inv(Thet);

  //Normalization
  H      = H*Thet;
  tmp    = kron(eye<mat>(p,p),ThetI)*A(span(0,m*p-1),span(0,m*p-1))*kron(eye<mat>(p,p),Thet);
  A(span(0,m-1),span(0,m*p-1))   = tmp(span(0,m-1),span(0,m*p-1));
  Q(span(0,m-1),span(0,m-1))     = ThetI*Q(span(0,m-1),span(0,m-1))*trans(ThetI);
  mat X  = Z.cols(0,m-1)*trans(ThetI);

  mat Ys = X*trans(H);

  List Out;
  Out["A"]    = A;
  Out["Q"]    = Q;
  Out["H"]    = H;
  Out["R"]    = R;
  Out["Lik"]  = Lik;
  Out["X"]    = X;
  Out["itc"]  = itc;
  Out["Ys"]   = Ys;

  return(Out);
}








