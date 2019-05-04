// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.h"
using namespace arma;
using namespace Rcpp;


//Return the appropriate mixed frequency helper matrix
// [[Rcpp::export]]
arma::sp_mat J_MF(arma::uword days, //number of high frequency periods in the low frequency period
                  arma::uword m,    //number of factors
                  arma::uword ld,   //type --- either level or difference
                  arma::uword sA){  //total number of columns (i.e. number of factors or size of A matrix)
  mat jm;
  if(days == 1){
    jm = eye<mat>(m,m);
  }else{
    if(ld == 0){
      mat weight = ones(1,days)/days;
      jm     = kron(weight,eye<mat>(m,m));
    }
    else if(ld == 1){
      mat weight(1,2*days-1);
      weight.cols(0,days-1) = trans(regspace(1,days))/days;
      weight.cols(days,2*days-2) = trans(regspace(days-1,1))/days;
      jm   = kron(weight,eye<mat>(m,m));
    }
  }
  sp_mat Jm(m,sA);
  Jm(span(0,m-1),span(0,jm.n_cols-1)) =  MakeSparse(jm);
  return(Jm);
}

//Principal Components
// [[Rcpp::export]]
List PrinComp(arma::mat Y,     // Observations Y
              arma::uword m){   // number of components
  
  //Preliminaries
  uword k = Y.n_cols;
  
  //Covariance matrix (may not be PSD due to missing data)
  vec Yj, Yi, Yt;
  mat Sig(k,k,fill::zeros);
  uvec ind;
  double tmp;
  for(uword j=0; j<k; j++){ //loop over variables
    Yj   = Y.col(j);
    for(uword i = 0; i<=j; i++){
      Yi  = Y.col(i);
      Yt  = Yj%Yi; //element wise multiplication
      ind = find_finite(Yt);
      tmp = sum(Yt(ind))/ind.n_elem;
      Sig(j,i) = tmp;
      Sig(i,j) = tmp;
    }
  }
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sig);
  eigvec = fliplr(eigvec); //put in decending order
  mat loadings   = eigvec.cols(0,m-1);
  mat components = Y*loadings;
  
  List Out;
  
  Out["Sig"]        = Sig;
  Out["loadings"]   = loadings;
  Out["components"] = components;
  
  return(Out);
}

//Bayesian linear regression --- does NOT accept missing variables
// [[Rcpp::export]]
List BReg(arma::mat X,   // RHS variables
          arma::mat Y,   // LHS variables
          bool Int,      // Estimate intercept term?
          arma::mat Bp,  // prior for B
          double lam,    // prior tightness
          double nu,     //prior "deg of freedom"
          arma::uword reps = 1000,
          arma::uword burn = 1000){
  
  uword k    = Y.n_cols;
  uword m    = X.n_cols;
  uword T    = Y.n_rows;
  mat Lam    = lam*eye<mat>(m,m);
  mat tmp;
  
  if(Int){ //if we should estimate an intercept term
    m        = m+1;
    Lam      = lam*eye<mat>(m,m);
    Lam(0,0) = 0;
    tmp      = zeros<mat>(k,1);
    Bp       = join_horiz(tmp, Bp);
    tmp      = ones<mat>(T,1);
    X        = join_horiz(tmp,X);
  }
  
  //declairing variables
  cube Bstore(k,m,reps), Qstore(k,k,reps);
  mat v_1, V_1, Mu, B, Beta, scale, q;
  vec mu;
  
  //Burn Loop
  
  for(uword rep = 0; rep<burn; rep++){
    
    Rcpp::checkUserInterrupt();
    
    v_1   = trans(X)*X+Lam;
    v_1   = (trans(v_1)+v_1)/2;
    v_1   = inv_sympd(v_1);
    Mu    = v_1*(trans(X)*Y+Lam*trans(Bp));
    scale = eye(k,k)+trans(Y-X*Mu)*(Y-X*Mu)+trans(Mu-trans(Bp))*Lam*(Mu-trans(Bp)); // eye(k) is the prior scale parameter for the IW distribution and eye(k)+junk the posterior.
    scale = (scale+trans(scale))/2;
    q     = rinvwish(1,nu+T,scale); // Draw for q
    mu    = vectorise(trans(Mu));   // vectorized posterior mean for beta
    V_1   = kron(v_1,q);            // covariance of vectorized parameters beta
    Beta  = mvrnrm(1,mu,V_1);       //Draw for B
    B     = reshape(Beta,k,m);      //recovering original dimensions
  }
  
  //Sampling Loop
  
  for(uword rep = 0; rep<reps; rep++){
    
    Rcpp::checkUserInterrupt();
    
    v_1   = trans(X)*X+Lam;
    v_1   = (trans(v_1)+v_1)/2;
    v_1   = inv_sympd(v_1);
    Mu    = v_1*(trans(X)*Y+Lam*trans(Bp));
    scale = eye(k,k)+trans(Y-X*Mu)*(Y-X*Mu)+trans(Mu-trans(Bp))*Lam*(Mu-trans(Bp)); // eye(k) is the prior scale parameter for the IW distribution and eye(k)+junk the posterior.
    scale = (scale+trans(scale))/2;
    q     = rinvwish(1,nu+T,scale); // Draw for q
    mu    = vectorise(trans(Mu));   // vectorized posterior mean for beta
    V_1   = kron(v_1,q);            // covariance of vectorized parameters beta
    Beta  = mvrnrm(1,mu,V_1);       //Draw for B
    B     = reshape(Beta,k,m);      //recovering original dimensions
    Qstore.slice(rep) = q;
    Bstore.slice(rep) = B;
  }
  
  //For B
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<m;cl++){
      B(rw,cl) = as_scalar(median(vectorise(Bstore.tube(rw,cl))));
    }
  }
  //For q
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<k;cl++){
      q(rw,cl) = as_scalar(median(vectorise(Qstore.tube(rw,cl))));
    }
  }
  
  List Out;
  Out["B"]   = B;
  Out["q"]   = q;
  Out["Bstore"]   = Bstore;
  Out["Qstore"]   = Qstore;
  
  return(Out);
  
}


//Bayesian linear regression with diagonal covariance to shocks. Accepts missing obs.
// [[Rcpp::export]]
List BReg_diag(arma::mat X,   // RHS variables
               arma::mat Y,   // LHS variables
               bool Int,      //Estimate intercept terms?
               arma::mat Bp,  // prior for B
               double lam,    // prior tightness
               arma::vec nu,  //prior "deg of freedom"
               arma::uword reps = 1000, //MCMC sampling iterations
               arma::uword burn = 1000){ //MCMC burn in iterations
  
  uword k    = Y.n_cols;
  uword m    = X.n_cols;
  uword T    = Y.n_rows;
  mat Lam    = lam*eye<mat>(m,m);
  mat tmp;
  
  if(Int){ //if we should estimate an intercept term
    m        = m+1;
    Lam      = lam*eye<mat>(m,m);
    Lam(0,0) = 0;
    tmp      = zeros<mat>(k,1);
    Bp       = join_horiz(tmp, Bp);
    tmp      = ones<mat>(T,1);
    X        = join_horiz(tmp,X);
  }
  
  //declairing variables
  cube Bstore(k,m,reps);
  mat v_1, Beta, Qstore(k,reps), xx;
  vec mu, q(k,fill::zeros);
  double scl;
  uvec ind_rm, ind;
  vec y, x;
  mat B(k, m, fill::zeros);
  uvec indX = find_nonfinite(X.col(0));
  if(m>1){
    for(uword j = 1; j<m; j++){
      x = X.col(j);
      indX = unique( join_cols(indX, find_nonfinite(x))); //index of missing X values
    }
  }
  
  //Burn Loop
  
  for(uword rep = 0; rep<burn; rep++){
    Rcpp::checkUserInterrupt();
    
    for(uword j=0; j<k; j++){
      y       = Y.col(j);
      ind     = unique(join_cols(indX,find_nonfinite(y))); //index of elements to remove
      xx      = X;
      
      //this seems a tedious way to shed non-contiguous indexes
      for(uword n = ind.n_elem; n>0; n--){
        xx.shed_row(ind(n-1));
        y.shed_row(ind(n-1));
      }
      
      v_1   = trans(xx)*xx+Lam;
      v_1   = (trans(v_1)+v_1)/2;
      v_1   = inv_sympd(v_1);
      mu    = v_1*(trans(xx)*y+Lam*trans(Bp.row(j)));
      scl   = as_scalar(trans(y-xx*mu)*(y-xx*mu)+trans(mu-trans(Bp.row(j)))*Lam*(mu-trans(Bp.row(j)))); // prior variance is zero... a little odd but it works
      q(j)  = invchisq(nu(j)+y.n_rows,scl); //Draw for r
      Beta  = mvrnrm(1, mu, v_1*q(j));
      B.row(j) = trans(Beta.col(0));
    }
    
  }
  
  // Sampling loop
  
  for(uword rep = 0; rep<reps; rep++){
    
    Rcpp::checkUserInterrupt();
    
    for(uword j=0; j<k; j++){
      y       = Y.col(j);
      ind     = unique(join_cols(indX,find_nonfinite(y))); //index of elements to remove
      xx      = X;
      //this seems a tedious way to shed non-contiguous indexes
      for(uword n = ind.n_elem; n>0; n--){
        xx.shed_row(ind(n-1));
        y.shed_row(ind(n-1));
      }
      v_1   = trans(xx)*xx+Lam;
      v_1   = (trans(v_1)+v_1)/2;
      v_1   = inv_sympd(v_1);
      mu    = v_1*(trans(xx)*y+Lam*trans(Bp.row(j)));
      scl   = as_scalar(trans(y-xx*mu)*(y-xx*mu)+trans(mu-trans(Bp.row(j)))*Lam*(mu-trans(Bp.row(j)))); // prior variance is zero... a little odd but it works
      q(j)  = invchisq(nu(j)+y.n_rows,scl); //Draw for r
      Beta  = mvrnrm(1, mu, v_1*q(j));
      B.row(j) = trans(Beta.col(0));
    }
    Qstore.col(rep)   = q;
    Bstore.slice(rep) = B;
    
  }
  
  
  //For B
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<m;cl++){
      B(rw,cl) = as_scalar(median(vectorise(Bstore.tube(rw,cl))));
    }
  }
  //For q
  
  for(uword rw=0;rw<k;rw++){
    q(rw) = as_scalar(median(Qstore.row(rw)));
  }
  
  
  List Out;
  Out["B"]   = B;
  Out["q"]   = q;
  Out["Bstore"]   = Bstore;
  Out["Qstore"]   = Qstore;
  
  return(Out);
}


// Disturbance smoother. Output is a list.
// [[Rcpp::export]]
List DSmooth(      arma::mat B,     // companion form of transition matrix
                   arma::sp_mat Jb, // helper matrix for transition equation
                   arma::mat q,     // covariance matrix of shocks to states
                   arma::mat H,     // measurement equation
                   arma::mat R,     // covariance matrix of shocks to observables; Y are observations
                   arma::mat Y,     //data
                   arma::uvec freq,  // frequency of each series (# low freq. periods in one obs)
                   arma::uvec LD){   // 0 if level, 1 if one diff.
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword p  = Jb.n_cols/m; //number of lags (must agree with lev/diff structure of data)
  uword k  = H.n_rows; //number of observables
  uword sA = m*p; //size of companion matrix A
  
  // For frequencies that do not change rows of HJ will be fixed.
  sp_mat HJ(k,sA);
  mat    hj;
  for(uword j = 0; j<k; j++){
    hj  = H(j,span::all)*J_MF(freq(j), m, LD(j), sA);
    HJ  = sprow(HJ,hj,j); //replace row j of HJ with vector hj
  }
  
  //Making the A matrix
  sp_mat BJb    = MakeSparse(B*Jb);
  sp_mat tmp_sp(sA-m,m);
  tmp_sp = join_horiz(speye<sp_mat>(sA-m,sA-m), tmp_sp);
  sp_mat A      = join_vert(BJb, tmp_sp);
  
  //Making the Q matrix
  mat qq(sA,sA,fill::zeros);
  qq(span(0,m-1),span(0,m-1)) = q;
  sp_mat Q(qq);
  
  //Using the long run variance
  mat P0, P1, S, C;
  double max_diff = 1;
  tmp_sp = Q;
  sp_mat P_new(sA,sA);
  sp_mat P_old(sA,sA);
  sp_mat P_diff;
  for(uword j = 0; j<1000; j++){
    P_new = tmp_sp + P_new;
    tmp_sp = A*tmp_sp*trans(A);
    P_diff = abs(P_new - P_old);
    max_diff = P_diff.max();
    P_old = P_new;
    if(max_diff>1e-5){
      break;
    }
  }
  P0  = P_new; //long run variance
  P1  = P0;
  
  //Declairing variables for the filter
  field<mat> Kstr(T); //store Kalman gain
  field<vec> PEstr(T); //store prediction error
  field<mat> Hstr(T), Sstr(T); //store H and S^-1
  mat VarY, ZP(T+1,sA,fill::zeros), Z(T,sA,fill::zeros), Zs(T,sA,fill::zeros), Lik, K, Rn, Si, tmp_mat;
  sp_mat Hn;
  vec PE, Yt, Yn, Yp;
  uvec ind, indM;
  double tmp;
  double tmpp;
  Lik << 0;
  vec Zp(sA, fill::zeros); //initialize to zero --- more or less arbitrary due to difuse variance
  
  mat zippo(1,1,fill::zeros);
  mat zippo_sA(sA,1);
  
  // -------- Filtering --------------------
  for(uword t=0; t<T; t++) {
    //Allowing for missing Y values
    Yt     = trans(Y.row(t));
    ind    = find_finite(Yt);
    Yn     = Yt(ind);
    // if nothing is observed
    if(Yn.is_empty()){
      Z.row(t) = trans(Zp);
      P0       = P1;
      Hstr(t)  = trans(zippo_sA);
      Sstr(t)  = zippo;
      PEstr(t) = zippo;
      Kstr(t)  = zippo_sA;
    } else{
      //if variables are observed
      Hn        = sp_rows(HJ,ind); //rows of HJ corresponding to observations
      Hstr(t)   = Hn; //Store for smoothing
      Rn        = R.rows(ind);  //rows of R corresponding to observations
      Rn        = Rn.cols(ind); //cols of R corresponding to observations
      Yp        = Hn*Zp; //prediction step for Y
      S         = Hn*P1*trans(Hn)+Rn; //variance of Yp
      S         = symmatu((S+trans(S))/2); //enforce pos. semi. def.
      Si        = inv_sympd(S); //invert S
      Sstr(t)   = Si; //sotre Si for smoothing
      K         = P1*trans(Hn)*Si; //Kalman gain
      PE        = Yn-Yp; // prediction error
      PEstr(t)  = PE; //store prediction error
      Kstr(t)   = K;  //store Kalman gain
      Z.row(t)  = trans(Zp+K*PE); //updating step for Z
      P0        = P1-P1*trans(Hn)*Si*Hn*P1; // variance Z(t+1)|Y(1:t+1)
      P0        = symmatu((P0+trans(P0))/2); //enforce pos semi def
      log_det(tmp,tmpp,S); //calculate log determinant of S for the likelihood
      Lik    = -.5*tmp-.5*trans(PE)*Si*PE+Lik; //calculate log likelihood
    }
    // Prediction for next period
    Zp  = A*trans(Z.row(t)); //prediction for Z(t+1) +itcZ
    ZP.row(t+1) = trans(Zp);
    P1     = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
    P1     = symmatu((P1+trans(P1))/2); //enforce pos semi def
  }
  ZP.shed_row(T);
  
  //Smoothing following Durbin Koopman 2001/2012
  mat r(T+1,sA,fill::zeros);
  mat L;
  
  //r is 1 indexed while all other variables are zero indexed
  for(uword t=T; t>0; t--) {
    L     = (A-A*Kstr(t-1)*Hstr(t-1));
    r.row(t-1) = trans(PEstr(t-1))*Sstr(t-1)*Hstr(t-1) + r.row(t)*L;
  }
  
  Zs.row(0)   = r.row(0)*P_new;
  
  //Forward again
  for(uword t = 0; t<T-1; t++){
    Zs.row(t+1)   = Zs.row(t)*trans(A) + r.row(t+1)*Q; //smoothed values of Z
  }
  
  mat Ys = Zs*trans(HJ); //fitted values of Y
  
  List Out;
  Out["Ys"]   = Ys;
  Out["Lik"]  = Lik;
  Out["Zz"]   = Z;
  Out["Z"]    = Zs;
  Out["Zp"]   = ZP;
  Out["Kstr"] = Kstr;
  Out["PEstr"]= PEstr;
  Out["r"]    = r;
  
  return(Out);
}

//Disturbance smoothing --- output is only smoothed factors for simulations
// [[Rcpp::export]]
arma::mat DSMF(           arma::mat B,     // companion form of transition matrix
                          arma::sp_mat Jb, // helper matrix for transition equation
                          arma::mat q,     // covariance matrix of shocks to states
                          arma::mat H,     // measurement equation
                          arma::mat R,     // covariance matrix of shocks to observables; Y are observations
                          arma::mat Y,     // data
                          arma::uvec freq, //frequency of each seres
                          arma::uvec LD){  // 0 for levels, 1 for first difference
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword p  = Jb.n_cols/m; //number of lags (must agree with lev/diff structure of data)
  uword k  = H.n_rows; //number of observables
  uword sA = m*p; //size of companion matrix A
  
  // For frequencies that do not change rows of HJ will be fixed.
  sp_mat HJ(k,sA);
  mat    hj;
  for(uword j = 0; j<k; j++){
    hj  = H(j,span::all)*J_MF(freq(j), m, LD(j), sA);
    HJ  = sprow(HJ,hj,j); //replace row j of HJ with vector hj
  }
  
  //Making the A matrix
  sp_mat BJb    = MakeSparse(B*Jb);
  sp_mat tmp_sp(sA-m,m);
  tmp_sp = join_horiz(speye<sp_mat>(sA-m,sA-m), tmp_sp);
  sp_mat A      = join_vert(BJb, tmp_sp);
  
  //Making the Q matrix
  mat qq(sA,sA,fill::zeros);
  qq(span(0,m-1),span(0,m-1)) = q;
  sp_mat Q(qq);
  
  //Using the long run variance
  mat P0, P1, S, C;
  double max_diff = 1;
  tmp_sp = Q;
  sp_mat P_new(sA,sA);
  sp_mat P_old(sA,sA);
  sp_mat P_diff;
  for(uword j = 0; j<1000; j++){
    P_new = tmp_sp + P_new;
    tmp_sp = A*tmp_sp*trans(A);
    P_diff = abs(P_new - P_old);
    max_diff = P_diff.max();
    P_old = P_new;
    if(max_diff>1e-5){
      break;
    }
  }
  P0  = P_new; //long run variance
  P1  = P0;
  
  //Declairing variables for the filter
  //mat P11 = P1; //output long run variancce for testing.
  field<mat> Kstr(T); //store Kalman gain
  field<vec> PEstr(T); //store prediction error
  field<mat> Hstr(T), Sstr(T); //store H and S^-1
  mat VarY, Z(T,sA,fill::zeros), Zs(T,sA,fill::zeros), Lik, K, Rn, Mn, Si, tmp_mat;
  sp_mat Hn;
  vec Z1, PE, Yt, Yn, Yp;
  uvec ind, indM;
  Lik << 0;
  double tmp;
  double tmpp;
  vec Zp(sA,fill::zeros); //initial factor values (arbitrary as variance difuse)
  
  mat zippo(1,1,fill::zeros);
  mat zippo_sA(sA,1,fill::zeros);
  
  // -------- Filtering --------------------
  for(uword t=0; t<T; t++) {
    //Allowing for missing Y values
    Yt     = trans(Y.row(t));
    ind    = find_finite(Yt);
    Yn     = Yt(ind);
    // if nothing is observed
    if(Yn.is_empty()){
      Z.row(t) = trans(Zp);
      P0       = P1;
      Hstr(t)  = trans(zippo_sA);
      Sstr(t)  = zippo;
      PEstr(t) = zippo;
      Kstr(t)  = zippo_sA;
    } else{
      //if variables are observed
      Hn        = sp_rows(HJ,ind);
      Hstr(t)   = Hn;
      Rn        = R.rows(ind);
      Rn        = Rn.cols(ind);
      Yp        = Hn*Zp; //prediction step for Y
      S         = Hn*P1*trans(Hn)+Rn; //variance of Yp
      S         = symmatu((S+trans(S))/2);
      Si        = inv_sympd(S);
      Sstr(t)   = Si;
      K         = P1*trans(Hn)*Si; //Kalman Gain
      PE        = Yn-Yp; // prediction error
      PEstr(t)  = PE;
      Kstr(t)   = K;
      Z.row(t)  = trans(Zp+K*PE); //updating step for Z
      P0        = P1-P1*trans(Hn)*Si*Hn*P1; // variance Z(t+1)|Y(1:t+1)
      P0        = symmatu((P0+trans(P0))/2);
      log_det(tmp,tmpp,S);
      Lik    = -.5*tmp-.5*trans(PE)*Si*PE+Lik;
      // Prediction for next period
      Zp     = A*trans(Z.row(t)); //prediction for Z(t+1) +itcZ
      //ZP.row(t) = trans(Zp);
      P1     = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
      P1     = symmatu((P1+trans(P1))/2);
    }
  }
  
  
  //Smoothing
  mat r(T+1,sA,fill::zeros);
  mat L;
  
  //t is 1 indexed, all other vars are 0 indexed
  for(uword t=T; t>0; t--) {
    L     = (A-A*Kstr(t-1)*Hstr(t-1));
    r.row(t-1) = trans(PEstr(t-1))*Sstr(t-1)*Hstr(t-1) + r.row(t)*L;
  }
  
  Zs.row(0)   = r.row(0)*P_new;
  
  //Forward again
  for(uword t = 0; t<T-1; t++){
    Zs.row(t+1)   = Zs.row(t)*trans(A) + r.row(t+1)*Q;
  }
  return(Zs);
}


//Forward recursion using draws for eps and e
// [[Rcpp::export]]
arma::field<arma::mat> FSimMF(    arma::mat B,     // companion form of transition matrix
                                  arma::sp_mat Jb, // helper matrix for transition equation
                                  arma::mat q,     // covariance matrix of shocks to states
                                  arma::mat H,     // measurement equation
                                  arma::mat R,     // covariance matrix of shocks to observables; Y are observations
                                  arma::mat Y,     // data
                                  arma::uvec freq, // frequency
                                  arma::uvec LD){  // level 0, or diff 1
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword sA = Jb.n_cols; //size of companion matrix A
  uword p  = sA/m; //number of lags
  uword k  = H.n_rows; //number of observables
  
  // For frequencies that do not change rows of HJ will be fixed.
  sp_mat HJ(k,sA);
  mat    hj;
  for(uword j = 0; j<k; j++){
    hj  = H(j,span::all)*J_MF(freq(j), m, LD(j), sA);
    HJ  = sprow(HJ,hj,j); //replace row j of HJ with vector hj
  }
  
  //Draw Eps (for observations) and E (for factors)
  vec mu_Eps(k,fill::zeros);
  vec mu_E(m,fill::zeros);
  mat Eps = trans(mvrnrm(T,mu_Eps,R));
  mat E   = trans(mvrnrm(T,mu_E,q));
  
  mat Q  = kron(eye<mat>(p,p), q);
  vec Z0(sA,fill::zeros);
  mat z0 = mvrnrm(1,Z0,Q); // Difuse initial conditions so not so important
  
  //Declairing variables for the forward recursion
  mat Z(T+1,sA), Yd(T,k);
  sp_mat Hn;
  vec x, yt, yd(k), eps;
  uvec ind;
  Z.row(0) = trans(z0);
  
  //Forward Recursion
  for(uword t=0; t<T; t++) {
    yt        = trans(Y.row(t));
    ind       = find_finite(yt); //identify missing values that they are replicated in simulated data
    Hn        = sp_rows(HJ,ind);
    yd.fill(datum::nan);
    eps       = trans(Eps.row(t));
    yd(ind)   = Hn*trans(Z.row(t)) + eps(ind);
    Yd.row(t) = trans(yd);
    //next period factors
    x      = B*Jb*trans(Z.row(t)) + trans(E.row(t)); //prediction for Z(t+1)
    Z(t+1,span(0,m-1))  = trans(x);
    if(p>1){
      Z(t+1,span(m,sA-1)) = Z(t,span(0,sA-m-1));
    }
  }
  Z.shed_row(T); //we don't use predictions for period T (zero indexed)
  
  
  field<mat> Out(3);
  Out(0)   = Z;
  Out(1)   = Yd;
  Out(2)   = Eps;
  
  return(Out);
}

// [[Rcpp::export]]
arma::field<arma::mat> Identify(arma::mat H,
                                arma::mat q){
  mat Sig_H;
  Sig_H = (trans(H)*H)/(H.n_rows);
  Sig_H = (trans(Sig_H)+Sig_H)/2;
  mat H_vec, Q_vec;
  vec H_val, Q_val;
  eig_sym(H_val, H_vec, Sig_H);
  mat PhiI = trans(H_vec)*diagmat(sqrt(H_val))*trans(H_vec);
  mat Q = PhiI*q*trans(PhiI);
  Q  = (Q + trans(Q))/2;
  eig_sym(Q_val, Q_vec, Q);
  Q_vec = Q_vec*diagmat(sign(Q_vec)); //enforce that diagonals must be positive
  mat Phi = H_vec*diagmat(1/sqrt(H_val))*trans(H_vec)*Q_vec;
  PhiI = trans(Q_vec)*PhiI;
  
  arma::field<arma::mat> Out(2);
  Out(0) = Phi;
  Out(1) = PhiI;
  
  return(Out);
}

