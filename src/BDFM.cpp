// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// A note on dimensions: these programs are written so that the long axis of the data is in rows.
// That means for the complete data set time is in rows.
// When the program uses an individual observation (i.e. one period) it takes the transpose of the data so that companion matrices and the like keep their standard form.
// That is, for a single period series are also indexed by rows.

// v 2019.02.15

// //cumsum for mixed frequency data

// arma::vec mf_cumsum(arma::uword fq,
//                     arma::vec y){
//   for(uword j = fq; j<y.n_elem; j++){
//     y(j) = y(j) + y(j-fq);
//   }
//   return(y);
// }

//Quick regression omitting missing values
// [[Rcpp::export]]
arma::mat QuickReg(arma::mat X,
                   arma::mat Y){
  uvec ind_rm, ind;
  vec y, x;
  mat B(X.n_cols, Y.n_cols, fill::zeros), xx;
  uvec indX = find_nonfinite(X.col(0));
  for(uword k = 1; k<X.n_cols; k++){
    x = X.col(k);
    indX = unique( join_cols(indX, find_nonfinite(x))); //index of missing X values
  }
  for(uword j=0; j<Y.n_cols; j++){
    y       = Y.col(j);
    ind     = unique(join_cols(indX,find_nonfinite(y))); //index of elements to remove
    xx      = X;
    //this seems a tedious way to shed non-contiguous indexes
    for(uword k = ind.n_elem; k>0; k--){
      xx.shed_row(ind(k-1));
      y.shed_row(ind(k-1));
    }
    B.col(j) = solve(trans(xx)*xx, trans(xx)*y);
  }
  return(B);
}

//Quick univariate regression omitting missing values
// [[Rcpp::export]]
List UVreg(arma::vec x,
           arma::vec y,
           arma::uword rm_outlier = 0){
  
  double xi, B, sig2, sd;
  vec u;
  uvec ind = find_nonfinite(x);
  ind = unique( join_cols(ind, find_nonfinite(y)));
  for(uword k = ind.n_elem; k>0; k--){
    x.shed_row(ind(k-1));
    y.shed_row(ind(k-1));
  }
  xi = 1/as_scalar(trans(x)*x);
  B = xi*as_scalar(trans(x)*y);
  u = y-B*x;
  sig2 = as_scalar(trans(u)*u)/(y.n_elem - 1);
  for(uword j = 0; j < rm_outlier; j++){
    ind = find(abs(u)>5*sqrt(sig2));
    for(uword k = ind.n_elem; k>0; k--){
      x.shed_row(ind(k-1));
      y.shed_row(ind(k-1));
    }
    xi = 1/as_scalar(trans(x)*x);
    B = xi*as_scalar(trans(x)*y);
    u = y-B*x;
    sig2 = as_scalar(trans(u)*u)/(y.n_elem - 1);
  }
  sd = sqrt((sig2)*xi);
  
  List Out;
  Out["B"]  = B;
  Out["sd"] = sd;
  Out["u"] = u;
  return(Out);
}

// //Quick univariate regression omitting missing values
// // [[Rcpp::export]]
// List UVreg(arma::vec x,
//            arma::vec y,
//            bool itc = true){
//   uvec ind = find_nonfinite(x);
//   ind = unique( join_cols(ind, find_nonfinite(y)));
//   for(uword k = ind.n_elem; k>0; k--){
//     x.shed_row(ind(k-1));
//     y.shed_row(ind(k-1));
//   }
//   mat X;
//   if(itc){
//     X = join_rows(ones<colvec>(x.n_elem),x);
//   }else{
//     X = x;
//   }
//   mat Xi = inv_sympd(trans(X)*X);
//   mat B = Xi*(trans(X)*y);
//   vec u = y-X*B;
//   mat V = (as_scalar(trans(u)*u)/(y.n_elem - B.n_elem))*Xi;
//   vec sd = sqrt(V.diag());
//   List Out;
//   Out["B"]  = B;
//   Out["sd"] = sd;
//   Out["u2"] = square(u);
//   return(Out);
// }

// // [[Rcpp::export]]
// arma::uword MonthDays(double year,
//                       double month){
//   double days;
//   if((month == 1) || (month == 3) || (month == 5) || (month == 7) || (month == 8) || (month == 10) || (month == 12)){
//     days = 31;
//   }
//   else if((month == 4) || (month == 6) || (month == 9) || (month == 11) ){
//     days = 30;
//   }
//   else if(round((year-1940)/4) == ((year-1940)/4) ){
//     days = 29;
//   }
//   else{
//     days = 28;
//   }
//   return(days);
// }

// //return last day for the given month
// // [[Rcpp::export]]
// arma::mat end_of_month(arma::mat Dates){
//   for(uword t=0; t<Dates.n_rows; t++){
//     Dates(t,2) = MonthDays(Dates(t,0), Dates(t,1));
//   }
//   return(Dates);
// }

arma::sp_mat MakeSparse(arma::mat A){
  uword n_rows   = A.n_rows;
  uword n_cols   = A.n_cols;
  uvec ind       = find(A);
  umat locations = ind2sub(size(A),ind);
  vec  values    = A(ind);
  sp_mat C(locations,values,n_rows,n_cols);
  return(C);
}

arma::sp_mat sp_rows(arma::sp_mat A,
                     arma::uvec r   ){
  uword n_rows   = A.n_rows;
  //  uword n_cols   = A.n_cols;
  uword n_r      = r.size();
  uvec  tmp      = regspace<uvec>(0,n_rows-1);
  tmp      = tmp.elem(r);
  umat  location = join_vert(trans(regspace<uvec>(0,n_r-1)),trans(tmp));
  sp_mat J(location,ones<vec>(n_r),n_r,n_rows);
  sp_mat C       = J*A;
  return(C);
}

arma::sp_mat sp_cols(arma::sp_mat A,
                     arma::uvec r   ){
  //  uword n_rows   = A.n_rows;
  uword n_cols   = A.n_cols;
  uword n_r      = r.size();
  uvec  tmp      = regspace<uvec>(0,n_cols-1);
  tmp            = tmp.elem(r);
  umat  location = join_vert(trans(tmp),trans(regspace<uvec>(0,n_r-1)));
  sp_mat J(location,ones<vec>(n_r),n_cols,n_r);
  sp_mat C       = A*J;
  return(C);
}


//Replace row r of sparse matrix A with the (sparse) vector a.
//Should be reasonably fast with Armadillo 8 or newer
arma::sp_mat sprow(arma::sp_mat A,
                   arma::mat a,
                   arma::uword r   ){
  //This intitally used find(a) to inentify non-zero elements of a, but that
  //did not replace elements that are non-zero in A and zero in a
  uword n_cols     = A.n_cols;
  if(n_cols>a.n_elem){
    a = join_horiz(a, zeros<mat>(1,n_cols-a.n_elem));
  }
  for(uword n      = 0; n < n_cols; n++){
    A(r,n)         = a(n);
  }
  return(A);
}

//Create the companion form of the transition matrix B
// [[Rcpp::export]]
arma::mat comp_form(arma::mat B){
  uword r = B.n_rows;
  uword c = B.n_cols;
  mat A   = join_vert(B, join_horiz(eye<mat>(c-r,c-r), zeros<mat>(c-r,r)));
  return(A);
}

//mvrnrm and rinvwish by Francis DiTraglia

// [[Rcpp::export]]
arma::mat mvrnrm(int n, arma::vec mu, arma::mat Sigma){
  /*-------------------------------------------------------
# Generate draws from a multivariate normal distribution
#--------------------------------------------------------
#  n        number of samples
#  mu       mean vector
#  Sigma    covariance matrix
#-------------------------------------------------------*/
  RNGScope scope;
  int p = Sigma.n_cols;
  mat X = reshape(vec(rnorm(p * n)), p, n);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  X = eigvec * diagmat(sqrt(eigval)) * X;
  X.each_col() += mu;
  return(X);
}

/*-------------------------------------------------------
# Generate Draws from an Inverse Wishart Distribution
# via the Bartlett Decomposition
#--------------------------------------------------------
# NOTE: output is identical to riwish from MCMCpack
#       provided the same random seed is used
#--------------------------------------------------------
#   n     number of samples
#   S     scale matrix
#   v     degrees of freedom
#-------------------------------------------------------*/
// [[Rcpp::export]]
arma::cube rinvwish(int n, int v, arma::mat S){
  RNGScope scope;
  int p = S.n_rows;
  mat L = chol(inv_sympd(S), "lower");
  cube sims(p, p, n, fill::zeros);
  for(int j = 0; j < n; j++){
    mat A(p,p, fill::zeros);
    for(int i = 0; i < p; i++){
      int df = v - (i + 1) + 1; //zero-indexing
      A(i,i) = sqrt(R::rchisq(df));
    }
    for(int row = 1; row < p; row++){
      for(int col = 0; col < row; col++){
        A(row, col) = R::rnorm(0,1);
      }
    }
    mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
    sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}

// [[Rcpp::export]]
double invchisq(double nu, double scale){
  /*-------------------------------------------------------
# Generate draws from a scaled inverse chi squared distribution
#--------------------------------------------------------
#  nu       "degrees of freedom"
#  scale    scale parameter
#-------------------------------------------------------*/
   vec    x = randn<vec>(nu)/sqrt(scale);
   double s = 1/sum(square(x));
   return(s);
}

//Stack times series data in VAR format
// [[Rcpp::export]]
arma:: mat stack_obs(arma::mat nn, arma::uword p, arma::uword r = 0){
  uword rr = nn.n_rows;
  uword mn = nn.n_cols;
  if(r == 0){
    r = rr-p+1;
  }
  if(rr-p+1 != r){
    stop("Length of input nn and length of data r do not agree.");
  }
  mat N(r,mn*p, fill::zeros);
  uword indx = 0;
  for(uword j = 1; j<=p; j++){
    N.cols(indx,indx+mn-1) = nn.rows(p-j,rr-j);
    indx = indx+mn;
  }
  return(N);
}

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



// -------------------------------------------------------------------
// ---------------- Uniform Frequency Programs -----------------------
//--------------------------------------------------------------------

// Uniform frequency disturbance smoother. Output is a list.
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

  // specifying difuse initial values
  mat P0, P1, S, C;
  P0  = 100000*eye<mat>(sA,sA); //difuse initial factor variance
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

  Zs.row(0)   = r.row(0)*100000*eye<mat>(sA,sA);

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

//Disturbance smoothing for uniform frequency models --- output is only smoothed factors for simulations
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

  // specifying difuse initial values
  mat P0, P1, S, C;
  P0  = 100000*eye<mat>(sA,sA);
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

  Zs.row(0)   = r.row(0)*100000*eye<mat>(sA,sA);

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
  sp_mat Hn, Mn;
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
                  bool store_Y = false, //Store distribution of Y?
                  arma::uword store_idx = 0, // index to store distribution of predicted values
                  arma::uword reps = 1000, //repetitions
                  arma::uword burn = 500,
                  bool verbose = false){ //burn in periods

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
  mat  Ystore;
  vec  Y_median;
  if(store_Y){
    Ystore = zeros<mat>(Y.n_rows, reps);
    Y_median = zeros<vec>(Y.n_rows);
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
      Jh        = J_MF(freq(store_idx), m, LD(store_idx), sA);
      //Ystore(t,rep) =  as_scalar(Zsim.row(t)*trans(Jh)*trans(H.row(0)));
      Ystore.col(rep) = Zsim*trans(Jh)*trans(H.row(store_idx));
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
    for(uword rw=0;rw<Y.n_rows;rw++){
      Y_median(rw) = median(Ystore.row(rw));
    }
  }else{
    Ystore = zeros<mat>(0,0);
    Y_median = zeros<vec>(0);
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
  Out["Y_median"] = Y_median;

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

// [[Rcpp::export]]
List KSeas(arma::mat B,
               double q,
               arma::mat M, //seas. adjustment loadings (the solution)
               double r,
               arma::mat Y, //Data must be entered as uniform frequency --- deal with this in R code
               arma::mat N){ //seasonal adjusment factors



  uword T  = Y.n_rows;
  uword p  = B.n_cols;
  uword m  = N.n_cols;

  // Helper matrix J

  sp_mat   HJ(1,p+1);
  HJ(0,0) = 1;
  B       = join_horiz(B,zeros(1,1)); //Adding an extra lag for WE estimation
  sp_mat   A(comp_form(B));
  sp_mat   Q(p+1,p+1);
  Q(0,0)  = q;
  mat R(1,1);
  R(0,0)  = r;

  // Remove seasonal factors to smooth (solution is identical to smoothing with exogenous factors)
  mat Ydm   = Y - N*trans(M);
  List Smth = Ksmoother(A, Q, HJ, R, Ydm);

  mat Z     = Smth["Z"];
  cube Ps   = Smth["Ps"];
  mat Lik   = Smth["Lik"];

  mat xx, Zx, axz, XZ, azz, ZZ, axx, tmp;

  //For B and qB

  xx  = Z.col(0);
  Zx  = Z.cols(1,p);
  axz = sum(Ps(span(0,0),span(1,p),span::all),2);
  XZ  = trans(Zx)*xx + trans(axz);
  azz = sum(Ps(span(1,p),span(1,p),span::all),2);
  ZZ  = trans(Zx)*Zx + azz;
  B = trans(solve(ZZ,XZ));

  axx = sum(Ps(span(0,0),span(0,0),span::all),2);
  q = as_scalar((trans(xx-Zx*trans(B))*(xx-Zx*trans(B)) + axx - axz*trans(B) - B*trans(axz) + B*azz*trans(B) )/T);


  // For H, M and R

  //declairing variables
  uword n_elm;
  uvec ind;
  vec yy, y, beta, axv;
  double h;
  mat x, XX, Yr;


  xx    =  join_horiz(Z.col(0), N);
  ind   = find_finite(Y);
  y     = Y(ind);
  x     = xx.rows(ind);
  n_elm = y.size();
  axx   = zeros<mat>(m+1,m+1);
  axv   = vectorise(Ps.tube(0,0));
  axx(0,0) = sum(axv(ind));
  XX        = trans(x)*x+axx;
  beta      = solve(XX,trans(x)*y);
  h         = beta(0);
  M         = trans(beta(span(1,m)));
  r         = as_scalar(trans(y-x*beta)*(y-x*beta)  + h*axx(0,0)*h )/n_elm;

  //Normalization

  q         = h*q*h;

  // mat Y_SA  = h*Z.col(0);
  // mat Y_hat = h*Z.col(0) + N*trans(M);

  List Out;

  // Out["A"] = A;
  // Out["Q"] = Q;
  // Out["R"] = R;
  // Out["HJ"] = HJ;
  // Out["Ydm"] = Ydm;
  Out["B"]    = B;
  Out["M"]    = M;
  Out["q"]    = q;
  Out["r"]    = r;
  Out["Lik"]  = Lik;
  // Out["Y_SA"] = Y_SA;
  // Out["Y_hat"]= Y_hat;
  Out["h"]    = h;

  return(Out);
}






