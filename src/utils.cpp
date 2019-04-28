// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

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
