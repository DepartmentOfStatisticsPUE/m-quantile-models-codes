// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec psi_Q(arma::colvec u, double q, double k){
  double sm = median(abs(u)) / 0.6745;
  int n = u.size();
  Rcpp::NumericVector tmp = as<NumericVector>(wrap(k/abs(u/sm)));
  arma::colvec w = Rcpp::pmin(1.0, tmp);
  arma::colvec ww = 2 *(1-q) * w;
  for (int i = 0; i < n; i++){
   if(u[i] > 0)ww[i] = 2 * q * w[i];
  }
  return ww%u;
}

// [[Rcpp::export]]
arma::colvec psi_Q_der(arma::colvec u, double q, double k){
  double sm = median(abs(u)) / 0.6745;
  int n = u.size();
  arma::colvec w = abs(u) / sm;
  for (int i = 0; i < n; i++){
  if( w[i] <= k) w[i] = 1;
  else w[i] = 0;
  }
  arma::colvec ww = 2 * (1-q) * w;
  for (int i = 0; i < n; i++){
   if(u[i] > 0)ww[i] = 2 * q * w[i];
  }
  return ww;
}



// [[Rcpp::export]]
arma::mat singVal(arma::mat m) {

    arma::mat u;
    arma::vec s;
    arma::mat v;

    arma::svd(u,s,v,m);
    
    arma::mat diag_s(m.n_rows, m.n_rows);
    diag_s.fill(0.0); 
    diag_s.diag() = sqrt(s);
    arma::mat sqrtM = u * diag_s * v.t();    
    
    return sqrtM;
}


// // [[Rcpp::export]]
//Rcpp::List makeVMat2(arma::mat Z, double sigma_v, double sigma_e, int m) {
//  double n = Z.n_rows;
//    
//    arma::mat In = arma::eye<arma::mat>(n, n);
//    arma::mat Im = arma::eye<arma::mat>(m, m);
//    arma::mat R = sigma_e * In;
//    arma::mat G = sigma_v * Im;
//    arma::mat V = Z * G * Z.t() + R;
//    arma::mat Vtmp = arma::inv(trimatu(chol(V)));
//    arma::mat Vinv = Vtmp * Vtmp.t();
//    return Rcpp::List::create(Rcpp::Named("V", V),
//    Rcpp::Named("Vinv", Vinv));
//}

// [[Rcpp::export]]
Rcpp::List makeVMat(arma::mat Z, double sigma_v, double sigma_e, int m) {
    double n = Z.n_rows;
    arma::colvec ni = arma::ones(m);
    
    for(int i = 0; i < m; ++i) {
      ni(i) = sum(Z.col(i));
    }
    
    arma::mat V = arma::eye<arma::mat>(n, n);
    arma::mat Vinv = arma::eye<arma::mat>(n, n);
    
    double firsti = 0;
    double firstj = 0;
    double lasti = 0;
    double lastj = 0;
    
    for(int i = 0; i < m; ++i) {
    
    
    arma::colvec Zi = arma::ones(ni(i)) ;
    arma::mat ZZi = Zi * Zi.t();
    
    arma::mat Ini = arma::eye<arma::mat>(ni(i), ni(i));
    arma::mat R = sigma_e * Ini;
    arma::mat Rinv = 1/sigma_e * Ini;
        
    arma::mat test = R+ZZi*sigma_v;
    arma::mat test2 = Rinv-ZZi*(sigma_v/(sigma_v*ni(i)+sigma_e))*(1/sigma_e);
    
    lasti = firsti+ni(i)-1;
    lastj = firstj+ni(i)-1;
    
    V(arma::span(firsti,lasti),arma::span(firstj,lastj)) = test;
    Vinv(arma::span(firsti,lasti),arma::span(firstj,lastj)) = test2;
      
    firsti = lasti+1;
    firstj = lastj+1;
     
    }
    
    return Rcpp::List::create(Rcpp::Named("V", V),
    Rcpp::Named("Vinv", Vinv));
}    


// [[Rcpp::export]]
arma::colvec RandEff(arma::colvec y, arma::mat x, arma::mat Z, double sigma_v, double sigma_e, double m, arma::colvec beta, double tol, int maxit, double k, double k_v, double qtl, double k_val) {
  
  double n = Z.n_rows;
  
  arma::mat In = arma::eye<arma::mat>(n, n);
  arma::mat Im = arma::eye<arma::mat>(m, m);
  arma::mat R = sigma_e * In;
  arma::mat G = sigma_v * Im;
    
  // Variance-Covarianze matrix
  Rcpp::List Vlist = makeVMat(Z, sigma_v, sigma_e, m);
  arma::mat V = Vlist(0);
  arma::mat Vinv = Vlist(1);

  // sqrt of U + inverse
  arma::mat sqrtU(V.n_rows, V.n_rows);
  sqrtU.fill(0.0);
  sqrtU.diag() = sqrt(V.diag());
  
  arma::mat sqrtUinv(V.n_rows, V.n_rows);
  sqrtUinv.fill(0.0);
  sqrtUinv.diag() = 1/sqrtU.diag();
  
  // sqrt of R
  arma::mat sqrtR(n, n); 
  sqrtR.fill(0.0); 
  sqrtR.diag() = sqrt(R.diag());
  
  // sqrt of G
  arma::mat sqrtG(m, m); 
  sqrtG.fill(0.0); 
  sqrtG.diag() = sqrt(G.diag());
  
  // sqrt of R inverse
  arma::mat sqrtRinv(n, n); 
  sqrtRinv.fill(0.0); 
  sqrtRinv.diag() = 1/sqrt(R.diag());
  
  // sqrt of G inverse
  arma::mat sqrtGinv(m, m); 
  sqrtGinv.fill(0.0); 
  sqrtGinv.diag() = 1/sqrt(G.diag());
  
  arma::colvec xbeta = (y - x * beta);
  arma::colvec resid = sqrtUinv * xbeta;
  arma::colvec vv = G * Z.t() * Vinv * xbeta;
  
  arma::mat Part1 = Z.t() * sqrtRinv;
  arma::mat E_psi_der1 = k_val * In;
  arma::mat E_psi_der2 = k_val * Im;
  
  // Algorithm
  
  double diff = 1000;
  int i = 0;
  
  while (diff > tol)
  {
    i++;
    arma::colvec v_robust = vv;
    arma::colvec Zv_robust = Z * v_robust;
    arma::colvec resid1 = sqrtRinv * (xbeta - Zv_robust);
    arma::colvec resid1_rob = psi_Q(resid1, qtl, k_v);
    
    arma::colvec v1 = sqrtGinv * v_robust;
    arma::colvec v1_rob = psi_Q(v1, qtl, k_v);
    
    arma::colvec H = Part1 * resid1_rob - sqrtGinv * v1_rob;
    arma::mat HH = Part1 * E_psi_der1 * sqrtRinv * Z + sqrtGinv * E_psi_der2 * sqrtGinv;
    
    vv = v_robust + inv(HH) * H;
    diff = sum(pow(vv-v_robust, 2));
    
    if (i > maxit) break;
  }
  
  return vv;
}


// [[Rcpp::export]]
Rcpp::List EstFunc_ML2(arma::colvec y, arma::mat x, arma::mat Z, double sigma_v, double sigma_e, double m, arma::colvec beta, double tol, int maxit, double k, double K2_val, double qtl, double k_val) {
  
  double n = Z.n_rows;
  arma::mat ZZ = Z * Z.t();
  arma::mat In = arma::eye<arma::mat>(n, n);
  arma::mat Im = arma::eye<arma::mat>(m, m);
  arma::mat E_psi_der1 = k_val * In;
  arma::mat K2 = K2_val * In;
  arma::colvec sigma(2);
  sigma(0) = sigma_e;
  sigma(1) = sigma_v;
  
  // Algorithm
  
  double diff = 1000;
  int i = 0;
  
  while (diff > tol)
  {
    i++;
    
    arma::colvec sigma1 = sigma;
    
  //  STEP 1 Estimation of beta
    arma::colvec beta1 = beta;
    arma::colvec xbeta = x * beta1;
    arma::mat R = sigma1(0) * In;
    arma::mat G = sigma1(1) * Im;
    
    // Variance-Covarianze matrix
    Rcpp::List Vlist = makeVMat(Z, sigma1(1), sigma1(0), m);
    arma::mat V = Vlist(0);
    arma::mat Vinv = Vlist(1);
    
    // sqrt of U + inverse
    arma::mat sqrtU(V.n_rows, V.n_rows);
    sqrtU.fill(0.0);
    sqrtU.diag() = sqrt(V.diag());
  
    arma::mat sqrtUinv(V.n_rows, V.n_rows);
    sqrtUinv.fill(0.0);
    sqrtUinv.diag() = 1/sqrtU.diag();
    
    arma::colvec resid = sqrtUinv * (y-xbeta);
    arma::colvec resid_rob = psi_Q(resid, qtl, k);
    arma::colvec resid_rob_der = psi_Q_der(resid, qtl, k);
    arma::mat Resid_rob_deri(V.n_rows, V.n_rows);
    Resid_rob_deri.fill(0.0); 
    Resid_rob_deri.diag() = resid_rob_der;
    
    arma::colvec qq = Vinv * sqrtU * resid_rob;
    arma::colvec q1 = x.t() * qq;
    arma::mat KK5 = sqrtU * Vinv;
    arma::mat M1 = x.t() * sqrtUinv * E_psi_der1 * KK5 * x;
    beta = beta1 + inv(M1) * q1;
    
  //  STEP 2 Estimation of variance components
    arma::mat KK0 = resid_rob.t() * KK5;
    arma::colvec a(2);
    a(0) = as_scalar(KK0 * qq);
    a(1) = as_scalar(KK0 * ZZ * qq);
    
    arma::mat KK1 = K2 * Vinv;
    arma::mat KK2 = KK1 * Vinv;
    arma::mat KK3 = KK1 * ZZ * Vinv;
    
    arma::mat A(2,2);
    A(0,0) = trace(KK2);
    A(0,1) = trace(KK2 * ZZ);
    A(1,0) = trace(KK3);
    A(1,1) = trace(KK3 * ZZ);
    
    sigma = inv(A) * a;
    sigma = abs(sigma);
    diff = sum(pow(beta1-beta, 2))+sum(pow(sigma1-sigma, 2));
    if (i > maxit) break;
  }
  
     return Rcpp::List::create(Rcpp::Named("beta", beta),
 Rcpp::Named("sigma", sigma),Rcpp::Named("Iteration", i));
}
