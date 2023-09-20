#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
int mod(int x, int y){
  return x - floor(x/y)*y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat vec_cpp(arma::mat W0){
  return W0.as_col();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat vec2mat(arma::vec x, int nrow, int ncol) {
  arma::mat y(x);
  y.reshape(nrow, ncol);
  return y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
bool in_cpp(arma::mat x, double y){
  bool ret = true;
  arma::uvec idx = find(x == y);
  double len = idx.size();
  if(len == 0){
    ret = false;
  }
  return ret;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec tapply_sum_cpp(arma::vec x, arma::vec x_index, int ng2){
  arma::vec ret = zeros(ng2);
  arma::uvec idx;
  for(int i = 0; i < ng2; i++){
    idx = find(x_index == i);
    ret(i) = sum(x.elem(idx));
  }
  return(ret);
}

// [[Rcpp::depends(RcppArmadillo)]]
double tnormRcpp(double lo, double hi, double mu, double sig){
  
  double q1, q2, z;
  
  q1 = Rf_pnorm5(lo,mu,sig,1,0);
  q2 = Rf_pnorm5(hi,mu,sig,1,0);
  z = q1 + unif_rand()*(q2-q1);
  z = Rf_qnorm5(z, mu, sig, 1, 0);
  
  if(z > hi){
    z = lo;
  }
  
  if(z < lo){
    z = hi;
  }
  return(z);
}

// [[Rcpp::export]]
arma::rowvec kron(arma::vec A, arma::vec B, arma::uvec ret_ids){
  int n = A.n_elem;
  int m = B.n_elem;
  int r = ret_ids.n_elem;
  arma::vec ret(n*m);
  arma::rowvec row_ret(r);
  int count=0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      ret[count] = A[i]*B[j];
      count++;
    }
  }
  row_ret = ret.elem(ret_ids).t();
  return row_ret;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rmvnormRcpp(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  bool success = false;
  arma::mat S = sigma;
  int i = 1;
  
  arma::mat Y = randn(n, ncols);
  
  while(success == false && i < 5){
    
    success = chol(S, sigma);
    
    if(success == false){
      sigma += eye(ncols,ncols) * 1e-5;
    }
    
    i = i + 1;
  }
  
  if(success == false){
    //    throw std::range_error("sigma not positive definite");
    return arma::repmat(mu*0, 1, n).t();
  }
  
  return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}



// [[Rcpp::export]]
double dtnorm(double x, double lo, double hi, double mu, double sigma){
  double q1, q2, ret=0;
  bool greater = false, lesser = false, within = false;
  if(x >= lo){
    greater = true;
  }
  if(x <= hi){
    lesser = true;
  }
  if(greater & lesser){
    within = true;
  }
  q1 = Rf_pnorm5(lo,mu,sigma,1,0);
  q2 = Rf_pnorm5(hi,mu,sigma,1,0);
  if (within){
    ret =  R::dnorm(x, mu, sigma, FALSE)/(q2-q1);
  }
  return ret;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List myvec(int n) {
  Rcpp::List x(n);
  Rcpp::IntegerVector choices = {1, 2 ,3};
  for (int i = 0; i < n; ++i) {
    int nc = Rcpp::sample(choices, 1).at(0);
    int nr = Rcpp::sample(choices, 1).at(0);
    Rcpp::NumericVector entries = Rcpp::rbinom(nc * nr, 1, 0.5);
    x(i) = Rcpp::NumericMatrix(nc, nr, entries.begin());
  }
  return x;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List test_list(int n, const List& ls) {
  Rcpp::List x(n);
  Rcpp::IntegerVector choices = {1, 2 ,3};
  for (int i = 0; i < n; ++i) {
    x(i) = ls[i];
  }
  return List::create(_["x"] = x);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_star_group_new(const List& W_tilde, const List& W_star, const List& W_star_deg,
                          const List& W_star_nondeg, 
                          const List& z_ids_ls, const List& Y_ls, sp_mat H_t,
                          int ng_t, int ng_tp1, int J, arma::mat w_prop_sd, arma::mat s_gamma_vec,
                          arma::mat s_eta_vec, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                          const List& V_ls, const List& U_ls, const List& XV,
                          arma::mat wstar_ub, arma::mat wstar_lb,
                          arma::mat P, arma::mat A,   arma::vec b0_vec,
                          arma::vec b1_vec,
                          const List& obs_ids_ls, const List& order_ls,
                          int n_group, IntegerVector group_len_vec, const List& Wstar_pc) {
  arma::vec all_props = zeros(ng_t*J);
  arma::vec LH1 = zeros(ng_t*J), LH2 = zeros(ng_t*J), LH3= zeros(ng_t*J), LH4= zeros(ng_t*J);
  arma::mat Ytp1 = Y_ls[t+1];
  
  // = H_ls[t];
  arma::mat W_tp1= W_tilde[t+1], W_tp1_curr = W_tilde[t+1];
  arma::mat Wstar_pc_t = Wstar_pc[t], Wstar_pc_tp1 = Wstar_pc[t+1];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  arma::mat Utp1 = U_ls[t+1], Utp1_curr;
  arma::mat Vtp1 = V_ls[t+1], Vtp1_curr;
  arma::mat Wstar_tp1 = W_star[t+1], Wstar_tp1_curr;
  
  arma::mat XVtp1 = XV[t+1];
  
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  // arma::vec acc = zeros(ng_t*J);
  NumericVector acc(ng_t*J);
  
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A, mprop_tp1, m_tp1;
  int i;
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  double lh6 = 0;
  bool is_zero;
  arma::uvec neg_ids1, neg_ids2, z_ids_prop;
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop, Uprop_tp1 = Utp1, Vprop_tp1 = Vtp1;
  arma::mat Wstar_cens;
  int group_len;
  arma::vec lh6_vec = zeros(n_group), lh2_vec = zeros(n_group);
  arma::vec ff = zeros(2);
  arma::vec wstar_pc_prop = zeros(ng_t*J);
  
  for(int q = 0; q < n_group; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    lh6 = 0;
    lhr = 0;
    
    W_tp1_curr = W_tp1;
    Utp1_curr = Utp1;
    Vtp1_curr = Vtp1;
    IntegerVector ord = order_ls(q);
    group_len =  group_len_vec(q);
    
    Wstar_cens_prop = Wstar_t;
    for(int o = 0; o < group_len; o++){ // within the group, propose
      i = ord[o];
      wstar_t = Wstar_t[i];
      if(Wt[i] == 0){ // when W_tilde = deg. 0
        // wstar_t = Wstar_deg_t[i];
        // Wstar_t[i] = wstar_t;
        wstar_prop_t = tnormRcpp(-1, 1, wstar_t, 0.075);
        lh4 = lh4+ log(dtnorm(wstar_t, -1, 1, wstar_prop_t, 0.075))-
          log(dtnorm(wstar_prop_t, -1, 1, wstar_t, 0.075));
      } else{
        wstar_prop_t = tnormRcpp(wstar_lb[i], wstar_ub[i], wstar_t, w_prop_sd[i]);
        lh4 = lh4 + log(dtnorm(wstar_t, wstar_lb[i], wstar_ub[i], wstar_prop_t, w_prop_sd[i])) -
          log(dtnorm(wstar_prop_t, wstar_lb[i], wstar_ub[i], wstar_t, w_prop_sd[i]));
      }
       
      // growth
      lh1 = lh1+ R::dnorm(wstar_prop_t, m_t[i], s_gamma_vec[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_gamma_vec[i], TRUE);
      Wstar_cens_prop(i) = wstar_prop_t;
      wstar_pc_prop(i) = wstar_prop_t;
      all_props(i) = wstar_prop_t;
    }
    
    
    // project forward: prop

    neg_ids1 = find(Wstar_cens_prop < 0);
    Wstar_cens_prop.elem(neg_ids1).fill(0);
    Wprop_tp1 = vec2mat(H_t * vec_cpp(Wstar_cens_prop), ng_t, J);

    // project forward: curr
    Wstar_cens = Wstar_t;
    neg_ids2 = find(Wstar_cens < 0);
    Wstar_cens.elem(neg_ids2).fill(0);
    W_tp1 = vec2mat(H_t * vec_cpp(Wstar_cens), ng_t, J);
    //

    // cannot have 0 w with positive y
    for(int l = 0; l < ng_tp1*J; l++){
      if( ((Wprop_tp1[l] == 0) & (Ytp1[l] > 0)) |((W_tp1[l] == 0) & (Ytp1[l] > 0)) ){
        lh6 = - 100000;
      }
    }


    if(lh6 == 0){
      // Y_tp1 ~ LN(W_tile_tp1)
      for(int n = 0; n < ng_t*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1);
        if(is_zero == false){
          lh2 = lh2 + R::dnorm(log(Ytp1[n]), log(Wprop_tp1[n]), s_eta_vec[n], TRUE)-
            R::dnorm(log(Ytp1[n]), log(W_tp1[n]), s_eta_vec[n], TRUE) +
            Rf_pnorm5(b0_vec[n] + b1_vec[n]*Wprop_tp1[n],0, 1, 0, 1)-
            Rf_pnorm5(b0_vec[n] + b1_vec[n]*W_tp1[n],0, 1, 0, 1);
        } else{
          lh2 = lh2 + Rf_pnorm5(b0_vec[n] + b1_vec[n]*Wprop_tp1[n],0, 1, 1, 1)-
            Rf_pnorm5(b0_vec[n] + b1_vec[n]*W_tp1[n],0, 1, 1, 1);
        }
      }

      for(int nn = 0; nn < ng_t; nn++){
        Utp1.row(nn) = kron(W_tp1.row(nn).t(), W_tp1.row(nn).t(), U_ids);
        Vtp1.row(nn) = kron(W_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
        Uprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), Wprop_tp1.row(nn).t(), U_ids);
        Vprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
      }

      mprop_tp1 =  Wprop_tp1 + Vprop_tp1 * P + Uprop_tp1 * A;
      m_tp1 =  W_tp1 + Vtp1* P + Utp1*A;
      for(int l = 0; l < ng_t*J; l++){
        lh3 = lh3 + R::dnorm(Wstar_pc_tp1[l], mprop_tp1[l], s_gamma_vec[l], TRUE)-
          R::dnorm(Wstar_pc_tp1[l], m_tp1[l], s_gamma_vec[l], TRUE);
      }

      lhr = lh1 +lh2 + lh3+ lh4;
      u = runif(1)[0];
      if(lhr > log(u)){ // acce[t]
        for(int o = 0; o < group_len; o++){
          i = ord[o];
          ff[1] = wstar_pc_prop[i];
          Wstar_t(i) = max(ff);
          Wstar_pc_t(i) =  wstar_pc_prop[i];
          acc(i) = 1;
        }

        // Wstar_t = Wstar_cens_prop;
        // Wstar_t(i) = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        Utp1 = Uprop_tp1;
        Vtp1 = Vprop_tp1;
      }else{
        W_tp1 = W_tp1_curr; // previous iteration
        Utp1 = Utp1_curr;
        Vtp1 = Vtp1_curr;
      }
    } else{
      W_tp1 = W_tp1_curr; // previous iteration
      Utp1 = Utp1_curr;
      Vtp1 = Vtp1_curr;
    }
    
    
    
    
    lh6_vec(q) = lh6;
    lh2_vec(q) = lh2;
    // lh2_vec(i) = lhr;
  }
  return List::create(_["W_star"] = Wstar_t, _["W_tp1"] = W_tp1, _["U_tp1"] = Utp1, 
                      _["V_tp1"] = Vtp1, _["acc"] = acc, _["all_props"] = all_props,
                      _["lh6_vec"] = lh6_vec, _["lh2_vec"] = lh2_vec, _["lhr"] = lhr,
                        _["Wstar_pc_t"]=Wstar_pc_t);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_TTm1_star_group_new(const List& W_tilde, const List& W_star, const List& W_star_deg,
                               const List& W_star_nondeg,
                               const List& z_ids_ls, const List& Y_ls, sp_mat H_t,
                               int ng_t, int ng_tp1, int J, arma::mat w_prop_sd, arma::mat s_gamma_vec,
                               arma::mat s_eta_vec, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                               const List& V_ls, const List& U_ls, const List& XV,
                               arma::mat wstar_ub, arma::mat wstar_lb,
                               arma::mat P, arma::mat A, arma::vec b0_vec,
                               arma::vec b1_vec,
                               const List& order_ls,
                               int n_group, IntegerVector group_len_vec,
                               const List& W_star_pc) {
  arma::mat Yt = Y_ls[t];
  arma::mat Ytp1 = Y_ls[t+1];
  arma::mat W_tp1 = W_tilde[t+1], W_tp1_curr;
  arma::mat Wstar_pc_t = W_star_pc[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::vec acc = zeros(ng_t*J);
  // 
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A;
  int i;
  arma::vec  all_props = zeros(ng_t*J);
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  double lh6 = 0;
  bool is_zero;
  int group_len;
  
  arma::vec lh6_vec = zeros(n_group), lh2_vec = zeros(n_group), lhr_vec = zeros(n_group);
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop;
  arma::vec wstar_pc_prop = zeros(ng_t*J);
  arma::mat Wstar_cens;
  arma::uvec neg_ids1, neg_ids2 ;
  arma::vec ff = zeros(2);
  for(int q = 0; q < n_group; q++){
    W_tp1_curr = W_tp1;
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    lh6 = 0;
    
    W_tp1_curr = W_tp1;
    IntegerVector ord = order_ls(q);
    group_len =  group_len_vec(q);
    
    Wstar_cens_prop = Wstar_t;
    for(int o = 0; o < group_len; o++){
      i = ord[o];
      // when W_tilde = deg. 0
      wstar_t = Wstar_t[i];
      if(Wt[i] == 0){// Wtilde = deg 0
        wstar_prop_t = tnormRcpp(-1, 1, wstar_t, 0.075);
        lh4 = lh4 + log(dtnorm(wstar_t, -1, 1, wstar_prop_t, 0.075))-
          log(dtnorm(wstar_prop_t, -1, 1, wstar_t, 0.075));
      } else{
        wstar_prop_t = tnormRcpp(wstar_lb[i], wstar_ub[i], wstar_t, w_prop_sd[i]);
        lh4 = lh4 + log(dtnorm(wstar_t, wstar_lb[i], wstar_ub[i], wstar_prop_t, w_prop_sd[i])) -
          log(dtnorm(wstar_prop_t, wstar_lb[i], wstar_ub[i], wstar_t, w_prop_sd[i]));
      }
      
      // growth
      lh1 = lh1+ R::dnorm(wstar_prop_t, m_t[i], s_gamma_vec[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_gamma_vec[i], TRUE);
      Wstar_cens_prop[i] = wstar_prop_t;
      wstar_pc_prop[i] = wstar_prop_t;
      all_props[i] = wstar_prop_t;
    }
    
    // project forward
    
    neg_ids1 = find(Wstar_cens_prop < 0);
    Wstar_cens_prop.elem(neg_ids1).fill(0);
    Wprop_tp1 = vec2mat(H_t * vec_cpp(Wstar_cens_prop),ng_t,J);
    
    
    Wstar_cens = Wstar_t;
    neg_ids2 = find(Wstar_cens < 0);
    Wstar_cens.elem(neg_ids2).fill(0);
    W_tp1 = vec2mat(H_t * vec_cpp(Wstar_cens),ng_t, J);
    
    
    for(int l = 0; l < ng_tp1*J; l++){
      if( (Ytp1[l] > 0 & Wprop_tp1[l] ==  0) | (Ytp1[l] > 0 & W_tp1[l] == 0)){
        lh6 = - 100000;
      } 
    }
    
    if(lh6 == 0){
      for(int n = 0; n < ng_tp1*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1);
        if(is_zero == false){
          lh2 = lh2 + R::dnorm(log(Ytp1[n]),log( Wprop_tp1[n]), s_eta_vec[n], TRUE)-
            R::dnorm(log(Ytp1[n]), log(W_tp1[n]), s_eta_vec[n], TRUE)+
            Rf_pnorm5(b0_vec[n] + b1_vec[n]*Wprop_tp1[n],0, 1, 0, 1)- 
            Rf_pnorm5(b0_vec[n] + b1_vec[n]*W_tp1[n],0, 1, 0, 1);
        } else{
          lh2 = lh2 + Rf_pnorm5(b0_vec[n] + b1_vec[n]*Wprop_tp1[n],0, 1, 1, 1)- 
            Rf_pnorm5(b0_vec[n] + b1_vec[n]*W_tp1[n],0, 1, 1, 1);
        }
      }
      
      lhr = lh1 +lh2 + lh4 ;//+
      lhr_vec(q) = lhr;
      u = runif(1)[0];
      
      if(lhr > log(u)){
        for(int o = 0; o < group_len; o++){
          i = ord[o];
          ff(1) = wstar_pc_prop[i];
          Wstar_t(i) = max(ff);
          Wstar_pc_t(i) =  wstar_pc_prop[i];
          acc(i) = 1;
        }
        W_tp1 = Wprop_tp1;
      } else{
        W_tp1 = W_tp1_curr;
      }
    } else{
      W_tp1 = W_tp1_curr;
    }
    
    for(int n = 0; n < ng_t*J; n++){
      is_zero = in_cpp(z_tp1_ids,n+1);
      if(is_zero == false){
        lh2 = lh2 + R::dnorm(log(Ytp1[n]),log( Wprop_tp1[n]), s_eta_vec[n], TRUE)-
          R::dnorm(log(Ytp1[n]), log(W_tp1[n]), s_eta_vec[n], TRUE);
      }
    }
    lh6_vec(q) = lh6;
    lh2_vec(q) = lh2;
  }
  
  return List::create(_["W_star"] = Wstar_t, _["W_tp1"] = W_tp1, 
                      _["acc"] = acc, _["lh6"] = lh6_vec, _["lhr_vec"] = lhr_vec,
                        _["Wstar_pc_t"] = Wstar_pc_t);
}






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List updateH_new(const List& W_tilde, const List& W_star, const List& W_star_pc,
                 const List& Hj_ls, const List& Hj_prop_ls,
                 int j, const List& deg_id_ls,  const List& z_ids_ls, const List& Y_ls,
                 int ng, int J,  arma::mat s_gamma_vec,
                 arma::mat s_eta_vec, int TT, arma::uvec U_ids, arma::uvec V_ids,
                 const List& V_ls, const List& U_ls, const List& XV,
                 arma::mat P, arma::mat A, double lh_prop, double lh_prior){
  Rcpp::List W_prop_ls(TT), V_prop_ls(TT-1),U_prop_ls(TT-1) ;
  W_prop_ls(0)=W_tilde[0];
  arma::mat U_tp1_prop = U_ls[0];
  arma::mat V_tp1_prop = V_ls[0];
  V_prop_ls(0) = V_tp1_prop.t();
  U_prop_ls(0) = U_tp1_prop.t();
  bool is_zero, is_deg;
  double lhr = 0, lh2 = 0, lh3 = 0, lh4 = 0;
  bool break_flag = false;
  for(int t=0; t < TT-1 ; t++){
    arma::mat W_star_t = W_star[t];
    arma::mat   Wtilde_tp1_curr = W_tilde[t+1];
    arma::mat Wtilde_tp1_prop =  Wtilde_tp1_curr;
    arma::mat Hj_prop_t = Hj_prop_ls[t];
    Wtilde_tp1_prop.col(j) = Hj_prop_t *  W_star_t.col(j);
    W_prop_ls(t+1) = Wtilde_tp1_prop;
    
    arma::mat  Y_tp1 = Y_ls[t+1];
    arma::vec z_tp1_ids = z_ids_ls[t+1];
    arma::mat deg_tp1 = deg_id_ls[t+1];
    
    for(int nnn = 0; nnn < ng*J;nnn++){
      if( (deg_tp1[nnn] != 1) & (Wtilde_tp1_prop[nnn] == 0)){
        lh4 = -100000;
        break_flag = true;
      }
      if(break_flag) break;
    }
    
    if(break_flag) break;
    
    for(int n = 0; n < ng*J; n++){
      is_zero = in_cpp(z_tp1_ids,n+1);
      if(is_zero == false){
        lh2 = lh2 + R::dnorm(log(Y_tp1[n]),log( Wtilde_tp1_prop[n]), s_eta_vec[n], TRUE)-
          R::dnorm(log(Y_tp1[n]), log(Wtilde_tp1_curr[n]), s_eta_vec[n], TRUE);
      }
    }
  
    //
    
    if(t < (TT-2)){
      arma::mat W_star_tp1 = W_star[t+1];
      arma::mat W_star_pc_tp1 = W_star_pc[t+1];
      arma::mat U_tp1 = U_ls[t+1];
      arma::mat V_tp1 = V_ls[t+1];
      arma::mat U_tp1_prop = U_ls[t+1];
      arma::mat V_tp1_prop = V_ls[t+1];
      arma::mat XV_tp1 = XV[t+1];
      for(int nn = 0; nn < ng; nn++){
        U_tp1_prop.row(nn) = kron( Wtilde_tp1_prop.row(nn).t(),  Wtilde_tp1_prop.row(nn).t(), U_ids);
        V_tp1_prop.row(nn) = kron( Wtilde_tp1_prop.row(nn).t(), XV_tp1.col(nn), V_ids);
      }
      U_prop_ls(t+1) = U_tp1_prop.t();
      V_prop_ls(t+1) = V_tp1_prop.t();
      
      arma::mat m_tp1_prop = Wtilde_tp1_prop + V_tp1_prop*P + U_tp1_prop*A;
      arma::mat m_tp1 =   Wtilde_tp1_curr + V_tp1* P + U_tp1*A;
      for(int l = 0; l < ng*J; l++){
        lh3 = lh3 + R::dnorm(W_star_pc_tp1[l], m_tp1_prop[l], s_gamma_vec[l], TRUE)-
          R::dnorm(W_star_pc_tp1[l], m_tp1[l], s_gamma_vec[l], TRUE);
      }
    }
    
  }
  lhr = lh_prop + lh_prior + lh2 + lh3;
  return List::create( _["lhr"] = lhr, _["W_prop_ls"] = W_prop_ls,
                       _["V_prop_ls"] = V_prop_ls, _["U_prop_ls"] = U_prop_ls, 
                       _["lh2"]  = lh2,  _["lh3"] = lh3, _["lh4"] = lh4);
  
}

