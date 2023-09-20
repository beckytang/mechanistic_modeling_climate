.update_alpha_new <- function(pA, J, TT, X,  alpha, deg_id_ls,
                              a0_prior = matrix(0, nrow = 1, ncol = 10), a0_sd = 0.5){
  alpha_curr <- alpha
  for(i in 1:pA){
    for(j in 1:J){
      alpha_prop <- alpha_curr
      alpha_prop[i,j] <- rnorm(1, alpha_curr[i,j], 0.2)
      lhr <- lh1 <- 0
      for(t in 1:TT){
        is_deg <- 1*(deg_id_ls[[t]][,j] == 1)
        is_nondeg <- 1*(deg_id_ls[[t]][,j] != 1)
        
        XA_prop <- pnorm(X[[t]] %*% alpha_prop[,j])
        XA <- pnorm(X[[t]] %*% alpha_curr[,j])
        
        lh1 <- lh1 + sum(log( ((XA_prop^is_deg)) * ((1-XA_prop)^(is_nondeg))))-
          sum(log( ((XA^is_deg)) * ((1-XA)^(is_nondeg)) ))
      }
      lhr <- lhr + lh1
      if(i == 1){
        lhr <- lhr + dnorm(alpha_prop[i,j], a0_prior[i,j],0.5, log = T)-
          dnorm(alpha_curr[i,j], a0_prior[i,j], a0_sd, log = T)
      } else{
        lhr <- lhr + dnorm(alpha_prop[i,j],a0_sd, log = T)-
          dnorm(alpha_curr[i,j], 0, 5, log = T)
      }
      
      if(lhr > log(runif(1))){
        alpha_curr <- alpha_prop
      }
    }
  }
  return(alpha_curr)
}

.update_beta_new <- function(pB, J, TT, beta, W,  deg_id_ls, beta_prop_sd,
                             b0_prior, b0_sd){
  beta_curr <- beta
  acc_mat <- beta*0
  for(i in 1:pB){
    for(j in 1:J){
      beta_prop <- beta_curr
      beta_prop[i,j] <- rnorm(1, beta_curr[i,j], beta_prop_sd[i,j])
      lhr <- lh1 <- 0
      for(t in 1:TT){
        ng_t <- nrow(W[[t]])
        is_zero <- 1*(deg_id_ls[[t]][,j] == 2)
        is_pos <- 1*(deg_id_ls[[t]][,j] == 0)
        q_prop <- q_curr <- rep(NA, ng_t)
        q_curr <- pnorm(cbind(1, W[[t]][,j]) %*% beta_curr[,j])
        q_prop <- pnorm(cbind(1, W[[t]][,j]) %*% beta_prop[,j])
        
        lh1 <- lh1 + sum(log( (q_prop^(is_zero)) * ((1-q_prop)^(is_pos)) ))-
          sum(log( (q_curr^(is_zero)) * ((1-q_curr)^(is_pos)) ))
      }
      lhr <- lh1 + dnorm(beta_prop[i,j],b0_prior[i,j],b0_sd[i,j], log = T)-
        dnorm(beta_curr[i,j], b0_prior[i,j], b0_sd[i,j], log = T)
      
      
      if(lhr > log(runif(1))){
        beta_curr <- beta_prop
        acc_mat[i,j] <- 1
      }
    }
  }
  ret <- list(beta = beta_curr, acc = acc_mat)
  return(ret)
}


.update_A <- function(TT, W, W_star_pre_cens, U, V, S, A, P, zeroA_vec, nzA_vec, A_lb_vec, A_ub_vec, L, J){
  Avec <- vec(A)
  C <- 0
  c <- 0
  for (t in 1:(TT-1)){
    Ct_inv <-  solve(S) %x%  (U[[t]] %*% t(U[[t]]))
    C <- C + Ct_inv
    c <- c + Ct_inv %*% vec( solve((U[[t]] %*% t(U[[t]])) +solve(0.000001*diag(L)) ) %*% U[[t]] %*%
                               (W_star_pre_cens[[t]] - (W[[t]] + t(V[[t]]) %*% P)))
  }
  
  Vv <- solve(C)
  m <- Vv %*% c
  
  m_cond <- m[-zeroA_vec] -  Vv[-zeroA_vec, zeroA_vec] %*% solve(Vv[zeroA_vec, zeroA_vec]) %*% m[zeroA_vec]
  Vv_cond <- Vv[-zeroA_vec, -zeroA_vec] - Vv[-zeroA_vec, zeroA_vec] %*% solve(Vv[zeroA_vec, zeroA_vec]) %*% Vv[zeroA_vec, -zeroA_vec]
  
  #A_vec <- t(rmvn(1, m, Vv))
  Aprop <- t(rmvnormRcpp(1, t(m_cond), round(Vv_cond, 8)))
  ids <- which(Aprop >= A_lb_vec & Aprop <= A_ub_vec)
  Avec[nzA_vec[ids]] <- Aprop[ids]
  A_ret <- matrix(Avec, nrow = L, ncol = J, byrow = F)
  return(A_ret)
}

.update_P <- function(TT, W, W_star_pre_cens, U, V, S, P, A, zeroP_vec, nzP_vec, P_lb_vec, P_ub_vec, R, J,
                      P_prior_prec){
  Pvec <- vec(P)
  C <- 0
  c <- 0
  for (t in 1:(TT-1)){
    if(t  >= 1){
      Ct_inv <-  solve(S) %x%  (V[[t]] %*% t(V[[t]]))# + solve(0.0001*diag(R))
      C <- C + Ct_inv
      # c <- c + Ct_inv %*% vec( solve((V[[t]] %*% t(V[[t]])) + 0.00001*diag(R)) %*% V[[t]] %*%
      #                            (W_star[[t]] - (W[[t]] + t(U[[t]]) %*% A)))
      c <- c + Ct_inv %*% vec( solve((V[[t]] %*% t(V[[t]])  +solve(0.000001*diag(R)))) %*% V[[t]] %*%
                                 (W_star_pre_cens[[t]] - (W[[t]] + t(U[[t]]) %*% A)))
    }
  }
  
  
  Vv <- solve(C +  P_prior_prec)
  m <- Vv %*% c
  
  m_cond <- m[-zeroP_vec] -  Vv[-zeroP_vec, zeroP_vec] %*% solve(Vv[zeroP_vec, zeroP_vec]) %*% m[zeroP_vec]
  Vv_cond <- Vv[-zeroP_vec, -zeroP_vec] - Vv[-zeroP_vec, zeroP_vec] %*% solve(Vv[zeroP_vec, zeroP_vec]) %*% Vv[zeroP_vec, -zeroP_vec]
  
  Pprop <- t(rmvnormRcpp(1, t(m_cond), round(Vv_cond, 8)))
  ids <- which(Pprop >= P_lb_vec & Pprop <= P_ub_vec)
  Pvec[nzP_vec[ids]] <- Pprop[ids]
  P_ret <- matrix(Pvec, nrow = R, ncol = J, byrow = F)
  return(P_ret)
}
expit <- function(x){1/(1+exp(-x))}



fisheries_mod <- function(Y, n, J, TT, XV, X, Z, sst_mat, R, pA, pB, 
                              A_lb, A_ub, A_lb_vec, A_ub_vec, nzA_vec, zeroA_vec,
                              P_lb_vec, P_ub_vec, nzP_vec, zeroP_vec, s2_eta_init, s2_gamma_init,
                              t_opt_init, lambda_init, rho_init, delta_init, prop_sd_mat,
                              t_opt_lb_vec, t_opt_ub_vec, Dik_vec, d_vec, 
                              nbs_ls, A_mat, spec_vec, seq_length, w_mod_amount, rep_id, update_H,
                              G, burnin, thin_amt, s2_gamma_a = 1, s2_gamma_b = 1, s2_tau_a = 1, s2_tau_b = 1,
                              seed = 1){
  # initialize
  set.seed(seed)
  if(T){
    n_pos <-  apply(do.call(rbind, Y), 2, function(x){sum(x > 0)})
    n_pos_m1 <-  apply(do.call(rbind, Y[-1]), 2, function(x){sum(x > 0)})
    z_ids_ls <- lapply(Y, function(x){which(x == 0)})
    U_ids <- (1:J^2)[-rep_id]
    
    # s2_eta
    s2_eta <- s2_eta_init
    s_eta_vec <- sqrt(rep(s2_eta, each  = n))
    # s2_gamma
    s2 <- s2_gamma_init
    s_vec <- rep(sqrt(s2), each = n)
    S <- diag(s2)
    # A
    A <- matrix(0, nrow = L, ncol = J); 
    A[which(A_lb < 0)] <- -0.05
    A[which(A_ub > 0)] <- 0.05
    A[which(A_lb > -0.01 & A_lb < 0)] <- -0.01/4
    A[which(A_ub < 0.01 & A_ub > 0)] <- 0.01/4
    
    Astart <- A
    Avec <- matrixcalc::vec(A)
    
    # P
    P <- matrix(0, nrow = R, ncol = J); 
    q <- nrow(XV[[1]])
    P[cbind(seq(1, q*J, q), 1:J)] <- c(-0.05, 0.05, 0.05, -0.03, -0.04, 0.03, 0.03, -0.05)
    
    Pstart <- P
    Pvec <- matrixcalc::vec(P)
    P_prior_prec = solve(1 * diag(R*J))
    
    # alpha
    alpha <- matrix(0, nrow = pA, ncol = J)
    
    # beta
    beta <- matrix(0, nrow = pB, ncol = J)
    beta_acc <- matrix(0, nrow = pB, ncol = J)
    
    lambda_vec <- lambda_init
    t_opt_vec <- t_opt_init
    rho_vec <- rho_init
    delta_vec <- delta_init
    Dik_ls <- Dik_vec[d_vec]
    
    
    ###### bounds for w_stars and w_tildes
    W_star <- W_star_pre_cens <- list()
    wstar_ub_ls <- wstar_lb_ls <- list() 
    for(t in 1:(TT-1)){
      ng_t <- length(obs_ids_ls[[t]])
      ub_mat <- lb_mat <- matrix(NA, nrow = ng_t, ncol = J)
      
      match_ids_t <- which(obs_ids_ls[[t]] %in% obs_ids_ls[[t+1]])
      match_ids_tp1 <- which(obs_ids_ls[[t+1]] %in% obs_ids_ls[[t]])
      
      for(i in 1:ng_t){
        if(obs_ids_ls[[t]][i] %in% obs_ids_ls[[t+1]]){
          # curr_id <- 
          match_tp1 <- which(obs_ids_ls[[t+1]] == obs_ids_ls[[t]][i])
          ub_mat[i,] <- Y[[t+1]][match_tp1 ,]  * 1.3#
          # spiny dog fish are too big!
          ub_mat[i,J] <- Y[[t+1]][match_tp1 , J]  * 1.1#
          lb_mat[i,] <- Y[[t+1]][match_tp1 ,] / 1.5
          
          small_tp1 <- which(Y[[t+1]][match_tp1,] < 2)
          if(length(small_tp1) > 0){
            lb_mat[i, small_tp1] <- -5
          }
          zero_tp1 <- which(Y[[t+1]][match_tp1,] == 0)
          if(length(zero_tp1) > 0){ 
            ub_mat[i, zero_tp1] <- 20 ### CHANGE THIS
          }
        }
      }
      ub_mat[is.na(ub_mat)] <- 5 ### CHANGE THIS
      lb_mat[is.na(lb_mat)] <- -5
      
      wstar_ub_ls[[t]] <- round(ub_mat,2)
      wstar_lb_ls[[t]] <- round(lb_mat,2)
    }
    
    w_ub <- ceiling(max(unlist(Y))) + 2
    w_lb <- 0
    W_nondeg <- list()
    # degeneracy
    deg_id_ls <- list()
    for(t in 1:(TT)){
      # randomly choose which are degenerate
      ng_t <- length(obs_ids_ls[[t]])
      draw_deg_mat <- matrix(0,nrow = TT, ncol = ng_t*J)
      deg_temp <- matrix(0, nrow = ng_t, ncol = J)
      z_ids <- z_ids_ls[[t]]
      
      # begin by assuming most are chance 0s
      if(length(z_ids) > 0){
        deg_temp[z_ids] <- sample(c(1,2), size = length(z_ids), 
                                  replace = T)
      }
      deg_id_ls[[t]] <- deg_temp
    }
    
    # init H_ls
    # each element of H_ls is species-specific list
    H_ls <- list()
    for(j in 1:J){
      HH_ls <- list()
      for(t in 1:(TT-1)){
        xx2 <- abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = F) - t_opt_vec[j]) - 
          abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = T) - t_opt_vec[j])
        diag(xx2) <- abs(sst_mat[,t+1] - t_opt_vec[j]) -  abs(sst_mat[,t] - t_opt_vec[j])
        H_temp <- exp(-lambda_vec[j]*xx2 - delta_vec[j]*Dik_ls[[j]]- rho_vec[j]*A_mat)
        H_temp <- as.matrix(H_temp*nbs_ls[[j]])
        H_temp[(deg_id_ls[[t+1]][,j] == 1),] <- 0
        
        for(nn in 1:n){
          H_temp[,nn] <-  H_temp[,nn] / sum( H_temp[,nn])
        }
        HH_ls[[t]] <- H_temp
      }
      H_ls[[j]] <- HH_ls 
    }
    names(H_ls) <- spec_vec
    
    ##### initialize w and w_star
    W <- U <- V <- list()
    W_star_deg <- W_star <- W <- W_nondeg <-  list()
    W_star_nondeg <- list()
    for(t in 1:(TT-1)){
      if(t==1){
        W[[t]] <- W_nondeg[[t]] <- matrix(0, nrow = n, ncol = J)
        W1_init <- Y[[t]]
        W1_init[which(W1_init  == 0)] <-1
        W[[t]] <-W_nondeg[[t]] <-W1_init
      } else{
        Wstar_cens <- W_star[[t-1]]
        Wstar_cens[which(Wstar_cens <0)] <- 0
        # W[[t]] <- W_nondeg[[t]] <- matrix(bdiag(lapply(H_ls, function(x){x[[t-1]]})) %*% vec(Wstar_cens),
        #                                   ncol = J, nrow = n)
        temp_post_redist <- matrix(bdiag(lapply(H_ls, function(x){x[[t-1]]})) %*% vec(Wstar_cens),
                       ncol = J, nrow = n)
        prob_ids <- which(Y[[t]] > 0 & temp_post_redist == 0, arr.ind = T)

        if(length(prob_ids) > 0){
          flag <- T
          while(flag){
            for(p in 1:nrow(prob_ids)){
              s_col <- prob_ids[p,2]
              temp_post_redist[prob_ids[p,1],s_col ] <- 0.1
              
              borrow_row <-  which(temp_post_redist[,s_col] == max(temp_post_redist[,s_col]))[1]
              temp_post_redist[borrow_row, s_col] <- temp_post_redist[borrow_row, s_col] -0.1
            }
            prob_ids <- which(Y[[t]] > 0 & temp_post_redist == 0, arr.ind = T)
            if(length(prob_ids) == 0){ flag <- F}
          }
        }
        W[[t]] <- W_nondeg[[t]] <- temp_post_redist
        
      }
      W[[t]][which(deg_id_ls[[t]] == 1)] <- 0
      U[[t]] <-  apply(W[[t]], 1, function(x){(x %x% x)[-rep_id]})
      V[[t]] <- apply(matrix(1:n, ncol = 1), 1, function(x){W[[t]][x,] %x% XV[[t]][,x]})
      
      mm <- (W[[t]] + t(V[[t]]) %*% P + t(U[[t]]) %*% A) + rnorm(n*J, 0, sd =  s_vec)
      l_ids <- which(mm < wstar_lb_ls[[t]])
      u_ids <- which(mm > wstar_ub_ls[[t]])
      if(length(l_ids) > 0){
        mm[l_ids] <- wstar_lb_ls[[t]][l_ids] + 0.01
      }
      if(length(u_ids) > 0){
        mm[u_ids]<- wstar_ub_ls[[t]][u_ids] - 0.01
      }
      deg_ids <- deg_id_ls[[t]]
      temp <- matrix(0, nrow = ng_t, ncol = J)
      for(i in 1:(ng_t*J)){
        if(deg_ids[i] == 1){
          temp[i] <- 0
        } else{
          temp[i] <- max(0,mm[i])
        }
      }
      W_star[[t]] <- temp
      W_star_pre_cens[[t]] <- mm
      
      W_star_nondeg[[t]] <- mm
      W_star_deg[[t]] <- matrix(0, nrow = ng_t, ncol = J)
    }
    t <- TT
    Wstar_cens <- W_star[[t-1]]
    Wstar_cens[which(Wstar_cens <0)] <- 0
    W[[t]] <- W_nondeg[[t]] <-  matrix(bdiag(lapply(H_ls, function(x){x[[t-1]]})) %*% vec(Wstar_cens),
                                       ncol = J, nrow = n)
    W[[t]][which(deg_id_ls[[t]] == 1)] <- 0
    # done initializing
    
  }
  
  ### STORE
  thin_amt2 <- 10; store <- store2 <- 0
  S2_GAMMA <- matrix(NA,nrow = G/thin_amt,ncol = J )
  S2_ETA <-matrix(NA, nrow = G/thin_amt, ncol = J)
  ALPHA <- matrix(NA, nrow = G/thin_amt, ncol = pA*J)
  BETA <- matrix(NA, nrow = G/thin_amt, ncol = pB*J)
  A_POST <- matrix(NA, nrow = G/thin_amt, ncol = length(A))
  P_POST <- matrix(NA, nrow = G/thin_amt, ncol = length(P))
  W_POST <- W_STAR <-list()
  W_TTm1 <- matrix(NA, nrow = G/thin_amt, ncol = ng_t*J)
  LAMBDA <- T_OPT <- DELTA <- RHO <- matrix(NA,nrow = G/thin_amt,ncol = J )
  W_all_ls <- prop_neg_ls <- list()
  for(t in 1:TT){
    W_all_ls[[t]] <- matrix(0, nrow = ng_t, ncol = J)
    
  }
  for(t in 1:(TT-1)){
    prop_neg_ls[[t]] <- matrix(0, nrow = ng_t, ncol = J)
  }
  
  delta.acc <- lambda.acc <- t_opt.acc <- rho.acc <- rep(0, J)
  delta.acc_prev <- lambda.acc_prev <- t_opt.acc_prev <- rho.acc_prev <- rep(0, J)
  check_sd <- seq(1,G+burnin, seq_length)
  check_id <- 2
  w_acc_ls <- lapply(1:(TT-1), matrix, data = 0, nrow = n, ncol = J)
  redist_diff_ls <- growth_diff_ls <-redist_diff_rel_ls <- growth_diff_rel_ls<-growth_div_ls <- redist_div_ls<- lapply(1:(TT-1), matrix, data = 0, nrow = n, ncol = J)
  
  w_star_prop_sd_ls <-list()
  acc_nondeg <-  list() 
  for(t in 1:(TT-1)){
    ng_t <-  length(obs_ids_ls[[t]])
    acc_nondeg[[t]] <-   rep(0, ng_t * J)
    temp <- matrix(1, ncol = J, nrow = ng_t)
    w_star_prop_sd_ls[[t]] <- temp
  }
  
  mod_amount <- 200; w_count <- 0
  

  
  g<-1
  gg <- 0
  # SAMPLER
  for(g in 1:(G + burnin)){
    if(g %in%  check_sd[-1]){
      print(paste0("Iteration: ", g))
      print(paste0("Acc: ", round((lambda.acc - lambda.acc_prev) / (check_sd[check_id] - check_sd[check_id - 1] ),2),
                   ", ", round( (t_opt.acc- t_opt.acc_prev) / (check_sd[check_id] - check_sd[check_id - 1] ),2),
                   ", ", round((delta.acc - delta.acc_prev) / (check_sd[check_id] - check_sd[check_id - 1] ),2),
                   ", ", round( (rho.acc - rho.acc_prev) / (check_sd[check_id] - check_sd[check_id - 1] ),2)))
      acc_diff <- rbind(lambda.acc - lambda.acc_prev,t_opt.acc - t_opt.acc_prev,  delta.acc - delta.acc_prev, rho.acc - rho.acc_prev)
      g_ids <- which(acc_diff/ (check_sd[check_id] - check_sd[check_id - 1] ) > 0.5)
      if(length(g_ids) > 0){
        prop_sd_mat[g_ids] <- (prop_sd_mat* 1.1)[g_ids]#sqrt(prop_sd_mat^2 * 1.1)[g_ids]
      }
      l_ids <- which(acc_diff/ (check_sd[check_id] - check_sd[check_id - 1] ) < 0.15)

      if(length(l_ids) > 0){
        prop_sd_mat[l_ids] <- (prop_sd_mat * 0.9)[l_ids] #sqrt(prop_sd_mat^2 * 0.9)[l_ids]
      }
      check_id <- check_id + 1
      t_opt.acc_prev <- t_opt.acc; lambda.acc_prev <- lambda.acc; delta.acc_prev <- delta.acc; rho.acc_prev <- rho.acc
    }
    
    if(update_H){ # do it for each H_j
      for(j in sample(1:J)){
        Hj_prop_ls <-list()
        nbs <- nbs_ls[[j]];  Dik <- Dik_ls[[j]]
        t_opt <-  t_opt_vec[j]
        delta <- delta_vec[j]
        rho <- rho_vec[j]
        lambda_curr <- lambda_vec[j]
        lambda_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_curr, sd = prop_sd_mat[1,j])
        for(t in 1:(TT-1)){
          xx2 <- abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = F) - t_opt) - 
            abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = T) - t_opt)
          diag(xx2) <- abs(sst_mat[,t+1] - t_opt) -  abs(sst_mat[,t] - t_opt)
          
          P_ik_star <- exp(-lambda_prop*xx2 - delta*Dik- rho*A_mat)
          P_ik_star <- as.matrix(P_ik_star*nbs)
          P_ik_star[which(deg_id_ls[[t+1]][,j] == 1),] <- 0
          P_ik_star <- P_ik_star / matrix(colSums(P_ik_star), nrow = n, ncol = n, byrow = T)
          Hj_prop_ls[[t]] <- P_ik_star
        }
        lh_prop <-  log(dtruncnorm(lambda_curr , a = 0, b = Inf, mean = lambda_prop, sd = prop_sd_mat[1,j]))-
          log(dtruncnorm(lambda_prop, a = 0, b = Inf, mean = lambda_curr , sd =prop_sd_mat[1,j])) ;
        lh_prior <- dgamma(lambda_prop, 5, 1, log = T)- dgamma(lambda_curr , 5,1, log = T);
        
        ret_temp <- updateH_new(W_tilde = W, W_star = W_star, W_star_pc = W_star_pre_cens,
                                Hj_ls = H_ls[[j]],
                                Hj_prop_ls= Hj_prop_ls, j = j-1, deg_id_ls= deg_id_ls,  z_ids_ls= z_ids_ls,  Y_ls = Y,
                                ng = n, J = J,  s_gamma_vec = matrix(s_vec, ncol = 1),
                                s_eta_vec = matrix(s_eta_vec, ncol = 1), TT = TT, U_ids = U_ids-1, V_ids = (1:R)-1,
                                V_ls = lapply(V, t),   U_ls = lapply(U,t), XV= XV,
                                P = P, A= A,lh_prop = lh_prop, lh_prior = lh_prior)

        if(ret_temp$lh4 != -100000){
          if(ret_temp$lhr > log(runif(1))){
            lambda_vec[j] <- lambda_prop
            H_ls[[j]] <- Hj_prop_ls; W <- ret_temp$W_prop_ls; V <- ret_temp$V_prop_ls; U <- ret_temp$U_prop_ls
            lambda.acc[j] <- lambda.acc[j] + 1
          }
        }
        
        
        Hj_prop_ls <-list()
        delta <- delta_vec[j]
        rho <- rho_vec[j]
        lambda <- lambda_vec[j]
        t_opt_curr <-  t_opt_vec[j]
        t_opt_prop <- rtruncnorm(1, a = t_opt_lb_vec[j], b = t_opt_ub_vec[j], mean = t_opt_curr, sd = prop_sd_mat[2,j])
        for(t in 1:(TT-1)){
          xx2 <- abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = F) - t_opt_prop) - 
            abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = T) - t_opt_prop)
          diag(xx2) <- abs(sst_mat[,t+1] - t_opt_prop) -  abs(sst_mat[,t] - t_opt_prop)
          P_ik_star <- exp(-lambda*xx2 - delta*Dik- rho*A_mat)
          P_ik_star <- as.matrix(P_ik_star*nbs)
          P_ik_star[which(deg_id_ls[[t+1]][,j] == 1),] <- 0
          P_ik_star <- P_ik_star / matrix(colSums(P_ik_star), nrow = n, ncol = n, byrow = T)
          Hj_prop_ls[[t]] <- P_ik_star
        }
        
        lh_prop <-  log(dtruncnorm(t_opt_curr, a = t_opt_lb_vec[j], b = t_opt_ub_vec[j], mean = t_opt_prop, sd = prop_sd_mat[2,j]))-
          log(dtruncnorm(t_opt_prop, a = t_opt_lb_vec[j], b = t_opt_ub_vec[j], mean = t_opt_curr, sd =prop_sd_mat[2,j]))
        lh_prior <- dnorm(t_opt_prop, 12, 1, T)- dnorm(t_opt_curr, 12, 1,T)
        ret_temp <- updateH_new(W_tilde = W, W_star = W_star, W_star_pc = W_star_pre_cens,
                                Hj_ls = H_ls[[j]],
                                Hj_prop_ls= Hj_prop_ls, j = j-1, deg_id_ls= deg_id_ls,  z_ids_ls= z_ids_ls,  Y_ls = Y,
                                ng = n, J = J,  s_gamma_vec = matrix(s_vec, ncol = 1),
                                s_eta_vec = matrix(s_eta_vec, ncol = 1), TT = TT , U_ids = U_ids-1, V_ids = (1:R)-1,
                                V_ls = lapply(V, t),   U_ls = lapply(U,t), XV= XV,
                                P = P, A= A,lh_prop = lh_prop, lh_prior = lh_prior)
        if(ret_temp$lh4 != -100000){
          if(ret_temp$lhr > log(runif(1))){
            t_opt_vec[j] <- t_opt_prop
            H_ls[[j]] <- Hj_prop_ls; W <- ret_temp$W_prop_ls; V <- ret_temp$V_prop_ls; U <- ret_temp$U_prop_ls
            t_opt.acc[j] <- t_opt.acc[j] + 1
          } 
        }
        
        
        Hj_prop_ls <-list()
        rho <- rho_vec[j]
        t_opt <- t_opt_vec[j]
        lambda <- lambda_vec[j]
        delta_curr <- delta_vec[j]
        delta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = delta_curr, sd = prop_sd_mat[3,j])
        for(t in 1:(TT-1)){
          xx2 <- abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = F) - t_opt) - 
            abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = T) - t_opt)
          diag(xx2) <- abs(sst_mat[,t+1] - t_opt) -  abs(sst_mat[,t] - t_opt)
          P_ik_star <- exp(-lambda*xx2 - delta_prop*Dik- rho*A_mat)
          P_ik_star <- as.matrix(P_ik_star*nbs)
          P_ik_star[which(deg_id_ls[[t+1]][,j] == 1),] <- 0
          P_ik_star <- P_ik_star / matrix(colSums(P_ik_star), nrow = n, ncol = n, byrow = T)
          Hj_prop_ls[[t]] <- P_ik_star
        }
        
        lh_prop <- log(dtruncnorm(delta_curr, a = 0, b = Inf, mean =delta_prop, sd = prop_sd_mat[3,j]))-
          log(dtruncnorm(delta_prop, a = 0, b = Inf, mean = delta_curr, sd = prop_sd_mat[3,j]))
        lh_prior <- dgamma(delta_prop, 1,1, log = T)-dgamma(delta_curr, 1,1, log = T)
        ret_temp <- updateH_new(W_tilde = W, W_star = W_star, W_star_pc = W_star_pre_cens,
                                Hj_ls = H_ls[[j]],
                                Hj_prop_ls= Hj_prop_ls, j = j-1, deg_id_ls= deg_id_ls,  z_ids_ls= z_ids_ls,  Y_ls = Y,
                                ng = n, J = J,  s_gamma_vec = matrix(s_vec, ncol = 1),
                                s_eta_vec = matrix(s_eta_vec, ncol = 1), TT = TT , U_ids = U_ids-1, V_ids = (1:R)-1,
                                V_ls = lapply(V, t),   U_ls = lapply(U,t), XV= XV,
                                P = P, A= A,lh_prop = lh_prop, lh_prior = lh_prior)
        if(ret_temp$lh4 != -100000){
          if(ret_temp$lhr > log(runif(1))){
            delta_vec[j] <- delta_prop
            H_ls[[j]] <- Hj_prop_ls; W <- ret_temp$W_prop_ls; V <- ret_temp$V_prop_ls; U <- ret_temp$U_prop_ls
            delta.acc[j] <- delta.acc[j] + 1
          }
        }
        
        Hj_prop_ls <-list()
        delta <- delta_vec[j]
        t_opt <- t_opt_vec[j]
        lambda <- lambda_vec[j]
        rho_curr <- rho_vec[j]
        rho_prop <- rtruncnorm(1, a = 0, b = Inf, mean = rho_curr, sd = prop_sd_mat[4,j])
        for(t in 1:(TT-1)){
          xx2 <- abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = F) - t_opt) - 
            abs(matrix(sst_mat[,t+1], nrow = n, ncol = n, byrow = T) - t_opt)
          diag(xx2) <- abs(sst_mat[,t+1] - t_opt) -  abs(sst_mat[,t] - t_opt)
          
          P_ik_star <- exp(-lambda*xx2 - delta*Dik- rho_prop*A_mat)
          P_ik_star <- as.matrix(P_ik_star*nbs)
          P_ik_star[which(deg_id_ls[[t+1]][,j] == 1),] <- 0
          P_ik_star <- P_ik_star / matrix(colSums(P_ik_star), nrow = n, ncol = n, byrow = T)
          Hj_prop_ls[[t]] <- P_ik_star
        }
        
        lh_prop <- log(dtruncnorm(rho_curr, a = 0, b = Inf, mean = rho_prop, sd = prop_sd_mat[4,j]))-
          log(dtruncnorm(rho_prop, a = 0, b = Inf, mean = rho_curr, sd = prop_sd_mat[4,j]))
        lh_prior <- dgamma(rho_prop, 1,1, log = T)-dgamma(rho_curr, 1,1, log = T)
        ret_temp <- updateH_new(W_tilde = W, W_star = W_star, W_star_pc = W_star_pre_cens,
                                Hj_ls = H_ls[[j]],
                                Hj_prop_ls= Hj_prop_ls, j = j-1, deg_id_ls= deg_id_ls,  z_ids_ls= z_ids_ls,  Y_ls = Y,
                                ng = n, J = J,  s_gamma_vec = matrix(s_vec, ncol = 1),
                                s_eta_vec = matrix(s_eta_vec, ncol = 1), TT = TT , U_ids = U_ids-1, V_ids = (1:R)-1,
                                V_ls = lapply(V, t),   U_ls = lapply(U,t), XV= XV,
                                P = P, A= A,lh_prop = lh_prop, lh_prior = lh_prior)
        if(ret_temp$lh4 != -100000){
          if(ret_temp$lhr > log(runif(1))){
            rho_vec[j] <- rho_prop
            H_ls[[j]] <- Hj_prop_ls; W <- ret_temp$W_prop_ls; V <- ret_temp$V_prop_ls; U <- ret_temp$U_prop_ls
            rho.acc[j] <- rho.acc[j] + 1
          }
        }
        
        
      }  
    }
    
    
    P <- .update_P(TT, W,  W_star_pre_cens, U, V, S, P, A, zeroP_vec, nzP_vec, P_lb_vec, P_ub_vec, R, J,
                   P_prior_prec)
    A <- .update_A(TT, W, W_star_pre_cens, U, V, S, A, P, zeroA_vec, nzA_vec, A_lb_vec, A_ub_vec, L, J)
    
    
    
    
    
    # update s2_gamma's
    for(j in 1:J){
      qa <- length(unlist(obs_ids_ls))/2
      to_add <- 0
      temp <- rep(NA, TT-1)
      for(t in 1:(TT-1)){
        temp[t] <-  crossprod(W_star_pre_cens[[t]][,j] - (W[[t]] + t(V[[t]]) %*% P + t(U[[t]]) %*% A)[,j])
      }
      to_add <- sum(temp)
      ra <- to_add/2
      s[j] <-  1/rgamma(1,qa + s2_gamma_a[j], ra + s2_gamma_b[j]) # this is s2, technically!
    }
    s_vec <- sqrt(rep(s, each = n))
    S <- diag(s)
    
    
    ## update s2_eta
    for(j in 1:J){
      to_add <- 0
      temp <- rep(NA, TT)
      for(t in 1:TT){
        temp[t] <- crossprod((log(Y[[t]][,j])- log(W[[t]][,j]))[which(deg_id_ls[[t]][,j] ==0)])
      }
      to_add <- sum(temp)
      s2_eta[j] <- 1/rgamma(1, s2_tau_a[j] + 0.5*n_pos_m1[j], s2_tau_b[j] + 0.5*to_add)
    }
    s_eta_vec <- sqrt(rep(s2_eta, each = n))
    
    
    ### update alpha
    alpha <-  alpha <- .update_alpha_new(pA, J, TT, X, alpha, deg_id_ls,
                                         a0_prior = matrix(rep(-1,J), nrow = 1))
    
    ### update beta
    beta_ret <-  .update_beta_new(pB, J, TT, beta, W, deg_id_ls,
                                  beta_prop_sd = matrix(c(rep(0.2, J),
                                                          rep(0.25,J)), ncol = J, byrow = T),
                                  # b0_prior = beta_true,
                                  b0_prior = matrix(0, nrow = pB, ncol = J, byrow = T),
                                  b0_sd = matrix(5, nrow = pB,ncol = J, byrow = T))
    beta <- beta_ret$beta
    beta_acc <- beta_acc + beta_ret$acc
    
    
    
    
    ## update w_stars
    if(T){
      for(t in 1:(TT-2)){
        H_big <- bdiag(lapply(H_ls, function(x){x[[t]]}))
        dd <- sample(1:(n*J) - 1)
        max <- J
        xx <- seq_along(dd)
        d1 <- split(dd, ceiling(xx/max))
        ret <- sampW_star_group_new(W_tilde = W, W_star = W_star, W_star_deg = W_star_deg,
                                    W_star_nondeg = W_star_nondeg, z_ids_ls = z_ids_ls,
                                    Y_ls = Y, H_t = H_big, 
                                    ng_t = n, ng_tp1 = n, J = J, w_prop_sd = w_star_prop_sd_ls[[t]],
                                    s_gamma_vec = matrix(s_vec, ncol = 1), s_eta_vec = matrix(s_eta_vec, ncol = 1),
                                    t = t-1, TT = TT-1, U_ids = U_ids-1, V_ids = (1:R)-1,
                                    V_ls = lapply(V, t), U_ls = lapply(U,t), XV = XV,
                                    wstar_lb = wstar_lb_ls[[t]], wstar_ub = wstar_ub_ls[[t]],
                                    P = P, A = A,  b0_vec = rep(beta[1,], each = n),
                                    b1_vec = rep(beta[2,], each = n),
                                    obs_ids_ls = obs_ids_ls, order_ls = d1, n_group = length(d1),
                                    group_len_vec = unlist(lapply(d1, length)),
                                    Wstar_pc = W_star_pre_cens)
        
        w_acc_ls[[t]] <- w_acc_ls[[t]] + ret$acc
        W_star[[t]] <- ret$W_star
        W_star_pre_cens[[t]] <- ret$Wstar_pc_t;
        acc_nondeg[[t]] <- acc_nondeg[[t]] + ret$acc
        deg_temp <- matrix(0, nrow = ng_t, ncol = J)
        deg0_ids <- which(ret$W_tp1 == 0 & Y[[t+1]] == 0)
        chance0_ids <- which(ret$W_tp1 > 0 & Y[[t+1]] == 0)
        if(length(deg0_ids) > 0){
          deg_temp[deg0_ids] <- 1
        }
        if(length(chance0_ids) > 0){
          deg_temp[chance0_ids] <- 2
        }
        deg_id_ls[[t+1]] <- deg_temp
        W[[t+1]] <- ret$W_tp1
        U[[t+1]] <- t(ret$U_tp1)
        V[[t+1]] <- t(ret$V_tp1)
      }
      t <- TT-1
      H_big <- bdiag(lapply(H_ls, function(x){x[[t]]}))
      dd <- sample(1:(n*J) - 1)
      max <- J
      xx <- seq_along(dd)
      d1 <- split(dd, ceiling(xx/max))[1:10]
      ret <- sampW_TTm1_star_group_new(W_tilde = W, W_star = W_star, W_star_deg = W_star_deg,
                                       W_star_nondeg = W_star_nondeg,  z_ids_ls = z_ids_ls,
                                       Y_ls = Y, H_t = H_big, ng_t = n, ng_tp1 = n,
                                       J = J, w_prop_sd = w_star_prop_sd_ls[[t]],
                                       s_gamma_vec = matrix(s_vec, ncol = 1), s_eta_vec = matrix(s_eta_vec, ncol = 1),
                                       t = t-1, TT = TT-1, U_ids = U_ids-1, V_ids = (1:R)-1,
                                       V_ls = lapply(V, t), U_ls = lapply(U,t), XV = XV,
                                       wstar_lb = wstar_lb_ls[[t]], wstar_ub = wstar_ub_ls[[t]],
                                       P = P, A = A, b0_vec = rep(beta[1,], each = n),
                                       b1_vec = rep(beta[2,], each = n),
                                       order_ls = d1, n_group = length(d1),
                                       group_len_vec = unlist(lapply(d1, length) ),
                                       W_star_pc = W_star_pre_cens)
      W_star[[t]] <- ret$W_star;
      W_star_pre_cens[[t]] <- ret$Wstar_pc_t;
      acc_nondeg[[t]]<- acc_nondeg[[t]] + ret$acc
      deg_temp <- matrix(0, nrow = ng_t, ncol = J)
      deg0_ids <- which(ret$W_tp1 == 0 & Y[[t+1]] == 0)
      chance0_ids <- which(ret$W_tp1 > 0 & Y[[t+1]] == 0)
      if(length(deg0_ids) > 0){
        deg_temp[deg0_ids] <- 1
      }
      if(length(chance0_ids) > 0){
        deg_temp[chance0_ids] <- 2
      }
      deg_id_ls[[t+1]] <- deg_temp
      
      W[[t+1]] <- ret$W_tp1;
    
    }
    
    if(g > burnin){
      for(t in 1:(TT-1)){
        prop_neg_ls[[t]] <-  prop_neg_ls[[t]] + 1*(W_star_pre_cens[[t]] < 0)
      }
    }
    
    ### store
    if((g > burnin) & (g %% thin_amt == 0)){
      gg <- gg + 1
      A_POST[gg,]<- A
      P_POST[gg,] <- P
      LAMBDA[gg,] <- lambda_vec
      T_OPT[gg,] <- t_opt_vec
      DELTA[gg,] <- delta_vec
      RHO[gg,] <- rho_vec
      S2_GAMMA[gg,] <- s
      S2_ETA[gg,] <- s2_eta
      ALPHA[gg,] <- alpha
      BETA[gg,] <- beta
      W_TTm1[gg,] <- W_nondeg[[TT]] # fills in by species
      for(t in 1:TT){
        W_all_ls[[t]] <- W_all_ls[[t]]+ W[[t]]
      }
      
      if(g %% w_mod_amount != 0){
        for(t in 1:(TT-1)){
          div_temp <- matrix(1, nrow = n, ncol = J)
          temp <-  (W_star[[t]] - W[[t]])/ W[[t]]
          temp[which(deg_id_ls[[t]] == 1)] <- div_temp[which(deg_id_ls[[t]] == 1)] <- 0
          temp[which(!is.finite(temp))] <-div_temp[which(!is.finite(temp))] <- 0
          growth_diff_rel_ls[[t]] <-  growth_diff_rel_ls[[t]] + temp
          growth_diff_ls[[t]] <-  growth_diff_ls[[t]] + temp*W[[t]]
          growth_div_ls[[t]] <-   growth_div_ls[[t]] + div_temp
          
          
          
          div_temp <- matrix(1, nrow = n, ncol = J)
          temp <-  (W[[t+1]] -W_star[[t]]) / abs( W_star[[t]])
          temp[which(!is.finite(temp))] <-div_temp[which(!is.finite(temp))] <- 0
          redist_diff_rel_ls[[t]] <- redist_diff_rel_ls[[t]] + temp
          redist_diff_ls[[t]] <- redist_diff_ls[[t]] + temp*abs( W_star[[t]])
          redist_div_ls[[t]] <- redist_div_ls[[t]] + div_temp
          
        }
      }
      
    }
  }
  
  W_all_mean <- lapply(W_all_ls,function(x){x/gg})
  redist_diff_ls2 <- growth_diff_ls2 <- redist_diff_rel_ls2 <- growth_diff_rel_ls2 <- list()
  for(t in 1:(TT-1)){
    growth_diff_ls2[[t]] <-  growth_diff_ls[[t]] / ( growth_div_ls[[t]] + 1)
    redist_diff_ls2[[t]] <- redist_diff_ls[[t]]/(redist_div_ls[[t]] + 1)
    growth_diff_rel_ls2[[t]] <-  growth_diff_rel_ls[[t]] / ( growth_div_ls[[t]] + 1)
    redist_diff_rel_ls2[[t]] <- redist_diff_rel_ls[[t]]/(redist_div_ls[[t]] + 1)
  }
  H_acc_tab <- rbind(t_opt.acc, lambda.acc, delta.acc, rho.acc)/g
  rownames(H_acc_tab) <- c("t_opt", "lambda", "delta", "rho")
  colnames(H_acc_tab) <- spec_vec
  w_acc_ls <- do.call(rbind,lapply(w_acc_ls, function(x){x[1,]/g}))
  prop_neg_ls <- lapply(prop_neg_ls, function(x){x/(G)})
  ret_ls <- list(A_POST = A_POST, P_POST = P_POST, 
                 LAMBDA = LAMBDA, T_OPT = T_OPT, DELTA = DELTA,
                 RHO = RHO, S2_GAMMA = S2_GAMMA, S2_ETA = S2_ETA,
                 ALPHA = ALPHA, BETA = BETA, W_TTm1 = W_TTm1,
                 H_acc_tab = H_acc_tab, w_acc_ls = w_acc_ls,
                 redist_diff_ls2 = redist_diff_ls2,
                 growth_diff_ls2 = growth_diff_ls2,
                 redist_diff_rel_ls2 = redist_diff_rel_ls2,
                 growth_diff_rel_ls2 = growth_diff_rel_ls2,
                 W_all_mean = W_all_mean
                 )
  
  return(ret_ls)
}


.distmat <- function(x1,y1,x2,y2){
  xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
  yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
  t(sqrt(xd + yd)) 
}

