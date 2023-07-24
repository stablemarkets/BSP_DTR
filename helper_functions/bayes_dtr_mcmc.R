### mcmc helper functions
bayes_dtr_mcmc = function(d, iter, burnin, thin, int_size=NULL, 
                          adapt_start, adapt_end){
  K = 4
  
  ### initial jumping dist covariances...to be tuned during adaptation period.
  jump_cov_tau1 = .001*diag(7)
  jump_cov_tau2 = jump_cov_tau3 = .001*diag(11)
  jump_cov_T1  = .001*diag(7)
  jump_cov_T2 = jump_cov_T3 = jump_cov_T4 = .001*diag(11)
  
  #### pre-process data for MCMC ####
  n_iter = length(seq(burnin, iter, thin))
  n_adapt = length(seq(adapt_start, adapt_end,1))
  adapt_ng = adapt_start:adapt_end
  
  N = nrow(d)
  
  #### Create model matrices ####
  ##### Decision point 1 ####
  ## X3, X4, V1 are binary.
  sum_V1 = sum(d$V1)
  sum_X3 = sum(d$X3)
  sum_X4 = sum(d$X4)
  
  x_L1 = x_X1 = x_X2 = model.matrix( as.formula( ~ 1), d )
  x_tau1 = model.matrix( as.formula( ~ -1 + A1 + L1 + V1 + X1 + X2 + X3 + X4 ), d )
  x_T1 =   model.matrix( as.formula( ~ -1 + A1 + L1 + V1 + X1 + X2 + X3 + X4 ), d )
  
  ##### Decision point 2 ####
  d2 = d[d$kappa>=2, ]
  d2$U1gt4 = 1*(d2$U1>4)
  N2 = nrow( d2 )
  x_V2 = model.matrix( as.formula( ~ V1 + A1 ), d2 )
  x_L2 = model.matrix( as.formula( ~ L1 + A1 ), d2 )
  x_tau2 = model.matrix( as.formula(~ -1 + A2 + A1 + L2 + L1 + V2 + V1 + X1 + X2 + X3 + X4 + U1gt4 ), d2 )
  x_T2 = model.matrix( as.formula(~ -1 + A2 + A1 + L2 + L1 + V2 + V1 + X1 + X2 + X3 + X4 + U1gt4), d2 )
  
  ##### Decision point 3 ####
  d3 = d[d$kappa>=3, ]
  d3$U2gt4 = 1*(d3$U2>4)
  N3 = nrow( d3 )
  x_L3 = model.matrix( as.formula( ~ L2 + A2), d3 )
  x_V3 = model.matrix( as.formula( ~ V2 + A2 ), d3 )
  x_tau3 = model.matrix( as.formula(~ -1 + A3 + A2 + L3 + L2 + V3 + V2 + X1 + X2 + X3 + X4 + U2gt4 ), d3 )
  x_T3 = model.matrix( as.formula(~ -1 + A3 + A2 + L3 + L2 + V3 + V2 + X1 + X2 + X3 + X4 + U2gt4), d3 )

  ##### Decision point 4 ####
  d4 = d[d$kappa>=4, ]
  d4$U3gt4 = 1*(d4$U3>4)
  N4 = nrow( d4 )
  x_L4 = model.matrix( as.formula( ~ L3 + A3 ), d4 )
  x_V4 = model.matrix( as.formula( ~ V3 + A3 ), d4 )
  x_T4 = model.matrix( as.formula( ~ -1 + A4 + A3 + L4 + L3 + V4 + V3 + X1 + X2 + X3 + X4 + U3gt4 ), d4 )
  
  ##### Dimensionality of matrices  ####
  ptau1 = ncol(x_tau1); pT1 = ncol(x_T1); 
  
  pL1 = ncol(x_L1)

  ## same dimensionality in decision points >1  
  pVk=ncol(x_V2);
  pLk=ncol(x_L2);
  ptauk=ncol(x_tau2); pTk=ncol(x_T2);
  
  # Create intervals for baseline hazards for waiting times 
  sp1 = split_time(d$U1)
  sp2 = split_time(d2$U2)
  sp3 = split_time(d3$U3)
  sp4 = split_time(d4$U4)
  
  #### Shell for storing MCMC draws ####
  res = list( X1_mu = matrix(NA, nrow=n_iter, ncol=1), 
              X1_psi = matrix(NA, nrow=n_iter, ncol=1),
              X2_mu = matrix(NA, nrow=n_iter, ncol=1), 
              X2_psi = matrix(NA, nrow=n_iter, ncol=1),
              X3_mu = matrix(NA, nrow=n_iter, ncol=1),
              X4_mu = matrix(NA, nrow=n_iter, ncol=1),
              
              L1_mu = matrix(NA, nrow=n_iter, ncol=pL1),
              psi1 = numeric(length = n_iter),
              Lk_mu = array(dim=c(pLk, K-1, n_iter) ),
              psik = array(dim=c(1, K-1, n_iter) ),
              
              V1_mu = matrix(NA, nrow=n_iter, ncol=1),
              Vk_mu = array(dim=c(pVk, K-1, n_iter) ),
              
              hazT1 = matrix(NA, nrow=n_iter, ncol=pT1),
              hazTk = array(dim=c(pTk, K-1, n_iter) ),
              
              haztau1 = matrix(NA, nrow=n_iter, ncol=ptau1),
              haztauk = array(dim=c(ptauk, K-2, n_iter) ), 
              
              lambdaT1 = matrix(NA, nrow=n_iter, ncol = sp1$tau ),
              lambdaT2 = matrix(NA, nrow=n_iter, ncol = sp2$tau ),
              lambdaT3 = matrix(NA, nrow=n_iter, ncol = sp3$tau ),
              lambdaT4 = matrix(NA, nrow=n_iter, ncol = sp4$tau ), 
              
              lambda1 = matrix(NA, nrow=n_iter, ncol = sp1$tau ),
              lambda2 = matrix(NA, nrow=n_iter, ncol = sp2$tau ),
              lambda3 = matrix(NA, nrow=n_iter, ncol = sp3$tau ), 
              
              sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4)
  
  #### shell for storing adaptation
  adapt = list(hazT1 = matrix(NA, nrow=n_adapt, ncol=pT1),
               hazTk = array(dim=c(pTk, K-1, n_adapt) ),
               haztau1 = matrix(NA, nrow=n_adapt, ncol=ptau1),
               haztauk = array(dim=c(ptauk, K-2, n_adapt) ))
  
  #### Initialize Parameters ####
  lambda1 = lambdaT1 = rep(1, sp1$tau)
  lambda2 = lambdaT2 = rep(1, sp2$tau)
  lambda3 = lambdaT3 = rep(1, sp3$tau)
  lambdaT4 = rep(1, sp4$tau)
  
  
  X1_psi = X2_psi = psi1 = psi2 = psi3 = psi4 = 10 
  
  haztau1 = rep(0, ptau1)
  hazT1 = rep(0, pT1)
  
  haztau2 = haztau3 = rep(0, ptauk)
  hazT2 = hazT3 = hazT4 = rep(0, pTk)
  
  V2_mu = V3_mu = V4_mu = rep(0, pVk)

  #### Set Prior Hyperparameters  ####
  prior_mean_X1 = matrix(0, nrow = 1);
  prior_var_X1 = matrix(rep(9, 1), nrow=1 )
  
  prior_mean_X2 = matrix(0, nrow = 1) 
  prior_var_X2 = matrix(rep(9, 1), nrow=1 )
  
  prior_mean_L1 = matrix(0, nrow = pL1)
  prior_var_L1 = matrix(rep(9, pL1), nrow=1 )
  
  prior_mean_Lk = matrix(0, nrow = pLk)
  prior_var_Lk = rep(4, pLk)
  
  prior_mean_Vk = matrix(0, nrow = pVk)
  prior_var_Vk = rep(4, pVk)

  gX1 = bX1 = 3
  gX2 = bX2 = 3
  
  g1 = b1 = 3
  gk = bk = 3
  
  #### Begin MCMC Loop ####
  for(i in 1:iter ){
    
    ###### update Baseline confounder models X1-X4        ######
    X1_mu = rcond_post_beta(d$X1, x_X1, X1_psi, prior_mean_X1, prior_var_X1)
    X1_psi = rcond_post_psi(n = N, beta = X1_mu, xm = x_X1, y = d$X1, gX1, bX1)
    
    X2_mu = rcond_post_beta(d$X2, x_X2, X2_psi, prior_mean_X2, prior_var_X2)
    X2_psi = rcond_post_psi(n = N, beta = X2_mu, xm = x_X2, y = d$X2, gX2, bX2)
    
    X3_mu = update_beta( N, sum_X3 )
    X4_mu = update_beta( N, sum_X4 )

    ##### update time-varying binary variable V parameters ######
    V1_mu = update_beta(N, sum_V1)
    V2_mu = update_logit( V2_mu, x_V2, pVk, d2$V2)
    V3_mu = update_logit( V3_mu, x_V3, pVk, d3$V3)
    V4_mu = update_logit( V4_mu, x_V4, pVk, d4$V4)
    
    ###### update time-vorying continuous variable L parameters ######
    L1_mu = rcond_post_beta(d$L1, x_L1, psi1, prior_mean_L1, prior_var_L1)
    psi1 = rcond_post_psi(n = N, beta = L1_mu, xm = x_L1, y = d$L1, g1, b1)
    
    L2_mu = rcond_post_beta(d2$L2, x_L2, psi2, prior_mean_Lk, prior_var_Lk)
    psi2 = rcond_post_psi(n = N2, beta = L2_mu, xm = x_L2, y = d2$L2, gk, bk)
    
    L3_mu = rcond_post_beta(d3$L3, x_L3, psi3, prior_mean_Lk, prior_var_Lk)
    psi3 = rcond_post_psi(n = N3, beta = L3_mu, xm = x_L3, y = d3$L3, gk, bk)

    L4_mu = rcond_post_beta(d4$L4, x_L4, psi4, prior_mean_Lk, prior_var_Lk)
    psi4 = rcond_post_psi(n = N4, beta = L4_mu, xm = x_L4, y = d4$L4, gk, bk)
    
    ###### update Treatment Hazard parameters  ######
    haztau1_update = update_haz_gp(haztau1, lambda1, x_tau1, ptau1,N,
                                   d$delta_tau1, d$U1, sp1$ts, sp1$te, sp1$tau, 
                                   theta_jump_v = jump_cov_tau1 )
    haztau1 = haztau1_update$theta 
    lambda1 = haztau1_update$lambda
    
    haztau2_update = update_haz_gp(haztau2, lambda2, x_tau2, ptauk, N2,
                                   d2$delta_tau2, d2$U2, sp2$ts, sp2$te,sp2$tau, 
                                   theta_jump_v = jump_cov_tau2)
    haztau2 = haztau2_update$theta 
    lambda2 = haztau2_update$lambda
    
    haztau3_update = update_haz_gp(haztau3, lambda3, x_tau3, ptauk,N3,
                                   d3$delta_tau3, d3$U3, sp3$ts, sp3$te,sp3$tau, 
                                   theta_jump_v = jump_cov_tau3)
    haztau3 = haztau3_update$theta
    lambda3 = haztau3_update$lambda
    
    ###### update Death Hazard parameters  ######
    hazT1_update = update_haz_gp(hazT1, lambdaT1, x_T1, pT1, N,
                                 d$delta_T1, d$U1, sp1$ts, sp1$te, sp1$tau,
                                 theta_jump_v = jump_cov_T1 )
    hazT1 = hazT1_update$theta 
    lambdaT1 = hazT1_update$lambda
    
    hazT2_update = update_haz_gp(hazT2, lambdaT2, x_T2, pTk, N2,
                                 d2$delta_T2, d2$U2, sp2$ts, sp2$te, sp2$tau,
                                 theta_jump_v = jump_cov_T2)
    hazT2 = hazT2_update$theta 
    lambdaT2 = hazT2_update$lambda
    
    hazT3_update = update_haz_gp(hazT3, lambdaT3, x_T3, pTk, N3,
                                 d3$delta_T3, d3$U3, sp3$ts, sp3$te, sp3$tau,
                                 theta_jump_v = jump_cov_T3)
    hazT3 = hazT3_update$theta 
    lambdaT3 = hazT3_update$lambda
    
    hazT4_update = update_haz_gp(hazT4, lambdaT4, x_T4, pTk, N4,
                                 d4$delta_T4, d4$U4, sp4$ts, sp4$te, sp4$tau,
                                 theta_jump_v = jump_cov_T4)
    hazT4 = hazT4_update$theta 
    lambdaT4 = hazT4_update$lambda
    
    #### Store post-burning draws ####
    if( i %% 1000==0 ){
      stage = ifelse(i<burnin,' (burnin)',' (sampling)')
      cat(paste0('iteration ',i,stage,'\n'))
    }
    
    if(adapt_start <= i & i <= adapt_end ){
      idx = i - adapt_start + 1
      adapt$haztau1[idx, ] = haztau1 
      adapt$haztauk[,1,idx] = haztau2
      adapt$haztauk[,2,idx] = haztau3
      adapt$hazT1[idx, ] = hazT1
      adapt$hazTk[, 1, idx] = hazT2
      adapt$hazTk[, 2, idx] = hazT3
      adapt$hazTk[, 3, idx] = hazT4
    }
    
    if(i == adapt_end){
      jump_cov_T1 = cov( adapt$hazT1 )
      jump_cov_T2 = cov( t(adapt$hazTk[, 1, ]) )
      jump_cov_T3 = cov( t(adapt$hazTk[, 2, ]) )
      jump_cov_T4 = cov( t(adapt$hazTk[, 3, ]) )
      
      jump_cov_tau1 = cov( adapt$haztau1 )
      jump_cov_tau2 = cov( t(adapt$haztauk[, 1, ]) )
      jump_cov_tau3 = cov( t(adapt$haztauk[, 2, ]) )
    }
    
    if( (i >= burnin) & (i %% thin==0) ){
      
      ## store post-burnin draws
      idx = (i - burnin)/thin + 1
      
      res$X1_mu[idx, ] = X1_mu 
      res$X1_psi[idx, ] = X1_psi
      res$X2_mu[idx, ] = X2_mu 
      res$X2_psi[idx, ] = X2_psi
      res$X3_mu[idx, ] = X3_mu
      res$X4_mu[idx, ] = X4_mu
      
      res$V1_mu[idx , ] = V1_mu
      res$Vk_mu[, 1, idx ] = V2_mu
      res$Vk_mu[, 2, idx ] = V3_mu
      res$Vk_mu[, 3, idx ] = V4_mu
      
      
      res$L1_mu[idx, ] = L1_mu
      res$psi1[idx] = psi1
      res$haztau1[idx, ] = haztau1 
      res$hazT1[idx, ] = hazT1
      
      res$Lk_mu[, 1, idx] = L2_mu
      res$Lk_mu[, 2, idx] = L3_mu
      res$Lk_mu[, 3, idx] = L4_mu
        
      res$psik[, 1, idx] = psi2
      res$psik[, 2, idx] = psi3
      res$psik[, 3, idx] = psi4
        
      res$hazTk[, 1, idx] = hazT2
      res$hazTk[, 2, idx] = hazT3
      res$hazTk[, 3, idx] = hazT4
      
      res$haztauk[,1,idx] = haztau2
      res$haztauk[,2,idx] = haztau3
      
      res$lambdaT1[idx, ] = lambdaT1
      res$lambdaT2[idx, ] = lambdaT2
      res$lambdaT3[idx, ] = lambdaT3
      res$lambdaT4[idx, ] = lambdaT4

      res$lambda1[idx, ] = lambda1
      res$lambda2[idx, ] = lambda2
      res$lambda3[idx, ] = lambda3
    }
      
  }
  return(res)
}



