

gcomp_bayes = function(mcmc_res, M, rule, tau_vec){
  
  simulate_time = function(M, eta, lambda, ts, te, tau, int_size){
    u = rexp(M, 1)*exp(-1*eta)
    ext_lambda = c(lambda, lambda[tau])
    ext_ts = c(ts, ts[tau])
    cslambda = cumsum(lambda*int_size)
    hlwr = c(0, cslambda)
    hupr = c(cslambda, Inf)
    
    kvec = findInterval(u, hupr) + 1
    
    draw = (u - hlwr[kvec]) / ext_lambda[kvec] + ext_ts[kvec]
    
    return(draw)
  }
  
  ## assuming equally spaced intervals at all time points
  int_size1 = mcmc_res$sp1$te[1] -  mcmc_res$sp1$ts[1]
  int_size2 = mcmc_res$sp2$te[1] -  mcmc_res$sp2$ts[1]
  int_size3 = mcmc_res$sp3$te[1] -  mcmc_res$sp3$ts[1]
  int_size4 = mcmc_res$sp4$te[1] -  mcmc_res$sp4$ts[1]
  
  mcmc_iter = length(mcmc_res$psi1)
  
  prop_gt_tau = function( tau ) mean( st > tau )
  
  post_scurve = matrix(NA, nrow=length(tau_vec), ncol=mcmc_iter)
  
  for( j in 1:mcmc_iter ){
    ######## -----------         time point 1  -------------------------########
    ## simulate baseline confounders
    X1 = rnorm(M, mcmc_res$X1_mu[j, 1], sqrt(mcmc_res$X1_psi[j,1] ) )
    X2 = rnorm(M, mcmc_res$X2_mu[j, 1], sqrt(mcmc_res$X2_psi[j,1] ) )
    X3 = rbern(M, mcmc_res$X3_mu[j, 1 ])
    X4 = rbern(M, mcmc_res$X4_mu[j, 1 ])
    
    X_bsl = cbind(X1, X2, X3, X4)
    
    L1 = rnorm(M, mcmc_res$L1_mu[j, 1], sqrt( mcmc_res$psi1[j] ) )    
    V1 = rbern(M, mcmc_res$V1_mu[j, 1 ])
      
    ## simulate time to death from trt 1, T1
    a1 = rule(L1)
    x_1 = cbind(a1, L1, V1, X_bsl )
    eta_T1 = x_1 %*% t(mcmc_res$hazT1[ j , , drop=F] )
    eta_tau1 =  x_1 %*% t(mcmc_res$haztau1[j, , drop=F])
    
    T1 = simulate_time( M, eta_T1, mcmc_res$lambdaT1[j,], mcmc_res$sp1$ts, 
                        mcmc_res$sp1$te, mcmc_res$sp1$tau, int_size1 )
    tau1 = simulate_time( M, eta_tau1, mcmc_res$lambda1[j,], mcmc_res$sp1$ts, 
                          mcmc_res$sp1$te, mcmc_res$sp1$tau, int_size1 )
    
    ######## -----------         time point 2  -------------------------########
    
    mu_L2 = cbind(1, L1, a1) %*%  mcmc_res$Lk_mu[, 1 , j, drop=F ]
    L2 = rnorm(M, mu_L2, sqrt( mcmc_res$psik[ , 1, j ] ) )
    
    mu_V2 = invlogit( cbind(1, V1, a1) %*%  mcmc_res$Vk_mu[, 1 , j, drop=F ] )
    V2 = rbern(M, mu_V2 )
    
    a2 = rule(L2)
    
    x_2 = cbind( a2, a1, L2, L1, V2, V1, X_bsl, 1*(tau1>4) )
    
    eta_T2 = x_2 %*% mcmc_res$hazTk[,1, j, drop=F]
    eta_tau2 = x_2 %*% mcmc_res$haztauk[,1,j,drop=F]
    T2 = simulate_time( M, eta_T2, mcmc_res$lambdaT2[j,], mcmc_res$sp2$ts, 
                        mcmc_res$sp2$te, mcmc_res$sp2$tau, int_size2 )
    tau2 = simulate_time( M, eta_tau2, mcmc_res$lambda2[j,], mcmc_res$sp2$ts, 
                          mcmc_res$sp2$te, mcmc_res$sp2$tau, int_size2 )
    
    ######## -----------         time point 3  -------------------------########

    mu_L3 = cbind(1, L2, a2) %*%  mcmc_res$Lk_mu[, 2 , j, drop=F ]
    L3 = rnorm(M, mu_L3, sqrt( mcmc_res$psik[ , 2, j ] ) )
    
    mu_V3 = invlogit( cbind(1, V2, a2) %*%  mcmc_res$Vk_mu[, 2 , j, drop=F ] )
    V3 = rbern(M, mu_V3 )
    
    a3 = rule(L3)
    
    x_3 = cbind( a3, a2, L3, L2, V3, V2, X_bsl, 1*(tau2>4) )
    
    eta_T3 = x_3 %*% mcmc_res$hazTk[,2, j, drop=F]
    eta_tau3 = x_3 %*% mcmc_res$haztauk[,2,j,drop=F]
    T3 = simulate_time( M, eta_T3, mcmc_res$lambdaT3[j,], mcmc_res$sp3$ts, 
                        mcmc_res$sp3$te, mcmc_res$sp3$tau, int_size3 )
    tau3 = simulate_time( M, eta_tau3, mcmc_res$lambda3[j,], mcmc_res$sp3$ts, 
                          mcmc_res$sp3$te, mcmc_res$sp3$tau, int_size3 )
    
    ######## -----------         time point 4  -------------------------########
    
    mu_L4 = cbind(1, L3, a3) %*%  mcmc_res$Lk_mu[, 3 , j, drop=F ]
    L4 = rnorm(M, mu_L4, sqrt( mcmc_res$psik[ , 3, j ] ) )
    
    mu_V4 = invlogit( cbind(1, V3, a3) %*%  mcmc_res$Vk_mu[, 3 , j, drop=F ] )
    V4 = rbern(M, mu_V4 )
    
    a4 = rule(L4)
    
    x_4 = cbind( a4, a3, L4, L3, V4, V3, X_bsl, 1*(tau3>4) )
    
    eta_T4 = x_4 %*% mcmc_res$hazTk[,3, j, drop=F]
    T4 = simulate_time( M, eta_T4, mcmc_res$lambdaT4[j,], mcmc_res$sp4$ts, 
                        mcmc_res$sp4$te, mcmc_res$sp4$tau, int_size4 )
    
    ## add up times
    st = ifelse(T1 < tau1, T1, tau1)
    atrisk = T1>tau1
    st[ atrisk ] = st[  atrisk ] + ifelse(T2[atrisk] < tau2[atrisk], T2[atrisk], tau2[atrisk])
    atrisk = T1>tau1 & T2>tau2 
    st[ atrisk ] = st[  atrisk ] + ifelse(T3[atrisk] < tau3[atrisk], T3[atrisk], tau3[atrisk])
    atrisk = T1>tau1 & T2>tau2 & T3>tau3
    st[ atrisk ] = st[  atrisk ] + T4[atrisk]
    
    post_scurve[ , j] = sapply(tau_vec, prop_gt_tau )
    
    if( j %% 1000==0 ){
      cat(paste0('iteration ',j,'\n'))
    }
  }
  return(post_scurve)
}