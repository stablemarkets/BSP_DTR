
split_time = function(obs_time, int_size=.5){
  tvec = seq(0, max(obs_time, na.rm=T)+.001,length.out = 30 )
  ts = tvec[1:(length(tvec)-1)]
  te = tvec[2:length(tvec)]
  tau = length(ts)
  return(list(ts=ts, te=te, tau=tau))
}

H0 = function(t, rate){ return(-1*pexp(t, rate, F, T)) }

## as opposed to log_lik_i, this computes log-likelihood evaluation for 
## entire sample.
log_lik = function(theta_x, lambda, xm, y, delta, n, ts, te, tau){
  dt = te-ts
  log_lik  = 0
  
  for( t in 1:tau){
    inds1 = y > te[t]
    inds2 = y > ts[t] & y<=te[t]
    
    w1 = sum( dt[t]*exp( xm[inds1,, drop=F] %*% theta_x ) )
    w2 = sum( (y[inds2] - ts[t])*exp( xm[inds2, , drop=F] %*% theta_x ) )
    
    log_lik = log_lik  - 1*(w1 + w2)*lambda[t]
  }
  
  log_lik = log_lik + sum( (xm %*% theta_x) * delta )
  
  return(log_lik)
}

## vector subject likelihood evaluations - useful for imputing subject-specific
## covariates
log_lik_i = function(eta, lambda, y, delta, n, ts, te, tau){
  
  dt = te-ts
  log_lik  = rep(0, n)
  
  for( t in 1:tau){
    inds1 = as.numeric(y > te[t])
    inds2 = as.numeric(y > ts[t] & y<=te[t])
    
    w1 = dt[t]*exp( eta )*inds1
    w2 = (y - ts[t])*exp( eta )*inds2
    
    log_lik = log_lik  - 1*(w1 + w2)*lambda[t]
  }
  
  log_lik = log_lik + (eta) * delta
  return(log_lik)
}


update_haz_bsl = function(lambda, ts, te, tau, delta, y, theta, xm, bsc){
  
  dt = te - ts
  
  b = rep(bsc, tau)
  
  ## divide by dt to center prior around hazard rate, not increment in cumhaz
  a = b*( H0(te, 1/mean(y) ) - H0(ts, 1/mean(y) ) ) /dt
  
  for( t in 1:tau ){
    
    ## update lambda
    inds1 = y > te[t]
    inds2 = y> ts[t] & y<=te[t]
    
    w1 = sum( dt[t]*exp( xm[inds1,, drop=F] %*% theta ) )
    w2 = sum( (y[inds2] - ts[t])*exp( xm[inds2, , drop=F] %*% theta ) )
    nd = sum( inds2 & delta!=0 )
    
    rate = b[t] + w1 + w2
    shape = a[t] + nd
    
    ## add tiny number to prevent numerically 0 draws
    lambda[t] = rgamma(n = 1, shape = shape, rate = rate) + 1e-100
  }
  
  # update eps via conjugacy
  #eps = rgamma(1, shape=a0 + tau, rate = (b0 + sum(cv)) )
  
  return(lambda)
}

update_haz_theta = function(theta, lambda,
                            xm, y, delta, n, ts, te, tau, theta_jump_v, 
                            theta_prior_var, p){
  
  # theta_star = rnorm(p, theta, sqrt(theta_jump_v) )
  
  theta_star = as.numeric(mvtnorm::rmvnorm(1, theta, sigma = (5.76/p)*theta_jump_v ) )
  
  log_post_star = log_lik(theta_star, lambda, xm, y, delta, n, ts, te, tau) + 
    sum(dnorm(theta_star, 0,sqrt(theta_prior_var), T))
  
  log_post_curr = log_lik(theta, lambda, xm, y, delta, n, ts, te, tau) + 
    sum(dnorm(theta, 0,sqrt(theta_prior_var), T))
  
  a = min(c(1, exp(log_post_star - log_post_curr)))
  
  if( runif(1)<a  ){ theta = theta_star }
  return(theta)
}

update_haz_gp = function(theta, lambda, x, p, n, event, time, ts, te, tau, 
                      bsc=.001, theta_jump_v = .001, theta_prior_var = 9){
  
  lambda = update_haz_bsl(lambda, ts, te, tau, event, time, theta, x, bsc)
  
  theta = update_haz_theta(theta, lambda, x, time, event, n, ts, te, tau, 
                          theta_jump_v, theta_prior_var, p)
  
  return(list(lambda=lambda, theta=theta))
}


update_beta = function( n_y, sum_y){
  ## draw from postioer under flat beta(1,1) prior
  draw = rbeta(n = 1, shape1 = 1 + sum_y, shape2 = 1 + n_y - sum_y )
  return(draw)
}

update_logit = function(beta_cur, x, p, y, prop_sd=.1){
  
  beta_star = rnorm(p, beta_cur, prop_sd)
  
  
  lik_star = sum(dbern(y, invlogit( x %*% beta_star ), T ))
  pr_star = sum(dnorm(beta_star, mean = 0, sd = 3, log = T))
  
  lik_cur = sum(dbern(y, invlogit( x %*% beta_cur ), T ))
  pr_cur = sum(dnorm(beta_cur, mean = 0, sd = 3, log = T))
  
  ratio = min(c(1, exp( lik_star + pr_star - lik_cur - pr_cur )  ))
  
  U = 1*(runif(1) < ratio)
  draw = U*beta_star + (1-U)*beta_cur
  if(sum(is.na(draw))>1) browser()
  return(draw)
}

rcond_post_beta = function(y, xm, psi, beta_prior_mean, beta_prior_var){
  
  mu_beta = beta_prior_mean
  #v_beta = diag(beta_prior_var)
  v_beta_inv = diag(1/beta_prior_var)
  xtx = t(xm)%*%xm
  
  post_cov = chol2inv(chol( v_beta_inv + (1/psi)*xtx))
  
  post_mean = post_cov %*% (v_beta_inv %*% mu_beta + (1/psi)*t(xm)%*%y )
  
  draw = t(mvtnorm::rmvnorm(n = 1, mean = post_mean, sigma = post_cov))
  return(draw)  
}
rcond_post_psi = function(n, beta, y, xm, g1, b1){
  shape_k = g1 + n/2
  rate_k = .5*(  sum( (y - xm %*% beta)^2 ) )  + b1
  draw = invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(draw)
}



