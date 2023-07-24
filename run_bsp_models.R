## This code assumes the working directory is the top-level directory
##   which contains the README.md file. If it is not, use setwd()
##   to set the working directory.
set.seed(1)

## load dependencies R packages
library(LaplacesDemon)
library(flexsurv)
library(invgamma)
library(latex2exp)

## load MCMC and gcomputation helper functions 
source("helper_functions/bayes_dtr_mcmc.R")
source("helper_functions/mcmc_helpers.R")
source("helper_functions/gcomp_bayes.R")

####---------------------------- Load sample data --------------------------####
load("data/sample_data.Rdata")

## data description
## each row corresponds to a single patient; maximum of 4 treatment courses.
#### kappa: number of treatment courses lived through kappa in {1,2,3,4}
#### X1...X4: time-fixed covariates measured before trt decision in course 1.
#             X1-X2 are continuous, X3-X4 are binary.
#### Lk: for k=1,2,3,4. time-varying continuous covariate assessed before trt 
#        decision in course k=1,2,3,4.
#### Vk: for k=1,2,3,4. time-varying binary covariate assessed before trt 
#        decision in course k=1,2,3,4.
#### Ak: for k=1,2,3,4. time-varying treatment decision in course k=1,2,3,4.
#### Uk: for k=1,2,3,4. obs. waiting time from start of course k to next event.
#### delta_Tk: indicator, =1 if Uk is time to death.
#### delta_tauk: indicator, =1 if Uk is time to course k+1
#                note if delta_tauk=delta_Tk=0 indicators Uk is censoring time.

## for those with a maximum of Kappa=j, Lk+1, Vk+1, Ak+1, Uk+1, delta_Tk+1, 
## and delta_tauk+1 are NA for all k>j. E.g. if someone only lived through 
## kappa =2 treatment courses, the data for course 3 and 4 are coded as NA.

## bayes_dtr_mcmc: function runs method and outputs posterior draws of all 
## unknown parameters, as described in the paper.

####----------------------------- Posterior Sampling -----------------------####
st = Sys.time()
post_draws = bayes_dtr_mcmc(d = d,
                            ## total number of mcmc/burnin iterations
                            ##   e.g. will yield 10000-5000 = 5000 draws.
                            iter = 10000,  burnin = 5000,
                            thin = 1, ## thinning interval if need to save RAM
                            ## when to start/end tuning metropolis-hastings 
                            ##   proposal covariance - should always be 
                            ##   within burnin period.
                            adapt_start=2000, adapt_end = 5000)
et = Sys.time()
run_time = et - st
print(run_time)

####----------------------------- G-Computation ----------------------------####

## example dynamic treatment rules that takes as input a time-varying covariate,
##   Lk and treats if Lk>0 (rule 1) or treats if Lk<0 (rule 2).
##   the goal is to compare survival under the rules.
rule1 = function(Lk) 1*(Lk > 0)
rule2 = function(Lk) 1*(Lk < 0)

## times at which to evaluate survival curve, here we pick 100 equally spaced
## poitns between t=1 and t=20
tau_vec = seq(1, 20, length.out=100) 

## gcomp_bayes() computes marginal survival curve under rule.
##   given posterior draws in post_draws as described in paper.
##   it runs M g-comp simulations per each of the iter-burnin posterior draws.
##   it will print updates every 1,000 iterations. 
## output will be a matrix with as many rows as the length of tau_vec and as
##   many columns as there are posterior draws. Each row contains posterior draws
### of the probability of surviving beyond the time point corresponding to that row.

st = Sys.time()
post_scurve1 = gcomp_bayes(mcmc_res = post_draws,
                          ## number of g-comp simulations per posterior draw
                          M = 5000, 
                          rule = rule1, 
                          tau_vec = tau_vec)

post_scurve2 = gcomp_bayes(mcmc_res = post_draws,
                           ## number of g-comp simulations per posterior draw
                           M = 5000, 
                           rule = rule2, 
                           tau_vec = tau_vec)
et = Sys.time()
run_time = et - st
print(run_time)

####--------------------------- Plot Survival Curve ------------------------####

## compute posterior mean / 95% credible intervals for rule 1
post_mean1 = rowMeans(post_scurve1)
pci1 = apply(post_scurve1, 1, quantile, probs=c(.025, .975)) 

## compute posterior mean / 95% credible intervals for rule 2
post_mean2 = rowMeans(post_scurve2)
pci2 = apply(post_scurve2, 1, quantile, probs=c(.025, .975)) 

setwd("output/")

png("survivalplots.png", width = 500, height = 500)

par(mfrow=c(1,1))

plot(tau_vec, post_mean1, col='white', type='l', 
     lty=1, ylim=c(.295,1.01),
     main=TeX('Marginal survival curves, rule $r$'),
     ylab = TeX('$\\Psi^r(t)=P(T^r > t)$'),axes=T,
     xlab=TeX('Time since trt. sequence start, $t$'), 
     cex.main=1, cex.lab=1, cex.axis=1)

## plot curve for rule 1 (blue)
polygon(c(tau_vec,rev(tau_vec)),c(pci1[1,],rev(pci1[2,])),
        col = rgb(0, 0, .9, .2), border = FALSE)
lines(tau_vec, rowMeans(post_scurve1), col='blue', lwd=1.5)

## plot curve for rule 2 (red)
polygon(c(tau_vec,rev(tau_vec)),c(pci2[1,],rev(pci2[2,])),
        col = rgb(1, 0, .9, .2), border = FALSE)
lines(tau_vec, rowMeans(post_scurve2), col='red', lwd=1.5)

legend('bottomleft',
       legend=c('Rule 1 - GP post. mean & 95% CI',
                'Rule 2 - GP post. mean & 95% CI'), 
       fill = c( rgb(0, 0, .9, .2), 
                 rgb(1, 0, .9, .2)), 
       cex = 1,
       border = 1, bty='n')

dev.off()


