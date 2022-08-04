rm(list=ls())

## Libraries
library(tidyverse)

## Load in the Posterior Draws
# load("results/Results_anchor2_K7-11-12_KMtypesupplied_22221.RData"); rep_id = 1; modK="7-11-12"
# load("results/Results_anchor2_K7-11-12_KMtypesupplied_22222.RData"); rep_id = 2; modK="7-11-12"
# load("results/Results_anchor2_K7-11-12_KMtypesupplied_22223.RData"); rep_id = 3; modK="7-11-12"
# load("results/Results_anchor2_K7-11-12_KMtypesupplied_22224.RData"); rep_id = 4; modK="7-11-12"

# load("results/Results_anchor2_K5-7-12_KMtypesupplied_22221_cont1.RData"); rep_id = 1; modK="5-7-12"
# load("results/Results_anchor2_K5-7-12_KMtypesupplied_22222_cont1.RData"); rep_id = 2; modK="5-7-12"
# load("results/Results_anchor2_K5-7-12_KMtypesupplied_22223_cont1.RData"); rep_id = 3; modK="5-7-12"
load("results/Results_anchor2_K5-7-12_KMtypesupplied_22224.RData"); rep_id = 4; modK="5-7-12"

allrdafiles = list.files("results/")
indx_files_use = grep("anchor2_K5-7-12_KMtypesupplied_2222\\d", allrdafiles)
# indx_files_use = grep("anchor2_K5-7-12_KMtypesupplied_2222\\d", allrdafiles)
(rdafiles_use = allrdafiles[indx_files_use])

(nchains = length(rdafiles_use))
source("scripts/0_combine_chains.R")
rep_id = 1234


### at observation level

get_resid <- function(rn) {
  n = length(river.list[[rn]]$df$NO3)
  niter = nrow(int.draws[it_indx,])
  
  mean_at_obs = matrix( rep(int.draws[it_indx, rn], times=n), ncol=n ) + 
    tcrossprod(time.draws[ cbind(it_indx, tclust.draws[it_indx, rn]) ], river.list[[rn]]$df$time) +
    (matrix( rep(amp.draws[ cbind(it_indx, aclust.draws[it_indx, rn])], times=n), ncol=n ) * 
       cos(t(outer(2*pi*river.list[[rn]]$df$time, phase.draws[ cbind(it_indx, pclust.draws[it_indx, rn]) ], FUN="-")))
    )
  
  e =  (matrix(rep(log(river.list[[rn]]$df$NO3), each=niter), nrow=niter) - mean_at_obs) # niter by n matrix
  
  e_standard = e / rep(sqrt(sig2.draws[it_indx, rn]), times=n) # niter by n matrix
  
  e_whiten = forwardsolve(river.list[[rn]]$Rchol, t(e_standard)) %>% t() # niter by n matrix
  
  logdens = dnorm(e_whiten, mean=0, sd=1, log=TRUE)
  logdenst = dt(e_whiten, df=3, log=TRUE) # HACK
  
  return(list(e=e, ew=e_whiten, llik=logdens, llikt=logdenst))
}


n.draws
(niter_tot = nrow(int.draws))
n.draws == niter_tot/nchains
str(llike.site.draws)

nkeep = 1000
it_indx = floor(seq(from=1, to=niter_tot, length=nkeep)); head(it_indx); tail(it_indx)
chain_id = floor((it_indx - 1) / n.draws) + 1
table(chain_id)


resid_all = list()
for (rn in 1:length(river.list)) {
  resid_all[[rn]] = get_resid(rn)
  cat(rn, "of", length(river.list), "\r")
}

set.seed(1)
# indxs = lapply(1:length(river.list), function(rn) sample(1:ncol(resid_all[[rn]]$e), 10, replace=FALSE))
indxs = lapply(1:length(river.list), function(rn) 1:ncol(resid_all[[rn]]$e))
str(indxs)

# rn = 20
# str(resid_all[[rn]])
# hist(resid_all[[rn]]$llik[,20]) # obs (across it)
# hist(resid_all[[rn]]$llik[1,]) # it (across obs)

llik_obs = lapply(1:length(resid_all), function(rn) resid_all[[rn]]$llik[,indxs[[rn]]]) %>% do.call(cbind, .)
# llik_obs = lapply(1:length(resid_all), function(rn) resid_all[[rn]]$llik) %>% do.call(cbind, .)
str(llik_obs)

library("loo")
# indx = sample(1:ncol(llik_obs), 20e3, replace=FALSE)
indx = 1:ncol(llik_obs)
str(indx)

# r_eff = relative_eff(exp(llik_obs[,indx]), chain_id=rep(1,nrow(llik_obs)), cores=2) # single chain
r_eff = relative_eff(exp(llik_obs[,indx]), chain_id=chain_id, cores=4) # multiple chains
summary(r_eff)
loo_fit = loo(llik_obs[,indx], r_eff=r_eff, cores=4)
print(loo_fit)

## WAIC
waic(llik_obs[,indx])


## DIC (assumes normality of posterior, which is not true in mixture)
# (Dbar = -2*mean(rowSums(llik_obs)))
# (pD = 2*mean(apply(llik_obs, 1, var)))
# pD + Dbar

