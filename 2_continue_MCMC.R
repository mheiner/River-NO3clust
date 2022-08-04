
## Libraries that I need
library(tidyverse)
library(parallel)
library(truncnorm)
source("scripts/0_fast_dmvnorm.R")
source("scripts/0_functions.R")
source("scripts/0_slice_sampling.R")

# args = "Results_anchor2_K7-11-12_KMtypesupplied_22221.RData"
args = commandArgs(trailingOnly = TRUE)
state = as.character(args[1])
load(state)

if ("it_prev" %in% ls()) {
  it_prev = it_prev + it
} else {
  it_prev = it
}
it = 0

seed = c(seed, seed[length(seed)] + 1111)
set.seed(seed[length(seed)])

n.draws = 2000
thin = 10
burn = 0

kp <- 0
n.it <- burn+thin*n.draws
kp.seq <- seq(burn+thin, n.it, by=thin)

## Matrices to Hold Draws
tclust.beta.draws <- aclust.beta.draws <-
  pclust.beta.x.draws <-
  pclust.beta.y.draws <- matrix(NA, nrow=n.draws, ncol=ncol(X.clust)-ncol(V.land)+ncol(X.land))
colnames(tclust.beta.draws) <- colnames(aclust.beta.draws) <-
  colnames(pclust.beta.x.draws) <- colnames(pclust.beta.y.draws) <-
  c(colnames(X.clust)[1:(ncol(X.clust)-3)], colnames(X.land))

time.draws <- matrix(NA, nrow=n.draws, ncol=K.time)
amp.draws <- matrix(NA, nrow=n.draws, ncol=K.amp)
phase.draws <- matrix(NA, nrow=n.draws, ncol=K.phase)

tclust.draws <- aclust.draws <- pclust.draws <- 
  matrix(NA, nrow=n.draws, ncol=length(time.clust))

time.cuts.draws <- matrix(0, nrow=n.draws, ncol=length(time.cuts))
amp.cuts.draws <- matrix(0, nrow=n.draws, ncol=length(amp.cuts))

llike.draws <- rep(0, n.draws)

# these are big; optionally comment them out in save() at the end
int.draws <- sig2.draws <- llike.site.draws <- matrix(0.0, nrow=n.draws, ncol=length(river.list)) # R makes copies instead of references


## out file
if(!("n_cv" %in% ls()) || n_cv == 0) {
  file.name <- paste0("results/Results_anchor", n.fixed, "_K", modK, "_KMtype", KMtype, "_", seed[1], "_cont", length(seed)-1, ".RData")
} else {
  file.name <- paste0("results/Results_anchor", n.fixed, "_K", modK, "_KMtype", KMtype, "_CV", n_cv, "_", seed[1], "_cont", length(seed)-1, ".RData")
}

## Run MCMC Loop
cat("continuing from iteration:", it_prev, "\n")
cat("resetting counter\n\n")
timestamp()
cat("it", 0, "of", n.it, "\n\n")
for(it in 1:n.it){
  
  ## Update station-specific intercepts and variances
  river.list <- lapply(river.list, update.intcpt.sig2, n0=sig2_n0, sig20=sig20)
  
  
  ## Update clustered time effects
  sum.x2 <- sapply(river.list, function(alist){
    (alist$df$time%*%alist$Rinv%*%alist$df$time)/alist$sig2})
  sum.xy <- sapply(river.list, function(alist){
    sumxy <- sum(alist$df$time%*%alist$Rinv%*%(log(alist$df$NO3)-alist$intcpt-
                         amp[amp.clust[alist$rn]]*cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])))
    return(sumxy/alist$sig2)
  })
  sum.x2 <- aggregate(x=sum.x2, by=list(clust=time.clust), FUN=sum)$x
  sum.xy <- aggregate(x=sum.xy, by=list(clust=time.clust), FUN=sum)$x
  post.var <- 1/(sum.x2+(1/time.pri.var))
  post.mean <- post.var*(sum.xy+time.pri.mean/time.pri.var)
  
  time.eff_cand <- numeric(K.time)
  
  trend_tries = 1
  trend_unsorted = TRUE
  while (trend_unsorted && trend_tries <= 500) {
    time.eff_cand[-kindx_trend0] = rnorm(K.time-1, post.mean[-kindx_trend0], sqrt(post.var[-kindx_trend0]))
    trend_unsorted = is.unsorted(time.eff_cand, strictly=TRUE)
    trend_tries = trend_tries + 1
  }
  
  if (trend_unsorted) { # rejection didn't work, just go one-at-a-time
    cat("Trends: rejection method didn't work. Resorting to full conditionals. Iteration:", it, "\n")
    time.eff[1] = rtruncnorm(1, mean=post.mean[1], sd=sqrt(post.var[1]), a=-Inf, b=time.eff[2])
    for (k in setdiff(2:(K.time-1), kindx_trend0)) {
      time.eff[k] = rtruncnorm(1, mean=post.mean[k], sd=sqrt(post.var[k]), a=time.eff[k-1], b=time.eff[k+1])
    }
    time.eff[K.time] = rtruncnorm(1, mean=post.mean[K.time], sd=sqrt(post.var[K.time]), a=time.eff[K.time-1], b=Inf)
  } else { # rejection did work, use it
    time.eff = time.eff_cand
  }


  ## Update clustered amplitude effects
  sum.x2 <- sapply(river.list, function(alist){
    cosx <- cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])
    (cosx%*%alist$Rinv%*%cosx)/alist$sig2})
  sum.xy <- sapply(river.list, function(alist){
    cosx <- cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])
    sumxy <- sum(cosx%*%alist$Rinv%*%(log(alist$df$NO3)-alist$intcpt-
                                                 time.eff[time.clust[alist$rn]]*alist$df$time))
    return(sumxy/alist$sig2)
  })
  sum.x2 <- aggregate(x=sum.x2, by=list(clust=amp.clust), FUN=sum)$x
  sum.xy <- aggregate(x=sum.xy, by=list(clust=amp.clust), FUN=sum)$x
  post.var <- 1/(sum.x2+(1/amp.pri.var))
  post.mean <- post.var*(sum.xy+amp.pri.mean/amp.pri.var)

  amp_cand <- numeric(K.amp)

  amp_tries = 1
  amp_unsorted = TRUE
  while (amp_unsorted && amp_tries <= 500) {
    amp_cand = rtruncnorm(K.amp, post.mean, sqrt(post.var), a=0, b=Inf)
    amp_unsorted = is.unsorted(amp_cand, strictly=TRUE)
    amp_tries = amp_tries + 1
  }

  if (amp_unsorted) { # rejection didn't work, just go one-at-a-time
    cat("Amplitudes: rejection method didn't work. Resorting to full conditionals. Iteration:", it, "\n")
    amp[1] = rtruncnorm(1, mean=post.mean[1], sd=sqrt(post.var[1]), a=0.0, b=amp[2])
    for (k in 2:(K.amp-1)) {
      amp[k] = rtruncnorm(1, mean=post.mean[k], sd=sqrt(post.var[k]), a=amp[k-1], b=amp[k+1])
    }
    amp[K.amp] = rtruncnorm(1, mean=post.mean[K.amp], sd=sqrt(post.var[K.amp]), a=amp[K.amp-1], b=Inf)
  } else { # rejection did work, use it
    amp = amp_cand
  }

  
  ## Update clustered phase effects via slice sampling
  phase = mclapply(1:K.phase, function(k) {
    
    indx_k = which(phase.clust == k)
    
    if (length(indx_k) >= 1) {
      
      ## define function to evaluate llik for phase cluster k
      llik_phase = function(ph) {
        llik_rn = mclapply(indx_k, function(rn) {
          
          dev = log(river.list[[rn]]$df$NO3) - river.list[[rn]]$intcpt - 
            time.eff[time.clust[rn]] * river.list[[rn]]$df$time - 
            amp[amp.clust[rn]] * cos(2.0*pi*river.list[[rn]]$df$time - ph)
          
          res <- forwardsolve(river.list[[rn]]$Rchol, dev)
          
          return(-0.5 * sum(res^2) / river.list[[rn]]$sig2)
        }, mc.cores=n.cores) %>% do.call(c, .)
        
        return( sum(llik_rn) )
      }
      
      ## define the slice
      ly = llik_phase(phase[k]) + log(runif(1, 0.0, 1.0))
      
      ## slice sample phase[k]
      phase[k] = slice_shrink(phase[k], ly, llik_phase, L=phase.cuts[k], R=phase.cuts[k+1]) 
      
    } else {
      phase[k] = runif(1, min=phase.cuts[k], max=phase.cuts[k+1])
    } 
  },
  mc.cores=min(n.cores, K.phase)
  ) %>% do.call(c,.)
  
  
  ## Update time cluster membership indicator
  time.clust <- mclapply(river.list, function(alist){

    xb = X.clust[alist$rn,]%*%tclust.beta
    resid0 = log(alist$df$NO3) - alist$intcpt -
      amp[amp.clust[alist$rn]]*cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])
    
    ## Check llike on all possible combinations
    llike <- sapply(1:K.time, function(k){
      tclust.prob <- pnorm(time.cuts[k+1], xb, sd=1) -
                     pnorm(time.cuts[k], xb, sd=1)
      resid <- resid0 - alist$df$time*time.eff[k]
      ret <- -0.5 * t(resid) %*% alist$Rinv %*% resid / alist$sig2 + log(tclust.prob)
      return(ret)
    })
    log.sum <- max(llike) + log(sum(exp(llike-max(llike))))
    
    ## Normalize to probability scale
    post.probs <- exp(llike-log.sum)
    
    ## Sample which cluster to assign to
    return(sample(1:K.time, size=1, prob=post.probs))
    
  }, mc.cores=n.cores) %>% do.call(c, .)
  
  if (n.fixed > 0) {
    time.clust[fixed.time.clust] <- col(fixed.time.clust)
  }
  
  
  ## Update amplitude and phase cluster indicators separately
  amp.clust <- mclapply(river.list, function(alist){

    XcBa = X.clust[alist$rn,]%*%aclust.beta
    resid0 <- log(alist$df$NO3) - alist$intcpt - alist$df$time*time.eff[time.clust[alist$rn]]
    cos0 <- cos(2*pi*alist$df$time - phase[phase.clust[alist$rn]])
    
    ## Check llike on all possible assignments
    llike <- sapply(1:K.amp, function(k){
      aclust.prob <- pnorm(amp.cuts[k+1], XcBa, sd=1) -
                     pnorm(amp.cuts[k], XcBa, sd=1)
      
      resid = resid0 - amp[k]*cos0
      
      ret <- -0.5 * t(resid)%*%alist$Rinv%*%resid / alist$sig2 +
            log(aclust.prob)
      return(ret)
    })
    log.sum <- max(llike) + log(sum(exp(llike-max(llike))))
    
    ## Normalize to probability scale
    post.probs <- exp(llike-log.sum)
    
    ## Sample which cluster to assign to
    return(sample(1:K.amp, size=1, prob=post.probs))
    
  }, mc.cores=n.cores) %>% do.call(c, .)
  
  if (n.fixed > 0) {
    amp.clust[fixed.amp.clust] <- col(fixed.amp.clust)
  }


  phase.clust <- mclapply(river.list, function(alist){

    XcBp_x = X.clust[alist$rn,]%*%pclust.beta.x
    XcBp_y = X.clust[alist$rn,]%*%pclust.beta.y
    resid0 = log(alist$df$NO3) - alist$intcpt - alist$df$time*time.eff[time.clust[alist$rn]]
    amps0 = amp[amp.clust[alist$rn]]
    
    ## Check llike on all possible assignments
    llike <- sapply(1:K.phase, function(k){
      phase.prob <- pprojnorm(c(XcBp_x, XcBp_y),
                              lwr=phase.cuts[k], upr=phase.cuts[k+1])
      resid <- resid0 - amps0*cos(2*pi*alist$df$time - phase[k])
      ret <- -0.5 * t(resid)%*%alist$Rinv%*%resid / alist$sig2 +
             log(phase.prob)
      return(ret)
    })
    log.sum <- max(llike) + log(sum(exp(llike-max(llike))))
    
    ## Normalize to probability scale
    post.probs <- exp(llike-log.sum)
    
    ## Sample which cluster to assign to
    return(sample(1:K.phase, size=1, prob=post.probs))
    
  }, mc.cores=n.cores) %>% do.call(c, .)

  
  
  
  ### Update cut points for membership to clusters
  
  ## time trend cut points
  tclust_means = X.clust %*% tclust.beta

  if (K.time > 2) {
    
    llik_tcut = function(cutpoints_to_update) {
      if (is.unsorted(c(time.cuts[2], cutpoints_to_update), strictly = TRUE)) {
        out = -Inf
      } else {
        cp_now = c(-Inf, time.cuts[2], cutpoints_to_update, Inf)
        probs = pnorm(cp_now[time.clust+1], tclust_means, 1.0) - pnorm(cp_now[time.clust], tclust_means, 1.0)
        out = sum(log(probs)) + sum(dnorm(cutpoints_to_update, mean=time.cuts[2], sd=15.0, log=TRUE)) # wide proper prior
        if (is.infinite(out)) stop("Error: distance between trend cut points too small.")
      }
      return(out)
    }
    
    ly = llik_tcut(time.cuts[3:K.time]) + log(runif(1, 0.0, 1.0))
    
    time.cuts[3:K.time] = slice_shrink_hyperrect(time.cuts[3:K.time], ly, llik_tcut, w=rep(0.7, K.time - 2))
  }
  
  rm(llik_tcut, ly, tclust_means)
  

  ### amplitude cluster cut points
  aclust_means = X.clust %*% aclust.beta

  if (K.amp > 2) {
    
    llik_acut = function(cutpoints_to_update) {
      if (is.unsorted(c(amp.cuts[2], cutpoints_to_update), strictly = TRUE)) {
        out = -Inf
      } else {
        cp_now = c(-Inf, amp.cuts[2], cutpoints_to_update, Inf)
        probs = pnorm(cp_now[amp.clust+1], aclust_means, 1.0) - pnorm(cp_now[amp.clust], aclust_means, 1.0)
        out = sum(log(probs)) + sum(dnorm(cutpoints_to_update, mean=amp.cuts[2], sd=15.0, log=TRUE)) # wide proper prior
        if (is.infinite(out)) stop("Error: distance between amplitude cut points too small.")
      }
      return(out)
    }
    
    ly = llik_acut(amp.cuts[3:K.amp]) + log(runif(1, 0.0, 1.0))
    
    amp.cuts[3:K.amp] = slice_shrink_hyperrect(amp.cuts[3:K.amp], ly, llik_acut, w=rep(0.7, K.amp - 2))
  }
  
  rm(llik_acut, ly, aclust_means)
  
  ## Update latent membership continuous RVs
  u.time <- rtruncnorm(length(time.clust), a=time.cuts[time.clust],
                       b=time.cuts[time.clust+1], mean=X.clust%*%tclust.beta,
                       sd=1)
  u.amp <- rtruncnorm(length(amp.clust), a=amp.cuts[amp.clust],
                      b=amp.cuts[amp.clust+1], mean=X.clust%*%aclust.beta,
                      sd=1)
  u.phase <- rprojnorm(n=length(phase.clust), 
                       mn=cbind(X.clust%*%pclust.beta.x, X.clust%*%pclust.beta.y),
                       lwr=phase.cuts[phase.clust],
                       upr=phase.cuts[phase.clust+1])
  r.phase <- sapply(1:length(u.phase), sample.radius)

  xy.phase <- cbind(r.phase*cos(u.phase), r.phase*sin(u.phase))
  
  ## Update cluster membership coefficients
  bta.post.mn <- bta.post.var%*%(t(X.clust)%*%u.time)
  tclust.beta <- bta.post.mn + bta.post.var.chol%*%rnorm(nrow(bta.post.mn))
  
  bta.post.mn <- bta.post.var%*%(t(X.clust)%*%u.amp)
  aclust.beta <- bta.post.mn + bta.post.var.chol%*%rnorm(nrow(bta.post.mn))
  
  bta.post.mn <- bta.post.var%*%(t(X.clust)%*%xy.phase[,1])
  pclust.beta.x <- bta.post.mn + bta.post.var.chol%*%rnorm(nrow(bta.post.mn))

  bta.post.mn <- bta.post.var%*%(t(X.clust)%*%xy.phase[,2])
  pclust.beta.y <- bta.post.mn + bta.post.var.chol%*%rnorm(nrow(bta.post.mn))

  
  ## Keep draws
  if(it %in% kp.seq){
    kp <- kp + 1
    time.draws[kp,] <- time.eff
    amp.draws[kp,] <- amp
    phase.draws[kp,] <- phase

    tclust.beta.draws[kp,] <- c(tclust.beta[1:(length(tclust.beta)-3)],
                                V.land%*%tclust.beta[(length(tclust.beta)-2):length(tclust.beta)])
    aclust.beta.draws[kp,] <- c(aclust.beta[1:(length(aclust.beta)-3)],
                                V.land%*%aclust.beta[(length(aclust.beta)-2):length(aclust.beta)])
    pclust.beta.x.draws[kp,] <- c(pclust.beta.x[1:(length(pclust.beta.x)-3)],
                                  V.land%*%pclust.beta.x[(length(pclust.beta.x)-2):length(pclust.beta.x)])
    pclust.beta.y.draws[kp,] <- c(pclust.beta.y[1:(length(pclust.beta.y)-3)],
                                  V.land%*%pclust.beta.y[(length(pclust.beta.y)-2):length(pclust.beta.y)])
    tclust.draws[kp,] <- time.clust
    aclust.draws[kp,] <- amp.clust
    pclust.draws[kp,] <- phase.clust
    time.cuts.draws[kp,] <- time.cuts
    amp.cuts.draws[kp,] <- amp.cuts
    
    int.draws[kp,] <- sapply(river.list, function(x) x$intcpt)
    sig2.draws[kp,] <- sapply(river.list, function(x) x$sig2)
    
    llike.site.draws[kp,] <- sapply(river.list, get.llike)
  }
  
  if (it %% 1000 == 0) {
    timestamp()
    cat("it", it, "of", n.it, "\n\n")
    save.image(file=file.name) # save results as we go...
  }
}

llike.draws <- rowSums(llike.site.draws)

## Save Draws at the End
save.image(file=file.name)

quit(save="no")