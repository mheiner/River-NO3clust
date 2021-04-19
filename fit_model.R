##
## Code for clustering time, amplitude and phase
## for France time series
##

## Libraries that I need
library(tidyverse)
library(parallel)
library(truncnorm)
source("scripts/0_AMCMCUpdate.R")
source("scripts/0_fast_dmvnorm.R")

args = as.numeric(commandArgs(trailingOnly = TRUE))

seed = args[1] # 32321 used for paper
set.seed(seed)

## Number of Cores
n.cores <- min(90, detectCores())

## Model Settings
K.time <- args[2]
K.amp <- args[3]

print(paste("seed =", seed))
print(paste("Ktime =", K.time))
print(paste("Kamp =", K.amp))

K.phase <- 12
modK = paste(K.time, K.amp, K.phase, sep="-")

## Priors
time.pri.var <- 1000
time.pri.mean <- 0
int.pri.mean <- 0
int.pri.var <- 1000
amp.pri.mean <- 1
amp.pri.var <- 1000

## MCMC Settings
n.draws <- 100 # 2000 for paper
thin <- 1 # 12 for paper
burn <- 20 # 24 for paper
kp <- 0
n.it <- burn+thin*n.draws
kp.seq <- seq(burn+thin, n.it, by=thin)

## Functions that I need
source("scripts/0_functions.R")

## Read in the Data
source("scripts/0_readData.R")

## Initialize clustering, fixed estimates
source("scripts/0_init.R")


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



## AMCMC for cuts
time.cuts.amcmc <- list(mn=matrix(0, nrow=length(time.trans.cuts), ncol=1),
                        var=matrix(0, nrow=length(time.trans.cuts), ncol=length(time.trans.cuts)))
amp.cuts.amcmc <- list(mn=matrix(0, nrow=length(amp.trans.cuts), ncol=1),
                        var=matrix(0, nrow=length(amp.trans.cuts), ncol=length(amp.trans.cuts)))
amcmc.it <- 250

## Run MCMC Loop
pb <- txtProgressBar(min = 0, max = n.it, style = 3)
for(it in 1:n.it){
  
  ## Update station-specific intercepts and variances
  river.list <- lapply(river.list, update.intcpt.sig2)
  
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
  
  time.eff <- numeric(K.time)
  time.eff[1] <- rnorm(1, mean=post.mean[1], sd=sqrt(post.var[1]))
  for (k in 2:K.time) {
    time.eff[k] <- rtruncnorm(1, mean=post.mean[k], sd=sqrt(post.var[k]), a=time.eff[k-1], b=Inf)
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

  amp <- numeric(K.amp)
  amp[1] <- rtruncnorm(1, mean=post.mean[1], sd=sqrt(post.var[1]), a=0, b=Inf)
  for (k in 2:K.amp) {
    amp[k] <- rtruncnorm(1, mean=post.mean[k], sd=sqrt(post.var[k]), a=amp[k-1], b=Inf)
  }
  
  ## Update clustered phase effects via Griddy gibbs
  for(k in 1:K.phase){
    if(sum(phase.clust==k) >= 1){
      p.llike <- mclapply(which(phase.clust==k), function(rn){
        res <- forwardsolve(river.list[[rn]]$Rchol, mapply(FUN=function(p){
          log(river.list[[rn]]$df$NO3) - river.list[[rn]]$intcpt - 
            river.list[[rn]]$df$time*time.eff[time.clust[rn]] - 
            amp[amp.clust[rn]]*cos(2*pi*river.list[[rn]]$df$time-p)},
          p=seq(phase.cuts[k], phase.cuts[k+1], length=100)))
      return(-0.5*colSums(res^2)/river.list[[rn]]$sig2)
      }, mc.cores=n.cores) %>% do.call(cbind, .)
      if( is.character(p.llike) ) {
        print("p.llike is a character vector:")
        save(p.llike, file="pllike_error_dump.rda")
        print(p.llike[grep("[a-z]", p.llike)[1]])
        p.llike = matrix(as.numeric(p.llike), dim(p.llike))
      }
      post.probs <- rowSums(p.llike)
      post.probs <- post.probs - log.sum(post.probs)
      phase[k] <- sample(seq(phase.cuts[k], phase.cuts[k+1], length=100),
                         1, prob=exp(post.probs)) 
    } else {
      phase[k] <- sample(seq(phase.cuts[k], phase.cuts[k+1], length=100),1)
    }
  }
  
  
  ## Update time cluster membership indicator
  time.clust <- mclapply(river.list, function(alist){
  
    ## Check llike on all possible combinations
    llike <- sapply(1:K.time, function(x){
      tclust.prob <- pnorm(time.cuts[x+1], X.clust[alist$rn,]%*%tclust.beta, sd=1) -
        pnorm(time.cuts[x], X.clust[alist$rn,]%*%tclust.beta, sd=1)
      resid <- log(alist$df$NO3)-alist$intcpt-alist$df$time*time.eff[x] -
        amp[amp.clust[alist$rn]]*cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])
      ret <- -0.5*t(resid)%*%alist$Rinv%*%resid/alist$sig2 + log(tclust.prob)
      return(ret)
    })
    log.sum <- max(llike) + log(sum(exp(llike-max(llike))))
    
    ## Normalize to probability scale
    post.probs <- exp(llike-log.sum)
    
    ## Sample which cluster to assign to
    return(sample(1:K.time, size=1, prob=post.probs))
    
  }, mc.cores=n.cores) %>% do.call(c, .)
  time.clust[fixed.time.clust] <- col(fixed.time.clust)
  
  ## Update amplitude and phase cluster indicator
  allclusts <- mclapply(river.list, function(alist){
    ## Expand on cluster indicators
    all.clusters <- expand.grid(1:K.amp, 1:K.phase)

    ## Check llike on all possible combinations
    llike <- apply(all.clusters, 1, function(x){
      aclust.prob <- pnorm(amp.cuts[x[1]+1], X.clust[alist$rn,]%*%aclust.beta, sd=1) -
        pnorm(amp.cuts[x[1]], X.clust[alist$rn,]%*%aclust.beta, sd=1)
      phase.prob <- pprojnorm(c(X.clust[alist$rn,]%*%pclust.beta.x,
                                X.clust[alist$rn,]%*%pclust.beta.y),
                              lwr=phase.cuts[x[2]], upr=phase.cuts[x[2]+1])
      resid <- log(alist$df$NO3)-alist$intcpt-alist$df$time*time.eff[time.clust[alist$rn]] -
        amp[x[1]]*cos(2*pi*alist$df$time-phase[x[2]])
      ret <- -0.5*t(resid)%*%alist$Rinv%*%resid/alist$sig2 +
        log(aclust.prob) + log(phase.prob)
      return(ret)
    })
    log.sum <- max(llike) + log(sum(exp(llike-max(llike))))

    ## Normalize to probability scale
    post.probs <- exp(llike-log.sum)

    ## Sample which cluster to assign to
    return(all.clusters[sample(1:nrow(all.clusters), size=1, prob=post.probs),])

  }, mc.cores=n.cores) %>% do.call(rbind, .)
  amp.clust <- allclusts[,1]
  amp.clust[fixed.amp.clust] <- col(fixed.amp.clust)
  phase.clust <- allclusts[,2]
  
  ## Update cutpoints for membership to clusters
  if(it>amcmc.it){
    prop.var <- (2.4^2/length(time.trans.cuts))*(0.0001*diag(length(time.trans.cuts))+
                                                   time.cuts.amcmc$var)
  } else {
    prop.var <- 0.0001*diag(length(time.trans.cuts))
  }
  prop.trans.cuts <- time.trans.cuts + t(chol(prop.var))%*%rnorm(nrow(prop.var))
  prop.cuts <- c(-Inf, time.cuts[2], time.cuts[2]+cumsum(exp(prop.trans.cuts)), Inf)
  prop.probs <- pnorm(prop.cuts[time.clust+1], mean=X.clust%*%tclust.beta, sd=1) -
    pnorm(prop.cuts[time.clust], mean=X.clust%*%tclust.beta, sd=1)
  cur.probs <- pnorm(time.cuts[time.clust+1], mean=X.clust%*%tclust.beta, sd=1) -
    pnorm(time.cuts[time.clust], mean=X.clust%*%tclust.beta, sd=1)
  MH <- sum(log(prop.probs)-log(cur.probs))
  if(log(runif(1, 0, 1)) < MH){
    time.trans.cuts <- prop.trans.cuts
    time.cuts <- prop.cuts
  }
  time.cuts.amcmc <- AMCMC.update(time.trans.cuts, time.cuts.amcmc$mn,
                                  time.cuts.amcmc$var, it)

  if(it>amcmc.it){
    prop.var <- (2.4^2/length(amp.trans.cuts))*(0.0001*diag(length(amp.trans.cuts))+
                                                  amp.cuts.amcmc$var)
  } else {
    prop.var <- 0.0001*diag(length(amp.trans.cuts))
  }
  prop.trans.cuts <- amp.trans.cuts + t(chol(prop.var))%*%rnorm(nrow(prop.var))
  prop.cuts <- c(-Inf, amp.cuts[2], amp.cuts[2]+cumsum(exp(prop.trans.cuts)), Inf)
  prop.probs <- pnorm(prop.cuts[amp.clust+1], mean=X.clust%*%aclust.beta, sd=1) -
    pnorm(prop.cuts[amp.clust], mean=X.clust%*%aclust.beta, sd=1)
  cur.probs <- pnorm(amp.cuts[amp.clust+1], mean=X.clust%*%aclust.beta, sd=1) -
    pnorm(amp.cuts[amp.clust], mean=X.clust%*%aclust.beta, sd=1)
  MH <- sum(log(prop.probs)-log(cur.probs))
  if(log(runif(1, 0, 1)) < MH){
    amp.trans.cuts <- amp.trans.cuts
    amp.cuts <- prop.cuts
  }
  amp.cuts.amcmc <- AMCMC.update(amp.trans.cuts, amp.cuts.amcmc$mn,
                                 amp.cuts.amcmc$var, it)



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
    
  ## Increment progress bar
  setTxtProgressBar(pb, it)
}
close(pb)

llike.draws <- rowSums(llike.site.draws)

## Save Draws at the End
sysDate <- format(Sys.Date(), "%m%d%Y")
file.name <- paste0("./results/GGConstrainedResults_K", modK, "_", sysDate, ".RData")
save(file=file.name, list=c("time.draws", "amp.draws", "phase.draws",
                                "tclust.beta.draws", "aclust.beta.draws",
                                "pclust.beta.x.draws", "pclust.beta.y.draws",
                                "tclust.draws", "aclust.draws", "pclust.draws",
                            "time.cuts.draws", "amp.cuts.draws", "phase.cuts",
                            "llike.draws",
                            "llike.site.draws", "int.draws", "sig2.draws"
                            ))
