log.sum <- function(x){
  the.max <- max(x)
  lsum <- the.max+log(sum(exp(x-the.max)))
  return(lsum)
}

update.intcpt.sig2 <- function(alist){
  
  ## Center the observations
  y <- log(alist$df$NO3)-alist$df$time*time.eff[time.clust[alist$rn]] -
    amp[amp.clust[alist$rn]]*cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])
  ybar <- mean(y)
  n <- length(y)
  
  ## Calculate posterior mean using non-informative N(pri.mn, pri.var) prior
  sig2.n <- alist$sig2*alist$sumR/(n^2)
  post.var <- 1/((1/int.pri.var)+(1/sig2.n))
  post.mn <- post.var*((int.pri.mean/int.pri.var)+(ybar/sig2.n))
  
  ## Sample new intercept
  alist$intcpt <- rnorm(1, post.mn, sqrt(post.var))
  
  ## Sample new sig2
  y.cntr <- y-alist$intcpt
  ss <- t(y.cntr)%*%alist$Rinv%*%y.cntr
  astar <- 2.1+n/2
  bstar <- 0.01+0.5*ss
  alist$sig2 <- 1/rgamma(1, shape=astar, rate=bstar)
  
  ## Return updated list
  return(alist)
  
}

sample.radius <- function(ind){
  mu <- rbind(c(pclust.beta.x), c(pclust.beta.y))%*%X.clust[ind,]
  r.rng <- c(0,sqrt(sum(mu^2))+4)
  # r.rng <- c(0,4)
  b <- cos(u.phase[ind])*sum(X.clust[ind,]*pclust.beta.x) +
    sin(u.phase[ind])*sum(X.clust[ind,]*pclust.beta.y)
  rseq <- seq(r.rng[1], r.rng[2], length=1000)
  wgt <- rseq*exp(-0.5*rseq^2+b*rseq)
  return(sample(rseq, size=1, prob=wgt))
}

get.llike <- function(alist){
  
  ## Center the observations
  y <- log(alist$df$NO3) - alist$intcpt - alist$df$time*time.eff[time.clust[alist$rn]] -
    amp[amp.clust[alist$rn]]*cos(2*pi*alist$df$time-phase[phase.clust[alist$rn]])
  n <- length(y)
  
  ## Calculate the llike
  ll <- -n/2*log(alist$sig2) - 0.5*t(y)%*%alist$Rinv%*%y/alist$sig2 + 0.5*alist$Rinv_logdet$modulus
  
  ## Return log-likelihood
  return(ll)
  
}


dprojnorm <- function(angle, mn=c(0,0), log=FALSE){
  u <- cbind(cos(angle), sin(angle))
  if(length(mn)==2){
    ut.mu <- rowSums(scale(u, center=FALSE, scale=1/mn))
    mn <- matrix(mn, nrow=1)
  } else if(length(mn)>2) {
    ut.mu <- rowSums(u*mn)
  }
  pdfout <- -log(2*pi) - 0.5*rowSums(mn^2) +
    log(1+ut.mu*pnorm(ut.mu)/dnorm(ut.mu))
  if(log==TRUE){
    return(pdfout)
  } else {
    return(exp(pdfout))
  }
}

pprojnorm <- function(mn, lwr=-pi, upr=pi, prec=5000){
  ang <- atan2(rnorm(prec, mn[2], 1), rnorm(prec, mn[1], 1))
  return(mean(ang >= lwr & ang <= upr))
}

rprojnorm <- function(n, mn, lwr=-pi, upr=pi, len=30){
  # lseq <- seq(-pi, pi, length=len)
  # samp.vals <- sapply(1:nrow(mn), function(x){
  #     tlseq <- lseq[lseq>=lwr[x] & lseq<=upr[x]]
  #     wgt <- dprojnorm(tlseq, mn[x,])
  #     wgt <- wgt/sum(wgt)
  #     return(sample(tlseq, size=1, prob=wgt, replace=TRUE))
  # })
  samp.vals <- mclapply(1:nrow(mn), function(ind){
    xyvals <- expand.grid(x=seq(mn[ind,1]-4, mn[ind,1]+4, length=len),
                          y=seq(mn[ind,2]-4, mn[ind,2]+4, length=len))
    avals <- atan2(xyvals[,2], xyvals[,1])
    kp <- which(avals >= lwr[ind] & avals <= upr[ind])
    
    if (length(kp) == 1) {
      return(avals[kp]) # sample function misinterprets a vector of length one
    } else if (length(kp) > 1) {
      lwgt <- dnorm(xyvals[kp,1], mn[ind,1], 1, log=TRUE) + dnorm(xyvals[kp,2], mn[ind,2], 1, log=TRUE)
      wgt <- exp(lwgt - max(lwgt))
      return(sample(avals[kp], size=1, replace=FALSE, prob=wgt))
    } else {
      return(runif(1, lwr[ind], upr[ind]))
    }
    
  }, mc.cores=n.cores)
  out = unlist(samp.vals)
  if (is.character(out)) {
    print("rprojnorm returned character values:")
    print(out[1])
    save(out, file="rprojnorm_error_dump.rda")
    out = as.numeric(out)
  }
  if (any(is.na(out))) {
    print("rprojnorm returned missing/NaN values. Values in mn:")
    print(summary(mn))
  } else if(any(is.infinite(out))) {
    print("rprojnorm returned Inf values. Values in mn:")
    print(summary(mn))
  }
  return(out)
}

to.phase.cuts <- function(tpc){
  pc <- c(0, plogis(tpc), 1)
  probs <- sapply(1:length(pc), function(x){
    if(x==1){
      return(pc[x])
    } else{
      return(prod(1-pc[1:(x-1)])*pc[x])
    }
  })
  return(-pi+2*pi*cumsum(probs))
}