source("scripts/0_fast_dmvnorm.R")

## Estimate time, amplitude and phase for each river
river.list <- split(rivers, f=rivers$CDSTATIONM)
get.coefs <- function(a.dframe){
  river.coef <- coef(lm(log(NO3)~time+cos(2*pi*time)+sin(2*pi*time), 
                        data=a.dframe))[-1]
  beta.time <- river.coef[1]
  amp <- sqrt(sum(river.coef[-1]^2))
  phase <- atan2(river.coef[3], river.coef[2])
  return(c(beta.time, amp, phase))
}
river.coefs <- sapply(river.list, get.coefs) %>% t()

## Cluster based on estimates and define cluster centers (coefficient start vals)
time.eff <- kmeans(river.coefs[,1], centers=K.time, iter.max=100)
time.clust <- match(time.eff$cluster, order(time.eff$centers))
time.eff <- sort(time.eff$centers)

amp <- kmeans(river.coefs[,2], centers=K.amp, iter.max=100)
amp.clust <- match(amp$cluster, order(amp$centers))
amp <- sort(amp$centers)


phase.clust <- cut(river.coefs[,3], breaks=seq(-pi, pi, length=K.phase+1), 
                   include.lowest=TRUE) %>% as.numeric()
phase <- aggregate(river.coefs[,3], by=list(pc=phase.clust),
                   FUN=mean)$x

## Fix a few rivers to each cluster
n.fixed <- 15
fixed.time.clust <- fields::rdist(river.coefs[,1], time.eff) %>% 
  apply(., 2, function(x){order(x)[1:min(n.fixed, min(table(time.clust)))]})
if (!is.matrix(fixed.time.clust)) {
  fixed.time.clust <- as.matrix(fixed.time.clust) %>% t()
}

fixed.amp.clust <- fields::rdist(river.coefs[,2], amp) %>% 
  apply(., 2, function(x){order(x)[1:min(n.fixed, min(table(amp.clust))), drop=FALSE]})

if (!is.matrix(fixed.amp.clust)) {
  fixed.amp.clust <- as.matrix(fixed.amp.clust) %>% t()
}

## Define cut points based on percentage
prop <- cumsum(table(time.clust))/length(time.clust)
time.cuts <- c(-Inf, qnorm(prop, 0, 1))
time.trans.cuts <- log(diff(time.cuts)[2:(K.time-1)])
prop <- cumsum(table(amp.clust))/length(amp.clust)
amp.cuts <- c(-Inf, qnorm(prop, 0, 1))
amp.trans.cuts <- log(diff(amp.cuts)[2:(K.time-1)])
phase.cuts <- seq(-pi, pi, length=K.phase+1)

## Initiate Matrices in each river
init.regression <- function(rn){
  ## Set data.frame
  dframe <- river.list[[rn]]
  dframe$phase <- phase[phase.clust[rn]]
  model <- log(NO3)~time+cos(2*pi*time-phase)
  
  ## Set up X-matrix for this regression
  X <- model.matrix(model, data=dframe)
  
  ## Choose AR1 parameter and set up correlation matrix
  ## Choose AR1 parameter and set up correlation matrix
  res <- matrix(resid(lm(model, data=dframe)), nrow=1)
  lwr <- 1/(max(fields::rdist(river.list[[rn]]$time)) %>%
              fields::Matern.cor.to.range(., nu=0.5, cor.target=0.95))
  upr <- 1/(fields::rdist(river.list[[rn]]$time) %>%
              apply(., 2, function(x){min(x[x>0])}) %>%
              min() %>% fields::Matern.cor.to.range(., nu=0.5, cor.target=0.05))
  phi <- optim(0, fn=function(p){
    R <- exp(-p*fields::rdist(river.list[[rn]]$time))
    sig2.est <- as.numeric(res%*%chol2inv(chol(R))%*%t(res)/length(res))
    neg.llike <- -1*(fast.dmvnorm(res, rep(0, ncol(res)), sig2.est*R, log=TRUE))
    return(neg.llike)
  }, method="L-BFGS-B", lower=lwr, upper=upr)$par
  R <- exp(-phi*fields::rdist(river.list[[rn]]$time))
  Rchol <- t(chol(R))
  Rinv <- chol2inv(chol(R))
  sumR <- sum(R)
  Rinv_logdet <- determinant(Rinv, logarithm=TRUE)
  
  ## Initialize station intercept
  intcpt <- with(dframe, {
    mean(log(NO3)-time*time.eff[time.clust[rn]] -
           amp[amp.clust[rn]]*cos(2*pi*time-phase[phase.clust[rn]]))
  })
  
  ## Initialize variance parameter
  # sig2 <- sigma(gls.fit)^2
  sig2 <- sigma(lm(model, data=dframe))
  
  ## Return list object
  return(list(df=dframe %>% dplyr::select(-phase),
              X=X,
              Rinv=Rinv,
              Rchol=Rchol,
              Rinv_logdet=Rinv_logdet,
              sig2=sig2,
              intcpt=intcpt,
              rn=rn,
              sumR=sumR))
  
}
river.list <- lapply(1:length(river.list), init.regression)




## Regression matrices for cluster probs
X.clust <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(Elevation=mean(elevation+2), Department=unique(LBDEPARTEM), Area=mean(area),
            bioGeoRegion=unique(bioGeoRegion), IDPR=mean(IDPR),
            annual_specific_runoff=mean(annual_specific_runoff),
            p_sedimentay=mean(p_sedimentay),
            p_waterbody=mean(p_waterbody),
            p_wetland=mean(100 - p_waterbody - p_not_wetland)) %>% # mutate(., bioGeoRegion=as.factor(bioGeoRegion)) %>% mutate(., bioGeoRegion=relevel(bioGeoRegion, ref=2)) %>%
  model.matrix(~ # scale(log(Area)) + 
                 scale(log(Elevation)) + Department + bioGeoRegion + scale(IDPR) +
                 scale(annual_specific_runoff) + 
                 # scale(log(annual_specific_runoff)) +
                 scale(p_sedimentay) +
                 scale(log(p_waterbody+0.01)) + scale(log(p_wetland+0.01)), data=.)
X.land <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(p_agricole_tot=mean(p_agricole_tot), p_forest=mean(p_forest),
            p_urban=mean(p_urban), p_other_land_uses=mean(p_other_land_uses)) %>%
  dplyr::select(-CDSTATIONM) %>% scale(.)
V.land <- eigen(cor(X.land))$vectors[,1:(ncol(X.land)-1)]
X.clust <- cbind(X.clust, X.land%*%V.land)
colnames(X.clust)[(ncol(X.clust)-2):ncol(X.clust)] <- paste0("V", 1:3)
bta.post.var <- (t(X.clust)%*%X.clust+(1/10)*diag(ncol(X.clust))) %>%
  chol(.) %>% chol2inv(.)
bta.post.var.chol <- chol(bta.post.var) %>% t(.)
tclust.beta <- aclust.beta <-
  pclust.beta.x <- pclust.beta.y <- rep(0, ncol(X.clust))

print(paste("dim of X.clust:", paste(dim(X.clust), collapse=", ")))

