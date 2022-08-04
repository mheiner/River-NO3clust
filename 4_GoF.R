#### Goodness of fit and posterior exploration
#### also includes code to generate Figures 4-7

rm(list=ls())

## Libraries
library(tidyverse)
library("RColorBrewer")

## Load in the Posterior Draws
load("results/Results_anchor2_K5-7-12_KMtypesupplied_22224.RData"); rep_id = 4; modK="5-7-12" # USED in paper (after two rounds of chain continuation)

## may combine multiple chains
allrdafiles = list.files("results/")
indx_files_use = grep("anchor2_K5-7-12_KMtypesupplied_2222\\d", allrdafiles)
(rdafiles_use = allrdafiles[indx_files_use])
(nchains = length(rdafiles_use))
source("scripts/0_combine_chains.R")
rep_id = 1234


(Ks = as.numeric(unlist(strsplit(modK, "-"))))
(K.time = Ks[1])
(K.amp = Ks[2])
(K.phase = Ks[3])
(n_anchor = n.fixed)



## Read in the Original Data
source("scripts/0_readData.R")
summary(rivers$time)


## get initial values and fixed quantities (e.g., error correlation matrix, cluster anchors)
# source("scripts/0_init_trend0fix.R")


aclust.probs <- apply(aclust.draws, 2, function(x){
  table(factor(x, levels=1:ncol(amp.draws))) %>% prop.table()
}) %>% t()
colnames(aclust.probs) <- paste0("AP", 1:ncol(amp.draws))
tclust.probs <- apply(tclust.draws, 2, function(x){
  table(factor(x, levels=1:ncol(time.draws))) %>% prop.table()
}) %>% t()
colnames(tclust.probs) <- paste0("TP", 1:ncol(time.draws))
pclust.probs <- apply(pclust.draws, 2, function(x){
  table(factor(x, levels=1:ncol(phase.draws))) %>% prop.table()
}) %>% t()
colnames(pclust.probs) <- paste0("PP", 1:ncol(phase.draws))
probs.df <- data.frame(tclust.probs, aclust.probs, pclust.probs)

tclass.probs = data.frame(class=apply(tclust.probs, 1, which.max), maxp=apply(tclust.probs, 1, max))
head(tclass.probs)
boxplot(maxp ~ class, data=tclass.probs, ylim=c(0,1))

aclass.probs = data.frame(class=apply(aclust.probs, 1, which.max), maxp=apply(aclust.probs, 1, max))
head(aclass.probs)
boxplot(maxp ~ class, data=aclass.probs, ylim=c(0,1))

pclass.probs = data.frame(class=apply(pclust.probs, 1, which.max), maxp=apply(pclust.probs, 1, max))
head(pclass.probs)
boxplot(maxp ~ class, data=pclass.probs, ylim=c(0,1))


dim(int.draws)
(N = ncol(int.draws))

N == length(river.list)

(niter = nrow(int.draws))
dim(time.draws)
dim(amp.draws)
dim(phase.draws)
dim(tclust.draws)
unique(tclust.draws[1,])

dim(aclust.draws)
unique(aclust.draws[1,])

dim(pclust.draws)
unique(pclust.draws[1,])

## regulatory thresholds for NO3 (references in paper)
lthreshold_WHO = log(50)
lthreshold_eutroU = log(19)
lthreshold_eutroL = log(4.3)

ntt = 200
tt = seq(0.0, 6.0, length=ntt)

get_lNO3mean <- function(rn) {
  n = length(river.list[[rn]]$df$NO3)
  niter = nrow(int.draws)
  
  mean_at_obs = matrix( rep(int.draws[,rn], times=n), ncol=n ) + 
    tcrossprod(time.draws[ cbind(1:niter, tclust.draws[,rn]) ], river.list[[rn]]$df$time) +
    (matrix( rep(amp.draws[ cbind(1:niter, aclust.draws[,rn])], times=n), ncol=n ) * 
       cos(t(outer(2*pi*river.list[[rn]]$df$time, phase.draws[ cbind(1:niter, pclust.draws[,rn]) ], FUN="-")))
    )
  
  # mean_at_obs_i95 = apply(mean_at_obs, 2, function(x) quantile(x, c(0.025, 0.975)))
  
  e = matrix(rep(log(river.list[[rn]]$df$NO3), each=niter), nrow=niter) - mean_at_obs # niter by n matrix
  e_mean = colMeans(e)
  
  e_whiten = forwardsolve(river.list[[rn]]$Rchol, t(e)) %>% t() # niter by n matrix
  ew_mean = colMeans(e_whiten)
  
  mean_at_grid = matrix( rep(int.draws[,rn], times=ntt), ncol=ntt ) + 
    tcrossprod(time.draws[ cbind(1:niter, tclust.draws[,rn]) ], tt) +
    (matrix( rep(amp.draws[ cbind(1:niter, aclust.draws[,rn])], times=ntt), ncol=ntt ) * 
       cos(t(outer(2*pi*tt, phase.draws[ cbind(1:niter, pclust.draws[,rn]) ], FUN="-")))
    )
  
  # mean_at_grid_i95 = apply(mean_at_grid, 2, function(x) quantile(x, c(0.025, 0.975)))
  # R2 = 1.0 - mean(sig2.draws[,rn]) / var(log(river.list[[rn]]$NO3))
  R2 = 1.0 - mean(rowSums(e_whiten^2)) / sum( (log(river.list[[rn]]$df$NO3) - mean(log(river.list[[rn]]$df$NO3)))^2 )
  
  start_exceedsU = int.draws[,rn] > lthreshold_eutroU
  start_belowL = int.draws[,rn] < lthreshold_eutroL
  
  end_exceedsU = (int.draws[,rn] + 6.0*time.draws[ cbind(1:niter, tclust.draws[,rn]) ]) > lthreshold_eutroU
  end_belowL = (int.draws[,rn] + 6.0*time.draws[ cbind(1:niter, tclust.draws[,rn]) ]) < lthreshold_eutroL
  
  end_peak_exceedsU = (int.draws[,rn] + 6.0*time.draws[ cbind(1:niter, tclust.draws[,rn]) ] + amp.draws[ cbind(1:niter, aclust.draws[,rn]) ]) > lthreshold_eutroU
  end_peak_belowL = (int.draws[,rn] + 6.0*time.draws[ cbind(1:niter, tclust.draws[,rn]) ] + amp.draws[ cbind(1:niter, aclust.draws[,rn]) ]) < lthreshold_eutroL
  
  return(list(mean_at_obs=colMeans(mean_at_obs), #mean_at_obs_i95=mean_at_obs_i95, 
              # e=e, 
              e_mean=e_mean, 
              # ew=e_whiten, 
              ew_mean=ew_mean,
              R2=R2,
              mean_at_grid=colMeans(mean_at_grid), #, mean_at_grid_i95=mean_at_grid_i95
              p_start_exceedsU = mean(start_exceedsU),
              p_start_belowL = mean(start_belowL),
              p_end_exceedsU = mean(end_exceedsU),
              p_end_belowL = mean(end_belowL),
              p_end_peak_exceedsU = mean(end_peak_exceedsU),
              p_end_peak_belowL = mean(end_peak_belowL)
  ))
}

rivermean.list <- lapply(1:length(river.list), get_lNO3mean)


## kurtosis
library("e1071")
kurt = sapply(rivermean.list, function(x) kurtosis(x$e_mean))
hist(kurt, breaks=50)
mean(kurt > 6) # level of excess kurtosis of student's t with 5 df
mean(kurt > 3) # level of excess kurtosis of Laplace distribution
hist(kurt[kurt<6], breaks=50)



## serial correlation in residuals
library("EnvStats") # for testing serial correlation in residuals w/ serialCorrelationTest()

rn = 1
rivermean.list[[rn]]$e_mean
dim(river.list[[rn]]$Rinv)
aa = acf(rivermean.list[[rn]]$e_mean)
aaw = acf(rivermean.list[[rn]]$ew_mean)
serialCorrelationTest(rivermean.list[[rn]]$e_mean)$p.value

for (rn in 1:length(river.list)) {
  rivermean.list[[rn]]$e_mean_test = serialCorrelationTest( rivermean.list[[rn]]$e_mean / sqrt(mean(sig2.draws[,rn])), test="rank.von.Neumann")
  rivermean.list[[rn]]$ew_mean_test = serialCorrelationTest( rivermean.list[[rn]]$ew_mean / sqrt(mean(sig2.draws[,rn])), test="rank.von.Neumann")
  # rivermean.list[[rn]]$R2w = 1.0 - sum(forwardsolve(river.list[[rn]]$Rchol, rivermean.list[[rn]]$e_mean)^2) / sum( (log(river.list[[rn]]$df$NO3) - mean(log(river.list[[rn]]$df$NO3)))^2 )
  cat(rn, " of ", N, "\r")
}

e_tests = sapply(1:N, function(rn) c(rivermean.list[[rn]]$e_mean_test$p.value, rivermean.list[[rn]]$ew_mean_test$p.value)) %>% t()
dim(e_tests)
head(e_tests)
mean(e_tests[,1] < 0.05)
mean(e_tests[,2] < 0.05)

e_acf = sapply(1:N, function(rn) c(max(rivermean.list[[rn]]$e_mean_test$estimate), max(rivermean.list[[rn]]$ew_mean_test$estimate))) %>% t()
hist(e_acf[,2] / e_acf[,1])
hist(apply(abs(e_acf), 1, diff))
hist(e_acf[,1]) # original mean residuals
hist(e_acf[,2]) # whitened mean residuals

rn = 4
(rn = which.max(apply(abs(e_acf), 1, diff))) # also try which.max :/
e_acf[rn,]
plot.ts(cbind(rivermean.list[[rn]]$e_mean, rivermean.list[[rn]]$e_mean_w))
acf(rivermean.list[[rn]]$e_mean)
acf(rivermean.list[[rn]]$ew_mean)




## R squared
dim(sig2.draws)
sig2_pm = colMeans(sig2.draws)
hist(sig2_pm)
hist(sig2.draws[,2150], breaks=50)

length(river.list) == length(sig2_pm)
length(river.list) == N
R2 = sapply(rivermean.list, function(x) x$R2)
R2_meanfit = sapply(1:N, function(rn) 1.0 - sum(rivermean.list[[rn]]$e_mean^2) / sum( (log(river.list[[rn]]$df$NO3) - mean(log(river.list[[rn]]$df$NO3)))^2 ))

hist(R2[R2>0], xlim=c(0,1))
hist(R2_meanfit[R2_meanfit>0], xlim=c(0,1))

sum(R2 < 0)
mean(R2 < 0)
mean(R2 < .2)

## optionally save and load for comparison
# R2_57 = sapply(rivermean.list, function(x) x$R2); save(R2_57, file="results/R2_57.rda")
# R2_711 = sapply(rivermean.list, function(x) x$R2); save(R2_711, file="results/R2_711.rda")
# load("results/R2_57.rda"); load("results/R2_711.rda")

R2rat75 = R2_711 / R2_57
hist(R2rat75)
summary(R2rat75)
summary(R2rat75 * 100 - 100)

mean( abs(R2rat75 - 1.0) > 0.8)
mean( abs(R2rat75 - 1.0) < 0.1)

hist(R2rat75[abs(R2rat75 - 1.0) < 1.5]*100 - 100, xlim=c(-200,200), breaks=50)

pdf(file="plots/R2_711to57.pdf", height=4, width=5)
hist(R2rat75[abs(R2rat75 - 1.0) < 0.5]*100 - 100, 
     xlim=c(-30, 30), breaks=50, main="7-11 from 5-7", xlab="Percent change", border=FALSE)
dev.off()






R2 = R2_57
R2 = sapply(rivermean.list, function(x) x$R2)


pdf(file=paste0("plots/hist_R2_", modK, "_anchor", n_anchor, "_seed", seed[1], ".pdf"), height=3.5, width=4)
hist(R2, breaks=30)
hist(R2[R2>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="", breaks=20)
hist(R2[R2<0.0])
dev.off()

pdf(file=paste0("plots/hist_R2_by_mod.pdf"), height=4.5, width=5)
hist(R2_57[R2_57>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="5-7", breaks=20)
hist(R2_711[R2_711>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="7-11")
dev.off()




indx_negR2 = which(R2 < 0)
length(indx_negR2)
ordR2 = order(R2)
R2[ordR2[1]]
R2[ordR2[N]]



hist(log(rivers$NO3))
nobs = sapply(river.list, function(alist) nrow(alist$df)) # number of obervations at each station
head(nobs)
which.max(nobs)


plot(tclass.probs$maxp, R2)
summary(lm(R2 ~ tclass.probs$maxp))

plot(aclass.probs$maxp, R2)
summary(lm(R2 ~ aclass.probs$maxp))

plot(pclass.probs$maxp, R2)
summary(lm(R2 ~ pclass.probs$maxp))




### Fit by cluster
## fix settings interactively
k = 1
head(probs.df)

rindx = which(probs.df[,paste0("TP",k)] > 0.9)
rindx = which(probs.df[,paste0("AP",k)] > 0.9)
rindx = which(probs.df[,paste0("PP",k)] > 0.9)

length(rindx)

hist(R2[rindx])

(rn = rindx[1])

# shorten rindx if necessary
rindx = rindx[1:20]
rindx = sample(rindx, 50)

mesg_plt = "time"
mesg_plt = "amp"
mesg_plt = "phase"

pdf(file=paste0("plots/fit_by_clust_", modK, "_anchor", n_anchor, "_7yr_seed_", seed[1], "_", mesg_plt, k, ".pdf"), width=5, height=4)
hist(R2[rindx])
for (rn in rindx) {
  plot(log(NO3) ~ time, data=river.list[[rn]]$df, 
       main=paste0(river.list[[rn]]$df$CDSTATIONM[1], "\n", river.list[[rn]]$df$LBDEPARTEM[1]), 
       xlim=c(0,6), xlab="year", type="o"
       # ylim=range(log(rivers$NO3))
       # ylim=c(-1, 5)
  )
  lines(tt,  rivermean.list[[rn]]$mean_at_grid, col="red", lwd=1.5)
  mtext(side=3, line=1, at=5.5, paste("R-squared", round(rivermean.list[[rn]]$R2, 2)))
}
dev.off()



rn = 150 # select a river network (or a sequence, below)
rn = 4230
rn = 1150

rn = 4514
pp = 0.90
(rn = ordR2[1]); R2[ordR2[1]]
(rn = ordR2[1 + ii]); R2[ordR2[1 + ii]]
(rn = ordR2[floor(N/2)]); R2[ordR2[floor(N/2)]]
(rn = ordR2[N - ii]); R2[ordR2[N - ii]]
(rn = ordR2[N]); R2[ordR2[N]]
(rn = ordR2[floor(N*pp)]); R2[ordR2[floor(N*pp)]]


pp_seq = c(1e-5, 0.25, 0.5, 0.75, 1.0)
pp_seq = sort(c(pp_seq, c(0.01, 0.26, 0.51, 0.76, 0.99)))
(rn_seq = ordR2[ceiling(N*pp_seq)]); R2[rn_seq]

for( i in 1:length(rn_seq) ) {
  pdf(file=paste0("plots/fit_resid_model_", modK, "_anchor", n_anchor, "_seed", seed[1], "_", sprintf("%03d", round(100*pp_seq[i])), "ile_id", sprintf("%04d.pdf", rn_seq[i])), 
      width=6, height=3.5)
  
  plot(log(NO3) ~ time, data=river.list[[rn_seq[i]]]$df, 
       main=paste0(river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1]), 
       xlim=c(0,6), xlab="year",
       # ylim=range(log(rivers$NO3))
       ylim=c(-1, 5)
  )
  lines(log(NO3) ~ time, data=river.list[[rn_seq[i]]]$df, type="o")
  lines(tt,  rivermean.list[[rn_seq[i]]]$mean_at_grid, col="red", lwd=1.5)
  mtext(side=3, line=1, at=5.5, paste("R-squared", round(rivermean.list[[rn_seq[i]]]$R2, 2)))
  
  acf(rivermean.list[[rn_seq[i]]]$ew_mean, 
      main=paste0("posterior mean residual, ", river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1])) 
  pacf(rivermean.list[[rn_seq[i]]]$ew_mean, 
       main=paste0("posterior mean residual, ", river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1]))
  
  par(mfrow=c(1,2))
  hist(rivermean.list[[rn_seq[i]]]$ew_mean, xlab="e",
       main=paste("Residuals:", river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1]))
  qqnorm(rivermean.list[[rn_seq[i]]]$ew_mean); qqline(rivermean.list[[rn_seq[i]]]$ew_mean)
  
  dev.off()
  cat(i, "\r")
}






### Maps with estimates (Figure 5)

R2df = data.frame(CDSTATIONM=sapply(river.list, function(x) x$df$CDSTATIONM[1]),
                  lon=sapply(river.list, function(x) x$df$lon[1]),
                  lat=sapply(river.list, function(x) x$df$lat[1]),
                  R2=sapply(rivermean.list, function(x) x$R2),
                  intcpt=colMeans(int.draws))
R2df %>% arrange(desc(R2)) %>% head(n=10)

france = map_data("world", region = "France")
france.plot = ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)
locs <- rivers %>% group_by(CDSTATIONM) %>% 
  summarise(lat=unique(lat)[1], lon=unique(lon)[1], alt=unique(elevation)[1])
locs$tclass <- apply(tclust.draws, 2, function(x){
  xtab <- table(factor(x, levels=1:ncol(time.draws)))
  return(sort(colMeans(time.draws))[which.max(xtab)])
})
locs$aclass <- apply(aclust.draws, 2, function(x){
  xtab <- table(factor(x, levels=1:ncol(amp.draws)))
  return(sort(colMeans(amp.draws))[which.max(xtab)])
})
locs$pclass <- apply(pclust.draws, 2, function(x){
  xtab <- table(factor(x, levels=1:ncol(phase.draws)))
  return(sort(colMeans(phase.draws))[which.max(xtab)])
})

int.pm = colMeans(int.draws)
time.pm = sapply(1:N, function(i) time.draws[cbind(1:niter, tclust.draws[,i])]) %>% colMeans
amp.pm = sapply(1:N, function(i) amp.draws[cbind(1:niter, aclust.draws[,i])]) %>% colMeans
phase.pm =colMeans(phase.draws)[apply(pclust.probs, 1, which.max)]

locs$int.pm = int.pm
locs$time.pm = time.pm
locs$amp.pm = amp.pm
locs$phase.pm = phase.pm

R2df = left_join(R2df, locs[,-grep("lat|lon", colnames(locs))], by="CDSTATIONM", keep=FALSE)
head(R2df)

all(R2df$intcpt == R2df$int.pm)



## Map setup
library("RColorBrewer")
france.plot.dep = ggplot() + geom_polygon(data=map_data("france"), fill="white", colour = "gray40", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)
(rdylbu7 = brewer.pal(7, "RdYlBu"))
(blylrd7 = rev(brewer.pal(7, "RdYlBu")))
rdbu11 = brewer.pal(11, "RdBu")
(burd11 = rev(rdbu11))

france.plot.clean = france.plot + 
  coord_cartesian( xlim = c(-5, 7.8), ylim=c(42, 52)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(-0.03,0,-0.01,0), "null")
  )

france.plot.dep.clean = france.plot.dep + 
  coord_cartesian( xlim = c(-5, 7.8), ylim=c(42, 52) ) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(-0.03,0,-0.01,0), "null")
  )

col_lims = list()
col_lims[["5-7-12"]] = list(int=c(-1.2, 4.8), time=c(-0.092, 0.092), amp=c(0.0, 1.35))
col_lims[["7-11-12"]] = list(int=c(-1.2, 4.9), time=c(-0.13, 0.13), amp=c(0.0, 1.65))

cpt = seq(col_lims[[modK]][["time"]][1], col_lims[[modK]][["time"]][2], length=12) # 12 for the number of colors on pallette

if(K.time == 7) { cpt[7] = 0.0110; cpt[8] = 0.025}

## factor variables for trend and amplitude class (not included in paper)
R2df$tf = cut(R2df$tclass, breaks = cpt)
cpt
sort(unique(as.numeric(R2df$tf)))
sort(unique(R2df$tclass))
levels(R2df$tf)[sort(unique(as.numeric(R2df$tf)))] <- paste(round(sort(unique(R2df$tclass)), 3))

levels(R2df$tf)

cpa = c(0.0, 0.03, 0.09, 0.15, 0.22, 0.33, 0.44, 0.55, 0.80, 1.0, 1.5,  1.75)

R2df$af = cut(R2df$aclass, breaks=cpa)
cpa
sort(unique(R2df$aclass))
sort(unique(as.numeric(R2df$af)))
R2df$af = cut(R2df$aclass, breaks = cpa)
levels(R2df$af)[sort(unique(as.numeric(R2df$af)))] <- paste(round(sort(unique(R2df$aclass)), 3))
levels(R2df$af)



## maps with R-squared and parameter estimates
pdf(file=paste0("plots/map_R2_", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), width=5, height=4)

plt = france.plot.clean + geom_point(data=R2df, 
                                     aes(x=lon, y=lat, color=R2), size=0.75) +
  scale_color_distiller(palette="Spectral")
print(plt)

plt = france.plot.clean + geom_point(data=R2df %>% filter(R2 >= 0.0), 
                                     aes(x=lon, y=lat, color=R2), size=0.75) +
  scale_color_distiller(palette="Spectral")
print(plt)

plt = france.plot.clean + geom_point(data=R2df %>% filter(R2 < 0.0), 
                                     aes(x=lon, y=lat, color=R2), size=0.75) +
  scale_color_distiller(palette="Spectral")
print(plt)



plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=(round(intcpt,3))), size=.75) +
  labs(color="Intercept") + 
  scale_color_gradient2(limits=col_lims[[modK]][["int"]], midpoint=2.3, low=blylrd7[1], mid=blylrd7[4], high=blylrd7[7])
print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=(round(intcpt,3))), size=.75) +
  labs(color="Intercept:\nNO3 mg/L") + 
  scale_color_gradient2(limits=col_lims[[modK]][["int"]], midpoint=2.3, 
                        low=blylrd7[1], mid=blylrd7[4], high=blylrd7[7],
                        breaks=log(c(1,2,5,10,20,50)), labels=c(1,2,5,10,20,50))
print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=tf), size=.75) +
  labs(color="Trend") + 
  scale_color_manual(values=burd11[sort(unique(as.numeric(R2df$tf)))]) + guides(colour = guide_legend(reverse=T))
print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=(round(time.pm,3)), alpha=pmax(R2,0)), size=.75) +
  labs(color="Trend", alpha="R-squared") + 
  scale_color_gradient2(limits=col_lims[[modK]][["time"]], midpoint=0, low=burd11[3], mid=burd11[6], high=burd11[9])
print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=(round(100*time.pm,3)), alpha=pmax(R2,0)), size=.75) +
  labs(color="Trend:\n% annual\nchange\nNO3", alpha="R-squared") + 
  scale_color_gradient2(limits=100*col_lims[[modK]][["time"]], midpoint=0, low=burd11[3], mid=burd11[6], high=burd11[9])
print(plt)


library("scales")
gnor11 = seq_gradient_pal("forestgreen", "darkorange1")(seq(0,1, length=11))

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=af, alpha=pmax(R2,0)), size=.75) +
  labs(color="Amplitude", alpha="R-squared") + 
  scale_color_manual(values=gnor11[sort(unique(as.numeric(R2df$af)))]) + guides(colour = guide_legend(reverse=T), alpha = guide_legend(reverse=T))
print(plt)

plt = france.plot.dep.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=af, alpha=pmax(R2,0)), size=.75) +
  labs(color="Amplitude", alpha="R-squared") + 
  scale_color_manual(values=gnor11[sort(unique(as.numeric(R2df$af)))]) + guides(colour = guide_legend(reverse=T), alpha = guide_legend(reverse=T))
print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=amp.pm, alpha=pmax(R2,0)), size=.75) +
  labs(color="Amplitude", alpha="R-squared") + 
  scale_color_gradient2(limits=c(0.0, 1.5), midpoint=0.75, 
                        low=gnor11[1], mid=gnor11[6], high=gnor11[11])

print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=2*amp.pm, alpha=pmax(R2,0)), size=.75) +
  labs(color="Amplitude:\nhigh / low\nratio\nNO3", alpha="R-squared") + 
  scale_color_gradient2(limits=c(0.0, 2*1.5), midpoint=2*0.75, 
                        low=gnor11[1], mid=gnor11[6], high=gnor11[11],
                        breaks=log(c(1,2,5,10,20)), labels=c(1,2,5,10,20))
print(plt)


plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=round(pclass,3), alpha=pmax(R2,0)), size=.75)  +
  labs(color="Phase", alpha="R-squared") + 
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10)) + guides(alpha = guide_legend(reverse=T))

print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=round(pclass,3), alpha=pmax(R2,0)), size=.75)  +
  labs(color="Phase:\npeak NO3", alpha="R-squared") + 
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10),
                        breaks=seq(-pi, pi, length=13), 
                        labels=c("Jul", "", "", "Oct", "", "", "Jan", "", "", "Apr", "", "", "Jul")) + 
  guides(alpha = guide_legend(reverse=T))
print(plt)


plt = france.plot.dep.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=round(pclass,3), alpha=pmax(R2,0)), size=.75)  +
  labs(color="Phase", alpha="R-squared") + 
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10)) + guides(alpha = guide_legend(reverse=T))

print(plt)


plt = france.plot.clean + geom_point(data=cbind(R2df, probs.df[,c('TP1', 'TP2')]),
                                     aes(x=lon, y=lat, color=TP1+TP2, alpha=pmax(0,R2)), size=0.75) +
  labs(alpha="R-squared") +
  scale_color_distiller(palette="Spectral", limits=c(0,1))
print(plt)

plt = france.plot.clean + geom_point(data=cbind(R2df, probs.df['TP3']),
                                     aes(x=lon, y=lat, color=TP3, alpha=pmax(0,R2)), size=0.75) +
  labs(alpha="R-squared") +
  scale_color_distiller(palette="Spectral", limits=c(0,1))
print(plt)

plt = france.plot.clean + geom_point(data=cbind(R2df, probs.df[,c('TP4', 'TP5')]),
                                     aes(x=lon, y=lat, color=TP4+TP5, alpha=pmax(0,R2)), size=0.75) +
  labs(alpha="R-squared") +
  scale_color_distiller(palette="Spectral", limits=c(0,1))
print(plt)

dev.off()








## R squared vs estimates
plot(int.pm, R2)
plot(int.pm[R2>0], R2[R2>0], xlab="log(NO3)", ylab="R squared", main="Intercept\nby location", ylim=c(0,1))
plot(exp(int.pm[R2>0]), R2[R2>0], xlab="NO3", ylab="R squared", main="Intercept\nby location", ylim=c(0,1))
summary(lm(R2[R2>0] ~ int.pm[R2>0]))

plot(time.pm, R2)
plot(time.pm[R2>0], R2[R2>0], xlab="trend", ylab="R squared", main="Time trend\nby location", ylim=c(0,1))
plot(100*exp(time.pm[R2>0]) - 100, R2[R2>0], xlab="Percent annual change in peak NO3", ylab="R squared", main="Time trend\nby location", ylim=c(0,1))
summary(lm(R2[R2>0] ~ time.pm[R2>0]))

plot(amp.pm, R2)
plot(amp.pm[R2>0], R2[R2>0], xlab="amp", ylab="R squared", main="Aplitude\nby location", ylim=c(0,1)) # not surprising
plot(log(amp.pm[R2>0]), R2[R2>0], xlab="log(amp)", ylab="R squared", main="Aplitude\nby location", ylim=c(0,1)) # not surprising
plot(exp(2*amp.pm[R2>0]), R2[R2>0], xlab="Peak-to-trough ratio, NO3", ylab="R squared", main="Aplitude\nby location", ylim=c(0,1)) # not surprising
summary(lm(R2[R2>0] ~ amp.pm[R2>0]))

plot(phase.pm, R2)
boxplot(R2 ~ phase.pm, axes=F)
boxplot(R2[R2>0] ~ phase.pm[R2>0], axes=F, ylim=c(0,1))
axis(side=1, at=1:13, labels=c("Jul", "", "Sep", "", "Nov", "", "Jan", "", "Mar", "", "May", "", "July")); axis(side=2)
# axis(side=1, labels=seq(-5, 6, by=1), at=1:12); axis(side=2)
plot(phase.pm[R2>0], R2[R2>0], xlab="phase", ylab="R squared", main="Phase\nby location", ylim=c(0,1), xlim=c(-pi, pi))
plot(phase.pm[R2>0], R2[R2>0], xlab="Time of peak NO3", ylab="R squared", main="Phase\nby location", axes=F, ylim=c(0,1), xlim=c(-pi, pi))
axis(side=2)
axis(side=1, at=seq(-pi, pi, length=13), labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul"))

pdf(paste0("plots/Rsqared_by_estimates_", modK, "_anchor", n_anchor, "_seed", seed[1], "_7yr.pdf"), height=4, width=6)
plot(int.pm[R2>0], R2[R2>0], xlab="log(NO3)", ylab="R squared", main="Intercept\nby location", ylim=c(0,1), pch=".")
plot(100*exp(time.pm[R2>0]) - 100, R2[R2>0], xlab="Percent annual change in NO3", ylab="R squared", main="Time trend\nby location", ylim=c(0,1), pch=".")
plot(exp(2*amp.pm[R2>0]), R2[R2>0], xlab="Peak-to-trough ratio, NO3", ylab="R squared", main="Aplitude\nby location", ylim=c(0,1), pch=".")
plot(jitter(phase.pm[R2>0]), R2[R2>0], xlab="Time of peak NO3", ylab="R squared", main="Phase\nby location", axes=F, ylim=c(0,1), xlim=c(-pi, pi), pch=".")
axis(side=2)
axis(side=1, at=seq(-pi, pi, length=13), labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul"))
boxplot(R2[R2>0] ~ phase.pm[R2>0], xlab="Time of peak NO3", ylab="R squared", main="Phase\nby location", axes=F, ylim=c(0,1))
axis(side=2)
axis(side=1, at=1:13, labels=c("Jul", "", "Sep", "", "Nov", "", "Jan", "", "Mar", "", "May", "", "July")); axis(side=2)
dev.off()

## generally, rivers that have means far from the stable clusters don't fit as well (they're straddling two less-ideal fits)








### Relationships among estimates (Figure 6)

K_phase = 12
## Data frame with river summaries
pclass <- apply(pclust.draws, 2, function(x){
  xtab <- table(factor(x, levels=1:K_phase))
  # return(sort(colMeans(phase.draws))[which.max(xtab)])
  which.max(xtab)
})
pclass
table(pclass)

head(probs.df)

threshold_class = function(pend_L, pend_H) {
  if (pend_L > 0.9) {
    out = "Below"
  } else if (pend_H > 0.9) {
    out = "Above"
  } else {
    out = "Uncertain"
  }
  out
}

dat_station = lapply(1:N, function(i) data.frame(i=i, CDSTATIONM=river.list[[i]]$df$CDSTATIONM[1], 
                                                 lon=river.list[[i]]$df$lon[1], lat=river.list[[i]]$df$lat[1],
                                                 elevation=river.list[[i]]$df$elevation[1],
                                                 IDPR = river.list[[i]]$df$IDPR[1],
                                                 Agriculture = river.list[[i]]$df$p_agricole_tot[1],
                                                 Runoff = river.list[[i]]$df$annual_specific_runoff[1],
                                                 PctWater = river.list[[i]]$df$p_waterbody[1],
                                                 LBDEPARTEM=river.list[[i]]$df$LBDEPARTEM[1],
                                                 Rsquared=rivermean.list[[i]]$R2,
                                                 int=int.pm[i], trend=time.pm[i], amp=amp.pm[i], phase=phase.pm[i],
                                                 itime = 100*exp(time.pm[i]) - 100,
                                                 iamp = exp(2*amp.pm[i]),
                                                 tclass=which.max(table(factor(tclust.draws[,i], levels=1:max(tclust.draws[,i])))),
                                                 aclass=which.max(table(factor(aclust.draws[,i], levels=1:max(aclust.draws[,i])))),
                                                 pclass=which.max(table(factor(pclust.draws[,i], levels=1:max(pclust.draws[,i])))),
                                                 threshold_class = threshold_class(rivermean.list[[i]]$p_end_belowL, rivermean.list[[i]]$p_end_exceedsU),
                                                 threshold_peak_class = threshold_class(rivermean.list[[i]]$p_end_peak_belowL, rivermean.list[[i]]$p_end_peak_exceedsU) )
) %>% do.call(rbind, .)
head(dat_station)
write.csv(file=paste0("results/FranceRivers_estimates", modK, "_anchor", n_anchor, "_seed", seed[1], "7yr.csv"), dat_station, row.names=FALSE)


## interpretable time and amplitude
itime.pm = 100*exp(time.pm) - 100
iamp.pm = exp(2*amp.pm)


library("ggplot2")
library("ggExtra")
library("grid")
library("gridExtra")

### Relationships among estimates
p = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=itime, y=int, color=threshold_class)) + 
  geom_hline(yintercept=lthreshold_WHO, color="gray50", lty=2) + 
  geom_hline(yintercept=lthreshold_eutroL, color="gray50", lty=1) +
  geom_hline(yintercept=lthreshold_eutroU, color="gray50", lty=1) + 
  geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + 
  labs(x="Percent annual change in NO3", y="Intercept\nNO3 mg/L", color="Threshold, 2016") +
  scale_y_continuous(breaks=log(c(1,2,5,10,20,50)), labels=c(1,2,5,10,20,50)) +
  scale_color_manual(values=c(brewer.pal(3, "Set1")[1:2], "black")) +
  guides(colour = guide_legend(override.aes = list(size=2)))
p1 = ggMarginal(p + theme(legend.position="bottom", legend.key.height=unit(0.022,"npc")), type="histogram")
p1 = arrangeGrob(p1, top=textGrob("a)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p1)

p = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=2*amp, y=int, color=threshold_peak_class)) + 
  geom_vline(xintercept=0, color="gray50") + 
  geom_hline(yintercept=lthreshold_WHO, color="gray50", lty=2) + 
  geom_hline(yintercept=lthreshold_eutroL, color="gray50", lty=1) +
  geom_hline(yintercept=lthreshold_eutroU, color="gray50", lty=1) + 
  labs(x="High / low ratio (NO3, annual cycle)", y="Intercept\nNO3 mg/L", color="Threshold at peak, 2016") +
  scale_x_continuous(breaks=log(c(1,2,5,10,20)), labels=c(1,2,5,10,20)) +
  scale_y_continuous(breaks=log(c(1,2,5,10,20,50)), labels=c(1,2,5,10,20,50)) +
  geom_point(size=0.025) +
  scale_color_manual(values=c(brewer.pal(3, "Set1")[1:2], "black")) +
  guides(colour = guide_legend(override.aes = list(size=2)))
p2 = ggMarginal(p + theme(legend.position="bottom", legend.key.height=unit(0.022,"npc")), type="histogram")
p2 = arrangeGrob(p2, top=textGrob("b)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p2)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter((phase+pi)*13/(2*pi), factor = 1.00)), aes(x=jphase, y=int)) +
  geom_hline(yintercept=lthreshold_WHO, color="gray50", lty=2) + 
  geom_hline(yintercept=lthreshold_eutroL, color="gray50", lty=1) +
  geom_hline(yintercept=lthreshold_eutroU, color="gray50", lty=1) + 
  geom_boxplot(data=dat_station %>% filter(Rsquared > 0), aes(x=factor(phase)), outlier.size=0.5) +
  geom_point( size=0.025, aes(color=Agriculture) ) + # geom_jitter(width=0.02, inherit.aes=TRUE) + 
  labs(x="Time of peak NO3", y="Intercept\nNO3 mg/L", color="% Ag") +
  scale_x_discrete(labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul")) +
  scale_y_continuous(breaks=log(c(1,2,5,10,20,50)), labels=c(1,2,5,10,20,50)) + 
  scale_color_distiller(palette="Spectral")
p3 = ggMarginal(p + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram")

p0 = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=as.numeric(as.factor(phase)), y=int)) + geom_point( ) 
p03 = ggMarginal(p0 + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram", bins=13)
p3$grobs[p3$layout$name == "topMargPlot"] <- p03$grobs[p03$layout$name == "topMargPlot"]
p3 = arrangeGrob(p3, top=textGrob("c)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p3)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter(phase, 0.5)), aes(x=2*amp, y=itime, color=jphase)) + 
  geom_hline(yintercept=0, color="gray50") + geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + 
  labs(x="High / low ratio (NO3, annual cycle)", y="Percent annual change\nin NO3", color="Peak NO3") +
  scale_x_continuous(breaks=log(c(1,2,5,10,20)), labels=c(1,2,5,10,20)) + 
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10), 
                        breaks=seq(-pi, pi, length=13), labels=c("Jul", "", "", "Oct", "", "", "Jan", "", "", "Apr", "", "", "Jul"))
p4 = ggMarginal(p + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram")
p4 = arrangeGrob(p4, top=textGrob("d)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p4)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter((phase+pi)*13/(2*pi), factor = 1.00)), aes(x=jphase, y=itime)) +
  geom_hline(yintercept=0, color="gray50", lty=1) + 
  geom_boxplot(data=dat_station %>% filter(Rsquared > 0), aes(x=factor(phase)), outlier.size=0.5) +
  geom_point( size=0.025, aes(color=Agriculture) ) + # geom_jitter(width=0.02, inherit.aes=TRUE) + 
  labs(x="Time of peak NO3", y="Percent annual change\nin NO3", color="% Ag") +
  scale_x_discrete(labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul")) +
  scale_color_distiller(palette="Spectral")
p5 = ggMarginal(p + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram")

p0 = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=as.numeric(as.factor(phase)), y=itime)) + geom_point( ) 
p05 = ggMarginal(p0 + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram", bins=13)
p5$grobs[p5$layout$name == "topMargPlot"] <- p05$grobs[p05$layout$name == "topMargPlot"]
p5 = arrangeGrob(p5, top=textGrob("e)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p5)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter((phase+pi)*13/(2*pi), factor = 1.00)), aes(x=jphase, y=2*amp)) +
  geom_hline(yintercept=0, color="gray50", lty=2) + 
  geom_boxplot(data=dat_station %>% filter(Rsquared > 0), aes(x=factor(phase)), outlier.size=0.5) +
  geom_point( size=0.025, aes(color=Agriculture) ) + # geom_jitter(width=0.02, inherit.aes=TRUE) + 
  labs(x="Time of peak NO3", y="High / low ratio\n(NO3, annual cycle)", color="% Ag") +
  scale_x_discrete(labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul")) +
  scale_y_continuous(breaks=log(c(1,2,5,10,20)), labels=c(1,2,5,10,20)) + 
  scale_color_distiller(palette="Spectral")
p6 = ggMarginal(p + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram")

p0 = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=as.numeric(as.factor(phase)), y=2*amp)) + geom_point( ) 
p06 = ggMarginal(p0 + theme(legend.position="bottom", legend.key.height=unit(0.008,"npc")), type="histogram", bins=13)
p6$grobs[p6$layout$name == "topMargPlot"] <- p06$grobs[p06$layout$name == "topMargPlot"]
p6 = arrangeGrob(p6, top=textGrob("f)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p6)

pdf(paste0("plots/estimate_relationships_hist_", modK, "_anchor", n_anchor, "_seed", seed[1], "_color.pdf"), width=9, height=9.5)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix=matrix(1:6, ncol=2))
dev.off()


## statistics to go with the preceding plots
summary(dat_station$elevation)
hist(dat_station$elevation[dat_station$elevation < 500])
dat_station %>% filter(Rsquared > 0) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, elevation < median(elevation)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, elevation > median(elevation)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")

dat_station %>% filter(Rsquared > 0) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, elevation < quantile(elevation, .25)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, elevation > quantile(elevation, .75)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")

dat_station %>% filter(Rsquared > 0) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, IDPR < quantile(IDPR, .2)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, IDPR > quantile(IDPR, .8)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")

dat_station %>% filter(Rsquared > 0) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, Agriculture < quantile(Agriculture, .25)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")
dat_station %>% filter(Rsquared > 0, Agriculture > quantile(Agriculture, .75)) %>% select(c("int", "trend", "amp")) %>% cor(., method="spearman")

## measure associations across posterior samples
library("polycor")

dim(aclust.draws)
aclust.draws[1:10,1:5]

ii = 2
table(tclust.draws[ii,], aclust.draws[ii,])

rho_time_int = numeric(n.draws)
rho_time_amp = numeric(n.draws)
rho_amp_int = numeric(n.draws)

rholow_time_int = numeric(n.draws)
rholow_time_amp = numeric(n.draws)
rholow_amp_int = numeric(n.draws)

lowindx = which(dat_station$elevation < median(dat_station$elevation))

for (ii in 1:n.draws) {
  rho_time_amp[ii] = polychor(factor(tclust.draws[ii,], levels=1:K.time), factor(aclust.draws[ii,], levels=1:K.amp))
  rho_time_int[ii] = polyserial(int.draws[ii,], factor(tclust.draws[ii,], levels=1:K.time))
  rho_amp_int[ii] = polyserial(int.draws[ii,], factor(aclust.draws[ii,], levels=1:K.amp))
  
  rholow_time_amp[ii] = polychor(factor(tclust.draws[ii, lowindx], levels=1:K.time), factor(aclust.draws[ii, lowindx], levels=1:K.amp))
  rholow_time_int[ii] = polyserial(int.draws[ii, lowindx], factor(tclust.draws[ii, lowindx], levels=1:K.time))
  rholow_amp_int[ii] = polyserial(int.draws[ii, lowindx], factor(aclust.draws[ii, lowindx], levels=1:K.amp))
  
  cat(ii, "of", n.draws, "\r")
}

cor(int.pm, time.pm)
hist(rho_time_int); mean(rho_time_int); quantile(rho_time_int, c(0.025, 0.975))
hist(rholow_time_int); mean(rholow_time_int); quantile(rholow_time_int, c(0.025, 0.975))

cor(int.pm, amp.pm)
hist(rho_amp_int); mean(rho_amp_int); quantile(rho_amp_int, c(0.025, 0.95))
hist(rholow_amp_int); mean(rholow_amp_int); quantile(rholow_amp_int, c(0.025, 0.975))

cor(time.pm, amp.pm)
hist(rho_time_amp); mean(rho_time_amp); quantile(rho_time_amp, c(0.025, 0.975))
hist(rholow_time_amp); mean(rholow_time_amp); quantile(rholow_time_amp, c(0.025, 0.975))


library("lattice")
table(dat_station$tclass, dat_station$aclass)
prop.table(table(dat_station$tclass, dat_station$aclass), margin = 2)
levelplot(table(dat_station$tclass, dat_station$aclass))
chisq.test(dat_station$tclass, dat_station$aclass, simulate.p.value=TRUE)

table(dat_station$tclass, dat_station$pclass)
levelplot(table(dat_station$tclass, dat_station$pclass))
chisq.test(dat_station$tclass, dat_station$pclass, simulate.p.value=TRUE)

dat_station %>% filter(Rsquared > 0) %>% select(c("int", "trend", "amp")) %>% cor()






### Covariates
## Added-variable plots
library("car")

trend_df = data.frame(trend=time.pm, X.clust[,-1])
head(trend_df)
tlm = lm(trend ~ ., data=trend_df)
summary(tlm)
pdf(file=paste0("plots/avPlots_trend_", modK, "_anchor", n_anchor, "_seed", seed[1], ".pdf"), width=12, height=10)
avPlots(tlm, ask=FALSE)
dev.off()

amp_df = data.frame(amp=amp.pm, X.clust[,-1])
head(trend_df)
alm = lm(amp ~ ., data=amp_df)
summary(alm)
pdf(file=paste0("plots/avPlots_amplitude_", modK, "_anchor", n_anchor, "_seed", seed[1], ".pdf"), width=12, height=10)
avPlots(alm, ask=FALSE)
dev.off()

phx_df = data.frame(phx=cos(phase.pm), X.clust[,-1])
head(phx_df)
phxlm = lm(phx ~ ., data=phx_df)
summary(phxlm)
pdf(file=paste0("plots/avPlots_phaseX_", modK, "_anchor", n_anchor, "_seed", seed[1], ".pdf"), width=12, height=10)
avPlots(phxlm, ask=FALSE)
dev.off()

phy_df = data.frame(phy=sin(phase.pm), X.clust[,-1])
head(phy_df)
phylm = lm(phy ~ ., data=phy_df)
summary(phylm)
pdf(file=paste0("plots/avPlots_phaseY_", modK, "_anchor", n_anchor, "_seed", seed[1], ".pdf"), width=12, height=10)
avPlots(phylm, ask=FALSE)
dev.off()

