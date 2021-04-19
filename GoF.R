##
##

rm(list=ls())

## Libraries
library(tidyverse)
library("RColorBrewer")

## Load in the Posterior Draws
load("./results/GGConstrainedResults_K5-5-12_04082021.RData"); modK = "5-5-12"
load("./results/GGConstrainedResults_K7-7-12_04112021.RData"); modK = "7-7-12"

(Ks = as.numeric(unlist(strsplit(modK, "-"))))
(K.time = Ks[1])
(K.amp = Ks[2])
(K.phase = Ks[3])

## Read in the Original Data
source("scripts/0_readData.R")
summary(rivers$time)


## get initial values and fixed quantities (e.g., error correlation matrix, cluster anchors)
source("scripts/0_init.R")


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
  
  mean_at_obs_i95 = apply(mean_at_obs, 2, function(x) quantile(x, c(0.025, 0.975)))
  
  e = matrix(rep(log(river.list[[rn]]$df$NO3), each=niter), nrow=niter) - mean_at_obs
  e_mean = colMeans(e)
  
  mean_at_grid = matrix( rep(int.draws[,rn], times=ntt), ncol=ntt ) + 
    tcrossprod(time.draws[ cbind(1:niter, tclust.draws[,rn]) ], tt) +
    (matrix( rep(amp.draws[ cbind(1:niter, aclust.draws[,rn])], times=ntt), ncol=ntt ) * 
       cos(t(outer(2*pi*tt, phase.draws[ cbind(1:niter, pclust.draws[,rn]) ], FUN="-")))
    )

  mean_at_grid_i95 = apply(mean_at_grid, 2, function(x) quantile(x, c(0.025, 0.975)))
  # R2 = 1.0 - mean(sig2.draws[,rn]) / var(log(river.list[[rn]]$NO3))
  R2 = 1.0 - mean(rowSums(e^2)) / sum( (log(river.list[[rn]]$df$NO3) - mean(log(river.list[[rn]]$df$NO3)))^2 )

  return(list(mean_at_obs=mean_at_obs, mean_at_obs_i95=mean_at_obs_i95, 
              e=e, e_mean=e_mean, R2=R2,
              mean_at_grid=mean_at_grid, mean_at_grid_i95=mean_at_grid_i95))
}

rivermean.list <- lapply(1:length(river.list), get_lNO3mean)


## kurtosis
library("e1071")
kurt = sapply(rivermean.list, function(x) kurtosis(x$e_mean))
hist(kurt, breaks=50)
mean(kurt > 6) # level of excess kurtosis of student's t with 5 df
mean(kurt > 3) # level of excess kurtosis of Laplace distribution
hist(kurt[kurt<6], breaks=50)



## residuals
library("EnvStats") # for testing serial correlation in residuals w/ serialCorrelationTest()

rn = 1
rivermean.list[[rn]]$e_mean
dim(river.list[[rn]]$Rinv)
aa = acf(rivermean.list[[rn]]$e_mean)
serialCorrelationTest(rivermean.list[[rn]]$e_mean)$p.value

for (rn in 1:length(river.list)) {
  rivermean.list[[rn]]$e_mean_w = forwardsolve(river.list[[rn]]$Rchol, rivermean.list[[rn]]$e_mean / sqrt(mean(sig2.draws[,rn])) )
  rivermean.list[[rn]]$e_mean_test = serialCorrelationTest( rivermean.list[[rn]]$e_mean / sqrt(mean(sig2.draws[,rn])), test="rank.von.Neumann")
  rivermean.list[[rn]]$e_mean_wtest = serialCorrelationTest( rivermean.list[[rn]]$e_mean_w, test="rank.von.Neumann")
  rivermean.list[[rn]]$R2w = 1.0 - sum(forwardsolve(river.list[[rn]]$Rchol, rivermean.list[[rn]]$e_mean)^2) / sum( (log(river.list[[rn]]$df$NO3) - mean(log(river.list[[rn]]$df$NO3)))^2 )
  cat(rn, " of ", N, "\r")
}

e_tests = sapply(1:N, function(rn) c(rivermean.list[[rn]]$e_mean_test$p.value, rivermean.list[[rn]]$e_mean_wtest$p.value)) %>% t()
dim(e_tests)
head(e_tests)
mean(e_tests[,1] < 0.05)
mean(e_tests[,2] < 0.05)

e_acf = sapply(1:N, function(rn) c(max(rivermean.list[[rn]]$e_mean_test$estimate), max(rivermean.list[[rn]]$e_mean_wtest$estimate))) %>% t()
hist(e_acf[,2] / e_acf[,1])
hist(apply(abs(e_acf), 1, diff))
hist(e_acf[,1]) # original mean residuals
hist(e_acf[,2]) # whitened mean residuals

rn = 4
(rn = which.min(apply(abs(e_acf), 1, diff))) # also try which.max :/
e_acf[rn,]
plot.ts(cbind(rivermean.list[[rn]]$e_mean, rivermean.list[[rn]]$e_mean_w))
acf(rivermean.list[[rn]]$e_mean)
acf(rivermean.list[[rn]]$e_mean_w)




## R squared
dim(sig2.draws)
sig2_pm = colMeans(sig2.draws)
hist(sig2_pm)
hist(sig2.draws[,2150], breaks=50)

length(river.list) == length(sig2_pm)
length(river.list) == N
# R2 = 1.0 - sig2_pm / sapply(river.list, function(x) var(log(x$NO3)))
# R2 = 1.0 - sapply(rivermean.list, function(x) sum(x$e_mean^2) ) / sapply(river.list, function(x) sum( (log(x$NO3) - mean(log(x$NO3)))^2 ))
# for(rn in 1:N) rivermean.list[[rn]]$R2 = R2[rn]
R2 = sapply(rivermean.list, function(x) x$R2)
R2w = sapply(rivermean.list, function(x) x$R2w)
R2_meanfit = sapply(1:N, function(rn) 1.0 - sum(rivermean.list[[rn]]$e_mean^2) / sum( (log(river.list[[rn]]$df$NO3) - mean(log(river.list[[rn]]$df$NO3)))^2 ))

hist(R2[R2>0], xlim=c(0,1))
hist(R2_meanfit[R2_meanfit>0], xlim=c(0,1))

sum(R2 < 0)
mean(R2 < 0)
mean(R2 < .2)

hist(R2w[R2w>0], xlim=c(0,1))

sum(R2w < 0)
mean(R2w < 0)
mean(R2w < .2)

# R2_55 = sapply(rivermean.list, function(x) x$R2); save(R2_55, file="results/R2_55.rda")
# R2_77 = sapply(rivermean.list, function(x) x$R2); save(R2_77, file="results/R2_77.rda")
load("results/R2_55.rda"); load("results/R2_77.rda")

R2rat75 = R2_77 / R2_55
hist(R2rat75)
summary(R2rat75)
summary(R2rat75 * 100 - 100)

mean( abs(R2rat75 - 1.0) > 0.8)
mean( abs(R2rat75 - 1.0) < 0.1)

hist(R2rat75[abs(R2rat75 - 1.0) < 1.5]*100 - 100, xlim=c(-200,200), breaks=50)

pdf(file="plots/R2_77to55.pdf", height=4, width=5)
hist(R2rat75[abs(R2rat75 - 1.0) < 0.9]*100 - 100, 
     xlim=c(-80,80), breaks=50, main="7-7 from 5-5", xlab="Percent change", border=FALSE)
dev.off()


R2 = R2_55
R2 = sapply(rivermean.list, function(x) x$R2)




pdf(file=paste0("./plots/hist_R2_", modK, "_7yr.pdf"), height=3.5, width=4)
hist(R2, breaks=30)
hist(R2[R2>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="")
hist(R2[R2<0.0])
dev.off()

pdf(file=paste0("./plots/hist_R2_by_mod.pdf"), height=4.5, width=5)
hist(R2_55[R2_55>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="5-5")
hist(R2_77[R2_77>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="7-7")
# hist(R2_1010[R2_1010>0.0], border=FALSE, xlab="R squared", xlim=c(0,1), main="10-10")
dev.off()


indx_negR2 = which(R2 < 0)
length(indx_negR2)
ordR2 = order(R2)
R2[ordR2[1]]
R2[ordR2[N]]



hist(log(rivers$NO3))
nobs = sapply(river.list, nrow)
head(nobs)
which.max(nobs)


plot(tclass.probs$maxp, R2)
summary(lm(R2 ~ tclass.probs$maxp))

plot(aclass.probs$maxp, R2)
summary(lm(R2 ~ aclass.probs$maxp))

plot(pclass.probs$maxp, R2)
summary(lm(R2 ~ pclass.probs$maxp))


### Fit by cluster, uses probs.df from PosteriorAnalysis.R
k = 1
head(probs.df)

rindx = which(probs.df[,paste0("TP",k)] > 0.9)
rindx = which(probs.df[,paste0("AP",k)] > 0.9)
rindx = which(probs.df[,paste0("PP",k)] > 0.9)

length(rindx)

hist(R2[rindx])
hist(R2_77[rindx])
hist(R2_1010[rindx])

(rn = rindx[1])
rindx = rindx[1:20]
rindx = sample(rindx, 50)

mesg_plt = "time"
mesg_plt = "amp"
mesg_plt = "phase"

pdf(file=paste0("plots/fit_by_clust_", modK, "_7yr_", mesg_plt, k, ".pdf"), width=5, height=4)
hist(R2[rindx])
for (rn in rindx) {
  plot(log(NO3) ~ time, data=river.list[[rn]]$df, 
       main=paste0(river.list[[rn]]$df$CDSTATIONM[1], "\n", river.list[[rn]]$df$LBDEPARTEM[1]), 
       xlim=c(0,6), xlab="year", type="o"
       # ylim=range(log(rivers$NO3))
       # ylim=c(-1, 5)
  )
  lines(tt,  colMeans(rivermean.list[[rn]]$mean_at_grid), col="red", lwd=1.5)
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
# rn_seq = floor(seq(1, length(river.list), length=500))
(rn_seq = ordR2[ceiling(N*pp_seq)]); R2[rn_seq]

for( i in 1:length(rn_seq) ) {
  pdf(file=paste0("./plots/fit_resid_model_", modK, "_", sprintf("%03d", round(100*pp_seq[i])), "ile_id", sprintf("%04d.pdf", rn_seq[i])), width=6, height=3.5)
  
  plot(log(NO3) ~ time, data=river.list[[rn_seq[i]]]$df, 
       main=paste0(river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1]), 
       xlim=c(0,6), xlab="year",
       # ylim=range(log(rivers$NO3))
       ylim=c(-1, 5)
  )
  polygon(x=c(tt, rev(tt)), y=c(rivermean.list[[rn_seq[i]]]$mean_at_grid_i95[1,], rev(rivermean.list[[rn_seq[i]]]$mean_at_grid_i95[2,])), col="gray80", border=NA)
  lines(log(NO3) ~ time, data=river.list[[rn_seq[i]]]$df, type="o")
  # arrows(river.list[[rn_seq[i]]]$time, rivermean.list[[rn_seq[i]]]$mean_at_obs_i95[1,], river.list[[rn_seq[i]]]$time, rivermean.list[[rn_seq[i]]]$mean_at_obs_i95[2,], angle=0, code=0, col="red") # CIs
  lines(tt,  colMeans(rivermean.list[[rn_seq[i]]]$mean_at_grid), col="red", lwd=1.5)
  mtext(side=3, line=1, at=5.5, paste("R-squared", round(rivermean.list[[rn_seq[i]]]$R2, 2)))
  
  # nb: obs sometimes not equally spaced
  acf(rivermean.list[[rn_seq[i]]]$e_mean, 
      main=paste0("posterior mean residual, ", river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1])) 
  pacf(rivermean.list[[rn_seq[i]]]$e_mean, 
       main=paste0("posterior mean residual, ", river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1]))

  par(mfrow=c(1,2))
  hist(rivermean.list[[rn_seq[i]]]$e_mean, xlab="e",
       main=paste("Residuals:", river.list[[rn_seq[i]]]$CDSTATIONM[1], "\n", river.list[[rn_seq[i]]]$LBDEPARTEM[1]))
  qqnorm(rivermean.list[[rn_seq[i]]]$e_mean); qqline(rivermean.list[[rn_seq[i]]]$e_mean)
    
  dev.off()
  cat(i, "\r")
}

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
  summarise(lat=unique(lat)[1], lon=unique(lon)[1], alt=unique(altitude)[1])
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

R2df = left_join(R2df, locs[,-grep("lat|lon", colnames(locs))], by="CDSTATIONM", keep=FALSE)
head(R2df)


library("RColorBrewer")
france.plot.dep = ggplot() + geom_polygon(data=map_data("france"), fill="white", colour = "gray40", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)
(rdylbu7 = brewer.pal(7, "RdYlBu"))
(blylrd7 = rev(brewer.pal(7, "RdYlBu")))
rdbu11 = brewer.pal(11, "RdBu")
(burd11 = rev(rdbu11))

france.plot.clean = france.plot + # france.plot.dep too messy 
  # geom_point(data=X, aes(x=lon, y=lat, color=(bioGeoRegion)), size=.6) +
  # labs(color="bioGeoRegion") + 
  coord_cartesian( xlim = c(-5, 7.8), ylim=c(42, 52)) +
  # guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # plot.background=element_blank()
        plot.margin = unit(c(-0.03,0,-0.01,0), "null")
  )

france.plot.dep.clean = france.plot.dep + # france.plot.dep too messy 
  # geom_point(data=X, aes(x=lon, y=lat, color=(bioGeoRegion)), size=.6) +
  # labs(color="bioGeoRegion") + 
  coord_cartesian( xlim = c(-5, 7.8), ylim=c(42, 52) ) +
  # guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # plot.background=element_blank()
        plot.margin = unit(c(-0.03,0,-0.01,0), "null")
  )

col_lims = list()
col_lims[["5-5-12"]] = list(int=c(-1.2, 4.8), time=c(-0.1, 0.1), amp=c(0.0, 1.1))
col_lims[["7-7-12"]] = list(int=c(-1.2, 4.8), time=c(-0.13, 0.13), amp=c(0.0, 1.45))
col_lims[["10-10-12"]] = list(int=c(-1.2, 4.8), time=c(-0.15, 0.15), amp=c(0.0, 1.8))


cpt = seq(col_lims[[modK]][["time"]][1], col_lims[[modK]][["time"]][2], length=12)

R2df$tf = cut(R2df$tclass, breaks = cpt)
tflabs = 1:11
levels(R2df$tf)[sort(unique(as.numeric(R2df$tf)))] <- paste(round(sort(unique(R2df$tclass)), 3))

levels(R2df$tf)

R2df$af = cut(R2df$aclass, breaks = seq(col_lims[[modK]][["amp"]][1], col_lims[[modK]][["amp"]][2], length=12))
aflabs = 1:11
levels(R2df$af)[sort(unique(as.numeric(R2df$af)))] <- paste(round(sort(unique(R2df$aclass)), 3))
levels(R2df$af)



## maps with R-squared and parameter estimates
pdf(file=paste0("./plots/map_R2_", modK, "_7yr.pdf"), width=5, height=4)

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
  geom_point(data=R2df, aes(x=lon, y=lat, color=tf), size=.75) +
  labs(color="Trend") + 
  scale_color_manual(values=burd11[sort(unique(as.numeric(R2df$tf)))]) + guides(colour = guide_legend(reverse=T))
print(plt)

plt = france.plot.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=(round(tclass,3)), alpha=pmax(R2,0)), size=.75) +
  labs(color="Trend", alpha="R-squared") + 
  scale_color_gradient2(limits=col_lims[[modK]][["time"]], midpoint=0, low=blylrd7[1], mid=blylrd7[4], high=blylrd7[7])
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
  geom_point(data=R2df, aes(x=lon, y=lat, color=round(pclass,3), alpha=pmax(R2,0)), size=.75)  +
  labs(color="Phase", alpha="R-squared") + 
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10)) + guides(alpha = guide_legend(reverse=T))

print(plt)

plt = france.plot.dep.clean + 
  geom_point(data=R2df, aes(x=lon, y=lat, color=round(pclass,3), alpha=pmax(R2,0)), size=.75)  +
  labs(color="Phase", alpha="R-squared") + 
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10)) + guides(alpha = guide_legend(reverse=T))

print(plt)

dev.off()





## R squared vs estimates
int.pm = colMeans(int.draws)
time.pm = sapply(1:N, function(i) time.draws[cbind(1:niter, tclust.draws[,i])]) %>% colMeans
amp.pm = sapply(1:N, function(i) amp.draws[cbind(1:niter, aclust.draws[,i])]) %>% colMeans
# phase.pm = sapply(1:N, function(i) phase.draws[cbind(1:niter, pclust.draws[,i])]) %>% colMeans # not good for circular
phase.pm =colMeans(phase.draws)[apply(pclust.probs, 1, which.max)]

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

pdf(paste0("plots/Rsqared_by_estimates_", modK, "_7yr.pdf"), height=4, width=6)
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


K_phase = 12
## Data frame with river summaries
pclass <- apply(pclust.draws, 2, function(x){
  xtab <- table(factor(x, levels=1:K_phase))
  # return(sort(colMeans(phase.draws))[which.max(xtab)])
  which.max(xtab)
})
pclass
table(pclass)

head(probs.df) # need to get this from PosteriorAnalysis.R

dat_station = lapply(1:N, function(i) data.frame(i=i, CDSTATIONM=river.list[[i]]$df$CDSTATIONM[1], 
                                   lon=river.list[[i]]$df$lon[1], lat=river.list[[i]]$df$lat[1],
                                   elevation=river.list[[i]]$df$elevation[1],
                                   LBDEPARTEM=river.list[[i]]$df$LBDEPARTEM[1],
                                   Rsquared=R2[i],
                                   int=int.pm[i], trend=time.pm[i], amp=amp.pm[i], phase=phase.pm[i],
                                   itime = 100*exp(time.pm[i]) - 100,
                                   iamp = exp(2*amp.pm[i]),
                                   tclass=which.max(table(factor(tclust.draws[,i], levels=1:max(tclust.draws[,i])))),
                                   aclass=which.max(table(factor(aclust.draws[,i], levels=1:max(aclust.draws[,i])))),
                                   pclass=which.max(table(factor(pclust.draws[,i], levels=1:max(pclust.draws[,i])))))
                     ) %>% do.call(rbind, .)
head(dat_station)
write.csv(file=paste0("restults/FranceRivers_estimates", modK, "_7yr.csv"), dat_station, row.names=FALSE)


## interpretable time and amplitude
itime.pm = 100*exp(time.pm) - 100
iamp.pm = exp(2*amp.pm)


library("ggplot2")
library("ggExtra")
library("grid")
library("gridExtra")

p = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=itime, y=int)) + 
  geom_hline(yintercept=0, color="gray50") + geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + xlab("Percent annual change in NO3") + ylab("Intercept")
p1 = ggMarginal(p, type="histogram")
p1 = arrangeGrob(p1, top=textGrob("a)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p1)

p = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=amp, y=int)) + 
  geom_hline(yintercept=0, color="gray50") + geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + xlab("Amplitude of log-NO3") + ylab("Intercept")
p2 = ggMarginal(p, type="histogram")
p2 = arrangeGrob(p2, top=textGrob("b)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p2)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter(phase, factor = 0.5)), aes(x=jphase, y=int)) +
  geom_hline(yintercept=0, color="gray50") + # geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + # geom_jitter(width=0.02, inherit.aes=TRUE) + 
  xlab("Time of peak log-NO3") + ylab("Intercept") + 
  scale_x_continuous(breaks=seq(-pi, pi, length=13), labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul"))
p3 = ggMarginal(p, type="histogram")
p3 = arrangeGrob(p3, top=textGrob("c)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p3)

p = ggplot(dat_station %>% filter(Rsquared > 0), aes(x=amp, y=itime)) + 
  geom_hline(yintercept=0, color="gray50") + geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + xlab("Amplitude of log-NO3") + ylab("Percent annual change in NO3")
p4 = ggMarginal(p, type="histogram")
p4 = arrangeGrob(p4, top=textGrob("d)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p4)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter(phase, factor = 0.5)), aes(x=jphase, y=itime)) +
  geom_hline(yintercept=0, color="gray50") + # geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + # geom_jitter(width=0.02, inherit.aes=TRUE) + 
  xlab("Time of peak log-NO3") + ylab("Percent annual change in NO3")  + 
  scale_x_continuous(breaks=seq(-pi, pi, length=13), labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul"))
p5 = ggMarginal(p, type="histogram")
p5 = arrangeGrob(p5, top=textGrob("e)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p5)

p = ggplot(dat_station %>% filter(Rsquared > 0) %>% mutate(jphase = jitter(phase, factor = 0.5)), aes(x=jphase, y=amp)) +
  geom_hline(yintercept=0, color="gray50") + # geom_vline(xintercept=0, color="gray50") + 
  geom_point(size=0.025) + # geom_jitter(width=0.02, inherit.aes=TRUE) + 
  xlab("Time of peak log-NO3") + ylab("Amplitude of log-NO3")  + 
  scale_x_continuous(breaks=seq(-pi, pi, length=13), labels=c("Jul", "", "Sep", "", "Nov", "", "Jan 1", "", "Mar", "", "May", "", "Jul"))
p6 = ggMarginal(p, type="histogram")
p6 = arrangeGrob(p6, top=textGrob("f)", x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
print(p6)

pdf(paste0("plots/estimate_relationships_hist_", modK, "_7yr.pdf"), width=8, height=8)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix=matrix(1:6, ncol=2))
dev.off()

dat_station %>% filter(Rsquared > 0) %>% select(c("int", "trend", "amp")) %>% cor()


library("lattice")
table(dat_station$tclass, dat_station$aclass)
prop.table(table(dat_station$tclass, dat_station$aclass), margin = 2)
levelplot(table(dat_station$tclass, dat_station$aclass))
chisq.test(dat_station$tclass, dat_station$aclass, simulate.p.value=TRUE)

table(dat_station$tclass, dat_station$pclass)
levelplot(table(dat_station$tclass, dat_station$pclass))
chisq.test(dat_station$tclass, dat_station$pclass, simulate.p.value=TRUE)


france = map_data("world", region = "France")
france.plot = ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)

france.plot + 
  geom_point(data=dat_station, aes(x=lon, y=lat, color=as.factor(tclass)), size=.6) +
  labs(color="Time Effect")
france.plot + 
  geom_point(data=dat_station, aes(x=lon, y=lat, color=factor(round(aclass))), size=.6) +
  labs(color="Amplitude Effect")
france.plot + 
  geom_point(data=filter(dat_station, Rsquared > 0), aes(x=lon, y=lat, color=round(Rsquared,3)), size=.6)  +
  labs(color="R Squared") +
  scale_color_distiller(palette="Spectral")





hist(probs.df[,"TP5"])
dim(probs.df)
length(R2_55)

head(locs$pclass)
R2_lm = lm(R2_55 ~ locs$tclass + locs$aclass + as.factor(locs$pclass))
summary(R2_lm)
plot(R2_lm)
plot(locs$aclass, R2_55)



### Covariates
## intercept coincides with AIN dept
colnames(X.clust)
indx_AIN = which(rowSums(X.clust[,grep("Department", colnames(X.clust))]) == 0)
length(indx_AIN)


## these would be clearer with added variable plots

colnames(tclust.beta.draws) # it looks like the last three were back transformed from the principal components V.land estimates?
hist(time.pm[indx_AIN])

length(time.pm)
dim(X.clust)

plot(jitter(X.clust[,"DepartmentAISNE"]), time.pm)
plot(jitter(X.clust[,"DepartmentCOTES-D'ARMOR"]), time.pm)
plot(jitter(X.clust[,"DepartmentMAINE-ET-LOIRE"]), time.pm)

plot(X.clust[,"scale(log(Elevation))"], time.pm)

plot(jitter(X.clust[,"bioGeoRegionatlantic"]), time.pm)
plot(jitter(X.clust[,"bioGeoRegioncontinental"]), time.pm)
plot(jitter(X.clust[,"bioGeoRegionmediterranean"]), time.pm)

plot(X.clust[,"scale(IDPR)"], time.pm)
plot(X.clust[,"scale(annual_specific_runoff)"], time.pm)
plot(X.clust[,"scale(p_sedimentay)"], time.pm)
plot(X.clust[,"scale(log(p_waterbody + 0.01))"], time.pm)
plot(X.clust[,"V1"], time.pm)
plot(X.clust[,"V2"], time.pm)
plot(X.clust[,"V3"], time.pm)

## Amp
hist(amp.pm[indx_AIN])
mean(amp.pm[indx_AIN])

length(amp.pm)
dim(X.clust)

plot(jitter(X.clust[,"DepartmentARDECHE"]), amp.pm)
plot(jitter(X.clust[,"DepartmentEURE"]), amp.pm)
plot(jitter(X.clust[,"DepartmentVENDEE"]), amp.pm)

plot(X.clust[,"scale(log(Elevation))"], amp.pm)

plot(jitter(X.clust[,"bioGeoRegionatlantic"]), amp.pm)
plot(jitter(X.clust[,"bioGeoRegioncontinental"]), amp.pm)
plot(jitter(X.clust[,"bioGeoRegionmediterranean"]), amp.pm)

plot(X.clust[,"scale(IDPR)"], amp.pm)
plot(X.clust[,"scale(annual_specific_runoff)"], amp.pm)
plot(X.clust[,"scale(p_sedimentay)"], amp.pm)
plot(X.clust[,"scale(log(p_waterbody + 0.01))"], amp.pm)
plot(X.clust[,"V1"], amp.pm)
plot(X.clust[,"V2"], amp.pm)
plot(X.clust[,"V3"], amp.pm)


## Phase X
hist(cos(phase.pm)[indx_AIN])
mean(cos(phase.pm)[indx_AIN])

length(cos(phase.pm))
dim(X.clust)

plot(jitter(X.clust[,"DepartmentAISNE"]), cos(phase.pm))
plot(jitter(X.clust[,"DepartmentCOTES-D'ARMOR"]), cos(phase.pm))
plot(jitter(X.clust[,"DepartmentMAINE-ET-LOIRE"]), cos(phase.pm))

plot(X.clust[,"scale(log(Elevation))"], cos(phase.pm))

plot(jitter(X.clust[,"bioGeoRegionatlantic"]), cos(phase.pm))
plot(jitter(X.clust[,"bioGeoRegioncontinental"]), cos(phase.pm))
plot(jitter(X.clust[,"bioGeoRegionmediterranean"]), cos(phase.pm))

plot(X.clust[,"scale(IDPR)"], cos(phase.pm))
plot(X.clust[,"scale(annual_specific_runoff)"], cos(phase.pm))
plot(X.clust[,"scale(p_sedimentay)"], cos(phase.pm))
plot(X.clust[,"scale(log(p_waterbody + 0.01))"], cos(phase.pm))
plot(X.clust[,"V1"], cos(phase.pm))
plot(X.clust[,"V2"], cos(phase.pm))
plot(X.clust[,"V3"], cos(phase.pm))

## Phase Y
hist(sin(phase.pm)[indx_AIN])
mean(sin(phase.pm)[indx_AIN])

length(sin(phase.pm))
dim(X.clust)

plot(jitter(X.clust[,"DepartmentARIEGE"]), sin(phase.pm))
plot(jitter(X.clust[,"DepartmentHAUTE-LOIRE"]), sin(phase.pm))
plot(jitter(X.clust[,"DepartmentLOIRET"]), sin(phase.pm))

plot(X.clust[,"scale(log(Elevation))"], sin(phase.pm))

plot(jitter(X.clust[,"bioGeoRegionatlantic"]), sin(phase.pm))
plot(jitter(X.clust[,"bioGeoRegioncontinental"]), sin(phase.pm))
plot(jitter(X.clust[,"bioGeoRegionmediterranean"]), sin(phase.pm))

plot(X.clust[,"scale(IDPR)"], sin(phase.pm))
plot(X.clust[,"scale(annual_specific_runoff)"], sin(phase.pm))
plot(X.clust[,"scale(p_sedimentay)"], sin(phase.pm))
plot(X.clust[,"scale(log(p_waterbody + 0.01))"], sin(phase.pm))
plot(X.clust[,"V1"], sin(phase.pm))
plot(X.clust[,"V2"], sin(phase.pm))
plot(X.clust[,"V3"], sin(phase.pm))


