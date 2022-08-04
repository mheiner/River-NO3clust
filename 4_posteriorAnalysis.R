###
### Analysis of MCMC Draws from model fit
### Code to generate Figures 4-7 is found in GoF.R
###

rm(list=ls())

## Libraries
library(tidyverse)

## Load in the Posterior Draws
load("results/Results_anchor2_K5-7-12_KMtypesupplied_22224.RData"); rep_id = 4; modK="5-7-12"; n_anchor = 2 # USED in paper (after two rounds of chain continuation)

allrdafiles = list.files("results/")
indx_files_use = grep("anchor2_K5-7-12_KMtypesupplied_2222\\d", allrdafiles)
(rdafiles_use = allrdafiles[indx_files_use])
(nchains = length(rdafiles_use))
source("scripts/0_combine_chains.R")
rep_id = 1234



## Read in the Original Data
source("scripts/0_readData.R")
summary(rivers$time)


## Look at Trace plots of parameters
plt.df <- data.frame(time.draws, Draw=1:nrow(time.draws))
time.df <- pivot_longer(plt.df, cols=1:ncol(time.draws), names_to="Cluster", values_to="TimeEffect")
ggplot(time.df, aes(x=Draw, y=TimeEffect, color=Cluster)) + geom_line()

plt.df <- data.frame(amp.draws, Draw=1:nrow(amp.draws))
amp.df <- pivot_longer(plt.df, cols=1:ncol(amp.draws), names_to="Cluster", values_to="Amplitude")
ggplot(amp.df, aes(x=Draw, y=Amplitude, color=Cluster)) + geom_line()

plt.df <- data.frame(phase.draws, Draw=1:nrow(phase.draws))
phase.df <- pivot_longer(plt.df, cols=1:ncol(phase.draws), names_to="Cluster", values_to="Phase")
ggplot(phase.df, aes(x=Draw, y=Phase, color=Cluster)) + geom_line()


library("coda")
(efsz_time = effectiveSize(time.draws))
(efsz_amp = effectiveSize(amp.draws))
(efsz_ph = effectiveSize(phase.draws))

time.cut.diffs = apply(time.cuts.draws[,-c(1,6)], 1, diff) %>% t()
dim(time.cut.diffs)
plot(as.mcmc(time.cut.diffs))

amp.cut.diffs = apply(amp.cuts.draws[,-c(1,6)], 1, diff) %>% t()
dim(amp.cut.diffs)
plot(as.mcmc(amp.cut.diffs))

colMeans(time.cuts.draws)
diff(colMeans(time.cuts.draws))

colMeans(amp.cuts.draws)
diff(colMeans(amp.cuts.draws))




## Map of Cluster Probs
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

france = map_data("world", region = "France")
france.plot = ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)
locs <- rivers %>% group_by(CDSTATIONM) %>% 
  summarise(lat=unique(lat)[1], lon=unique(lon)[1], alt=unique(elevation)[1], dept=unique(LBDEPARTEM))

france.plot + geom_point(data=cbind(locs, probs.df['TP1']), 
                         aes(x=lon, y=lat, color=TP1, alpha=TP1), size=0.75) +
  scale_color_distiller(palette="Spectral")

france.plot + geom_point(data=cbind(locs, probs.df[c('TP1', 'TP2')]), 
                         aes(x=lon, y=lat, color=TP1+TP2, alpha=TP1+TP2), size=0.75) +
  scale_color_distiller(palette="Spectral")


mean((probs.df$TP1 + probs.df$TP2) > 0.95)
mean(probs.df$TP4 + probs.df$TP5 > 0.95)
mean(probs.df$TP5 > 0.95)
mean((probs.df$TP1 + probs.df$TP2 + probs.df$TP5) > 0.95)
mean(probs.df$TP3 < 0.05)
mean(probs.df$TP3 + probs.df$TP4 < 0.05)

### if using 7 trends
mean(probs.df$TP1 + probs.df$TP2 > 0.95)
mean(probs.df$TP6 + probs.df$TP7 > 0.95)
mean(probs.df$TP4 < 0.05)
###


france.plot + geom_point(data=cbind(locs, probs.df[c('TP4', 'TP5')]), 
                         aes(x=lon, y=lat, color=TP4+TP5, alpha=TP4+TP5), size=0.75) +
  scale_color_distiller(palette="Spectral")

france.plot + geom_point(data=cbind(locs, probs.df[c('TP5')]), 
                         aes(x=lon, y=lat, color=TP5, alpha=TP5), size=0.75) +
  scale_color_distiller(palette="Spectral")



france.plot + geom_point(data=cbind(locs, probs.df['AP6']), 
                         aes(x=lon, y=lat, color=AP6, alpha=AP6), size=0.75) +
  scale_color_distiller(palette="Spectral")

france.plot + geom_point(data=cbind(locs, probs.df['AP7']), 
                         aes(x=lon, y=lat, color=AP7, alpha=AP7), size=0.75) +
  scale_color_distiller(palette="Spectral")


tclass.probs = data.frame(class=apply(tclust.probs, 1, which.max), maxp=apply(tclust.probs, 1, max))
head(tclass.probs)
boxplot(maxp ~ class, data=tclass.probs, ylim=c(0,1))

aclass.probs = data.frame(class=apply(aclust.probs, 1, which.max), maxp=apply(aclust.probs, 1, max))
head(aclass.probs)
boxplot(maxp ~ class, data=aclass.probs, ylim=c(0,1))

pclass.probs = data.frame(class=apply(pclust.probs, 1, which.max), maxp=apply(pclust.probs, 1, max))
head(pclass.probs)
boxplot(maxp ~ class, data=pclass.probs, ylim=c(0,1))



head(pclust.probs)
pclust.noccup = apply(pclust.probs, 1, function(x) sum(x > 0.2))
hist(pclust.noccup)
table(pclust.noccup)
library("lattice")
levelplot(pclust.probs[pclust.noccup < 1,]) # these just tend to be more spread out -- nothing weird happening
levelplot(pclust.probs[pclust.noccup == 1,][1:20,])
levelplot(pclust.probs[pclust.noccup == 2,][1:20,])
levelplot(pclust.probs[pclust.noccup == 3,][1:20,])


head(tclust.probs)
tclust.noccup = apply(tclust.probs, 1, function(x) sum(x > 0.2))
table(tclust.noccup)
library("lattice")
levelplot(tclust.probs[tclust.noccup == 1,][111:120,])
levelplot(tclust.probs[tclust.noccup == 2,][111:120,])
levelplot(tclust.probs[tclust.noccup == 3,][50:83,])

# hypotheses on tclust
mean(tclust.probs[,3] < .05)
mean(rowSums(tclust.probs[,4:5]) > .95)
mean(rowSums(tclust.probs[,1:2]) > .95)



head(aclust.probs)
aclust.noccup = apply(aclust.probs, 1, function(x) sum(x > 0.2))
table(aclust.noccup)
library("lattice")
levelplot(aclust.probs[aclust.noccup == 1,][111:140,])
levelplot(aclust.probs[aclust.noccup == 2,][111:140,])
levelplot(aclust.probs[aclust.noccup == 3,][20:70,])






### Inspect MCMC relationships among cluster membership, sizes, and parameter estimates
tclust.size = apply(tclust.draws, 1, function(x){
  table(factor(x, levels=1:ncol(time.draws)))
}) %>% t()
head(tclust.size)
plot.ts(tclust.size[,4])
(pm_tclust.size = colMeans(tclust.size))

table(tclust.size[,5])

aclust.size = apply(aclust.draws, 1, function(x){
  table(factor(x, levels=1:ncol(amp.draws)))
}) %>% t()
head(aclust.size)
(pm_aclust.size = colMeans(aclust.size))

table(aclust.size[,4])

pclust.size = apply(pclust.draws, 1, function(x){
  table(factor(x, levels=1:ncol(phase.draws)))
}) %>% t()
head(pclust.size)
plot.ts(pclust.size[,5])
(pm_pclust.size = colMeans(pclust.size))

table(pclust.size[,2])


plot(as.mcmc(tclust.size))
plot(as.mcmc(aclust.size))
plot(as.mcmc(pclust.size))


plot.ts(amp.draws[,3])
(pm_amp = colMeans(amp.draws))

plot.ts(time.draws[,3])
(pm_time = colMeans(time.draws))

colMeans(phase.draws)

save(file=paste0("results/pm_time_amp_", modK, "_7yr.rda"), pm_amp, pm_time, pm_aclust.size, pm_tclust.size)


k = 4 # select which cluster
plot.ts(tclust.size[,k])
plot.ts(time.draws[,k])
plot(tclust.size[,k], time.draws[,k])

plot.ts(aclust.size[,k])
plot.ts(amp.draws[,k])
plot(aclust.size[,k], amp.draws[,k])

plot.ts(pclust.size[,k])
plot.ts(phase.draws[,k])
plot(pclust.size[,k], phase.draws[,k])

plot(pm_tclust.size, efsz_time)
plot(pm_aclust.size, efsz_amp)
plot(pm_pclust.size, efsz_ph)






## Analysis of Covariate effects on cluster probabilities
apply(tclust.beta.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>%
  write.csv(x=., file=paste0("results/TimeBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".csv"))

plt <- apply(tclust.beta.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5)  + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/TimeBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=9, height=18)

plt <- apply(tclust.beta.draws[,-grep("Department", colnames(tclust.beta.draws))], 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5)  + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/TimeBetaIntervals_small_", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=4, height=5)


apply(aclust.beta.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>%
  write.csv(x=., file=paste0("results/AmpBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".csv"))

plt <- apply(aclust.beta.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5) + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/AmpBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=9, height=18)

plt <- apply(aclust.beta.draws[,-grep("Department", colnames(aclust.beta.draws))], 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5) + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/AmpBetaIntervals_small_", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=4, height=5)


apply(pclust.beta.x.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>%
  write.csv(x=., file=paste0("results/PhaseXBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".csv"))

plt <- apply(pclust.beta.x.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5) + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/PhaseXBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=9, height=18)

plt <- apply(pclust.beta.x.draws[,-grep("Department", colnames(pclust.beta.x.draws))], 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5) + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/PhaseXBetaIntervals_small_", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=4, height=5)


apply(pclust.beta.y.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>%
  write.csv(x=., file=paste0("results/PhaseYBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".csv"))

plt <- apply(pclust.beta.y.draws, 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5) + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/PhaseYBetaIntervals", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=9, height=18)


plt <- apply(pclust.beta.y.draws[,-grep("Department", colnames(pclust.beta.y.draws))], 2, quantile, probs=c(0.025, 0.5, 0.975)) %>% t() %>% as.data.frame() %>% 
  cbind(var=rownames(.), .) %>% rename(., covariate=var, lower="2.5%", effect="50%", upper="97.5%") %>% 
  ggplot(aes(x=effect, y=covariate)) + geom_vline(aes(xintercept=0, color="red"), size=0.5) + 
  geom_linerange(aes(xmin=lower, xmax=upper)) + theme(legend.position="none")
ggsave(paste0("plots/PhaseYBetaIntervals_small_", modK, "_anchor", n_anchor, "_7yr_seed", seed[1], ".pdf"), plt, width=4, height=5)





### Department effects on clustering probabilities

library("RColorBrewer")
france.plot.dep = ggplot() + geom_polygon(data=map_data("france"), fill="white", colour = "gray30", aes(x = long, y = lat, group = group)) +
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


tclustBeta_med = apply(tclust.beta.draws, 2, median) %>% enframe() %>% mutate(name = gsub("Department", "", name)) %>% rename(dept=name, tclustBeta_med=value)
aclustBeta_med = apply(aclust.beta.draws, 2, median) %>% enframe() %>% mutate(name = gsub("Department", "", name)) %>% rename(dept=name, aclustBeta_med=value)
pxclustBeta_med = apply(pclust.beta.x.draws, 2, median) %>% enframe() %>% mutate(name = gsub("Department", "", name)) %>% rename(dept=name, pxclustBeta_med=value)
pyclustBeta_med = apply(pclust.beta.y.draws, 2, median) %>% enframe() %>% mutate(name = gsub("Department", "", name)) %>% rename(dept=name, pyclustBeta_med=value)

locs = left_join(locs, tclustBeta_med, by="dept")
locs = left_join(locs, aclustBeta_med, by="dept")
locs = left_join(locs, pxclustBeta_med, by="dept")
locs = left_join(locs, pyclustBeta_med, by="dept")

locs

pdf(file="plots/DeptEffectMaps.pdf", height=4, width=5)

plt = france.plot.dep.clean + 
  geom_point(data=locs, aes(x=lon, y=lat, color=(round(tclustBeta_med,3))), size=.5) +
  labs(color="Dept.\neffect") + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  annotate("text", x=min(locs$lon), y=max(locs$lat), label="Trend", hjust=0, size=7) 
print(plt)

plt = france.plot.dep.clean + 
  geom_point(data=locs, aes(x=lon, y=lat, color=(round(aclustBeta_med,3))), size=.5) +
  labs(color="Dept.\neffect") + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  annotate("text", x=min(locs$lon), y=max(locs$lat), label="Amplitude", hjust=0, size=7) 
print(plt)

plt = france.plot.dep.clean + 
  geom_point(data=locs, aes(x=lon, y=lat, color=(round(pxclustBeta_med,3))), size=.5) +
  labs(color="Dept.\neffect") + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  annotate("text", x=min(locs$lon), y=max(locs$lat), label="Phase D1\n(x-axis)", hjust=0, size=7) 
print(plt)

plt = france.plot.dep.clean + 
  geom_point(data=locs, aes(x=lon, y=lat, color=(round(pyclustBeta_med,3))), size=.5) +
  labs(color="Dept.\neffect") + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  annotate("text", x=min(locs$lon), y=max(locs$lat), label="Phase D2\n(y-axis)", hjust=0, size=7)
print(plt)

dev.off()
