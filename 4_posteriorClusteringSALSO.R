rm(list=ls())

## Libraries
library(tidyverse)

## Load in the Posterior Draws

# load("./GGConstrainedResults_K5-5-12_02082021.RData"); modK = "5-5-12"
# load("Results_anchor2_K5-7-12_KMtypesupplied_22221.RData"); rep_id = 1; modK="5-7-12" # production runs?
# load("Results_anchor2_K7-11-12_KMtypesupplied_22221.RData"); rep_id = 1; modK="7-11-12" # production runs?

allrdafiles = list.files("results/")
indx_files_use = grep("anchor2_K5-7-12_KMtypesupplied_2222\\d", allrdafiles)
# indx_files_use = grep("anchor2_K7-11-12_KMtypesupplied_2222\\d_cont3", allrdafiles)
(rdafiles_use = allrdafiles[indx_files_use])
(nchains = length(rdafiles_use))
source("scripts/0_combine_chains.R")
rep_id = 1234

n.draws
(niter_tot = nrow(int.draws))
n.draws == niter_tot/nchains
str(llike.site.draws)

nkeep = 1000
it_indx = floor(seq(from=1, to=niter_tot, length=nkeep)); head(it_indx); tail(it_indx)
chain_id = floor((it_indx - 1) / n.draws) + 1
table(chain_id)



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




ls()
dim(tclust.draws)

(K_time = ncol(time.draws))
apply(tclust.draws, 1, max)
tclust.draws[1:5, 1:20]

(pm_time = colMeans(time.draws))
trend_estimates = tclust.probs %*% pm_time %>% drop()

library("salso")

set.seed(1)
n_plot = 100
# indx_obs_use = 1:20
# indx_obs_use = 1:100
# indx_obs_use = 1:ncol(tclust.draws)
indx_obs_use = sort(sample(1:ncol(tclust.draws), n_plot, replace=FALSE))

trend_ord = order(trend_estimates[indx_obs_use])

tclust_psm = psm(tclust.draws[it_indx, indx_obs_use[trend_ord]])
pdf(file=paste0("plots/clust_heatmap_trend", modK, "_anch", n.fixed, "_seed", seed[1], rep_id, ".pdf"), height=8, width=7)
heatmap(tclust_psm, Rowv=NA, Colv=NA, scale="none", main=paste0("Trend\n", modK), labRow=NA, labCol=NA)
dev.off()
str(tclust_psm)
# plot(tclust_psm)
# tclust_psm[1:10, 1:10]

# indx_obs_use = 1:ncol(tclust.draws)
indx_obs_use = sort(sample(1:ncol(tclust.draws), n_plot, replace=FALSE))

set.seed(1)
trend_ord = order(trend_estimates[indx_obs_use])

# tclust_salso = salso(tclust.draws[it_indx, indx_obs_use[trend_ord]])
tclust_salso = salso(tclust.draws[it_indx, indx_obs_use[trend_ord]], loss=VI(a=0.7)) # higher a => fewer clusters

tclust_salso_summ = summary(tclust_salso)
plot(tclust_salso_summ)

tclust_salso_summ$sizes

tloss = partition.loss(tclust_salso_summ$estimate, tclust.draws[it_indx, indx_obs_use[trend_ord]], loss = VI(a=0.7))
str(tloss)
summary(tloss)


tclass = apply(tclust.probs[indx_obs_use[trend_ord],], 1, which.max)

length(unique(tclass))
unique(tclass)

tclassF = as.factor(tclass)
levels(tclassF)
tclassI = as.numeric(tclassF)

# indx = 1:n_plot
indx = 1:length(indx_obs_use)
# indx = indx_obs_use
length(indx)

compare_tclass = cbind(tclassI[indx], tclust_salso_summ$estimate[indx]) # ?? all lumpted together when using variation of information loss; more interesting with Binder loss

# table(tclust.draws[,100])

table(tclassI)
tclust_salso_summ$sizes

(ttab = table(tclassI=factor(tclassI[indx], levels=1:K_time), SALSO=factor(tclust_salso_summ$estimate[indx], levels=1:K_time)))
addmargins(ttab)
library("xtable")
xtable(addmargins(ttab))

time.draws.wrtEst = sapply(1:nrow(time.draws), function(ii) {
  tapply(X=time.draws[ii,][tclust.draws[ii, indx_obs_use[trend_ord]]], INDEX=tclust_salso_summ$estimate[indx], FUN = mean)
}) %>% t()
head(time.draws.wrtEst)

library("bayesplot")
pdf(file=paste0("plots/clustMeansSALSO_trend", modK, "_anch", n.fixed, "_seed", seed[1], rep_id, ".pdf"), height=3, width=5)
mcmc_areas(time.draws.wrtEst) + labs(title="Cluster Means with respect to Optimal Cluster\nTrend")
mcmc_areas(matrix(time.draws, ncol=K_time, dimnames=list(NULL, paste(1:K_time)))) + labs(title="Cluster Means\nTrend")
dev.off()
## relatively little change




### Amplitude
ls()
dim(aclust.draws)

(K_amp = ncol(amp.draws))
apply(aclust.draws, 1, max)
aclust.draws[1:5, 1:20]

(pm_amp = colMeans(amp.draws))
amp_estimates = aclust.probs %*% pm_amp %>% drop()

library("salso")

set.seed(1)
n_plot = 100
# indx_obs_use = 1:20
# indx_obs_use = 1:100
# indx_obs_use = 1:ncol(tclust.draws)
indx_obs_use = sort(sample(1:ncol(aclust.draws), n_plot, replace=FALSE))

amp_ord = order(amp_estimates[indx_obs_use])

aclust_psm = psm(aclust.draws[it_indx, indx_obs_use[amp_ord]])
pdf(file=paste0("plots/clust_heatmap_amp", modK, "_anch", n.fixed, "_seed", seed[1], rep_id, ".pdf"), height=8, width=7)
heatmap(aclust_psm, Rowv=NA, Colv=NA, scale="none", main=paste0("Amplitude\n", modK), labRow=NA, labCol=NA)
dev.off()
str(aclust_psm)
# plot(tclust_psm)
# tclust_psm[1:10, 1:10]

# indx_obs_use = 1:ncol(aclust.draws)
indx_obs_use = sort(sample(1:ncol(aclust.draws), n_plot, replace=FALSE))

set.seed(1)
amp_ord = order(amp_estimates[indx_obs_use])

# aclust_salso = salso(aclust.draws[,indx_obs_use[amp_ord]])
aclust_salso = salso(aclust.draws[it_indx, indx_obs_use[amp_ord]], loss=VI(a=0.7)) # higher a => fewer clusters

aclust_salso_summ = summary(aclust_salso)
plot(aclust_salso_summ)

aclust_salso_summ$sizes


aclass = apply(aclust.probs[indx_obs_use[amp_ord],], 1, which.max)

length(unique(aclass))
unique(aclass)

aclassF = as.factor(aclass)
levels(aclassF)
aclassI = as.numeric(aclassF)

# indx = 1:n_plot
indx = indx_obs_use

compare_aclass = cbind(aclassI[indx], aclust_salso_summ$estimate[indx]) # ?? all lumpted together when using variation of information loss; more interesting with Binder loss

# table(tclust.draws[,100])

table(aclassI)
aclust_salso_summ$sizes

(atab = table(aclassI=factor(aclassI[indx], levels=1:K_amp), SALSO=factor(aclust_salso_summ$estimate[indx], levels=1:K_amp)))
addmargins(atab)
library("xtable")
xtable(addmargins(atab))

amp.draws.wrtEst = sapply(1:nrow(amp.draws), function(ii) {
  tapply(X=amp.draws[ii,][aclust.draws[ii, indx_obs_use[amp_ord]]], INDEX=aclust_salso_summ$estimate[indx], FUN = mean)
}) %>% t()
head(amp.draws.wrtEst)

library("bayesplot")
pdf(file=paste0("plots/clustMeansSALSO_amp", modK, "_anch", n.fixed, "_seed", seed[1], rep_id, ".pdf"), height=3, width=5)
mcmc_areas(amp.draws.wrtEst) + labs(title="Cluster Means with respect to Optimal Cluster\nAmplitude")
mcmc_areas(matrix(amp.draws, ncol=K_amp, dimnames=list(NULL, paste(1:K_amp)))) + labs(title="Cluster Means\nAmplitude")
dev.off()
## relatively little change




### Phase
ls()
dim(pclust.draws)

(K_ph = ncol(phase.draws))
apply(pclust.draws, 1, max)
pclust.draws[1:5, 1:20]

(pm_ph = colMeans(phase.draws))
ph_estimates = pclust.probs %*% pm_ph %>% drop()

set.seed(1)
n_plot = 300
# indx_obs_use = 1:20
# indx_obs_use = 1:100
# indx_obs_use = 1:ncol(tclust.draws)
indx_obs_use = sort(sample(1:ncol(pclust.draws), n_plot, replace=FALSE))

ph_ord = order(ph_estimates[indx_obs_use])

pclust_psm = psm(pclust.draws[it_indx, indx_obs_use[ph_ord]])
pdf(file=paste0("plots/clust_heatmap_phase", modK, "_anch", n.fixed, "_seed", seed[1], rep_id, ".pdf"), height=8, width=7)
heatmap(pclust_psm, Rowv=NA, Colv=NA, scale="none", main=paste0("Phase\n", modK), labRow=NA, labCol=NA)
dev.off()
str(pclust_psm)
# plot(tclust_psm)
# tclust_psm[1:10, 1:10]


library("salso")

indx_obs_use = 1:ncol(pclust.draws)
indx_obs_use = sort(sample(1:ncol(pclust.draws), n_plot, replace=FALSE))

set.seed(1)
ph_ord = order(ph_estimates[indx_obs_use])

# pclust_salso = salso(pclust.draws[it_indx, indx_obs_use[ph_ord]])
pclust_salso = salso(pclust.draws[it_indx, indx_obs_use[ph_ord]], loss=VI(a=0.7)) # higher a => fewer clusters

pclust_salso_summ = summary(pclust_salso)
plot(pclust_salso_summ)

pclust_salso_summ$sizes


pclass = apply(pclust.probs[indx_obs_use[ph_ord],], 1, which.max)

length(unique(pclass))
unique(pclass)

pclassF = as.factor(pclass)
levels(pclassF)
pclassI = as.numeric(pclassF)

# indx = 1:n_plot
indx = indx_obs_use

compare_pclass = cbind(pclassI[indx], pclust_salso_summ$estimate[indx]) # ?? all lumpted together when using variation of information loss; more interesting with Binder loss

# table(tclust.draws[,100])

table(pclassI)
pclust_salso_summ$sizes

(ptab = table(pclassI=factor(pclassI[indx], levels=1:K_ph), SALSO=factor(pclust_salso_summ$estimate[indx], levels=1:K_ph)))
addmargins(ptab)

library("xtable")
xtable(addmargins(ptab))


ph.draws.wrtEst = sapply(1:nrow(phase.draws), function(ii) {
  tapply(X=phase.draws[ii,][pclust.draws[ii, indx_obs_use[ph_ord]]], INDEX=pclust_salso_summ$estimate[indx], FUN = mean)
}) %>% t()
head(ph.draws.wrtEst)

library("bayesplot")
mcmc_areas(ph.draws.wrtEst) + labs(title="Cluster Means with respect to Optimal Cluster\nPhase")
mcmc_areas(matrix(phase.draws, ncol=K_ph, dimnames=list(NULL, paste(1:K_ph)))) + labs(title="Cluster Means\nPhase")



## Combine phase clusters into seasons

head(pclassI)
pclassIS = factor(rep(NA, length(pclassI)), levels=1:4)
head(pclassIS)
pclassIS[which(pclassI %in% 1:3)] = 1
pclassIS[which(pclassI %in% 4:6)] = 2
pclassIS[which(pclassI %in% 7:9)] = 3
pclassIS[which(pclassI %in% 10:12)] = 4
table(pclassIS)

length(tclassI)
length(aclassI)
length(pclassIS)

table(pclassIS=factor(pclassIS[indx], levels=1:4), SALSO=factor(pclust_salso_summ$estimate[indx], levels=1:K_ph))



## Combine time, amp, phase clusters into single cluster index

sort(unique(tclass))
sort(unique(aclass))

tapClass = paste(tclassI, aclassI, pclassI, sep="")
tapClass = paste(tclassI, aclassI, pclassIS, sep="")
str(tapClass)
tapClassTab = table(tapClass)
as.data.frame(tapClassTab[order(tapClassTab, decreasing=TRUE)])

plot(tapClassTab[order(tapClassTab, decreasing = TRUE)])
str(tapClassTab)

plot(as.numeric(tapClassTab[order(tapClassTab, decreasing = TRUE)]))

n_topClass = 10
which_topClass = 1:3
which_topClass = 4:7
(topClass = dimnames(tapClassTab[order(tapClassTab, decreasing = TRUE)])[[1]][which_topClass])

locs
cbind(locs, tapClass)

france = map_data("world", region = "France")
france.plot = ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)

france.plot + 
  geom_point(data=cbind(locs, tapClass) %>% filter(tapClass %in% topClass), aes(x=lon, y=lat, color=tapClass), size=.6)

