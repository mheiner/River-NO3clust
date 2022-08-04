library("tidyverse")
library("coda")

rm(list=ls())

nchains = 2
nburn = 0

trends = list()
trend_sizes = list()

amps = list()
amp_sizes = list()

phs = list()
ph_sizes = list()

tcuts = list()
acuts = list()

landbetas_trend = list()
landbetas_amp = list()
landbetas_phx = list()
landbetas_phy = list()

ints = list()
sig2s = list()

lliks = list()


for (ii in 1:nchains) {

  ### select runs  
  # load(paste0("results/Results_anchor2_K5-7-12_KMtypesupplied_2222", ii, ".RData"))
  # load(paste0("results/Results_anchor2_K5-7-12_KMtypesupplied_2222", ii, "_cont1.RData"))
  load(paste0("results/Results_anchor2_K5-7-12_KMtypesupplied_2222", ii, ".RData")); modK = "5-7-12";

  # load(paste0("results/Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, ".RData"))
  # load(paste0("results/Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, "_cont1.RData"))
  # load(paste0("results/Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, "_cont2.RData")); modK = "7-11-12";
  # load(paste0("results/Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, "_cont3.RData")); modK = "7-11-12";
  
  # load(paste0("results/Results_anchor2_K5-7-12_KMtypesupplied_CV5_2222", ii, "_cont1.RData")); modK = "5-7-12"; # each chain (fit to different obs) looks fine
  # load(paste0("results/Results_anchor2_K7-11-12_KMtypesupplied_CV5_2222", ii, "_cont1.RData")); modK = "7-11-12"; # individual chains are still wandering for trends, amplitudes have an early shift

  ## only 2 chains each of the 15 anchor runs
  # load(paste0("results/Results_anchor15_K5-7-12_KMtypesupplied_2222", ii, "_cont1.RData")) # trend chains still moving, amp chains a little different
  # load(paste0("results/Results_anchor15_K7-11-12_KMtypesupplied_2222", ii, "_cont1.RData")) # struggling
  
  
  ###
  indx_trend0 = which(time.eff == 0.0)
  nkeep_now = nrow(time.draws)
  
  trends[[ii]] = time.draws[(nburn+1):nkeep_now, -indx_trend0] %>% as.mcmc()
  trend_sizes[[ii]] = apply(tclust.draws[(nburn+1):nkeep_now,], 1, function(x){
    table(factor(x, levels=1:ncol(time.draws)))
    }) %>% t() %>% as.mcmc()
  
  amps[[ii]] = amp.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  amp_sizes[[ii]] = apply(aclust.draws[(nburn+1):nkeep_now,], 1, function(x){
    table(factor(x, levels=1:ncol(amp.draws)))
    }) %>% t() %>% as.mcmc()
  
  phs[[ii]] = phase.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  ph_sizes[[ii]] = apply(pclust.draws[(nburn+1):nkeep_now,], 1, function(x){
    table(factor(x, levels=1:ncol(phase.draws)))
    }) %>% t() %>% as.mcmc()
  
  tcuts[[ii]] = time.cuts.draws[(nburn+1):nkeep_now, -c(1:2, length(time.cuts))] %>% as.mcmc()
  acuts[[ii]] = amp.cuts.draws[(nburn+1):nkeep_now, -c(1:2, length(amp.cuts))] %>% as.mcmc()
  
  landbetas_trend[[ii]] = tclust.beta.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  landbetas_amp[[ii]] = aclust.beta.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  landbetas_phx[[ii]] = pclust.beta.x.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  landbetas_phy[[ii]] = pclust.beta.y.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  
  ints[[ii]] = int.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  sig2s[[ii]] = sig2.draws[(nburn+1):nkeep_now,] %>% as.mcmc()
  
  lliks[[ii]] = llike.draws[(nburn+1):nkeep_now] %>% as.mcmc()
  
  cat(ii, "of", nchains, "\n")
}

it
it_prev

summary(trends)
summary(trends[[1]])

GR_trend = gelman.diag(trends)
GR_trendsize = gelman.diag(trend_sizes)
GR_amp = gelman.diag(amps)
GR_ampsize = gelman.diag(amp_sizes)
GR_ph = gelman.diag(phs)
GR_phsize = gelman.diag(ph_sizes)
GR_tcut = gelman.diag(tcuts)
GR_acut = gelman.diag(acuts)
GR_landbeta_trend = gelman.diag(landbetas_trend)
GR_landbeta_amp = gelman.diag(landbetas_amp)
GR_landbeta_phx = gelman.diag(landbetas_phx)
GR_landbeta_phy = gelman.diag(landbetas_phy)
GR_int = gelman.diag(ints)
GR_sig2s = gelman.diag(sig2s)
GR_landbeta_llik = gelman.diag(lliks)

GR_landbeta_trend = gelman.diag(lapply(landbetas_trend, function(xx) xx[,-105]))
GR_landbeta_amp = gelman.diag(lapply(landbetas_amp, function(xx) xx[,-105]))
GR_landbeta_phx = gelman.diag(lapply(landbetas_phx, function(xx) xx[,-105]))
GR_landbeta_phy = gelman.diag(lapply(landbetas_phy, function(xx) xx[,-105]))


summary(GR_int$psrf[,1])

summary(c(GR_landbeta_trend$psrf[,1], 
          GR_landbeta_amp$psrf[,1], 
          GR_landbeta_phx$psrf[,1], 
          GR_landbeta_phy$psrf[,1]))


# acfplot(trends[[1]])
# acfplot(trend_sizes[[1]])
# acfplot(amps[[1]])
# acfplot(amp_sizes[[1]])

ess_trend = effectiveSize(trends)
ess_amp = effectiveSize(amps)
ess_ph = effectiveSize(phs)
ess_trendsize = effectiveSize(trend_sizes)
ess_ampsize = effectiveSize(amp_sizes)
ess_phsize = effectiveSize(ph_sizes)
ess_tcut = effectiveSize(tcuts)
ess_acut = effectiveSize(acuts)
ess_tcut = effectiveSize(tcuts)



pdf(file=paste0("plots/traceplots_", modK, "_anch", n.fixed, "_seed", seed[1], ".pdf"), width=7, height=4)

# traceplot(trends, main=paste("trend"))
for(i in 1:(K.time-1)) {
  plot(as.matrix(trends[[1]][,i]), type="l", lty=2, ylab="trend", xlab="iteration", ylim=c(min(sapply(trends, function(x) x[,i])), max(sapply(trends, function(x) x[,i]))), 
       main=paste("Trend", ifelse(i >= indx_trend0, i+1, i), "\nGelman-Rubin:", round(GR_trend$psrf[i,1], 3), "  Effective size:", round(ess_trend[i])))
  for(j in 2:nchains) {
    lines(as.matrix(trends[[j]][,i]), lty=2, col=j)
  }
}
traceplot(trend_sizes, main="trend cluster size")

# traceplot(amps, main="amplitude")
for(i in 1:(K.amp)) {
  plot(as.matrix(amps[[1]][,i]), type="l", lty=2, ylab="amplitude", xlab="iteration", ylim=c(min(sapply(amps, function(x) x[,i])), max(sapply(amps, function(x) x[,i]))), 
       main=paste("Amplitude", i, "\nGelman-Rubin:", round(GR_amp$psrf[i,1], 3), "  Effective size:", round(ess_amp[i])))
  for(j in 2:nchains) {
    lines(as.matrix(amps[[j]][,i]), lty=2, col=j)
  }
}
traceplot(amp_sizes, main="amplitude cluster size")

# traceplot(phs, main="phase")
for(i in 1:(K.phase)) {
  plot(as.matrix(phs[[1]][,i]), type="l", lty=2, ylab="phase", xlab="iteration", ylim=c(min(sapply(phs, function(x) x[,i])), max(sapply(phs, function(x) x[,i]))), 
       main=paste("Phase:", month.name[c(7:12, 1:6)][i], "\nGelman-Rubin:", round(GR_ph$psrf[i,1], 3), "  Effective size:", round(ess_ph[i])))
  for(j in 2:nchains) {
    lines(as.matrix(phs[[j]][,i]), lty=2, col=j)
  }
}
traceplot(ph_sizes, main="phase cluster size")

# traceplot(tcuts, main="cut points: trend")
for(i in 1:ncol(tcuts[[1]])) {
  plot(as.matrix(tcuts[[1]][,i]), type="l", lty=2, ylab="cut point", xlab="iteration", ylim=c(min(sapply(tcuts, function(x) x[,i])), max(sapply(tcuts, function(x) x[,i]))), 
       main=paste("Trend cut point:", i + 2, "\nGelman-Rubin:", round(GR_tcut$psrf[i,1], 3), "  Effective size:", round(ess_tcut[i])))
  for(j in 2:nchains) {
    lines(as.matrix(tcuts[[j]][,i]), lty=2, col=j)
  }
}

# traceplot(acuts, main="cut points: amplitude")
for(i in 1:ncol(acuts[[1]])) {
  plot(as.matrix(acuts[[1]][,i]), type="l", lty=2, ylab="cut point", xlab="iteration", ylim=c(min(sapply(acuts, function(x) x[,i])), max(sapply(acuts, function(x) x[,i]))), 
       main=paste("Amplitude cut point:", i + 2, "\nGelman-Rubin:", round(GR_acut$psrf[i,1], 3), "  Effective size:", round(ess_acut[i])))
  for(j in 2:nchains) {
    lines(as.matrix(acuts[[j]][,i]), lty=2, col=j)
  }
}

dev.off()

pdf(file=paste0("plots/traceplots_landbetas_trend_", modK, "_anch", n.fixed, "_seed", seed[1], ".pdf"), width=7, height=4)
traceplot(landbetas_trend)
dev.off()

pdf(file=paste0("plots/traceplots_landbetas_amp_", modK, "_anch", n.fixed, "_seed", seed[1], ".pdf"), width=7, height=4)
traceplot(landbetas_amp)
dev.off()

pdf(file=paste0("plots/traceplots_landbetas_phx_", modK, "_anch", n.fixed, "_seed", seed[1], ".pdf"), width=7, height=4)
traceplot(landbetas_phx)
dev.off()

pdf(file=paste0("plots/traceplots_landbetas_phy_", modK, "_anch", n.fixed, "_seed", seed[1], ".pdf"), width=7, height=4)
traceplot(landbetas_phy)
dev.off()
