
# nchains = 4

trends = list()
tclusts = list()

amps = list()
aclusts = list()

phs = list()
pclusts = list()

tcuts = list()
acuts = list()

landbetas_trend = list()
landbetas_amp = list()
landbetas_phx = list()
landbetas_phy = list()

intcpts = list()
sig2s = list()

lliks = list()
lliks_site = list()


for (ii in 1:nchains) {
  
  ### select runs  
  # load(paste0("Results_anchor2_K5-7-12_KMtypesupplied_2222", ii, ".RData"))
  # load(paste0("Results_anchor2_K5-7-12_KMtypesupplied_2222", ii, "_cont1.RData"))
  # load(paste0("Results_anchor2_K5-7-12_KMtypesupplied_2222", ii, "_cont2.RData"))
  # load(paste0("Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, ".RData"))
  # load(paste0("Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, "_cont1.RData"))
  # load(paste0("Results_anchor2_K7-11-12_KMtypesupplied_2222", ii, "_cont2.RData"))
  
  load(paste0("results/", rdafiles_use[ii]))
  
  ###
  trends[[ii]] = time.draws
  tclusts[[ii]] = tclust.draws
  
  amps[[ii]] = amp.draws
  aclusts[[ii]] = aclust.draws
  
  phs[[ii]] = phase.draws
  pclusts[[ii]] = pclust.draws
  
  tcuts[[ii]] = time.cuts.draws
  acuts[[ii]] = amp.cuts.draws
  
  landbetas_trend[[ii]] = tclust.beta.draws
  landbetas_amp[[ii]] = aclust.beta.draws
  landbetas_phx[[ii]] = pclust.beta.x.draws
  landbetas_phy[[ii]] = pclust.beta.y.draws
  
  intcpts[[ii]] = int.draws
  sig2s[[ii]] = sig2.draws
  
  lliks[[ii]] = llike.draws
  lliks_site[[ii]] = llike.site.draws
  
  cat(ii, "of", nchains, "\n")
}

time.draws = do.call(rbind, trends)
tclust.draws = do.call(rbind, tclusts)
amp.draws = do.call(rbind, amps)
aclust.draws = do.call(rbind, aclusts)
phase.draws = do.call(rbind, phs)
pclust.draws = do.call(rbind, pclusts)
time.cuts.draws = do.call(rbind, tcuts)
amp.cuts.draws = do.call(rbind, acuts)
tclust.beta.draws = do.call(rbind, landbetas_trend)
aclust.beta.draws = do.call(rbind, landbetas_amp)
pclust.beta.x.draws = do.call(rbind, landbetas_phx)
pclust.beta.y.draws = do.call(rbind, landbetas_phy)
int.draws = do.call(rbind, intcpts)
sig2.draws = do.call(rbind, sig2s)
llike.draws = do.call(rbind, lliks)
llike.site.draws = do.call(rbind, lliks_site)
