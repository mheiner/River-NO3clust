rm(list=ls())

allrdafiles = list.files("results/")
indx_files_use = grep("anchor2_K5-7-12_KMtypesupplied_CV5_2222\\d", allrdafiles)
# indx_files_use = grep("anchor2_K7-11-12_KMtypesupplied_CV5_2222\\d", allrdafiles)
(rdafiles_use = allrdafiles[indx_files_use])
(nchains = length(rdafiles_use))


get_lppo <- function(rn) {
  n_obs = length(river.list[[rn]]$df$NO3)
  n_test = length(river.list.cv[[rn]]$df$NO3)
  niter = nrow(int.draws)
  
  mean_at_test = matrix( rep(int.draws[,rn], times=n_test), ncol=n_test ) + 
    tcrossprod(time.draws[ cbind(1:niter, tclust.draws[,rn]) ], river.list.cv[[rn]]$df$time) +
    (matrix( rep(amp.draws[ cbind(1:niter, aclust.draws[,rn])], times=n_test), ncol=n_test ) * 
       cos(t(outer(2*pi*river.list.cv[[rn]]$df$time, phase.draws[ cbind(1:niter, pclust.draws[,rn]) ], FUN="-")))
    )
  
  mean_at_obs = matrix( rep(int.draws[,rn], times=n_obs), ncol=n_obs ) + 
    tcrossprod(time.draws[ cbind(1:niter, tclust.draws[,rn]) ], river.list[[rn]]$df$time) +
    (matrix( rep(amp.draws[ cbind(1:niter, aclust.draws[,rn])], times=n_obs), ncol=n_obs ) * 
       cos(t(outer(2*pi*river.list[[rn]]$df$time, phase.draws[ cbind(1:niter, pclust.draws[,rn]) ], FUN="-")))
    )
  
  e_obs = matrix(rep(log(river.list[[rn]]$df$NO3), each=niter), nrow=niter) - mean_at_obs # niter by n_obs matrix
  
  Rxx_chol = river.list[[rn]]$Rchol
  Rxy = river.list.full[[rn]]$R[-indx_cv[[rn]], indx_cv[[rn]]]
  Ryy = river.list.cv[[rn]]$R
  
  Ry_given_x = Ryy - crossprod(forwardsolve(Rxx_chol, Rxy))
  Ry_given_x_chol = chol(Ry_given_x) %>% t()
  
  AA = forwardsolve(Rxx_chol, t(e_obs)) # n_obs by niter matrix
  BB = crossprod(Rxy, AA) %>% t() # niter by n_test
  
  mean_y = mean_at_test + BB # niter by n_test
  e_test = matrix(rep(log(river.list.cv[[rn]]$df$NO3), each=niter), nrow=niter) - mean_y # niter by n_test
  
  zz = forwardsolve(Ry_given_x_chol, t(e_test)) %>% t() / matrix( rep(sqrt(sig2.draws[,rn]), times=n_test), ncol=n_test ) # niter by n_test
  
  lpdens = dnorm(zz, mean=0, sd=1, log=TRUE)
  
  sum(colMeans(lpdens))
}

## for viewing one-at-a-time
# plot(river.list[[rn]]$df$time, log(river.list[[rn]]$df$NO3), type="o")
# points(river.list[[rn]]$df$time, colMeans(mean_at_obs), col="red", type="o")
# points(river.list.cv[[rn]]$df$time, colMeans(mean_at_test), col="red", pch=1)
# points(river.list.cv[[rn]]$df$time, colMeans(mean_y), col="red", pch=19)
# points(river.list.cv[[rn]]$df$time, log(river.list.cv[[rn]]$df$NO3), col="black", pch=4)
# colMeans(e_test)[order(river.list.cv[[rn]]$df$time)]
# colMeans(zz)[order(river.list.cv[[rn]]$df$time)]

lppo = list()

i = 1
for( i in 1:nchains ) {
  load(paste0("results/", rdafiles_use[i]))
  print(indx_cv[[1]])
  
  print(it)
  print(it_prev)
  
  Nst = length(river.list)
  lppo[[i]] = numeric(Nst)
  for (rn in 1:Nst) {
    lppo[[i]][rn] = get_lppo(rn)
    cat(rn, "of", Nst, "\r")
  }
  
  # hist(lppo)
  # hist(lppo[abs(lppo) < 50])
  
  print(summary(lppo[[i]]))
  print(median(lppo[[i]]))
  print(sum(lppo[[i]]))
  # mean(lppo)
  # median(lppo)
}

lppo_57 = lppo

## optionally go back and get the lppo for a competing model, save, and combine for comparison between models
lppo_711 = lppo

# save(lppo_57, lppo_711, file="results/cv_2222.rda")

load("results/cv_2222.rda")

i = 1
(nchains = length(lppo_57))
DD = matrix(NA, ncol=nchains, nrow=length(lppo_57[[1]]))

for( i in 1:nchains ) {
  DD[,i] = lppo_711[[i]] - lppo_57[[i]] # positive means higher density for 5-7 model
}

i = 1
hist(DD[,i])
hist(DD[which(abs(DD[,i])<2),i])
hist(DD[which(abs(DD[,i])<1),i])
hist(DD[which(abs(DD[,i])<0.5),i])
hist(DD[which(abs(DD[,i])<0.2),i])
# apply(DD, 2, t.test)
colSums(DD)
colMeans(DD)
summary(DD)
apply(DD, 2, median)

i = 1
mean(DD[,i] > 0)
sum(DD[,i] > abs(min(DD[,i])))
