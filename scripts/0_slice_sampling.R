
slice_shrink = function(x, ly, lf, L, R, maxiter=5000) {
  ## Algorithm in Fig 5 of Neal (2003)
  
  # ly = lf(x) + log(runif(1, 0.0, 1.0))
  
  U = runif(1, 0.0, 1.0)
  x1 = L + U * (R - L)
  i = 1
  
  while (lf(x1) < ly && i <= maxiter) {
    
    if (x1 < x) {
      L = x1
    } else {
      R = x1
    }
    
    U = runif(1, 0.0, 1.0)
    x1 = L + U * (R - L)
    
    i = i + 1
  }
  
  if(i > 5000) print("Slice shrink exceeded 5000 steps.")
  x1
}


slice_step_out = function(x, ly, lf, w, maxiter=5000) {
  ## Algorithm in Figs 3 and 5 of Neal (2003)
  
  # ly = lf(x) + log(runif(1, 0.0, 1.0))

  U = runif(1, 0.0, 1.0)
  L = x - w * U
  R = L + w
  i = 1
  
  while (ly < lf(L) && i <= maxiter) {
    L = L - w
    i = i + 1
  }
  
  while (ly < lf(R) && i <= maxiter) {
    R = R + w
    i = i + 1
  }

  if(i > maxiter) print("Slice step out exceeded 5000 steps.")
  slice_shrink(x, ly, lf, L, R)
}



slice_shrink_hyperrect = function(x, ly, lf, w, maxiter=10000) {
  ## Algorithm in Fig 8 of Neal (2003)
  
  # ly = lf(x) + log(runif(1, 0.0, 1.0))
  
  KK = length(x)
  
  ## randomly position hyperrectangle
  U = runif(KK, 0.0, 1.0)
  L = x - w*U
  R = L + w
  
  ## sample from hyperrectangle
  U = runif(KK, 0.0, 1.0)
  x1 = L + U * (R - L)
  i = 1
  
  while (lf(x1) < ly && i <= maxiter) {
    
    for (kk in 1:KK) {
      if (x1[kk] < x[kk]) {
        L[kk] = x1[kk]
      } else {
        R[kk] = x1[kk]
      }
    }
    
    U = runif(KK, 0.0, 1.0)
    x1 = L + U * (R - L)
    
    i = i + 1
  }
  
  if(i > 5000) print("Slice shrink exceeded 5000 steps.")
  x1
}



### test
# a = 5
# b = 2
# lf = function(x) dbeta(x, a, b, log=TRUE)
# 
# n = 10000
# theta = numeric(n)
# theta[1] = 0.5
# for (i in 2:n) {
#   ly = lf(theta[i-1]) + log(runif(1))
#   theta[i] = slice_step_out(theta[i-1], ly, lf, 0.1)
#   cat("i =", i, "of", n, "\r")
# }
# plot(theta, type="l")
# hist(theta, freq=FALSE)
# curve(exp(lf(x)), col="red", add=TRUE)

## multivariate test
# Rho = matrix(c(1.0, 0.8, 0.8, 1.0), nrow=2)
# sig = c(1.0, 2.0)
# Sig = diag(sig) %*% Rho %*% diag(sig)
# SigL = t(chol(Sig))
# mu = c(0.0, 3.0)
# 
# lf = function(x) {
#   aa = forwardsolve(SigL, (x - mu))
#   -0.5 * drop(crossprod(aa))
# } 
# 
# n = 10000
# theta = matrix(0, nrow=n, ncol=2)
# for (i in 2:n) {
#   ly = lf(theta[i-1,]) + log(runif(1))
#   theta[i,] = slice_shrink_hyperrect(theta[i-1,], ly, lf, w=c(2.5, 2.5))
#   cat("i =", i, "of", n, "\r")
# }
# YY = t(SigL %*% matrix(rnorm(n*2), nrow=2) + rep(mu, times=n))
# 
# plot(theta[,2], type='l')
# 
# plot(theta[,1], theta[,2])
# points(YY[,1], YY[,2], col="red", pch=".")
