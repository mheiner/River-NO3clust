fast.dmvnorm <- function(x, mu, Sigma,log=FALSE){
  
  low.tri <- t(chol(Sigma))
  u <- forwardsolve(low.tri,t(scale(x, center=mu, scale=FALSE)))
  pdfout <- -(nrow(Sigma)/2)*log(2*pi) - sum(log(diag(low.tri))) -0.5*colSums(u^2)
  if(log){
    return(pdfout)
  } else {
    return(exp(pdfout))
  }

}

######################################
## Example using the fast.dmvnorm() ##
######################################

## Set mu and Sigma
mu <- matrix(c(1,17), ncol=1, nrow=2)
Sig <- matrix(c(1, .9, .9, 1), nrow=2)

## Set Sample size
n <- 1000

## Simulate data from mvnorm(mu, Sigma)
simData <- scale(matrix(rnorm(2*n), ncol=2)%*%chol(Sig), center=c(-mu), scale=FALSE)
colMeans(simData)
var(simData)

## Evaluate Density of MVN using fast.dmvnorm()
system.time(dens <- fast.dmvnorm(simData, mu, Sig, log=FALSE))
system.time(dens2 <- mvtnorm::dmvnorm(simData, mu, Sig, log=FALSE))
cbind(dens, dens2)
