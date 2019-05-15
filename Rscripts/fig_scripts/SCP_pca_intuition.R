# Init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')


# Create data -------------------------------------------------------------

set.seed(1)

exps <- seq(1,50)
beta_0 <- rnorm(length(exps), 0, 3)
beta_1 <- rnorm(length(exps), 1, 0.2)

num_peps <- 1000
# mus <- runif(num_peps, 1, 60)
mus <- rlnorm(num_peps, meanlog=log(32), sdlog=log(1.3))
mu_error <- rnorm(length(mus), 0, 2)

dmat <- ((mus) %*% t(beta_1)) + # slope
  matrix(nrow=num_peps, ncol=length(exps), rep(beta_0, num_peps), byrow=T) + # intercept
  matrix(nrow=num_peps, ncol=length(exps), rnorm(num_peps*length(exps), 0, 2.5)) # noise

plot(dmat[,1], dmat[,2])


# covariance --------------------------------------------------------------

# center around 0
cmat <- dmat - apply(dmat, 1, mean)

# covariance matrix
sim <- t(cmat) %*% cmat

# eigenvalue decomposition
p <- svd(sim)

plot(p$u[,1], p$u[,2])

plot(seq(1,length(exps)), p$d / sum(p$d))


# PCs to betas ------------------------------------------------------------

par(mfrow=c(2,2))

plot(p$u[,1], beta_0, xlab='PC1', ylab='beta_0')
plot(p$u[,1], beta_1, xlab='PC1', ylab='beta_1')

plot(p$u[,2], beta_0, xlab='PC2', ylab='beta_0')
plot(p$u[,2], beta_1, xlab='PC2', ylab='beta_1')


# map loadings back to data, using 1st PC ---------------------------------

# use_pcs <- c(1, 2, 3)
use_pcs <- c(1)

# Xhat = XVV^T + mu
# where X = n x p centered data matrix
#       V = p x k matrix of eigenvectors after decomposition
#       mu = vector of row means
xhat <- (cmat %*% p$u[,use_pcs] %*% t(p$u[,use_pcs])) + apply(dmat, 1, mean)

# plot data, and projection using PC1 --------------------------------------

dev.off()

set.seed(22)

runs <- sample.int(50, size=8)

par(mfrow=c(2,2), pty='s', las=1, mar=c(3.2,2.2,1,1), mgp=c(2,0.5,0), tck=-0.02)

for(i in 1:(length(runs) / 2)) {
  a <- runs[((i-1)*2)+1]
  b <- runs[((i-1)*2)+2]
  
  plot(dmat[,a], dmat[,b], 
       xlab=paste0('Run ', a), 
       ylab=paste0('Run ', b), asp=1)
  
  points(xhat[,a], xhat[,b], col='red')
}

# dev.off()


# plot RTs vs. PC1 on a single plot ---------------------------------------

dev.off()

set.seed(22)

runs <- sample.int(50, size=8)

# par(mfrow=c(2,2), pty='s', las=1, mar=c(3.2,2.2,1,1), mgp=c(2,0.5,0), tck=-0.02)

plot(xhat[,1], dmat[,1])
points(xhat[,2], dmat[,2])
abline(a=0, b=1, col='red')

# for(i in 1:(length(runs) / 2)) {
#   a <- runs[((i-1)*2)+1]
#   b <- runs[((i-1)*2)+2]
#   
#   plot(dmat[,a], dmat[,b], 
#        xlab=paste0('Run ', a), 
#        ylab=paste0('Run ', b), asp=1)
#   
#   points(xhat[,a], xhat[,b], col='red')
# }
