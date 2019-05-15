# Init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')




# create data -------------------------------------------------------------

set.seed(20)

n <- 300


# a <- runif(n, 10, 55) + rnorm(n, 0, 2)
e1 <- rnorm(n, 0, 3)
e2 <- rnorm(n, 0, 3)

a <- rlnorm(n, meanlog=log(25), sdlog=log(1.3)) + e1
b <- ((a-e1) * 1) + e2

dev.off()
par(pty='s', mfrow=c(1,2), las=1, tck=-0.02, mgp=c(2,0.5,0), mar=c(3,3,1,1))

# B vs. A
fit <- lm(b ~ a-1)

plot(a, b, asp=1, xlim=c(5,60), ylim=c(5,60), xlab='A', ylab='B')
abline(a=fit$coefficients[1], b=fit$coefficients[1], col='red', lwd=2)
#abline(a=0, b=1, col='blue', lwd=2, lty=2)
text(par('usr')[1]+2, par('usr')[4]-2, adj=c(0, 1), cex=1.5,
     paste0('beta = ', formatC(fit$coefficients[1], digits=4)))

# A vs. B
fit <- lm(a ~ b-1)

plot(b, a, asp=1, xlim=c(5,60), ylim=c(5,60), xlab='B', ylab='A')
abline(a=fit$coefficients[1], b=fit$coefficients[1], col='red', lwd=2)
#abline(a=0, b=1, col='blue', lwd=2, lty=2)
text(par('usr')[2]-2, par('usr')[3]+2, adj=c(1, 0), cex=1.5,
     paste0('beta = ', formatC(fit$coefficients[1], digits=4)))

# sigmoidal ---------------------------------------------------------------


set.seed(1)

n <- 1000

a <- rlnorm(n, meanlog=log(30), sdlog=log(1.1)) + rnorm(n, 0, 1.5)
b <- (sigmoid(a, a=0.5, b=30) * max(a)) + rnorm(n, 0, 1.5)

dev.off()
plot(a, b, xlab='A', ylab='B')
