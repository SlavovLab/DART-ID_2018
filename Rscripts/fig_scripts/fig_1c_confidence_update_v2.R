library(tidyverse)
library(rmutil)
source('Rscripts/lib.R')

## Peptide Update ------

pdf(file='manuscript/Figs/confidence_update_v10.pdf', width=3.5, height=2.5)

layout(rbind(c(1, 2),
             c(1, 3)),
       widths=c(3.5, 1))

par(oma=c(2,0,0,0.75),
    mar=c(0.5,3,0.5,0.75),
    xaxs='i', yaxs='i', pty='m', cex.axis=1)

plot(0, 0, type='n',
     xlim=c(20, 25), ylim=c(-0.15, 1.95),
     xaxt='n', yaxt='n',
     #xlab="Retention Time (min)", ylab="Density",
     #main='Confidence Update',
     xlab=NA, ylab=NA)

mu <- 25.25
exp <- 0.95

denx <- seq(20, 25, length.out=1000)
dist_sd <- 0.275
#deny <- dnorm(denx, mean=mu*exp, sd=dist_sd)
deny <- dlaplace(denx, m=mu*exp, s=dist_sd)
polygon(c(denx, 25), c(deny, 0), 
        col=rgb(0, 0, 1, 0.3))
lines(denx[400:1000], deny[400:1000], col='blue', lwd=2)

#segments(x0=mu*exp, x1=mu*exp,
#         y0=0, y1=dnorm(mu*exp, mean=mu*exp, sd=dist_sd),
#         col='blue', lwd=2, lty=2)
abline(h=0, col='black', lwd=2)

# null distribution
den_null <- (dnorm(denx, mean=22, sd=5) * 1)
polygon(c(0, denx,25), c(0, den_null,0), 
        col=rgb(1, 0, 0, 0.3), border=NA)
lines(denx, den_null, col='red', lwd=2)

# points
x1 <- mu*exp-2
x2 <- mu*exp+0.3
points(c(x1, x2), rep(0,2), pch=c(21, 22), bg='black', col='black', lwd=2, cex=2)

#text(x=mu*exp,y=1.65, labels="Peptide X\nExperiment A", font=1, cex=1.2)

axis(side=1, tck=-0.02, mgp=c(0, 0.3, 0))
axis(side=2, tck=-0.02, las=1, mgp=c(0, 0.4, 0))

mtext('Retention time (min)', side=1, line=1.5, cex=1)
mtext('Probability density', side=2, line=1.85, cex=1, las=3)

par(lheight=0.85)
text(20.4, 1.7, 'Conditional RT\nLikelihood', adj=c(0, 0.5), cex=1)
# legend(20.25, 1.6, 
#        c(expression(atop('ID correct', delta=1)), 
#          expression('ID incorrect\n'*delta)), 
#        col=c('blue', 'red'), xjust=0, yjust=1,
#        pch=22, pt.bg=c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), pt.lwd=2, pt.cex=3,
#        lty=1, lwd=NA, bty='n', cex=1, y.intersp=1.5, inset=c(0,0), adj=c(0, 0.5))
text(21.25, 1.15, expression('ID Correct\n', (delta==1)), adj=c(0, 0.5))
text(21.25, 0.65, expression('ID Incorrect\n', (delta==0)), adj=c(0, 0.5))
points(rep(20.75, 2), c(1.25, 0.75), col=c('blue', 'red'), 
       pch=22, bg=c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), cex=3, lwd=2)

par(#oma=c(0, 0, 1, 0),
    mar=c(0.25, 1, 1.5, 0.5),
    pty='m')

plot(0, 0, type='n',
     xlim=c(0, 2.75), ylim=c(0, 3.5),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(-log10(c(5e-2, 0.3)), width=1, space=0.25, add=T,
        col=c('white', 'black'),
        xaxt='n', yaxt='n')
#arrows(x0=0.6, x1=2, y0=2, y1=1, col='black', length=0.075, code=2)

mtext('   PSM 1 [   ]', side=3, line=0.15, cex=0.75, font=2)
points(2.775, 3.925, pch=21, cex=1.25, xpd=T,bg='black', col='black', lwd=1.5)
mtext('ID Confidence', 2, line=0.35, cex=1, at=-0.5, las=3)

par(#oma=c(1,0,0,0),
    mar=c(0.5, 1, 1.25, 0.5),
    pty='m')

plot(0, 0, type='n',
     xlim=c(0, 2.75), ylim=c(0, 3.5),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(-log10(c(5e-2, 3e-3)), width=1, space=0.25, add=T,
        col=c('white', 'black'),
        xaxt='n', yaxt='n')
#arrows(x0=0.6, x1=1.2, y0=1.5, y1=2.5, col='black', length=0.075, code=2)

legend(x=-0.75, y=0, c('Spectra', 'DART-ID'), 
       pch=22, pt.cex=2, pt.bg=c('white', 'black'), cex=0.9, col='black',
       y.intersp=1, bty='n', xpd=NA)

mtext('   PSM 2 [   ]', side=3, line=0.15, cex=0.75, font=2)
points(2.775, 3.925, pch=22, cex=1.25, xpd=T, bg='black', col='black', lwd=1.5)

dev.off()


# for SCP 2019 ------------------------------------------------------------


# prediction - wide peaks -------------------------------------------------

# png(file='manuscript/Figs/SCP_2019_conf_update_prediction.png', width=3.5, height=2.5, units='in', res=400)
png(file='manuscript/Figs/SCP_2019_conf_update_alignment.png', width=3.5, height=2.5, units='in', res=400)

par(oma=c(2,0,0,0.75),
    mar=c(0.5,3,0.5,0.75),
    xaxs='i', yaxs='i', pty='m', cex.axis=1)

plot(0, 0, type='n',
     xlim=c(20, 25), 
#     ylim=c(-0.05, 0.4), # prediction
     ylim=c(-0.05, 1), # alignment
     xaxt='n', yaxt='n',
     xlab=NA, ylab=NA)

mu <- 25.25
exp <- 0.95

denx <- seq(10, 100, length.out=10000)

#dist_sd <- 2 # prediction
dist_sd <- 0.7 # alignment

deny <- dlaplace(denx, m=mu*exp, s=dist_sd)
polygon(c(denx, 25), c(deny, 0), 
        col=rgb(0, 0, 1, 0.3))
lines(denx[1:10000], deny[1:10000], col='blue', lwd=2)

#segments(x0=mu*exp, x1=mu*exp,
#         y0=0, y1=dnorm(mu*exp, mean=mu*exp, sd=dist_sd),
#         col='blue', lwd=2, lty=2)
abline(h=0, col='black', lwd=2)

# null distribution
den_null <- (dnorm(denx, mean=22, sd=5) * 1)
polygon(c(0, denx,25), c(0, den_null,0), 
        col=rgb(1, 0, 0, 0.3), border=NA)
lines(denx, den_null, col='red', lwd=2)

# points
# x1 <- mu*exp-2
# x2 <- mu*exp+0.3
# points(c(x1, x2), rep(0,2), pch=c(21, 22), bg='black', col='black', lwd=2, cex=2)

#text(x=mu*exp,y=1.65, labels="Peptide X\nExperiment A", font=1, cex=1.2)

axis(side=1, tck=-0.02, mgp=c(0, 0.3, 0))
axis(side=2, tck=-0.02, las=1, mgp=c(0, 0.4, 0))

mtext('Retention time (min)', side=1, line=1.5, cex=1)
mtext('Probability density', side=2, line=1.85, cex=1, las=3)

par(lheight=0.85)
print(par('usr'))
text(par('usr')[1] + (par('usr')[2] - par('usr')[1])*0.05, 
     par('usr')[4] - (par('usr')[4] - par('usr')[3])*0.05, 
     'Conditional RT Likelihood', adj=c(0, 1), cex=1)
# text(21.25, 
#      1.15, 
#      expression('ID Correct\n', (delta==1)), adj=c(0, 0.5))
# text(21.25, 0.65, expression('ID Incorrect\n', (delta==0)), adj=c(0, 0.5))
# points(rep(20.75, 2), c(1.25, 0.75), col=c('blue', 'red'), 
#        pch=22, bg=c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), cex=3, lwd=2)

dev.off()
