# init, load data ---------------------------------------------------------

library(tidyverse)
source('Rscripts/lib.R')

iprg <- read_tsv('/gd/bayesian_RT/Alignments/iPRG2015_20181228/ev_updated.txt')
tko <- read_tsv('/gd/bayesian_RT/Alignments/PXD011654_20190120/ev_updated.txt')


# alignment residuals -----------------------------------------------------

pdf(file='manuscript/Figs/residual_iprg.pdf', width=3.5, height=3)

par(mar=c(3,3,3,2), cex.axis=1, cex.main=1, mgp=c(1,3,1))

hist(iprg$residual[abs(iprg$residual) < 1], breaks=seq(-1,1,by=0.02), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
axis(1, at=seq(-1, 1, by=0.5), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 10, by=2), tck=-0.02, mgp=c(0,0.5,0), las=1)
mtext('Residual RT (min)', side=1, line=1.3)
mtext('Density', side=2, line=1.5)
mtext('iPRG 2015 (Label-free)', side=3, line=0.7, font=2)

dev.off()

pdf(file='manuscript/Figs/residual_tko.pdf', width=3.5, height=3)

par(mar=c(3,4,3,2), cex.axis=1, cex.main=1, mgp=c(1,3,1))

hist(tko$residual[abs(tko$residual) < 5], breaks=seq(-5,5,by=0.1), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
axis(1, at=seq(-5, 5, by=1), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 1.5, by=0.25), tck=-0.02, mgp=c(0,0.5,0), las=1)
mtext('Residual RT (min)', side=1, line=1.3)
mtext('Density', side=2, line=2.5)
mtext('TKO 2018 (TMT)', side=3, line=0.7, font=2)

dev.off()



# % increase in IDs -------------------------------------------------------

# calculate q-values
iprg <- iprg %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 1:nrow(iprg))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 1:nrow(iprg))[order(order(pep_updated))])

tko <- tko %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 1:nrow(tko))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 1:nrow(tko))[order(order(pep_updated))])

x <- seq(-4, -1, length.out=100)

# matrix - iprg, tko
dat <- matrix(nrow=length(x), ncol=2)

for(i in 1:length(x)) {
  # percent increase in IDs relative to spectral FDR
  dat[i,] <- c(
    (sum(iprg$qval_updated < (10 ** x[i])) / sum(iprg$qval < (10 ** x[i]))-1) * 100,
    (sum(tko$qval_updated < (10 ** x[i]))  / sum(tko$qval < (10 ** x[i]))-1) * 100
  )  
}

# plot --------------------------------------------------------------------

pdf(file='manuscript/Figs/fdr_increase_iprg_tko.pdf', width=3.5, height=2.5)

par(mar=c(2.5,3,0.4,1.2), cex.axis=1)

plot(0, 0, type='n', xlim=c(-4, -1), ylim=c(-5, 110),
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n', xaxs='i')

lines(x, rep(0, length(x)), col=rgb(0.5,0.5,0.5), lwd=2) # baseline
abline(v=-2, col='black', lty=2, lwd=2) # 1% FDR
lines(x, dat[,1], col='red', lwd=2) # iPRG
lines(x, dat[,2], col='blue', lwd=2) # TKO

axis(1, at=seq(-4, -1, by=1), labels=fancy_scientific(10 ** seq(-4, -1, by=1)),
     tck=-0.02, mgp=c(0, 0.3, 0))
p_labels <- seq(0, 110, by=10)
p_labels[seq(2,12,by=2)] <- NA
axis(2, at=seq(0, 110, by=10), labels=p_labels, tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext('FDR Threshold', side=1, line=1.35)
mtext('% Increase in PSMs', side=2, line=2)

legend('topright', c('iPRG 2015 (LF)', 'TKO 2018 (TMT)', 'Spectra'), col=c('red', 'blue', rgb(0.5,0.5,0.5)),
       lwd=3, lty=c(1, 1, 2), seg.len=0.8, cex=0.85, x.intersp=1, y.intersp=1, bg='white')

dev.off()


# Increase in PSMs - Scatter ----------------------------------------------

iprg_exp <- iprg %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 1:nrow(iprg))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) /  1:nrow(iprg))[order(order(pep_updated))]) %>%
  group_by(`Raw file`) %>%
  summarise(spectra=sum(qval < 1e-2),
            dart=sum(qval_updated < 1e-2))

tko_exp <- tko %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 1:nrow(tko))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) /  1:nrow(tko))[order(order(pep_updated))]) %>%
  group_by(`Raw file`) %>%
  summarise(spectra=sum(qval < 1e-2),
            dart=sum(qval_updated < 1e-2))


# plot --------------------------------------------------------------------

pdf(file='manuscript/Figs/peps_per_exp_scatter_iprg_tko.pdf', width=3.5, height=3.5)
layout(c(1,2))

par(mar=c(1.5,0,0,0), oma=c(0.8, 3.5, 1.8, 1.5), las=1, cex.axis=0.85)

plot(0, 0, type='n', xlim=c(20000,30000), ylim=c(20000, 30000),
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')

points(iprg_exp$spectra, iprg_exp$dart, pch=1, col='red')
abline(a=0, b=1, lwd=2, lty=1, col='black')
text(25000, 20500, 'iPRG 2015 (LF)', adj=c(0, 0), font=1, cex=1)
axis(1, at=seq(20000, 30000, by=2000), labels=formatC(seq(20, 30, by=2), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.1, 0))
axis(2, at=seq(20000, 30000, by=2000), labels=formatC(seq(20, 30, by=2), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.3, 0), las=1)

par(mar=c(1.25,0,0,0))

plot(0, 0, type='n', xlim=c(4500,14500), ylim=c(4500, 14500),
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
points(tko_exp$spectra, tko_exp$dart, pch=1, col='blue')
abline(a=0, b=1, lwd=2, lty=1, col='black')
text(9200, 5000, 'TKO 2018 (TMT)', adj=c(0, 0), font=1, cex=1)
axis(1, at=seq(4500, 14500, by=2500), labels=formatC(seq(4.5, 14.5, by=2.5), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.1, 0))
axis(2, at=seq(4500, 14500, by=2500), labels=formatC(seq(4.5, 14.5, by=2.5), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.3, 0), las=1)

mtext('Spectra', side=1, outer=T, line=-0.25)
mtext('DART-ID', side=2, outer=T, line=2, las=3)
mtext('(PSMs x 1000) per Experiment, 1% FDR      ', side=3, outer=T, line=0.4, font=2)

dev.off()

