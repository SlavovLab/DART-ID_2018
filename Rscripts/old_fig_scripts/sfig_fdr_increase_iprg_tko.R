# init, load data ---------------------------------------------------------

library(tidyverse)
source('Rscripts/lib.R')

iprg <- read_tsv('/gd/bayesian_RT/Alignments/iPRG2015_20181228/ev_updated.txt')
tko <- read_tsv('/gd/bayesian_RT/Alignments/PXD011654_20190120/ev_updated.txt')

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

legend('topright', c('Label-free', 'TMT-labelled', 'Spectra'), col=c('red', 'blue', rgb(0.5,0.5,0.5)),
       lwd=3, lty=c(1, 1, 2), seg.len=0.8, cex=0.85, x.intersp=1, y.intersp=1, bg='white')

dev.off()
