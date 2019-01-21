# init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

# load sqc w/ reverse seqs ------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_include_decoy_20190109/ev_updated.txt')
# load data with PSM FDR specified at 1% in MaxQuant
ev_001 <- read_tsv('/gd/bayesian_RT/dat_SQC_FDR_1/evidence.txt')

ev_001_f <- ev_001 %>%
  dplyr::select(c('Modified sequence', 'Raw file', 'Leading razor protein', 
                  'PEP', 'MS/MS scan number')) %>%
  # ceil PEPs to 1
  mutate_at(c('PEP'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev_001)))[order(order(PEP))])

ev_f <- ev %>%
  dplyr::select(c('Modified sequence', 'Raw file', 'Leading razor protein', 
                  'PEP', 'pep_updated', 'pep_new', 'residual', 'MS/MS scan number')) %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))])


# decoy sequences with aligned RTs? ---------------------------------------

sum(grepl('REV_', ev_f$`Leading razor protein`) & !is.na(ev_f$pep_new))

# frequency of decoy sequences

rev <- ev_f[grepl('REV_', ev_f$`Leading razor protein`),]

a <- as.data.frame(table(rev$`Leading razor protein`))
a <- a %>% arrange(desc(Freq))

dpep <- log10(rev$PEP / rev$pep_updated)
hist(dpep)

hist(-log10(rev$pep_updated[rev$pep_updated > 1e-4]), breaks=100)
hist(rev$residual[abs(rev$residual) < 3], breaks=100, 
     xlab='Residual RT', main='Decoys (Reverse Sequences)')

# number of decoys as func of FDR -----------------------------------------

rev_001 <- ev_001_f[grepl('REV_', ev_001_f$`Leading razor protein`),]
rev <- ev_f[grepl('REV_', ev_f$`Leading razor protein`),]

# only select REV sequences with at least one confident spectral PSM?
conf_rev_seqs <- rev %>% filter(PEP < 0.01) %>% pull(`Modified sequence`)

x <- seq(-3.5, -0.5, length.out=200)
fpvec <- vector(length=length(x))
fpvec_dart_1 <- vector(length=length(x))
fpvec_dart_2 <- vector(length=length(x))
for(i in 1:length(x)) {
  fpvec[i] <- sum(rev_001$qval <= 10**x[i], na.rm=T) / sum(ev_001_f$qval <= 10**x[i], na.rm=T)
  fpvec_dart_2[i] <- sum(rev$qval_updated <= 10**x[i], na.rm=T) / sum(ev_f$qval_updated <= 10**x[i], na.rm=T)
}

rev$pep_updated[!rev$`Modified sequence` %in% conf_rev_seqs] <- rev$PEP[!rev$`Modified sequence` %in% conf_rev_seqs]
# recalculate qvalues
rev$qval_updated=(cumsum(rev$pep_updated[order(rev$pep_updated)]) / 
                     seq(1, nrow(rev)))[order(order(rev$pep_updated))]

for(i in 1:length(x)) {
  fpvec_dart_1[i] <- sum(rev$qval_updated <= 10**x[i], na.rm=T) / sum(ev_f$qval_updated <= 10**x[i], na.rm=T)
}

pdf(file='manuscript/Figs/error_rate_line_v2.pdf', width=5, height=4)

par(mar=c(2.75,3.5,0.75,0.75), cex.axis=1)

plot(0, 0, type='n', xlim=c(-3, -1), ylim=c(-3.5, -0.5),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(x, log10(fpvec), col=cb[1], lwd=3)
lines(x, log10(fpvec_dart_1), col=cb[2], lwd=3)
#lines(x, log10(fpvec_dart_2), col=cb[2], lwd=3)
#abline(a=0, b=1, lty=2, col='black', lwd=2)
abline(v=-2, lty=2, col='black', lwd=2)

axis(1, at=seq(-3, -1, by=0.5), labels=c('0.1%', NA, '1%', NA, '10%'),
     tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(-3.5, -0.5, by=0.5), labels=sapply(c(NA, 1e-3, NA, 1e-2, NA, 1e-1, NA), fancy_scientific),
     tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext('PSM False Discovery Rate (q-value)', side=1, line=1.5, cex=1)
mtext('# Decoys / # All PSMs', side=2, line=2.25, cex=1)
# mtext('All Decoy Sequences', side=3, line=0.25, font=2, cex=1)

# legend('bottomright', c('Spectra', expression('DART-ID'[1+2]), expression('DART-ID'[1]), 'Ideal FDR'), 
#        col=c(cb[1], cb[2], paste0(cb[2], '66'), 'black'), lty=c(1, 1, 1, 2), 
#        lwd=c(3, 3, 3, 2), pt.cex=2, seg.len=1.75,
#        bty='n', x.intersp=1, y.intersp=1.25, inset=c(0.01, 0))
legend('bottomright', c('Spectra', expression('DART-ID'[1])),
       col=c(cb[1], cb[2]), lty=c(1, 1),
       lwd=c(3, 3), pt.cex=2, seg.len=1.75,
       bty='n', x.intersp=1, y.intersp=1.25, inset=c(0.01, 0))

dev.off()
