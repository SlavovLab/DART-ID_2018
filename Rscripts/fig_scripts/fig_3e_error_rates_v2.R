# init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

# Load Data ---------------------------------------------------------------

ev_dart <- read_tsv('/gd/bayesian_RT/Alignments/SQC_include_decoy_20190109/ev_updated.txt')

ev_dart.f <- ev_dart %>%
  dplyr::select(c('Modified sequence', 'Raw file', 'Leading razor protein', 
                  'PEP', 'pep_updated', 'pep_new', 'residual', 'MS/MS scan number')) %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev_dart)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / seq(1, nrow(ev_dart)))[order(order(pep_updated))])


fdr_levels <- c(0.001, 0.003, 0.01, 0.03, 0.1)

y_dart <- vector(length=length(fdr_levels))

for(i in 1:length(fdr_levels)) {
  thresh <- fdr_levels[i]
  
  remove_seqs <- ev_dart.f %>%
    group_by(`Modified sequence`) %>%
    summarise(p = min(PEP)) %>%
    filter(p > thresh) %>%
    pull(`Modified sequence`)
  
  e <- ev_dart.f %>% filter(!`Modified sequence` %in% remove_seqs)
  
  # recalculate qvalues
  e <- e %>% 
    mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(e)))[order(order(PEP))],
           qval_updated=(cumsum(pep_updated[order(pep_updated)]) /
                         seq(1, nrow(e)))[order(order(pep_updated))])
  
  y_dart[i] <- nrow(e %>% filter(qval_updated <= thresh) %>% filter(grepl('REV', `Leading razor protein`))) / 
    nrow(e %>% filter(qval_updated <= thresh))
}

plot(log10(fdr_levels), log10(y_dart), 
     xlab='Log10 PSM False Discovery Rate (q-value)', ylab='Log10 # Decoys / # All PSMs')
abline(a=0, b=1, col='red')

# MQ data -----------------------------------------------------------------

y_mq <- vector(length=length(fdr_levels))
fdr_levels <- c(0.001, 0.003, 0.01, 0.03, 0.1)
ev_fdr_files <- c('/gd/bayesian_RT/dat_SQC_FDR_1/evidence.txt',
                  '/gd/bayesian_RT/dat_SQC_FDR_1/evidence.txt',
                  '/gd/bayesian_RT/dat_SQC_FDR_1/evidence.txt',
                  '/gd/bayesian_RT/dat_SQC_FDR_3/evidence.txt',
                  '/gd/bayesian_RT/dat_SQC_FDR_10/evidence.txt')

for(i in 1:length(fdr_levels)) {

  ev_mq <- read_tsv(ev_fdr_files[i])

  ev_mq.f <- ev_mq %>%
    dplyr::select(c('Leading razor protein', 'PEP')) %>%
    # ceil PEPs to 1
    mutate_at(c('PEP'), funs(ifelse(. > 1, 1, .))) %>%
    # calculate q-values
    mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev_mq)))[order(order(PEP))])

  thresh <- fdr_levels[i]

  rev_mq <- ev_mq.f %>% filter(grepl('REV_', `Leading razor protein`))

  y_mq[i] <- sum(rev_mq$qval <= thresh, na.rm=T) / sum(ev_mq.f$qval <= thresh, na.rm=T)
}

plot(log10(fdr_levels), log10(y_mq),
     xlab='Log10 PSM False Discovery Rate (q-value)', ylab='Log10 # Decoys / # All PSMs')
abline(a=0, b=1, col='red')


# Plot --------------------------------------------------------------------

pdf(file='manuscript/Figs/error_rate_line_v4.pdf', width=1.75, height=1.5)

par(mar=c(0.5,2,0.1,0.1), oma=c(1,0,0,0), cex.axis=0.7, pty='s')

plot(0, 0, type='n', xlim=c(-3.3, -0.7), ylim=c(-3.3, -0.7),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')

# grid()
abline(v=log10(fdr_levels), col = "grey70", lty = "dotted", lwd=1)
abline(h=log10(fdr_levels), col = "grey70", lty = "dotted", lwd=1)

points(log10(fdr_levels), log10(y_dart), col=cb[2], pch=16, cex=0.8)
lines(log10(fdr_levels), log10(y_dart), col=cb[2], lwd=2)

points(log10(fdr_levels), log10(y_mq), col=cb[1], pch=16, cex=0.8)
lines(log10(fdr_levels), log10(y_mq), col=cb[1], lwd=2)

abline(a=0, b=1, col='black', lwd=1, lty=2)
axis(1, at=log10(fdr_levels), labels=c('0.1%', '0.3%', '1%', '3%', '10%'),
     tck=-0.02, mgp=c(0, -0.15, 0))
axis(2, at=log10(fdr_levels), labels=c('0.1%', NA, '1%', NA, '10%'),
     tck=-0.02, mgp=c(0, 0.2, 0), las=1)

mtext('PSM FDR (q-value)', side=1, line=0.5, cex=0.7, outer=F)
mtext('# Decoys / # All PSMs', side=2, line=1.6, cex=0.7)

legend('bottomright', c(expression('Spectra'), expression('DART-ID'[1])),
       col=c(cb[1], cb[2]), lty=c(1, 1),
       lwd=c(3, 3), pt.cex=2, seg.len=0.8, cex=0.7,
       bty='n', x.intersp=0.75, y.intersp=1, inset=c(-0.01, -0.05))


dev.off()

