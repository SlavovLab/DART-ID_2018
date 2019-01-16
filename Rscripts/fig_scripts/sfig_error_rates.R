# init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

# load data ---------------------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/FP18_human_yeast_20190104/ev_updated.txt')

#human <- read.fasta(file='/gd/MS/DBs_Params/FASTA/swissprot_human_20181210.fasta', seqtype='AA')
# use all proteins, including computationally posited ones
human <- read.fasta(file='/gd/MS/DBs_Params/FASTA/uniprot_human_all_20190108.fasta', seqtype='AA')

yeast_seqs <- ev %>% 
  filter(grepl('YEAST', `Leading razor protein`)) %>%
  filter(!grepl('REV__', `Leading razor protein`)) %>%
  filter(!is.na(pep_new)) %>%
  #dplyr::select(c('Modified sequence', 'Leading razor protein', 'PEP', 'pep_updated'))
  pull(Sequence)

a <- ev %>% 
  filter(grepl('YEAST', `Leading razor protein`)) %>%
  filter(!grepl('REV__', `Leading razor protein`)) %>%
  filter(!is.na(pep_new)) %>%
  dplyr::select(c('Modified sequence', 'Raw file', 'Leading razor protein', 'PEP', 'pep_new', 'Retention time'))

yeast_seqs <- unique(yeast_seqs)

b <- lapply(human, function(entry) {
  #any(!is.na(str_locate(entry, yeast_seqs)))
  which(!is.na(str_locate(entry, yeast_seqs)[,1]))
})

matches <- sort(unique(unlist(b)))


# load sqc w/ reverse seqs ------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_include_decoy_20190109/ev_updated.txt')

ev.f <- ev %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))])

# decoy sequences with aligned RTs? ---------------------------------------

sum(grepl('REV_', ev.f$`Leading razor protein`) & !is.na(ev.f$pep_new))

# frequency of decoy sequences

rev <- ev.f[grepl('REV_', ev.f$`Leading razor protein`),]

a <- as.data.frame(table(rev$`Leading razor protein`))
a <- a %>% arrange(desc(Freq))

dpep <- log10(rev$PEP / rev$pep_updated)
hist(dpep)

hist(-log10(rev$pep_updated[rev$pep_updated > 1e-4]), breaks=100)
hist(rev$residual[abs(rev$residual) < 3], breaks=100, 
     xlab='Residual RT', main='Decoys (Reverse Sequences)')

# number of decoys as func of FDR -----------------------------------------

# only select REV sequences with at least one confident spectral PSM?
conf_rev_seqs <- rev %>% filter(PEP < 0.01) %>% pull(`Modified sequence`)

revv <- rev

x <- seq(-3.5, -0.5, length.out=200)
fpvec <- vector(length=length(x))
fpvec_dart_1 <- vector(length=length(x))
fpvec_dart_2 <- vector(length=length(x))
for(i in 1:length(x)) {
  fpvec[i] <- sum(revv$qval <= 10**x[i], na.rm=T) / sum(ev.f$qval <= 10**x[i], na.rm=T)
  fpvec_dart_2[i] <- sum(revv$qval_updated <= 10**x[i], na.rm=T) / sum(ev.f$qval_updated <= 10**x[i], na.rm=T)
}

revv$pep_updated[!revv$`Modified sequence` %in% conf_rev_seqs] <- rev$PEP[!revv$`Modified sequence` %in% conf_rev_seqs]
# recalculate qvalues
revv$qval_updated=(cumsum(revv$pep_updated[order(revv$pep_updated)]) / 
                     seq(1, nrow(revv)))[order(order(revv$pep_updated))]

for(i in 1:length(x)) {
  fpvec_dart_1[i] <- sum(revv$qval_updated <= 10**x[i], na.rm=T) / sum(ev.f$qval_updated <= 10**x[i], na.rm=T)
}

pdf(file='manuscript/Figs/error_rate_line.pdf', width=5, height=4)

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
