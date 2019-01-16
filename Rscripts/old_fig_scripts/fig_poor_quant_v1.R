library(tidyverse)
library(ggridges)
library(pracma)
source('Rscripts/lib.R')

#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt')
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## load cormat data -----

source('Rscripts/validation_cormats.R')

# Precursor Area, PIF dists -----------------------------------------------

ev_a <- ev.fi %>% filter(qval < 0.01)
ev_b <- ev.fi %>% filter(qval > 0.01 & qval_updated < 0.01)

# count missing data
dcols <- grep('corrected', colnames(ev.fi))[5:10]
dmat <- data.matrix(ev_a %>% select(colnames(ev.fi)[dcols]))
missing_a <- apply(dmat==0, 1, sum)
dmat <- data.matrix(ev_b %>% select(colnames(ev.fi)[dcols]))
missing_b <- apply(dmat==0, 1, sum)

# boxplots ----------------------------------------------------------------

pdf(file='manuscript/Figs/poor_quant.pdf', width=7, height=2)

layout(rbind(c(1, 2, 3, 4)))

par(cex.axis=1, 
    mar=c(3.5,2,1,2.5),
    oma=c(0,1,0,0))

boxplot(list(log10(ev_a$Intensity), log10(ev_b$Intensity)),
        range=1.5, col=c(cb[1], cb[2]), ylim=c(4.75, 8),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
#axis(1, at=c(1, 2), labels=c('Spectra', 'DART-ID'), tck=-0.02, mgp=c(0, 0.2, 0))
axis(1, at=c(1, 2), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(c(1, 2), rep(4.3, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
axis(2, at=seq(4, 10), tck=-0.02, mgp=c(0, 0.5, 0), las=1)
mtext(expression('log'[10]*' Precursor Ion Area'), side=2, line=1.25, cex=0.85)

boxplot(list((ev_a$PIF), (ev_b$PIF)), 
        range=1.5, col=c(cb[1], cb[2]), ylim=c(0.55, 1),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
axis(1, at=c(1, 2), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0.3, 1, by=0.1), label=seq(0.3, 1, by=0.1)*100, tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2), rep(0.485, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Precursor Ion, %'), side=2, line=2, cex=0.85)

barplot(c(mean(ev_a$`Missed cleavages`), mean(ev_b$`Missed cleavages`)), width=1, space=0.5,
        range=1.5, col=c(cb[1], cb[2]), xlim=c(0.25, 3.15), ylim=c(0, 0.2),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
axis(1, at=c(-10, 1, 2.5, 10), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 0.25, by=0.05), tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2.5), rep(-0.02, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Missed Cleavages, %'), side=2, line=2.25, cex=0.85)

barplot(c(mean(missing_a) / 6, mean(missing_b) / 6), width=1, space=0.5,
        range=1.5, col=c(cb[1], cb[2]), xlim=c(0.25, 3.15), ylim=c(0, 0.15),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)
axis(1, at=c(-10, 1, 2.5, 10), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 1, by=0.05), labels=seq(0, 1, by=0.05)*100, 
     tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2.5), rep(-0.015, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Missing Data, %'), side=2, line=2, cex=0.85)

dev.off()
