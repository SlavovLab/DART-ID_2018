library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

## ---------
no_conf_ids <- apply(dmat, 1, sum) == 0
boxs <- list(Spectra=apply(dmat, 2, sum), 
             Percolator=apply(dmat_perc, 2, sum),
             #Percolator_prev=apply(dmat_perc[!no_conf_ids,], 2, sum),
             DART=apply(dmat_new, 2, sum),
             DART_prev=apply(dmat_new[!no_conf_ids,], 2, sum))

#save(boxs, file='dat/peptide_ids_20180816.rds')

## ---------

load('/gd/bayesian_RT/dat/peptide_ids_20180816.rds')

# horizontal boxplot ------------------------------------------------------

pdf(file='manuscript/Figs/peps_per_exp_v9.pdf', width=1.75, height=1.5)

par(mar=c(1,3,0,0.25),
    oma=c(0,0,1,0),
    pty='m', las=1, cex.axis=0.6)

boxplot(rev(boxs), horizontal=T,
        col=rev(c(cb[1], cb[3], cb[2], paste0(cb[2], '44'))), 
        xaxt='n', yaxt='n', ylim=c(50, 2650),
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0, 0))
axis(1, at=seq(0, 2500, by=500), labels=c(0, NA, 1000, NA, 2000, NA), tck=-0.02, 
     mgp=c(0, -0.1, 0))
#axis(2, at=seq(1,4), labels=NA, #labels=c(expression('DART-ID'[2]), expression('DART-ID'[1]), 
                              #expression('Percolator'[2]), expression('Percolator'[1]),
                              #'Percolator',
                              #'Spectra'),
#     tck=-0.02, mgp=c(0, 0.3, 0), las=1)
text(x=rep(-150, 4), y=seq(0.85,4,by=1),
     labels=c(expression('DART-ID'[1]), expression('DART-ID'[1+2]),'Percolator','Spectra'),
     xpd=T, srt=-30, cex=0.65, adj=c(0.95, 0.7))

mtext('        Peptides/Experiment', side=3, line=0.25, las=1, font=2, cex=0.85, outer=T)

dev.off()

# horizontal boxplot - HUPO ------------------------------------------------------

pdf(file='manuscript/Figs/peps_per_exp_v8_harrison.pdf', width=3, height=2.5)

par(mar=c(1.5,4.25,0,0.5),
    oma=c(0,0,1.6,0),
    pty='m', las=1, cex.axis=1)

boxplot(rev(boxs[1:3]), horizontal=T,
        col=rev(c(cb[1], cb[3], cb[2])), 
        xaxt='n', yaxt='n', ylim=c(0, 2700),
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0))
axis(1, at=seq(0, 2500, by=500), labels=c(0, NA, 1000, NA, 2000, NA), tck=-0.02, 
     mgp=c(0, 0.3, 0))
#axis(2, at=seq(1,4), labels=NA, #labels=c(expression('DART-ID'[2]), expression('DART-ID'[1]), 
#expression('Percolator'[2]), expression('Percolator'[1]),
#'Percolator',
#'Spectra'),
#     tck=-0.02, mgp=c(0, 0.3, 0), las=1)
text(x=rep(-180, 4), y=seq(0.95,4,by=1),
     labels=c(expression('DART-ID'),'Percolator','Spectra'),
     xpd=T, srt=-30, cex=1, adj=c(1, 0.5))

mtext(' Peptides per 60 min LC-MS/MS run', side=3, line=0.4, las=1, font=2, cex=1, outer=T)

dev.off()

## ---------

pdf(file='manuscript/Figs/peps_per_exp_v7.pdf', width=1.75, height=1.5)

par(mar=c(2,1.75,0.25,0.25),
    pty='m', las=1, cex.axis=0.65)

# plot(0, 0, type='n',
#      xlim=c(0, 3), ylim=c(0, 1000),
#      xlab=NA, ylab=NA,
#      xaxt='n', yaxt='n')

boxplot(boxs, 
        col=c(av[1], av[3], paste0(av[3], '66'), av[2], paste0(av[2], '66')), 
        xaxt='n', yaxt='n', ylim=c(75, 2525),
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0, 0.5))
axis(1, at=c(1, 2.5, 4.5), tck=-0.02, labels=NA,
     #labels=c('Spectra', 'Percolator', 'DART-ID'),
     srt=45)
axis(2, at=seq(0, 2500, by=500), labels=seq(0, 2500, by=500), tck=-0.02, 
     mgp=c(0, 0.1, 0), las=3)
text(c(1, 2.5, 4.5), par('usr')[3]-150, srt=45, adj=c(1, 0.5), xpd=T,
     labels=c('Spectra', 'Percolator', 'DART-ID'), cex=0.6)
#text(c())

mtext('Peptides ID\'d per Experiment           ', 2, line=1, las=3, cex=0.7)

dev.off()


# scatter -----------------------------------------------------------------

pdf(file='manuscript/Figs/peps_per_exp_scatter.pdf', width=5, height=5)

par(mar=c(3, 4, 3, 1), cex.axis=1)

dart_pars <- lm(boxs[['DART']] ~ boxs[['Spectra']])
perc_pars <- lm(boxs[['Percolator']] ~ boxs[['Spectra']])

plot(0, 0, type='n', xlim=c(0, 1700), ylim=c(0, 2500),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)
points(boxs[['Spectra']], boxs[['DART']], pch=16, cex=1.25, col=cb[2])
points(boxs[['Spectra']], boxs[['Percolator']], pch=16, cex=1.25, col=cb[3])

abline(a=0, b=1, col=cb[1], lwd=2) # Spectra
#abline(a=dart_pars$coefficients[1], b=dart_pars$coefficients[2], col=cb[2], lwd=2)
#abline(a=perc_pars$coefficients[1], b=perc_pars$coefficients[2], col=cb[3], lwd=2)

axis(1, at=seq(0, 2000, by=500), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 2500, by=500), tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext('# PSMs from Spectra', side=1, line=1.5, cex=1)
mtext('# Updated PSMs', side=2, line=2.75, cex=1)
mtext('Peptides Quantified Per SCoPE-MS set at 1% FDR', side=3, line=0.5, cex=1, font=2)

legend('bottomright', c('Spectra', 'DART-ID', 'Percolator'), 
       pch=c(NA, 16, 16), lty=c(1, NA, NA), lwd=c(3, NA, NA), seg.len=1, pt.cex=2,
       col=c(cb[1], cb[2], cb[3]), bty='n', inset=c(0.01, 0.01), cex=1,
       x.intersp=1, y.intersp=1.2)

dev.off()
