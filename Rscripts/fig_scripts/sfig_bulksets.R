# init, load data ---------------------------------------------------------

library(tidyverse)
source('Rscripts/lib.R')

iprg <- read_tsv('/gd/bayesian_RT/Alignments/iPRG2015_20181228/ev_updated.txt')
tko <- read_tsv('/gd/bayesian_RT/Alignments/PXD011654_20190120/ev_updated.txt')


# alignment residuals -----------------------------------------------------

pdf(file='manuscript/Figs/residual_iprg.pdf', width=3.5, height=3)

par(mar=c(3,3,3,2), cex.axis=1, cex.main=1, mgp=c(1,3,1))

hist(iprg$residual[abs(iprg$residual) < 1], breaks=seq(-1,1,by=0.02), freq=F,
     col='red', border=NA, xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
axis(1, at=seq(-1, 1, by=0.5), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 10, by=2), tck=-0.02, mgp=c(0,0.5,0), las=1)
mtext('Residual RT (min)', side=1, line=1.3)
mtext('Density', side=2, line=1.5)
mtext('Label-free', side=3, line=0.7, font=2)

dev.off()

pdf(file='manuscript/Figs/residual_tko.pdf', width=3.5, height=3)

par(mar=c(3,4,3,2), cex.axis=1, cex.main=1, mgp=c(1,3,1))

hist(tko$residual[abs(tko$residual) < 5], breaks=seq(-5,5,by=0.1), freq=F,
     col='blue', border=NA, xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
axis(1, at=seq(-5, 5, by=1), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 1.5, by=0.25), tck=-0.02, mgp=c(0,0.5,0), las=1)
mtext('Residual RT (min)', side=1, line=1.3)
mtext('Density', side=2, line=2.5)
mtext('TMT-labelled', side=3, line=0.7, font=2)

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


# dot plot exp ------------------------------------------------------------

fdr_thresholds <- c(10**-4, 10**-3, 10**-2, 10**-1)
df <- data.frame()
# averages
avgs <- matrix(nrow=length(fdr_thresholds), ncol=2)
for(i in 1:length(fdr_thresholds)) {
  fdr_threshold <- fdr_thresholds[i]
  
  df <- rbind(df, 
    iprg %>%
      group_by(`Raw file`) %>%
      summarise(p=((sum(qval_updated < fdr_threshold) / sum(qval < fdr_threshold)) - 1) * 100) %>%
      mutate(thresh=fdr_threshold, set='Label-free'),
    tko %>%
      group_by(`Raw file`) %>%
      summarise(p=((sum(qval_updated < fdr_threshold) / sum(qval < fdr_threshold)) - 1) * 100) %>%
      mutate(thresh=fdr_threshold, set='TMT-labelled')
  )
  
  avgs[i,] <- df %>% 
    filter(thresh == fdr_threshold) %>%
    group_by(set) %>%
    summarise(n=mean(p, na.rm=T)) %>% pull(n)
}

# dotplot -----------------------------------------------------------------

set.seed(1)

x_scatter_limit <- 0.2
iprg_exps <- length(unique(iprg$`Raw file`))
tko_exps <- length(unique(tko$`Raw file`))

pdf(file='manuscript/Figs/fdr_increase_iprg_tko_v2.pdf', width=3.5, height=2.5)

par(mar=c(2.75,3.25,0.5,0.5), las=1)

plot(0, 0, type='n', xlim=c(0.5, 4.5), ylim=c(0, 130),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)

for(i in 1:length(fdr_thresholds)) {
  fdr_threshold <- fdr_thresholds[i]
  
  y <- df %>% 
    filter(thresh == fdr_threshold) %>%
    arrange(set)
  
  x <- i + runif(nrow(y), -x_scatter_limit, x_scatter_limit)
  
  points(x, y$p, col=paste0( rep(c('#FF0000', '#0000FF'), c(iprg_exps, tko_exps)) , '55'), pch=16, cex=1)
}

# plot averages
lines(seq(1, 4), avgs[,1], col='red', lty=1, lwd=2)
lines(seq(1, 4), avgs[,2], col='blue', lty=1, lwd=2)
# abline(h=0, col=rgb(0,0,0,0.5), lty=2, lwd=1)

axis(1, at=seq(1, 4), labels=fancy_scientific(10**seq(-4, -1)),
     tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 120, by=20), tck=-0.02, mgp=c(0, 0.4, 0))

legend('topright', c('Label-free', 'TMT-labelled'), col=c('red', 'blue'), pch=16, pt.cex=1.25,
       lwd=3, lty=c(1, 1), seg.len=1.5, cex=0.9, x.intersp=1, y.intersp=1.2, bty='n')

mtext('FDR Threshold', side=1, line=1.25)
mtext('% Increase in PSMs', side=2, line=2, las=0)


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


# plot together -----------------------------------------------------------

# fancy scientific scales
# from: https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific_2 <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE, digits=1)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # make sure +0 just turns into 0
  l <- gsub("\\+0(\\d)", "\\1", l)
  # return this as an expression
  return(parse(text=l))
}

pdf(file='manuscript/Figs/peps_per_exp_scatter_iprg_tko_v2.pdf', width=3.5, height=3.5)

par(mar=c(3.25,4.25,2,0.5), cex.axis=0.85)

plot(0, 0, type='n', xlim=c(4e3, 31e3), ylim=c(4e3, 31e3),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)

points(iprg_exp$spectra, iprg_exp$dart, pch=16, col=rgb(1,0,0,0.5), cex=1.5)
points(tko_exp$spectra, tko_exp$dart, pch=16, col=rgb(0,0,1,0.5), cex=1.5)
abline(a=0, b=1, col=rgb(0,0,0,0.5), lty=2, lwd=1)

legend('bottomright', c('Label-free', 'TMT-labelled'), col=c('red', 'blue'),
       pch=16, pt.cex=1.25, cex=0.9, x.intersp=1, y.intersp=1.2, bty='n')

axis(1, at=seq(5e3, 30e3, by=5e3), labels=NA, tck=-0.02, mgp=c(0, 0.3, 0))
# manual x axis tick labels
text(x=seq(5e3, 30e3, by=5e3), y=2000, labels=seq(5e3, 30e3, by=5e3), 
     adj=c(1, 1), srt=30, xpd=T, cex=0.85)
axis(2, at=seq(5e3, 30e3, by=5e3), tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext('Spectra PSMs', side=1, line=2.2)
mtext('DART-ID PSMs', side=2, line=2.8)
mtext('PSMs per Experiment, 1% FDR', side=3, line=0.3, cex=1, font=2)

dev.off()


# plot --------------------------------------------------------------------

pdf(file='manuscript/Figs/peps_per_exp_scatter_iprg_tko.pdf', width=3.5, height=3.5)
layout(c(1,2))

par(mar=c(1.5,0,0,0), oma=c(0.8, 3.5, 1.8, 1.5), las=1, cex.axis=0.85)

plot(0, 0, type='n', xlim=c(20000,30000), ylim=c(20000, 30000),
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')

points(iprg_exp$spectra, iprg_exp$dart, pch=1, col='red')
abline(a=0, b=1, lwd=2, lty=1, col='black')
text(25000, 20500, 'Label-free', adj=c(0, 0), font=1, cex=1)
axis(1, at=seq(20000, 30000, by=2000), labels=formatC(seq(20, 30, by=2), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.1, 0))
axis(2, at=seq(20000, 30000, by=2000), labels=formatC(seq(20, 30, by=2), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.3, 0), las=1)

par(mar=c(1.25,0,0,0))

plot(0, 0, type='n', xlim=c(4500,14500), ylim=c(4500, 14500),
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
points(tko_exp$spectra, tko_exp$dart, pch=1, col='blue')
abline(a=0, b=1, lwd=2, lty=1, col='black')
text(9200, 5000, 'TMT-labelled', adj=c(0, 0), font=1, cex=1)
axis(1, at=seq(4500, 14500, by=2500), labels=formatC(seq(4.5, 14.5, by=2.5), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.1, 0))
axis(2, at=seq(4500, 14500, by=2500), labels=formatC(seq(4.5, 14.5, by=2.5), digits=1, format='f'), 
     tck=-0.02, mgp=c(0, 0.3, 0), las=1)

mtext('Spectra', side=1, outer=T, line=-0.25)
mtext('DART-ID', side=2, outer=T, line=2, las=3)
mtext('(PSMs x 1000) per Experiment, 1% FDR      ', side=3, outer=T, line=0.4, font=2)

dev.off()

