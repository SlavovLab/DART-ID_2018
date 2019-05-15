# init --------------------------------------------------------------------

library(tidyverse)
library(ggridges)
library(RColorBrewer)
source('Rscripts/lib.R')


# load data ---------------------------------------------------------------

load(file='/gd/bayesian_RT/dat/othermethods_20190422.rds')

# wrangle -----------------------------------------------------------------

a <- ev.f %>%
  select(id, PEP, pep_updated, ssrcalc_pep, biolccc_pep, elude_pep, irt_pep, MBR_pep) %>%
  gather(key='method', value='error', -id) %>%
  arrange(id) %>%
  mutate(method=factor(method, 
                    levels=rev(c('PEP', 'ssrcalc_pep', 'biolccc_pep', 'elude_pep', 'irt_pep', 'MBR_pep', 'pep_updated')),
                    labels=rev(c('Spectra', 'SSRCalc', 'BioLCCC', 'ELUDE', 'iRT', 'MaxQuant', 'DART-ID')))) %>%
  group_by(method) %>%
  mutate(med_error=median(log10(error), na.rm=T))

# plot --------------------------------------------------------------------

ggplot(a) +
  geom_density_ridges(aes(log10(error), y=method)) +
  scale_x_continuous(limits=c(-7.5, 0)) +
  labs(x='Log10 Error Probability', y=NULL) +
  theme(text=element_text(size=20))

# just spectra - for SCP slide --------------------------------------------

load('/gd/bayesian_RT/dat/pep_dists_20190405.rds')

df$pep[df$pep > 1] <- 1

thrown_out <- sum(df$method == 'Spectra' & df$pep > 0.01) / sum(df$method == 'Spectra') * 100

# ----------

png(file='manuscript/Figs/SCP_spectra_pep_dist.png', width=4, height=4, units='in', res=400)
par(las=1, mgp=c(2.5,1,0), cex.axis=1.25, cex.lab=1.4, cex.main=1.4)

hist(log10(df$pep[df$method == 'Spectra']), breaks=50, xlim=c(-10, 0),
     xaxt='n', yaxt='n', xlab='Error Probability', ylab='# Peptides x 10,000', main='SCoPE-MS: Identified Peptides')
abline(v=-2, col='red', lwd=2, lty=1)
rect(xleft=-2, ybottom=0, xright=0, ytop=1e6, col=rgb(1,0,0,0.2), border=NA)
text(-2.5, 1e5, paste0(formatC(thrown_out, digits=4), '% removed'), adj=c(1, 0.5), cex=1.4, col='red')
axis(1, at=seq(-10, 0, by=2), labels=fancy_scientific(10**seq(-10, 0, by=2)))
axis(2, at=seq(0, 1.5e5, by=5e4), labels=c(0, 5, 10, 15))

dev.off()
