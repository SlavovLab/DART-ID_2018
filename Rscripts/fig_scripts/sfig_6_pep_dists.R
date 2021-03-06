# Init --------------------------------------------------------------------

library(tidyverse)
library(ggridges)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180803_5exp_parametric_mixture_v2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180808_5exp_parametric/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')


df <- data.frame(pep=c(ev$PEP, ev$pep_updated),
                 method=rep(c('Spectra', 'DART-ID'), each=nrow(ev)))

df <- df[df$pep > 1e-15,]

save(df, file='/gd/bayesian_RT/dat/pep_dists_20190405.rds')

## ---------

load('/gd/bayesian_RT/dat/pep_dists_20190405.rds')


# plot --------------------------------------------------------------------

rng <- seq(-10, 0, by=2)

p <- ggplot(df) +
  geom_density_ridges(aes(x=log10(pep), y=method, group=method, fill=method), 
                      rel_min_height=0.01, scale=2) +
  scale_x_continuous(limits=c(-8, 0.3), breaks=rng, labels=fancy_scientific(10**rng)) +
  scale_y_discrete(expand=c(0.01, 0.01)) +
  scale_fill_manual(values=c(cb[2], cb[1]), guide=F) +
  labs(x='Posterior Error Probability (PEP)', y=NULL) +
  theme_ridges() + theme(
    #plot.margin = margin(t=0.1, r=0.1, b=0, l=0.1, unit='cm'),
    axis.text.y = element_text(size=18),
    axis.text.x = element_text(size=16),
    axis.title.x = element_text(size=18, hjust=0.5, vjust=0.5)
  )

ggsave(filename='manuscript/Figs/pep_dists_v2.pdf', plot=p, device='pdf', width=5, height=5)

