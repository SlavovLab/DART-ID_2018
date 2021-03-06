library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## New scatter ----------

#scatter <- 
ggplot(ev.f, aes(x=PEP, y=pep_new)) +
  stat_bin2d(bins=bins, drop=TRUE, geom='tile', aes(fill=..density..)) +
  geom_abline(slope=1, intercept=0, color='black', size=0.5) +
  #geom_segment(aes(x=5e-4, xend=1, y=1e-2, yend=1e-2),
  #             color='black', linetype='dotted', size=0.5) +
  geom_vline(xintercept=1e-2, linetype='dotted', color='black', size=0.5) +
  geom_hline(yintercept=1e-2, linetype='dotted', color='black', size=0.5) +
  scale_fill_gradient(low='white', high='red', 
                      breaks=c(0.002,0.004,0.006),
                      #labels=c(2, 4, 6)) +
                      labels=NULL) +
  scale_x_log10(expand=c(0.0,0), limits=c(1e-5, 1), 
                breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), labels=fancy_scientific) +
  #breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  scale_y_log10(expand=c(0.0,0), limits=c(1e-5, 1), 
                breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), labels=fancy_scientific) +
  #breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  labs(x='Spectra PEP', y='DART-ID PEP', fill='Density'
       ,title='Confidence shifts') +
  #guides(fill=guide_colorbar()) +
  theme_bert() + theme(
    plot.margin=unit(c(0.1,0.3,0.1,0.1), 'cm'),
    legend.position=c(0.03, 0.9),
    legend.justification=c(0,1),
    #legend.background=element_rect(fill='white'),
    legend.title=element_text(size=8, hjust=0.5),
    legend.text=element_text(size=8, hjust=0.5),
    legend.key.height=unit(0.35,'cm'),
    legend.key.width=unit(0.35,'cm'),
    plot.title=element_text(size=10, margin=margin(0,0,0.1,0,'cm'), hjust=0.5),
    aspect.ratio=1
  )

#return(scatter)