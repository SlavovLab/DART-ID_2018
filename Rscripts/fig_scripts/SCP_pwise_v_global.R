# init, load data ---------------------------------------------------------

library(tidyverse)
library(viridisLite)
library(RColorBrewer)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

# wrangle data ------------------------------------------------------------

ev.f <- ev %>%
  filter(pep_updated < 0.01) %>%
  filter(!is.na(exp_id)) %>%
  filter(!grepl('CON|REV', `Leading razor protein`)) %>%
  select(`Raw file`, `Modified sequence`, `Retention time`, `mu`, exp_id)


#FP18K
ref <- ev.f %>% 
  filter(exp_id == 105) %>%
  distinct(`Modified sequence`, .keep_all=TRUE) %>%
  arrange(`Modified sequence`)


a <- ev.f %>% 
  filter(!is.na(exp_id)) %>% distinct(`Raw file`, .keep_all=T) %>% 
  arrange(exp_id)


# plot --------------------------------------------------------------------


png(file='manuscript/Figs/SCP_pwise_v_global.png', width=6, height=3, units='in', res=400)
par(las=1, mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02)
layout(cbind(1,2))

# faux pairwise comparison. use FP18K as the reference 

# exps <- c(1,2,14,66,95,180,190,200)
exps <- c(1,2,14,180,190)
cols <- brewer.pal(9, 'Set1')
# cols <- av

plot(0,0,type='n', xlim=c(10, 60), ylim=c(10, 65),
     xlab='Reference - Observed RT (min)', ylab='Observed retention time (min)')

for(i in 1:length(exps)) {
  e <- exps[i]
  ev.a <- ev.f %>%
    filter(exp_id == e) %>%
    distinct(`Modified sequence`, .keep_all=TRUE) %>%
    arrange(`Modified sequence`)
  
  common_seqs <- intersect(ev.a$`Modified sequence`, ref$`Modified sequence`)
  
  points(ref[ref$`Modified sequence` %in% common_seqs,]$`Retention time`, 
         ev.a[ev.a$`Modified sequence` %in% common_seqs,]$`Retention time`, 
         col=paste0(cols[i], '66'), pch=16, cex=0.3)
}

# against mus

# exps <- c(1,2,14,66,95,180,190,200)
exps <- c(1,2,14,180,190)
cols <- brewer.pal(9, 'Set1')
# cols <- av

plot(0,0,type='n',
     xlim=c(10, 60), ylim=c(10, 65),
     xlab='Latent retention time (min)', ylab='Observed retention time (min)')

for(i in 1:length(exps)) {
  e <- exps[i]
  ev.a <- ev.f %>%
    filter(exp_id == e) %>%
    distinct(`Modified sequence`, .keep_all=TRUE) %>%
    arrange(`Modified sequence`)
  
  common_seqs <- intersect(ev.a$`Modified sequence`, ref$`Modified sequence`)
  
  points(ev.a[ev.a$`Modified sequence` %in% common_seqs,]$mu,
         ev.a[ev.a$`Modified sequence` %in% common_seqs,]$`Retention time`, 
         col=paste0(cols[i], '66'), pch=16, cex=0.3)
}

dev.off()

# intersect between all runs? ---------------------------------------------

exps <- c(150, 1,2,14,180,190)
ev.b <- ev.f %>%
  filter(exp_id %in% exps) %>%
  select(-`Raw file`) %>%
  distinct(`Modified sequence`, exp_id, .keep_all=TRUE) %>%
  spread(exp_id, `Retention time`)

