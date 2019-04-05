# Init --------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(pracma)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')
iprg <- read_tsv('/gd/bayesian_RT/Alignments/iPRG2015_20181228/ev_updated.txt')
tko <- read_tsv('/gd/bayesian_RT/Alignments/PXD011654_20190120/ev_updated.txt')

# Get boxplot data --------------------------------------------------------

# get experiment missingness, before/after DART, at peptide and protein level, for
# peptides/proteins quantified in 1) >20% experiments and 2) >50% experiments
bplot_data <- function(e) {
  
  # filter for sqc master sets only, and only keep correctly IDed proteins
  e.f <- e %>%
    filter(!grepl('SQC67[AB][16]|SQC67C1[3-9]|SQC67[CD]5|SQC68[DE]|IFN6[H-K]-Trg|SQC72D|SQC73[CD]|SQC74M|180416S_QC_SQC78A2',`Raw file`)) %>%
    filter(!grepl('REV__', `Leading razor protein`)) %>%
    filter(!grepl('CON__', Proteins)) %>%
    # filter(apply(ev[,data.cols]!=0, 1, sum) == length(data.cols)) %>%
    mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
      if(length(unlist(p)) == 1) return(p[1])
      else if(length(unlist(p)) == 3) return(p[2])
      else return(p[1])
    }))
  
  e.f <- e.f %>%
    # ceil PEPs to 1
    mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
    # calculate q-values
    mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(e.f)))[order(order(PEP))],
           qval_updated=(cumsum(pep_updated[order(pep_updated)]) / seq(1, nrow(e.f)))[order(order(pep_updated))])
  
  # get peptides, proteins in 1) >20% and 2) >50% experiments
  pep_exp_map <- e.f %>%
    group_by(`Raw file`, `Modified sequence`) %>%
    tally() %>%
    spread(`Raw file`, n)
  
  peps_per_exp <- apply(pep_exp_map[,-1], 1, sum, na.rm=T)
  peps_20 <- pep_exp_map$`Modified sequence`[( peps_per_exp / (ncol(pep_exp_map) - 1) ) > 0.2]
  peps_50 <- pep_exp_map$`Modified sequence`[( peps_per_exp / (ncol(pep_exp_map) - 1) ) > 0.5]
  
  # get peptides, proteins in 1) >20% and 2) >50% experiments
  prot_exp_map <- e.f %>%
    group_by(`Raw file`, Protein) %>%
    tally() %>%
    spread(`Raw file`, n)
  
  prots_per_exp <- apply(prot_exp_map[,-1], 1, sum, na.rm=T)
  prots_20 <- prot_exp_map$Protein[( prots_per_exp / (ncol(prot_exp_map) - 1) ) > 0.2]
  prots_50 <- prot_exp_map$Protein[( prots_per_exp / (ncol(prot_exp_map) - 1) ) > 0.5]
  
  # peptides - get "missingness" per experiment before/after dart
  # 1) >20%
  pep_exp_20 <- e.f %>%
    select(`Raw file`, `Modified sequence`, qval, qval_updated) %>%
    filter(`Modified sequence` %in% peps_20) %>%
    group_by(`Raw file`, `Modified sequence`) %>%
    summarise(qval=min(qval), qval_updated=min(qval_updated)) %>% ungroup() %>%
    mutate(spectra = qval < 0.01, dart = qval_updated < 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(spectra = sum(spectra, na.rm=T), dart = sum(dart, na.rm=T)) %>%
    mutate(spectra = spectra / length(peps_20), dart = dart / length(peps_20)) %>%
    mutate(level='Peptide', group='20%')
  
  # 2) >50%
  pep_exp_50 <- e.f %>%
    select(`Raw file`, `Modified sequence`, qval, qval_updated) %>%
    filter(`Modified sequence` %in% peps_50) %>%
    group_by(`Raw file`, `Modified sequence`) %>%
    summarise(qval=min(qval), qval_updated=min(qval_updated)) %>% ungroup() %>%
    mutate(spectra = qval < 0.01, dart = qval_updated < 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(spectra = sum(spectra, na.rm=T), dart = sum(dart, na.rm=T)) %>%
    mutate(spectra = spectra / length(peps_50), dart = dart / length(peps_50)) %>%
    mutate(level='Peptide', group='50%')
  
  pep_exp_all <- e.f %>%
    select(`Raw file`, `Modified sequence`, qval, qval_updated) %>%
    group_by(`Raw file`, `Modified sequence`) %>%
    summarise(qval=min(qval), qval_updated=min(qval_updated)) %>% ungroup() %>%
    mutate(spectra = qval < 0.01, dart = qval_updated < 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(spectra = sum(spectra, na.rm=T), dart = sum(dart, na.rm=T)) %>%
    mutate(spectra = spectra / length(peps_per_exp), dart = dart / length(peps_per_exp)) %>%
    mutate(level='Peptide', group='All')
  
  # proteins - get "missingness" per experiment before/after dart
  # 1) >20%
  prot_exp_20 <- e.f %>%
    select(`Raw file`, Protein, qval, qval_updated) %>%
    filter(Protein %in% prots_20) %>%
    group_by(`Raw file`, Protein) %>%
    summarise(qval=min(qval), qval_updated=min(qval_updated)) %>% ungroup() %>%
    mutate(spectra = qval < 0.01, dart = qval_updated < 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(spectra = sum(spectra, na.rm=T), dart = sum(dart, na.rm=T)) %>%
    mutate(spectra = spectra / length(prots_20), dart = dart / length(prots_20)) %>%
    mutate(level='Protein', group='20%')
  
  # 2) >50%
  prot_exp_50 <- e.f %>%
    select(`Raw file`, Protein, qval, qval_updated) %>%
    filter(Protein %in% prots_50) %>%
    group_by(`Raw file`, Protein) %>%
    summarise(qval=min(qval), qval_updated=min(qval_updated)) %>% ungroup() %>%
    mutate(spectra = qval < 0.01, dart = qval_updated < 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(spectra = sum(spectra, na.rm=T), dart = sum(dart, na.rm=T)) %>%
    mutate(spectra = spectra / length(prots_50), dart = dart / length(prots_50)) %>%
    mutate(level='Protein', group='50%')
  
  prot_exp_all <- e.f %>%
    select(`Raw file`, Protein, qval, qval_updated) %>%
    group_by(`Raw file`, Protein) %>%
    summarise(qval=min(qval), qval_updated=min(qval_updated)) %>% ungroup() %>%
    mutate(spectra = qval < 0.01, dart = qval_updated < 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(spectra = sum(spectra, na.rm=T), dart = sum(dart, na.rm=T)) %>%
    mutate(spectra = spectra / length(prots_per_exp), dart = dart / length(prots_per_exp)) %>%
    mutate(level='Protein', group='All')
  
  return(rbind(pep_exp_20, pep_exp_50, pep_exp_all, prot_exp_20, prot_exp_50, prot_exp_all))
  # return(rbind(pep_exp_all, prot_exp_all))
}

ev.f <- bplot_data(ev) %>% mutate(dataset='SCoPE-MS')
iprg.f <- bplot_data(iprg) %>% mutate(dataset='Label-free')
tko.f <- bplot_data(tko) %>% mutate(dataset='TMT-labelled')

df <- rbind(ev.f, iprg.f, tko.f)
df <- df %>% 
  gather('Method', 'x', -c(`Raw file`, level, group, dataset)) %>%
  filter(group == 'All') %>%
  select(-group) %>%
  mutate(x = 1 - x)

df$dataset <- fct_recode(df$dataset, `SCoPE-MS`='SCoPE-MS', `Bulk Label-free`='Label-free', `Bulk TMT-labelled`='TMT-labelled')
df$dataset <- fct_relevel(df$dataset, 'SCoPE-MS', 'Bulk Label-free', 'Bulk TMT-labelled')

df$level <- fct_recode(df$level, Peptides='Peptide', Proteins='Protein')
df$Method <- fct_recode(df$Method, Spectra='spectra', `DART-ID`='dart')
df$Method <- fct_relevel(df$Method, 'Spectra', 'DART-ID')


# PLOT --------------------------------------------------------------------

b_blue <- '#205AB0'
b_light_blue <- '#92D8FF'
b_red <- '#C94040'
b_light_red <- '#FFACAA'

p <- ggplot(df) +
  geom_boxplot(aes(x=dataset, y=x, color=Method, fill=Method), 
               outlier.color=rgb(0,0,0,0.5), outlier.shape=4, size=0.25, width=0.5,
               position=position_dodge2(padding=0.15)) +
  geom_vline(xintercept=c(1.5, 2.5), size=0.5, linetype='dashed', color=rgb(0,0,0,0.3)) +
  #geom_violin(bw=0.01, kernel='rectangular', size=0.25) +
  facet_grid(.~level) +
  scale_y_continuous(limits=c(0.3, 1), breaks=seq(0.3, 1, by=0.15))
  
dat <- ggplot_build(p)$data[[1]]
dat$level <- factor(dat$PANEL, labels=c('Peptides', 'Proteins'))

# p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), colour='black', size=0.25)

p <- p + 
  scale_fill_manual(values=c(b_light_blue, b_light_red)) +
  scale_color_manual(values=c(b_blue, b_red), guide=F) +
  labs(x=NULL, y='Missing Data per Run', fill=NULL) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size=12, angle=30, hjust=1, vjust=1),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14),
    strip.text.x = element_text(size=12),
    strip.text.y = element_text(size=14),
    legend.text = element_text(size=12),
    legend.key.size = unit(0.7, 'cm'),
    # panel.spacing.y = unit(0.4, 'cm'),
    legend.position = 'top',
    legend.margin = margin(b=0, unit='cm')
  )

# p

#ggsave('manuscript/Figs/missing_data_violin.pdf', plot=p, device='pdf', width=7, height=4, units='in')
ggsave('manuscript/Figs/missing_data_boxplot.pdf', plot=p, device='pdf', width=3.5, height=4.5, units='in')


# Significance Tests ------------------------------------------------------

for(l in levels(df$level)) {
  for(d in levels(df$dataset)) {
    print(paste0(l, ' - ', d))
    # paired, one-sided t-test, w/ equal variance:
    s <- t.test(df %>% filter(level == l & dataset == d) %>% filter(Method == 'Spectra') %>% pull(x),
                df %>% filter(level == l & dataset == d) %>% filter(Method == 'DART-ID') %>% pull(x),
                alternative='greater', paired=TRUE, var.equal=TRUE)
    print(s$p.value)
  }
}

# OUTPUT:
# [1] "Peptides - SCoPE-MS"
# [1] 6.267421e-91
# [1] "Peptides - Label-free"
# [1] 8.043409e-16
# [1] "Peptides - TMT-labelled"
# [1] 4.469861e-29
# [1] "Proteins - SCoPE-MS"
# [1] 1.792966e-100
# [1] "Proteins - Label-free"
# [1] 4.911075e-15
# [1] "Proteins - TMT-labelled"
# [1] 3.80131e-53

# all very statistically significant