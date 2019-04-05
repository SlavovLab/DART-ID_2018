# init, load data ---------------------------------------------------------

library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

# filtering ---------------------------------------------------------------

ev.f <- ev %>%
  # get UniProt accession ID
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })) %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))]) %>%
  # exclude CON, REV proteins
  filter(!grepl('CON__|REV__', Protein)) %>%
  # only select SQC master sets
  filter(grepl('SQC', `Raw file`)) %>%
  # filter out non J/U sets (QC40E,F, QC42A)
  filter(!grepl('SQC9', `Raw file`))  %>%
  # select 1% protein fdr
  filter(!is.na(prot_fdr)) %>%
  filter(prot_fdr < 0.01)

ev_a <- ev.f %>% filter(qval < 0.01)
ev_b <- ev.f %>% filter(qval > 0.01 & qval_updated < 0.01)
ev_c <- ev.f %>% filter(qval_updated < 0.01)


# peptides/PSMs per protein -----------------------------------------------

peps_per_prot <- function(e) {
  e %>%
    group_by(Protein) %>%
    summarise(peptides=length(unique(`Sequence`)),
              psms=log2(n()))
}

p_a <- peps_per_prot(ev_a) %>% mutate(group=1)
p_b <- peps_per_prot(ev_b) %>% mutate(group=2)
p_c <- peps_per_prot(ev_c) %>% mutate(group=3)

# for DART-ID new proteins (p_b), remove proteins from Spectra (p_a)
p_b <- p_b %>% filter(!Protein %in% p_a$Protein)

# combine data
p_all <- rbind(p_a, p_b, p_c)
# ceiling peptides/protein to 10
p_all <- p_all %>% mutate_at('peptides', funs(ifelse(. > 10, 10, .)))
# factorize group
p_all$group <- factor(p_all$group, labels=c('Spectra', 'DART-ID new proteins', 'DART-ID all proteins'))

# Plotting ----------------------------------------------------------------

# function to increase vertical spacing between legend keys
# @clauswilke
# trick from https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2/26971729
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.75, "npc"),
    height = grid::unit(0.75, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3

p <- 
ggplot(p_all) +
  geom_histogram(aes(peptides, fill=group), position='dodge', binwidth=1, 
                 color='black', size=0.25) +
  scale_x_continuous(breaks=seq(1, 10), labels=c(seq(1,9), '10+')) +
  scale_y_continuous(limits=c(0, 700), breaks=seq(0, 700, by=100), expand=c(0,0)) + 
  scale_fill_manual(values=c(cb[1], paste0(cb[2], '88'), cb[2])) +
  #guides(fill=guide_legend(keywidth=0.5, keyheight=0.5, default.unit='cm')) +
  labs(x='Peptide Sequences per Protein', y='Count', fill=NULL) +
  theme_bert() + 
  theme(
    plot.margin=unit(c(0.25,0.25,0.25,0.25), 'cm'),
    legend.position=c(0.65, 0.8),
    legend.key.size = unit(0.75,'cm'),
    legend.text=element_text(size=10),
    strip.text = element_text(size=12),
    axis.text=element_text(size=12)
  )

ggsave('manuscript/Figs/peptides_per_prot.pdf', plot=p, device='pdf',
       width=3.5, height=3, units='in')

p <- 
ggplot(p_all) +
  geom_histogram(aes(psms, fill=group), position='dodge', binwidth=1, 
                 color='black', size=0.25) +
  scale_x_continuous(breaks=seq(0, 12)) +
  scale_y_continuous(limits=c(0, 500), expand=c(0,0)) + 
  scale_fill_manual(values=c(cb[1], paste0(cb[2], '88'), cb[2])) +
  labs(x='Log2 PSMs per Protein', y='Count', fill=NULL) +
  theme_bert() + 
  theme(
    plot.margin=unit(c(0.25,0.25,0.25,0.25), 'cm'),
    legend.position=c(0.65, 0.8),
    legend.key.size = unit(0.75,'cm'),
    legend.text=element_text(size=10),
    strip.text = element_text(size=12),
    axis.text=element_text(size=12)
  )

ggsave('manuscript/Figs/psms_per_prot.pdf', plot=p, device='pdf',
       width=3.5, height=3, units='in')
