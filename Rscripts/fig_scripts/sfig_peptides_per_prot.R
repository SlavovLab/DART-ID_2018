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

# proteins with only 1 peptide? -------------------------------------------

peps_per_prot <- ev.f %>%
  # Spectra PSMs = 1
  # DART-ID PSMs = 2
  # all other are 0
  mutate(group=as.numeric(qval < 0.01)) %>%
  mutate(group=group + as.numeric((qval > 0.01 & qval_updated < 0.01) * 2)) %>%
  # remove group 0
  filter(group > 0) %>%
  mutate(group=factor(group, labels=c('Spectra', 'DART-ID'))) %>%
  group_by(Protein, group) %>%
  summarise(n=length(unique(`Modified sequence`)),
            m=n()) %>%
  # cap n at 10
  mutate_at('n', funs(ifelse(. > 10, 10, .))) %>%
  # log2 transform count
  mutate(m=log2(m)) %>%
  arrange(desc(n))

# remove proteins from DART-ID that are present already in Spectra
orig_prots <- peps_per_prot %>% filter(group == 'Spectra') %>% pull(Protein)
peps_per_prot <- peps_per_prot %>%
  filter(!(group == 'DART-ID' & Protein %in% orig_prots))


# Plotting ----------------------------------------------------------------

p <- 
ggplot(peps_per_prot) +
  geom_histogram(aes(n, fill=group), position='dodge', binwidth=1, 
                 color='black', size=0.25) +
  facet_wrap(group~., nrow=2, ncol=1) +
  scale_x_continuous(breaks=seq(1, 10), labels=c(seq(1,9), '10+')) +
  scale_y_continuous(limits=c(0, 700), expand=c(0,0)) + 
  scale_fill_manual(values=paste0(c(cb[1], cb[2]), 'FF'), guide=F) +
  labs(x='Peptide Sequences per Protein', y='Count', fill=NULL) +
  theme_bert() + 
  theme(
    strip.text = element_text(size=12)
  )

ggsave('manuscript/Figs/peptides_per_prot.pdf', plot=p, device='pdf',
       width=3.5, height=4, units='in')

p <- 
ggplot(peps_per_prot) +
  geom_histogram(aes(m, fill=group), position='dodge', binwidth=1, 
                 color='black', size=0.25) +
  facet_wrap(group~., nrow=2, ncol=1) +
  scale_x_continuous(breaks=seq(0, 12)) +
  scale_y_continuous(limits=c(0, 500), expand=c(0,0)) + 
  scale_fill_manual(values=paste0(c(cb[1], cb[2]), 'FF'), guide=F) +
  labs(x='Log2 PSMs per Protein', y='Count', fill=NULL) +
  theme_bert() + 
  theme(
    strip.text = element_text(size=12)
  )

ggsave('manuscript/Figs/psms_per_prot.pdf', plot=p, device='pdf',
       width=3.5, height=4, units='in')
