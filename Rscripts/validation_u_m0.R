# validate boosted PSMs with U/Macrophage Differentiation System ----------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

# Load MaxQuant Output ----------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/FP60_64_m0_20181121/ev_updated.txt')
#ev <- read_tsv('/gd/MS/SCoPE/mPOP/dat/FP60-64/evidence.txt')

# Load Olfactory Receptor Protein IDs -------------------------------------

# load olfactory receptor FASTA file
olfactory_fasta <- read.fasta('/gd/bayesian_RT/dat/uniprot-olfactory.fasta', seqtype='AA')
# extract protein IDs
olfactory_ids <- as.character(sapply(olfactory_fasta, function(l) { attr(l, 'name') }))
# extract uniprot ID from ID string
olfactory_ids <- str_match(olfactory_ids, '\\|([0-9A-Z\\-]+)\\|')[,2]


# Look for peptides belonging to these proteins ---------------------------

ev.f <- ev %>%
  # remove decoys and contaminants
  filter(!grepl('CON', Proteins)) %>%
  filter(!grepl('REV', `Leading razor protein`)) %>%
  
  # only take unique peptides
  filter(!grepl('\\;', Proteins)) %>%
  
  # extract uniprot IDs from leading razor protein
  mutate(Protein=gsub('^sp\\||\\|[A-Z0-9\\_]+$', '', `Leading razor protein`)) %>%
  
  # match olfactory receptor protein IDs
  filter(Protein %in% olfactory_ids) %>%
  
  dplyr::select(c('Modified sequence', 'Protein', 'Proteins', 'Raw file', 
                  'PEP', 'pep_updated', 'residual', 'Retention time', 'muij')) %>%
  
  # fix weird PEPs
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .)))

ev.sqc <- ev.f %>% 
  filter(grepl('SQC', `Raw file`)) %>%
  filter(!is.na(residual)) %>%
  mutate(dPEP=-log10(pep_updated/PEP))

ev.fp <- ev.f %>% 
  filter(grepl('FP', `Raw file`)) %>%
  filter(!is.na(residual)) %>%
  mutate(dPEP=-log10(pep_updated/PEP))

hist(ev.sqc$dPEP)
hist(ev.fp$dPEP)
