# init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

# load data ---------------------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

# false positives at different FDR levels ---------------------------------

# ceiling PEP to 1
ev.f <- ev %>%
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .)))

# get sequence of FDR thresholds, log-spaced
x <- logseq(1e-5, 1, 100)

# results matrix - first col = spectral false positives,
# second col = DART false positives
fp_mat <- matrix(nrow=length(x), ncol=2)

for(i in 1:length(x)) {
  thresh <- x[i]
  fp_mat[i,1] <- sum(ev.f$PEP[ev.f$PEP <= thresh])
  fp_mat[i,2] <- sum(ev.f$pep_updated[ev.f$PEP <= thresh])
}

plot(log10(x), fp_mat[,1], type='l', col='black')
lines(log10(x), fp_mat[,2], type='l', col='red')


# oh no -------------------------------------------------------------------

a <- ev %>%
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  group_by(`Raw file`) %>%
  summarise(fp_spectra=sum(PEP),
            fp_dart=sum(pep_updated))

plot(a$fp_spectra, a$fp_dart)
abline(a=0, b=1, col='red')



# load data ---------------------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/FP18_human_yeast_20190104/ev_updated.txt')

human <- read.fasta(file='/gd/MS/DBs_Params/FASTA/swissprot_human_20181210.fasta', seqtype='AA')

yeast_seqs <- ev %>% 
  filter(grepl('YEAST', `Leading razor protein`)) %>%
  filter(!grepl('REV__', `Leading razor protein`)) %>%
  filter(!is.na(pep_new)) %>%
  #dplyr::select(c('Modified sequence', 'Leading razor protein', 'PEP', 'pep_updated'))
  pull(Sequence)

a <- ev %>% 
  filter(grepl('YEAST', `Leading razor protein`)) %>%
  filter(!grepl('REV__', `Leading razor protein`)) %>%
  filter(!is.na(pep_new)) %>%
  dplyr::select(c('Modified sequence', 'Raw file', 'Leading razor protein', 'PEP', 'pep_new', 'Retention time'))

yeast_seqs <- unique(yeast_seqs)

b <- lapply(human, function(entry) {
  #any(!is.na(str_locate(entry, yeast_seqs)))
  which(!is.na(str_locate(entry, yeast_seqs)[,1]))
})

matches <- sort(unique(unlist(b)))
