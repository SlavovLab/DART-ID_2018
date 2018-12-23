library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## config -----------------------------------------------------------------------------------------

# for reference, all SQC experiments being used have the QC40F design:

# .____________________________________________________________________________.
# | 126C | 127N | 127C | 128N | 128C | 129N | 129C | 130N | 130C | 131N | 131C |
# .____________________________________________________________________________.
# | 100J | 100U | empty| empty|  1J  |  1U  |  1J  |  1U  |  1J  |  1U  | empty|
# .____________________________________________________________________________.

# MS_PREFIX <- 'G:/My Drive/MS/'
MS_PREFIX <- '/gd/MS/'

RI_COLUMN_EXPRESSION <- 'Reporter intensity corrected [0-9]+'
# RI_COLUMN_EXPRESSION <- 'Reporter intensity [0-9]+'

MIN_IDS_PER_EXP <- 300
FDR_THRESHOLD <- 0.01 # 0.1% FDR
PIF_THRESHOLD <- 0.8

# minimum correlation single cell channels must have to the carrier channel
SC_COR_THRESHOLD <- 0.3

# protein or peptide-level data?
DATA_COLLAPSE  = 'PEPTIDE' # = 'PROTEIN' #

# KNN Nearest Neighbors
KNN_N = 10

## load cormat data -----

ev.f <- ev %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))]) %>%
  # exclude CON, REV proteins
  filter(!grepl('CON__|REV__', `Leading razor protein`)) %>%
  # only select SQC master sets
  filter(grepl('SQC', `Raw file`)) %>%
  # filter out non J/U sets (QC40E,F, QC42A)
  filter(!grepl('SQC9', `Raw file`)) %>%
  # select 1% protein fdr
  filter(!is.na(prot_fdr)) %>%
  filter(prot_fdr < 0.01) %>%
  # remove duplicate file-sequence pairs
  # distinct(`Raw file`, `Modified sequence`, .keep_all=T) %>%
  # remove raw files w/ less than 300 observations
  # filter(!`Raw file` %in% names(table(`Raw file`))[ table(`Raw file`) < MIN_IDS_PER_EXP ]) %>%
  # get UniProt accession ID
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })) %>%
  dplyr::select('Raw file', 'Modified sequence', 'Protein', 'qval', 'qval_updated',
         matches(RI_COLUMN_EXPRESSION)) %>%
  # remove empty channels, carrier channels
  dplyr::select(-grep('Reporter intensity corrected [0-3]', colnames(ev.f)))

# Separate into Spectra and DART PSM sets ---------------------------------

ev_a <- ev.f %>%
  filter(qval < 0.01)

ev_b <- ev.f %>%
  filter(qval_updated < 0.01)


# Normalize RI data within each experiment --------------------------------

# do separate normalizations for sets A and B

ev_norm_a <- data.frame()
ev_norm_b <- data.frame()
exps <- sort(unique(ev.f$`Raw file`))

for(i in 1:length(exps)) {
  e <- exps[i]
  print(e)
  
  # get subset of ev for this experiment
  .ev_a <- ev_a %>% filter(`Raw file` == e)
  .ev_b <- ev_b %>% filter(`Raw file` == e)

  # column indicecs of RI intensities
  ri_cols <- grep(RI_COLUMN_EXPRESSION, colnames(.ev_a))
  
  # normalize data
  .ev_a <- normalize_ri_data_table(.ev_a, ri_cols)
  .ev_b <- normalize_ri_data_table(.ev_b, ri_cols)
  
  # attach to output DF
  ev_norm_a <- rbind(ev_norm_a, .ev_a)
  ev_norm_b <- rbind(ev_norm_b, .ev_b)
}
rm(.ev_a, .ev_b)


# spread data out into an expression matrix -------------------------------

ev_spread_a <- ev_norm_a %>%
  # collapse by protein. normalize within each experiment
  group_by(Protein, `Raw file`) %>%
  summarise_at(grep(RI_COLUMN_EXPRESSION, colnames(ev_norm_a)), funs(median)) %>%
  ungroup() %>%
  gather('Channel', 'Intensity', -c(`Raw file`, Protein)) %>%
  mutate(Channel=gsub('[^0-9]+\\s([0-9]+)', '\\1', Channel)) %>%
  mutate(id=paste0(`Raw file`, '_', Channel)) %>%
  select(-c(`Raw file`, Channel)) %>%
  # remove quantitation with values of NaN, Inf, or -Inf
  mutate(Intensity=ifelse(is.infinite(Intensity) | is.na(Intensity), NA, Intensity)) %>%
  spread('id', 'Intensity')

ev_spread_b <- ev_norm_b %>%
  # collapse by protein. normalize within each experiment
  group_by(Protein, `Raw file`) %>%
  summarise_at(grep(RI_COLUMN_EXPRESSION, colnames(ev_norm_b)), funs(median)) %>%
  ungroup() %>%
  gather('Channel', 'Intensity', -c(`Raw file`, Protein)) %>%
  mutate(Channel=gsub('[^0-9]+\\s([0-9]+)', '\\1', Channel)) %>%
  mutate(id=paste0(`Raw file`, '_', Channel)) %>%
  select(-c(`Raw file`, Channel)) %>%
  # remove quantitation with values of NaN, Inf, or -Inf
  mutate(Intensity=ifelse(is.infinite(Intensity) | is.na(Intensity), NA, Intensity)) %>%
  spread('id', 'Intensity')

# average J vs. U protein levels, take ratios -----------------------------

# its just the way the above "spread" matrices are structured
# that even columns (where protein column is column 1) are J cells
# and odd are U cells

N <- ncol(ev_spread_a)
P <- nrow(ev_spread_a)

# get the fold change, between the medians of relative J and U levels
fc_a <- apply(ev_spread_a[,seq(2,N,2)], 1, median, na.rm=T) / 
          apply(ev_spread_a[,seq(3,N,2)], 1, median, na.rm=T)
# get the significance of the change, based on a t-test run for each
sig_a <- vector(length=P)
for(i in 1:P) {
  sig_a[i] <- t.test(x=as.numeric(ev_spread_a[i,seq(2,N,2)]), 
                     y=as.numeric(ev_spread_a[i,seq(3,N,2)]),
                     alternative='two.sided',
                     na.action='na.omit')$p.value
}
# apply bonferroni correction on p-values
# sig_a <- sig_a * P
# apply FDR correction on p-values
sig_a <- (cumsum(sig_a[order(sig_a)]) / seq(1, length(sig_a)))[order(order(sig_a))]
# set 0 to minimum value
sig_a[sig_a == 0] <- .Machine$double.xmin
# put ceiling on significance at 20
sig_a[-log10(sig_a) > 20] <- 1e-20

plot(log2(fc_a), -log10(sig_a))
abline(a=2, b=0, col='red')

# count number of significant observations above 1 fold-change
sum(abs(log2(fc_a)) > 1 & -log10(sig_a) > 2)

N <- ncol(ev_spread_b)
P <- nrow(ev_spread_b)

fc_b <- apply(ev_spread_b[,seq(2,N,2)], 1, median, na.rm=T) / 
  apply(ev_spread_b[,seq(3,N,2)], 1, median, na.rm=T)

sig_b <- vector(length=P)
for(i in 1:P) {
  sig_b[i] <- t.test(x=as.numeric(ev_spread_b[i,seq(2,N,2)]), 
                     y=as.numeric(ev_spread_b[i,seq(3,N,2)]),
                     alternative='two.sided',
                     na.action='na.omit')$p.value
}
# sig_b <- sig_b * P
sig_b <- (cumsum(sig_b[order(sig_b)]) / seq(1, length(sig_b)))[order(order(sig_b))]
# set 0 to minimum value
sig_b[sig_b == 0] <- .Machine$double.xmin
# put ceiling on significance at 20
sig_b[-log10(sig_b) > 20] <- 1e-20


# make figure -------------------------------------------------------------

sig_points_b <- rep(1, length(fc_b))
sig_points_b[-log10(sig_b) <= 2] <- 2
sig_points_b[(abs(log2(fc_b)) > 1) & (-log10(sig_b) > 2)] <- 3

cols <- c(rgb(0,0,0,0.5), rgb(1,0,0,0.5), rgb(0,0,1,0.5))

plot(log2(fc_b), -log10(sig_b), 
     col=cols[as.numeric(sig_points_b)],
     pch=16, cex=1, 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)
abline(a=2, b=0, col='red')
abline(v=c(-1, 1), col='blue')

axis(1, at=seq(-4, 4, by=1), mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 20, by=5), labels=c(seq(0, 15, by=5), '>20'))

# count number of significant observations above 1 fold-change
sum(abs(log2(fc_b)) > 0.5 & -log10(sig_b) > 2)
