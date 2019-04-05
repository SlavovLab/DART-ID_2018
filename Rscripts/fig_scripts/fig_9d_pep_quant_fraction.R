# init, load data ---------------------------------------------------------

library(tidyverse)
library(viridisLite)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## config -----------------------------------------------------------------------------------------

# for reference, all SQC experiments being used have the QC40F design:

# .____________________________________________________________________________.
# | 126C | 127N | 127C | 128N | 128C | 129N | 129C | 130N | 130C | 131N | 131C |
# .____________________________________________________________________________.
# | 50J  | 50U  | empty| empty|  1J  |  1U  |  1J  |  1U  |  1J  |  1U  | empty|
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

## load cormat data ----------------------------------------------------

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
                matches(RI_COLUMN_EXPRESSION))

ev.f <- ev.f %>%
  # remove empty channels, carrier channels
  dplyr::select(-grep('Reporter intensity corrected [0-3]', colnames(ev.f)))

# Separate into Spectra and DART PSM sets ---------------------------------

ev_a <- ev.f %>%
  filter(qval < 0.01)

ev_b <- ev.f %>%
  filter(qval_updated < 0.01)

ev_c <- ev.f %>%
  filter(qval > 0.01 & qval_updated < 0.01)


# Normalize RI data within each experiment --------------------------------

# do separate normalizations for sets A and B

ev_norm_a <- data.frame()
ev_norm_b <- data.frame()
ev_norm_c <- data.frame()
exps <- sort(unique(ev.f$`Raw file`))

for(i in 1:length(exps)) {
  e <- exps[i]
  print(e)
  
  # get subset of ev for this experiment
  .ev_a <- ev_a %>% filter(`Raw file` == e)
  .ev_b <- ev_b %>% filter(`Raw file` == e)
  .ev_c <- ev_c %>% filter(`Raw file` == e)
  
  # column indicecs of RI intensities
  ri_cols <- grep(RI_COLUMN_EXPRESSION, colnames(.ev_a))
  
  # normalize data
  .ev_a <- normalize_ri_data_table(.ev_a, ri_cols)
  .ev_b <- normalize_ri_data_table(.ev_b, ri_cols)
  .ev_c <- normalize_ri_data_table(.ev_c, ri_cols)
  
  # attach to output DF
  ev_norm_a <- rbind(ev_norm_a, .ev_a)
  ev_norm_b <- rbind(ev_norm_b, .ev_b)
  ev_norm_c <- rbind(ev_norm_c, .ev_c)
}
rm(.ev_a, .ev_b, .ev_c)

# missing data ------------------------------------------------------------

ev_norm_a$TotalRI <- apply(data.matrix(ev_norm_a[,grep(RI_COLUMN_EXPRESSION, colnames(ev_norm_a))]), 1, sum, na.rm=T)
ev_spread_a <- ev_norm_a %>%
  group_by(`Modified sequence`, `Raw file`) %>%
  summarise(n=sum(TotalRI)) %>%
  ungroup() %>%
  gather('Channel', 'Intensity', -c(`Raw file`, `Modified sequence`)) %>%
  dplyr::select(-c(Channel)) %>%
  # remove quantitation with values of NaN, Inf, or -Inf
  mutate(Intensity=ifelse(is.infinite(Intensity) | is.na(Intensity), NA, Intensity)) %>%
  spread('Raw file', 'Intensity')

ev_norm_b$TotalRI <- apply(data.matrix(ev_norm_b[,grep(RI_COLUMN_EXPRESSION, colnames(ev_norm_b))]), 1, sum, na.rm=T)
ev_spread_b <- ev_norm_b %>%
  # collapse by protein. normalize within each experiment
  group_by(`Modified sequence`, `Raw file`) %>%
  summarise(n=sum(TotalRI)) %>%
  ungroup() %>%
  gather('Channel', 'Intensity', -c(`Raw file`, `Modified sequence`)) %>%
  dplyr::select(-c(Channel)) %>%
  # remove quantitation with values of NaN, Inf, or -Inf
  mutate(Intensity=ifelse(is.infinite(Intensity) | is.na(Intensity), NA, Intensity)) %>%
  spread('Raw file', 'Intensity')

common_seqs <- intersect(ev_spread_a$`Modified sequence`, ev_spread_b$`Modified sequence`)

missing_a <- apply(is.na(ev_spread_a %>% filter(`Modified sequence` %in% common_seqs) %>% arrange(`Modified sequence`) %>% dplyr::select(-'Modified sequence')), 1, sum) / (ncol(ev_spread_a)-1)
missing_b <- apply(is.na(ev_spread_b %>% filter(`Modified sequence` %in% common_seqs) %>% arrange(`Modified sequence`) %>% dplyr::select(-'Modified sequence')), 1, sum) / (ncol(ev_spread_b)-1)

# plot missing data per peptide -------------------------------------------------------

k <- 100
contour_cols <- viridis(k, alpha=1)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
dens <- get_density(1-missing_a, 1-missing_b, k)

pdf(file='manuscript/Figs/missing_data_per_pep.pdf', width=3.5, height=3)

par(mar=c(2.7, 3, 0.25, 1), oma=c(0.21, 0, 2.8, 0), cex.axis=1)

layout(t(c(1, 2)), widths=c(5, 1))

plot(0, 0, type='n', xaxt='n', yaxt='n', xlim=c(0, 1), ylim=c(0, 1), xlab=NA, ylab=NA, main=NA)
points(1-missing_a, 1-missing_b,
       pch=16, cex=0.5,
       col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))])
abline(a=0, b=1, col='red', lwd=2)

text(0.97, 0.15, paste0('', length(common_seqs), ' peptides',
                        '\n', '', ncol(ev_spread_a)-1, ' experiments'), 
     adj=c(1, 0.5), cex=1)

axis(1, at=seq(0, 1, by=0.2), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 1, by=0.2), tck=-0.02, mgp=c(0, 0.5, 0), las=1)

mtext('Spectra', side=1, line=1.5, cex=1)
mtext('DART-ID', side=2, line=2, cex=1)
mtext('Fraction of Experiments\nWhere Peptide is Quantified', side=3, line=0.15, cex=1, font=2, outer=T)

# colorbar

par(mar=c(2.7, 0.5, 1.25, 1.25), pty='m')
image(matrix(seq(0, 1, length.out=k), ncol=k), col=contour_cols,
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
mtext('Density', side=3, line=0.25)

dev.off()
