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
  mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev)))[order(order(PEP))],
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


# spread data out into an expression matrix -------------------------------

ev_spread_a <- ev_norm_a %>%
  # collapse by protein. normalize within each experiment
  group_by(Protein, `Raw file`) %>%
  summarise_at(grep(RI_COLUMN_EXPRESSION, colnames(ev_norm_a)), funs(median)) %>%
  ungroup() %>%
  gather('Channel', 'Intensity', -c(`Raw file`, Protein)) %>%
  mutate(Channel=gsub('[^0-9]+\\s([0-9]+)', '\\1', Channel)) %>%
  mutate(id=paste0(`Raw file`, '_', Channel)) %>%
  dplyr::select(-c(`Raw file`, Channel)) %>%
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
  dplyr::select(-c(`Raw file`, Channel)) %>%
  # remove quantitation with values of NaN, Inf, or -Inf
  mutate(Intensity=ifelse(is.infinite(Intensity) | is.na(Intensity), NA, Intensity)) %>%
  spread('id', 'Intensity')

ev_spread_c <- ev_norm_c %>%
  # collapse by protein. normalize within each experiment
  group_by(Protein, `Raw file`) %>%
  summarise_at(grep(RI_COLUMN_EXPRESSION, colnames(ev_norm_c)), funs(median)) %>%
  ungroup() %>%
  gather('Channel', 'Intensity', -c(`Raw file`, Protein)) %>%
  mutate(Channel=gsub('[^0-9]+\\s([0-9]+)', '\\1', Channel)) %>%
  mutate(id=paste0(`Raw file`, '_', Channel)) %>%
  dplyr::select(-c(`Raw file`, Channel)) %>%
  # remove quantitation with values of NaN, Inf, or -Inf
  mutate(Intensity=ifelse(is.infinite(Intensity) | is.na(Intensity), NA, Intensity)) %>%
  spread('id', 'Intensity')

# compare ratios of common proteins ---------------------------------------

common_prots <- intersect(ev_spread_a$Protein, ev_spread_c$Protein)

c_a <- ev_spread_a %>% filter(Protein %in% common_prots)
c_c <- ev_spread_c %>% filter(Protein %in% common_prots)

N <- ncol(c_a)
P <- nrow(c_a)

fc_a <- apply(c_a[,seq(2,N,2)], 1, mean, na.rm=T) / 
  apply(c_a[,seq(3,N,2)], 1, mean, na.rm=T)

N <- ncol(c_c)
P <- nrow(c_c)

fc_c <- apply(c_c[,seq(2,N,2)], 1, mean, na.rm=T) / 
  apply(c_c[,seq(3,N,2)], 1, mean, na.rm=T)


# ratio scatter plot ------------------------------------------------------

pdf(file='~/git/DART-ID_2018/manuscript/Figs/ratio_scatter_v5.pdf', width=5, height=5)

par(mar=c(3,3,3,1), pty='s', mgp=c(0, 0.3, 0), cex.axis=1)

plot(0, 0, type='n',
     xlim=c(-3.5, 3.5), ylim=c(-3.5, 3.5),
     xaxt='n', yaxt='n',
     xlab=NA, ylab=NA)
abline(a=0, b=1, col='red', lty=1, lwd=2)

k <- 60
contour_cols <- viridis(k, alpha=0.5)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
dens <- get_density(log2(fc_a), log2(fc_c), k)

points(log2(fc_a), log2(fc_c), 
       #col=rgb(0,0,0, 0.3), 
       col=contour_cols[findInterval(dens, seq(0, max(dens), length.out=k))],
       pch=16)

cor_text <- cor(log2(fc_a), log2(fc_c))
cor_text <- formatC(cor_text, digits=3)
text(x=-3, y=3, adj=c(0, 0.5),
     label=bquote(.(as.name('rho'))*.(' = ')*.(cor_text)),
     cex=1.2)

axis(1, at=seq(-3, 3), tck=-0.02, mgp=c(0, 0.5, 0))
axis(2, at=seq(-3, 3), tck=-0.02, las=1, mgp=c(0, 0.6, 0))

mtext('Spectra', side=1, line=1.5, cex=1)
mtext('DART-ID', side=2, line=1.5, cex=1)
mtext(parse(text='Log[2]~paste(T,"-",Cell)/Monocyte~Ratio~(Fold~Change)'))

dev.off()
