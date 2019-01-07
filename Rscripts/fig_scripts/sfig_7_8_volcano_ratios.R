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

# average J vs. U protein levels, take ratios -----------------------------

# its just the way the above "spread" matrices are structured
# that even columns (where protein column is column 1) are J cells
# and odd are U cells

N <- ncol(ev_spread_a)
P <- nrow(ev_spread_a)

# get the fold change, between the medians of relative J and U levels
fc_a <- apply(ev_spread_a[,seq(2,N,2)], 1, mean, na.rm=T) / 
  apply(ev_spread_a[,seq(3,N,2)], 1, mean, na.rm=T)
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
# put ceiling on significance at 10
sig_a[-log10(sig_a) > 10] <- 1e-10

plot(log2(fc_a), -log10(sig_a))
abline(a=2, b=0, col='red')

# count number of significant observations above 1 fold-change
sum(abs(log2(fc_a)) > 1 & -log10(sig_a) > 2)

N <- ncol(ev_spread_b)
P <- nrow(ev_spread_b)

fc_b <- apply(ev_spread_b[,seq(2,N,2)], 1, mean, na.rm=T) / 
  apply(ev_spread_b[,seq(3,N,2)], 1, mean, na.rm=T)

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
# put ceiling on significance at 10
sig_b[-log10(sig_b) > 10] <- 1e-10

# count number of significant observations above 1 fold-change
sum(abs(log2(fc_b)) > 1 & -log10(sig_b) > 2)

N <- ncol(ev_spread_c)
P <- nrow(ev_spread_c)

fc_c <- apply(ev_spread_c[,seq(2,N,2)], 1, mean, na.rm=T) / 
  apply(ev_spread_c[,seq(3,N,2)], 1, mean, na.rm=T)

sig_c <- vector(length=P)
for(i in 1:P) {
  sig_c[i] <- t.test(x=as.numeric(ev_spread_c[i,seq(2,N,2)]), 
                     y=as.numeric(ev_spread_c[i,seq(3,N,2)]),
                     alternative='two.sided',
                     na.action='na.omit')$p.value
}
# sig_c <- sig_c * P
sig_c <- (cumsum(sig_c[order(sig_c)]) / seq(1, length(sig_c)))[order(order(sig_c))]
# set 0 to minimum value
sig_c[sig_c == 0] <- .Machine$double.xmin
# put ceiling on significance at 10
sig_c[-log10(sig_c) > 10] <- 1e-10

# count number of significant observations above 1 fold-change
sum(abs(log2(fc_c)) > 1 & -log10(sig_c) > 2)

# make figure -------------------------------------------------------------

pdf(file='~/git/DART-ID_2018/manuscript/Figs/volcanoes.pdf', width=3.5, height=6)

layout(rbind(c(1), c(2)))

# Volcano A

par(mar=c(3, 3, 3, 1), mgp=c(0, 0.1, 0),
    cex.axis=1.1)

sig_points_a <- rep(1, length(fc_a))
sig_points_a[-log10(sig_a) <= 2] <- 2
sig_points_a[(abs(log2(fc_a)) > 1) & (-log10(sig_a) > 2)] <- 3

cols <- c(rgb(0,0,0,0.5), rgb(0,0,0,0.5), rgb(0,0,1,0.5))

plot(log2(fc_a), -log10(sig_a), 
     #col=cols[as.numeric(sig_points_b)],
     col=paste0(cb[1], '66'),
     pch=16, cex=1, 
     xlim=c(-3.5, 3.5), 
     #ylim=c(0, 10),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)
abline(v=0, col='black', lty=2, lwd=2)


axis(1, at=seq(-4, 4, by=1), mgp=c(0, 0.3, 0), tck=-0.02)
axis(2, at=seq(0, 10, by=2), labels=c(seq(0, 8, by=2), '>10'), 
     mgp=c(0, 0.5, 0), tck=-0.02, las=1)

mtext(parse(text='Log[2]~Ratio~(Fold~Change)'), side=1, line=1.5, cex=1)
mtext(parse(text='-Log[10]~~paste(q,"-",value)'), side=2, line=1.5, cex=1)
mtext('Spectra PSMs', side=3, line=0.5, cex=1, font=2)

# Volcano B

par(mar=c(3, 3, 3, 1), mgp=c(0, 0.1, 0),
    cex.axis=1.1)

sig_points_b <- rep(1, length(fc_b))
sig_points_b[-log10(sig_b) <= 2] <- 2
sig_points_b[(abs(log2(fc_b)) > 1) & (-log10(sig_b) > 2)] <- 3

cols <- c(rgb(0,0,0,0.5), rgb(0,0,0,0.5), rgb(0,0,1,0.5))

plot(log2(fc_b), -log10(sig_b), 
     #col=cols[as.numeric(sig_points_b)],
     col=paste0(cb[2], '66'),
     pch=16, cex=1, 
     xlim=c(-3.5, 3.5), 
     #ylim=c(0, 10),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)
abline(v=0, col='black', lty=2, lwd=2)


axis(1, at=seq(-4, 4, by=1), mgp=c(0, 0.3, 0), tck=-0.02)
axis(2, at=seq(0, 10, by=2), labels=c(seq(0, 8, by=2), '>10'), 
     mgp=c(0, 0.5, 0), tck=-0.02, las=1)

mtext(parse(text='Log[2]~Ratio~(Fold~Change)'), side=1, line=1.5, cex=1)
mtext(parse(text='-Log[10]~~paste(q,"-",value)'), side=2, line=1.5, cex=1)
mtext('DART-ID + Spectra PSMs', side=3, line=0.5, cex=1, font=2)

dev.off()


# Line Plot ---------------------------------------------------------------------

pdf(file='~/git/DART-ID_2018/manuscript/Figs/volcano_line.pdf', width=3.5, height=3)

par(mar=c(3, 4, 1, 1), mgp=c(0, 0.3, 0),
    cex.axis=1)

x <- seq(-10, 0, length.out=100)
nprots_a <- vector(length=length(x))
nprots_b <- vector(length=length(x))
for(i in 1:length(x)) {
  .x <- x[i]
  nprots_a[i] <- sum(sig_a <= (10 ** .x))
  nprots_b[i] <- sum(sig_b <= (10 ** .x))
}

plot(0, 0, type='n',
     xlim=c(-10, 0), ylim=c(0, 1400),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, main=NA)

lines(x, nprots_a, col=cb[1], lty=1, lwd=3)
lines(x, nprots_b, col=cb[2], lty=1, lwd=3)

xlabs <- fancy_scientific(10**seq(-10, 0, by=1))
xlabs[seq(2,10,by=2)] <- NA # skip every other label
#xlabs[1] <- parse(text='phantom(0)<10^-10')

axis(1, at=seq(-10, 0, by=1), labels=xlabs,
     mgp=c(0, 0.5, 0), tck=-0.02)
axis(2, at=seq(0, 1400, by=200),
     mgp=c(0, 0.5, 0), tck=-0.02, las=1)

mtext('False Discovery Rate, q-value', side=1, line=1.5, cex=1)
mtext('# Differentially Abundant Proteins', side=2, line=2.75, cex=1)

legend('topleft', c('Spectra PSMs', 'DART-ID + Spectra\nPSMs'),
       lwd=3, lty=1, col=c(cb[1], cb[2]), seg.len=1,
       bty='n', cex=1, x.intersp=0.6, y.intersp=1.2, inset=c(0.01, -0.04))


dev.off()


# metrics for the table ---------------------------------------------------

# total number of proteins
length(sig_a)
length(sig_b)

# number of proteins at 0.1% FDR
sum(sig_a <= 1e-3)
sum(sig_b <= 1e-3)

# number of proteins at 1% FDR
sum(sig_a <= 1e-2)
sum(sig_b <= 1e-2)

# % proteins at 0.1% FDR
sum(sig_a <= 1e-3) / length(sig_a) * 100
sum(sig_b <= 1e-3) / length(sig_b) * 100

# % proteins at 1% FDR
sum(sig_a <= 1e-2) / length(sig_a) * 100
sum(sig_b <= 1e-2) / length(sig_b) * 100


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
