library(tidyverse)
library(ggridges)
library(pracma)
source('Rscripts/lib.R')

#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt')
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## load cormat data -----

source('Rscripts/validation_cormats.R')

# Precursor Area, PIF dists -----------------------------------------------

ev_a <- ev.fi %>% filter(qval < 0.01)
ev_b <- ev.fi %>% filter(qval > 0.01 & qval_updated < 0.01)

# count missing data
dcols <- grep('corrected', colnames(ev.fi))[5:10]
dmat <- data.matrix(ev_a %>% select(colnames(ev.fi)[dcols]))
missing_a <- apply(dmat==0, 1, sum)
dmat <- data.matrix(ev_b %>% select(colnames(ev.fi)[dcols]))
missing_b <- apply(dmat==0, 1, sum)

# poor quant v2 -----------------------------------------------------------

a_intens <- log10(ev_a$Intensity[!is.na(ev_a$Intensity)])
b_intens <- log10(ev_b$Intensity[!is.na(ev_b$Intensity)])
# equal/unequal variance doesn't make a difference... but let's just do it right anyways
var.test(a_intens, b_intens, alternative='two.sided')
t.test(a_intens, b_intens, var.equal=F) # p < 2.2e-16. three stars.

a_pif <- ev_a$PIF[!is.na(ev_a$PIF)]
b_pif <- ev_b$PIF[!is.na(ev_b$PIF)]
var.test(a_pif, b_pif, alternative='two.sided')
t.test(a_pif, b_pif, var.equal=F) # p < 2.2e-16. three stars

a_mc <- ev_a$`Missed cleavages`
b_mc <- ev_b$`Missed cleavages`
var.test(a_mc, b_mc, alternative='two.sided')
t.test(a_mc, b_mc, var.equal=F) # p < 2.2e-16. three stars

var.test(missing_a, missing_b)
t.test(missing_a, missing_b, var.equal=F) # p < 2.2e-16, three stars

# plot --------------------------------------------------------------------

pdf(file='manuscript/Figs/poor_quant_v2.pdf', width=7, height=2)

layout(rbind(c(1, 2, 3, 4)))

par(cex.axis=1, 
    mar=c(3.5,2,1,2.5),
    oma=c(0,1,0,0))

boxplot(list(a_intens, b_intens),
        range=1.5, col=c(cb[1], cb[2]), ylim=c(4.75, 8.4),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)

# significance stars
lines(c(1, 2), c(8, 8), col='black')
# brackets
lines(c(1, 1), c(7.9, 8), col='black')
lines(c(2, 2), c(7.9, 8), col='black')
text(1.5, 8, '***', adj=c(0.5, 0.2), cex=2)

#axis(1, at=c(1, 2), labels=c('Spectra', 'DART-ID'), tck=-0.02, mgp=c(0, 0.2, 0))
axis(1, at=c(1, 2), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(c(1, 2), rep(4.3, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
axis(2, at=seq(4, 10), tck=-0.02, mgp=c(0, 0.5, 0), las=1)
mtext(expression('log'[10]*' Precursor Ion Area'), side=2, line=1.25, cex=0.85)



boxplot(list((ev_a$PIF), (ev_b$PIF)), 
        range=1.5, col=c(cb[1], cb[2]), ylim=c(0.55, 1.11),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)

# significance stars
lines(c(1, 2), c(1.05, 1.05), col='black')
# brackets
lines(c(1, 1), c(1.035, 1.05), col='black')
lines(c(2, 2), c(1.035, 1.05), col='black')
text(1.5, 1.05, '***', adj=c(0.5, 0.2), cex=2)

axis(1, at=c(1, 2), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0.3, 1, by=0.1), label=seq(0.3, 1, by=0.1)*100, tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2), rep(0.485, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Precursor Ion, %'), side=2, line=2, cex=0.85)



barplot(c(mean(ev_a$`Missed cleavages`), mean(ev_b$`Missed cleavages`)), width=1, space=0.5,
        range=1.5, col=c(cb[1], cb[2]), xlim=c(0.25, 3.15), ylim=c(0, 0.24),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)

# significance stars
lines(c(1, 2.5), c(0.21, 0.21), col='black')
# brackets
lines(c(1, 1), c(0.203, 0.21), col='black')
lines(c(2.5, 2.5), c(0.203, 0.21), col='black')
text(1.75, 0.21, '***', adj=c(0.5, 0.2), cex=2)

axis(1, at=c(-10, 1, 2.5, 10), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 0.25, by=0.05), tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2.5), rep(-0.02, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Missed Cleavages, %'), side=2, line=2.25, cex=0.85)



barplot(c(mean(missing_a) / 6, mean(missing_b) / 6), width=1, space=0.5,
        range=1.5, col=c(cb[1], cb[2]), xlim=c(0.25, 3.15), ylim=c(0, 0.16),
        outcex=0, outpch=4, outcol=rgb(0,0,0,0.1),
        xaxt='n', yaxt='n', xlab=NA, ylab=NA)

# significance stars
lines(c(1, 2.5), c(0.14, 0.14), col='black')
# brackets
lines(c(1, 1), c(0.135, 0.14), col='black')
lines(c(2.5, 2.5), c(0.135, 0.14), col='black')
text(1.75, 0.14, '***', adj=c(0.5, 0.2), cex=2)

axis(1, at=c(-10, 1, 2.5, 10), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 1, by=0.05), labels=seq(0, 1, by=0.05)*100, 
     tck=-0.02, mgp=c(0, 0.5, 0), las=1)
text(c(1, 2.5), rep(-0.015, 2), c('Spectra', 'DART-ID'), xpd=T, srt=30, cex=1.15, adj=c(1, 0))
mtext(expression('Missing Data, %'), side=2, line=2, cex=0.85)

dev.off()


# PCA plots ---------

exps <- sort(unique(ev.f$`Raw file`))
# take 81 experiments 40-50, 80-90, 100-110
#exps <- exps[10:149]
# or exclude experiments with a lot of missing values
exps <- exps[-c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,22,27,28,33,38,39,54,55,56,57,60,61,62,
                64,67,68,69,76,77,78,80,81,94,95,102,103,108,110,111,112,113,114,115,116,
                117,120,121,122,123,124,125,126,127,128,129,130,132,133,137,139,147,149,150)]

# filter for the list of experiments
ev.fa <- ev.f %>% filter(`Raw file` %in% exps)

ev.fb <- ev.fa %>% filter(qval_updated < 0.01)
ev.fa <- ev.fa %>% filter(qval < 0.01)

prots.a <- ev.fa %>% 
  group_by(Protein) %>% 
  summarise(l=length(unique(`Raw file`))) %>% 
  filter(l >= length(exps)*0.95) %>%
  pull(Protein)

prots.b <- ev.fb %>% 
  group_by(Protein) %>% 
  summarise(l=length(unique(`Raw file`))) %>% 
  filter(l >= length(exps)*0.95) %>%
  pull(Protein)

prots.b <- prots.b[!prots.b %in% prots.a]

cat(length(prots.a), length(prots.b))

## -------

# select only proteins in protein list 
ev.fa <- ev.fa %>% filter(Protein %in% prots.a)
ev.fb <- ev.fb %>% filter(Protein %in% prots.b)

# number of peptides (only used in figure title)
peps.a <- unique(ev.fa$`Modified sequence`)
peps.b <- unique(ev.fb$`Modified sequence`)

dcols <- grep('Reporter intensity corrected', colnames(ev.fa))

dmat.fa <- ev.fa %>%
  # collapse data by sequence, by mean
  group_by(`Modified sequence`, Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fa)[dcols], funs(mean)) %>%
  # collapse data by protein, by mean
  group_by(Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fa)[dcols], funs(mean)) %>%
  ungroup()

dmat.fb <- ev.fb %>%
  group_by(`Modified sequence`, Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fb)[dcols], funs(mean)) %>%
  group_by(Protein, `Raw file`) %>%
  summarise_at(colnames(ev.fb)[dcols], funs(mean)) %>%
  ungroup()

dmat.fa <- gather(dmat.fa, k, v, -c(Protein, `Raw file`))
dmat.fa <- spread(dmat.fa, 'Raw file', 'v')

dmat.fb <- gather(dmat.fb, k, v, -c(Protein, `Raw file`))
dmat.fb <- spread(dmat.fb, 'Raw file', 'v')

# impute
library(VIM)
dmat.fa[,-c(1:2)] <- VIM::kNN(dmat.fa[,-c(1:2)], imp_var=F)
dmat.fb[,-c(1:2)] <- VIM::kNN(dmat.fb[,-c(1:2)], imp_var=F)

## --------

# collapse experiments
cormat.a <- zeros(length(prots.a), length(exps)*length(dcols))
cormat.b <- zeros(length(prots.b), length(exps)*length(dcols))
for(i in 1:length(prots.a)) {
  cormat.a[i,] <- as.vector(data.matrix((dmat.fa %>% filter(Protein==prots.a[i]))[,-c(1,2)]))
}
for(i in 1:length(prots.b)) {
  cormat.b[i,] <- as.vector(data.matrix((dmat.fb %>% filter(Protein==prots.b[i]))[,-c(1,2)]))
}

# correlate cell types against each other
cormat.a <- cor(cormat.a)
cormat.b <- cor(cormat.b)

pars.a <- svd(cormat.a)
pars.b <- svd(cormat.b)



# plot PCA cells ----------------------------------------------------------

pdf(file='manuscript/Figs/pca.pdf', width=7, height=2.25)

par(oma=c(0, 2.5, 2.25, 2))

layout(rbind(c(1, 2)))

par(mar=c(2.4, 3.75, 0.5, 0.75), cex.axis=0.85)

plot(pars.a$u[,1], pars.a$u[,2], pch=16, cex=0.5,
     col=paste0(rep(c('#FF0000', '#0000FF'), nrow(pars.a$u)/2), '44'),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', xlim=c(-0.045, 0.045), ylim=c(-0.15, 0.17))

text(-0.019, 0.11, 'Jurkat', col='red', font=2, cex=1, adj=c(0, 0.5))
text(0.025, -0.1, 'U-937', col='blue', font=2, cex=1, adj=c(1, 0.5))

axis(1, at=seq(-0.1, 0.1, by=0.02), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(seq(-0.04, 0.04, by=0.02), -0.175, seq(-0.04, 0.04, by=0.02), 
     xpd=T, srt=30, adj=c(1, 1), cex=0.85)
axis(2, at=seq(-0.2, 0.2, by=0.05), las=1, tck=-0.02, mgp=c(0, 0.5, 0))

mtext(paste0('PC1 - ', formatC((pars.a$d[1] / sum(pars.a$d)) * 100, digits=3), '%'), side=1,
      cex=0.85, line=1.4)
mtext(paste0('PC2 - ', formatC((pars.a$d[2] / sum(pars.a$d)) * 100, digits=2), '%'), side=2,
      cex=0.85, line=2.5)
mtext('Only proteins from Spectra', side=3, cex=0.85, line=0.9)
mtext(paste0(length(prots.a), ' proteins | ', length(peps.a), ' peptides'), 
      side=3, cex=0.85, line=0.05)

#mtext(paste0('Separation of cellular proteomes - ',length(pars.a$d),' samples'), 
#      side=3, outer=T, cex=1, font=2, line=1.4)

mtext(paste0('Separation of cellular proteomes'), 
      side=3, outer=T, cex=1, font=2, line=1.4)


#par(mar=c(1.5, 1, 1, 1))

plot(pars.a$u[,1], pars.b$u[,2], pch=16, cex=0.5,
     col=paste0(rep(c('#FF0000', '#0000FF'), nrow(pars.b$u)/2), '44'),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', xlim=c(-0.045, 0.045), ylim=c(-0.15, 0.17))

text(-0.02, 0.11, 'Jurkat', col='red', font=2, cex=1, adj=c(0, 0.5))
text(0.025, -0.1, 'U-937', col='blue', font=2, cex=1, adj=c(1, 0.5))

axis(1, at=seq(-0.1, 0.1, by=0.02), labels=NA, tck=-0.02, mgp=c(0, 0.2, 0))
text(seq(-0.04, 0.04, by=0.02), -0.175, seq(-0.04, 0.04, by=0.02), 
     xpd=T, srt=30, adj=c(1, 1), cex=0.85)
axis(2, at=seq(-0.2, 0.2, by=0.05), las=1, tck=-0.02, mgp=c(0, 0.5, 0))

mtext(paste0('PC1 - ', formatC((pars.b$d[1] / sum(pars.b$d)) * 100, digits=3), '%'), side=1,
      cex=0.85, line=1.4)
mtext(paste0('PC2 - ', formatC((pars.b$d[2] / sum(pars.b$d)) * 100, digits=2), '%'), side=2,
      cex=0.85, line=2.5)
mtext('Only proteins from DART-ID', side=3, cex=0.85, line=0.9)
mtext(paste0(length(prots.b), ' proteins | ', length(peps.b), ' peptides'), 
      side=3, cex=0.85, line=0.05)

dev.off()


# proteins with only 1 peptide? -------------------------------------------

a <- ev_b %>%
  filter(!Protein %in% ev_a$Protein) %>%
  group_by(Protein) %>%
  summarise(n=length(unique(`Modified sequence`))) %>%
  arrange(desc(n))

sum(a$n == 1)
