# Init --------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(pracma)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## run protein quant ---------------

# increase in protein quant

# find out the indices of columns with reporter ion data
# should be 10 columns, but the code doesn't rely on this number
data.cols <- grep('Reporter intensity corrected', colnames(ev))
# ignore empty, carrier channels
# data.cols <- data.cols[c(5,6)]

# filter for sqc master sets only, and only keep correctly IDed proteins
ev.f <- ev %>%
  filter(grepl('SQC', `Raw file`)) %>%
  filter(!grepl('SQC67[AB][16]|SQC67C1[3-9]|SQC67[CD]5|SQC68[DE]|IFN6[H-K]-Trg|SQC72D|SQC73[CD]|SQC74M|180416S_QC_SQC78A2',`Raw file`)) %>%
  filter(!grepl('REV__', `Leading razor protein`)) %>%
  filter(!grepl('CON__', Proteins)) %>%
  # filter(apply(ev[,data.cols]!=0, 1, sum) == length(data.cols)) %>%
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  }))

ev.f <- ev.f %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / seq(1, nrow(ev.f)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / seq(1, nrow(ev.f)))[order(order(pep_updated))])

# only take rows w/ complete quantitation, excluding empty channels (3 and 4)
ev.f <- ev.f[apply(ev.f[,data.cols[-c(3,4)]] == 0, 1, sum) == 0,]

experiments <- sort(unique(ev.f$`Raw file`))
peptides <- unique(ev.f$`Modified sequence`)
proteins <- unique(ev.f$Protein)

# create expression matrix, but fill in each cell with that observation's q-value
dmat_peps_spec <- matrix(nrow=length(peptides), ncol=length(experiments))
dmat_peps_dart <- matrix(nrow=length(peptides), ncol=length(experiments))
dmat_prots_spec <- matrix(nrow=length(proteins), ncol=length(experiments))
dmat_prots_dart <- matrix(nrow=length(proteins), ncol=length(experiments))

for(i in 1:length(experiments)) {
  cat('\r', i, '/', length(experiments), '         ')
  flush.console()
  
  ev.a <- ev.f %>% 
    filter(`Raw file`==experiments[i])
  
  ev.a.pep <- ev.a %>%
    group_by(`Modified sequence`) %>%
    summarise(qval=mean(qval),
              qval_updated=mean(qval_updated))
  
  m <- match(peptides, ev.a.pep$`Modified sequence`)
  dmat_peps_spec[, i] <- ev.a.pep$qval[m]
  dmat_peps_dart[, i] <- ev.a.pep$qval_updated[m]
  
  ev.a.prot <- ev.a %>%
    group_by(Protein) %>%
    summarise(qval=mean(qval),
              qval_updated=mean(qval_updated))
  
  m <- match(proteins, ev.a.prot$Protein)
  dmat_prots_spec[, i] <- ev.a.prot$qval[m]
  dmat_prots_dart[, i] <- ev.a.prot$qval_updated[m]
}

# remove peptides + proteins that are observed in less than 10/200 experiments
# so our matrices are less sparse.
remove_peptides <- apply(!is.na(dmat_peps_spec), 1, sum, na.rm=T) < 10
remove_prots <- apply(!is.na(dmat_prots_spec), 1, sum, na.rm=T) < 10

dmat_peps_spec <- dmat_peps_spec[!remove_peptides,]
dmat_peps_dart <- dmat_peps_dart[!remove_peptides,]
dmat_prots_spec <- dmat_prots_spec[!remove_prots,]
dmat_prots_dart <- dmat_prots_dart[!remove_prots,]


# Create Heatmaps ---------------------------------------------------------

# TRUE if under threshold, FALSE if over or NA
conf_threshold <- 1e-2

# convert to image raster
conv_to_raster <- function(dmat) {
  img <- dmat < conf_threshold
  # set all NAs to FALSE
  img[is.na(img)] <- F
  
  img[img == T] <- rgb(1,0,0)
  img[img == F] <- rgb(1,1,1)
  
  as.raster(img)
}

img_peps_spec <- conv_to_raster(dmat_peps_spec)
img_peps_dart <- conv_to_raster(dmat_peps_dart)
img_prots_spec <- conv_to_raster(dmat_prots_spec)
img_prots_dart <- conv_to_raster(dmat_prots_dart)


# Peptide Maps ------------------------------------------------------------

plot_img <- function(img) {
  plot(0, 0, type='n', xlim=c(0, ncol(img)), ylim=c(0, nrow(img)),
       xaxt='n', yaxt='n', xaxs='i', yaxs='i', xlab=NA, ylab=NA)
  rasterImage(img[nrow(img):1,], 0, 0, ncol(img), nrow(img), interpolate=F)
}

# pdf(file='manuscript/Figs/peptide_missing_dat_maps.pdf', width=3.5, height=3)
# par(oma=c(2,2,1.5,0))
# layout(cbind(1,2))
# 
# par(mar=c(0.5,1,0.25,0.5))
# plot_img(img_peps_spec)
# 
# text(par('usr')[2]*0.05, par('usr')[4]*0.97,
#      'Spectra', adj=c(0, 1), cex=1, col='black', xpd=T, font=2)
# 
# axis(1, tck=-0.02, at=seq(0, ncol(img_peps_spec), by=40), mgp=c(0, 0.2, 0))
# axis(2, tck=-0.02, at=seq(0, nrow(img_peps_spec), by=1000), las=1,
#      labels=seq(0, nrow(img_peps_spec) / 1000, 1), mgp=c(0, 0.3, 0))
# 
# par(mar=c(0.5,0.5,0.25,1))
# plot_img(img_peps_dart)
# 
# text(par('usr')[2]*0.05, par('usr')[4]*0.97,
#      'DART-ID', adj=c(0, 1), cex=1, col='black', xpd=T, font=2)
# 
# axis(1, tck=-0.02, at=seq(0, ncol(img_peps_spec), by=40), mgp=c(0, 0.2, 0))
# axis(2, tck=-0.02, at=seq(0, nrow(img_peps_spec), by=1000), labels=NA)
# 
# mtext('Peptides Quantified - FDR 1%', side=3, line=0.2, outer=T, cex=1, font=1)
# mtext('Experiment #', side=1, line=0.75, outer=T, cex=1)
# mtext('Distinct Peptides x 1000', side=2, line=0.5, outer=T, cex=1)
# 
# dev.off()

pdf(file='manuscript/Figs/protein_missing_dat_maps.pdf', width=3.5, height=4.5)
par(oma=c(2,2.5,1.5,0))
layout(cbind(1,2))

par(mar=c(0.5,1,0.25,0.5))
plot_img(img_prots_spec)

text(par('usr')[2]*0.05, par('usr')[4]*0.97,
     'Spectra', adj=c(0, 1), cex=1, col='black', xpd=T, font=2)

axis(1, tck=-0.02, at=seq(0, ncol(img_prots_spec), by=40), mgp=c(0, 0.2, 0))
axis(2, tck=-0.02, at=seq(0, nrow(img_prots_spec), by=300), las=1,
     labels=seq(0, nrow(img_prots_spec), by=300), mgp=c(0, 0.3, 0))

par(mar=c(0.5,0.5,0.25,1))
plot_img(img_prots_dart)

text(par('usr')[2]*0.05, par('usr')[4]*0.97,
     'DART-ID', adj=c(0, 1), cex=1, col='black', xpd=T, font=2)

axis(1, tck=-0.02, at=seq(0, ncol(img_prots_dart), by=40), mgp=c(0, 0.2, 0))
axis(2, tck=-0.02, at=seq(0, nrow(img_prots_dart), by=400), labels=NA)

mtext('Proteins Quantified - FDR 1%', side=3, line=0.2, outer=T, cex=1, font=1)
mtext('SCoPE-MS Run #', side=1, line=0.75, outer=T, cex=1)
mtext('Distinct Proteins', side=2, line=1.5, outer=T, cex=1)

dev.off()


# for SCP 2019 presentation -----------------------------------------------

a <- conv_to_raster(t(dmat_prots_spec))

plot_img <- function(img) {
  plot(0, 0, type='n', xlim=c(0, ncol(img)), ylim=c(0, nrow(img)),
       xaxt='n', yaxt='n', xaxs='i', yaxs='i', xlab=NA, ylab=NA)
  rasterImage(img[nrow(img):1,], xleft=0, ybottom=0, xright=ncol(img), ytop=nrow(img), interpolate=F)
}

png(file='manuscript/Figs/scp_2019_missing_dat.png', width=4.5, height=3, units='in', res=400)
par(oma=c(2,2.5,1.5,0))

par(mar=c(0.5,1,0.25,0.5))
plot_img(a)

# text(par('usr')[2]*0.05, par('usr')[4]*0.97,
#      'Spectra', adj=c(0, 1), cex=1, col='black', xpd=T, font=2)

axis(2, tck=-0.02, at=seq(0, ncol(img_prots_spec), by=40), las=1, mgp=c(0, 0.3, 0))
axis(1, tck=-0.02, at=seq(0, nrow(img_prots_spec), by=400),
     labels=seq(0, nrow(img_prots_spec), by=400), mgp=c(0, 0.2, 0))

mtext('Proteins Quantified - FDR 1%', side=3, line=0.2, outer=T, cex=1, font=1)
mtext('Distinct Proteins', side=1, line=0.75, outer=T, cex=1)
mtext('SCoPE-MS Run #', side=2, line=1, outer=T, cex=1)

dev.off()
