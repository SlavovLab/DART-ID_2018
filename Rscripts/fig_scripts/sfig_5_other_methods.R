# init --------------------------------------------------------------------

library(tidyverse)
library(pracma)
library(RColorBrewer)
library(reshape2)
library(ggridges)
source('Rscripts/lib.R')

## load MaxQuant output ---------------------------------------------------------------------------

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ev_updated.txt')

# only take confident observations (from MaxQuant PEP, not DART-ID PEP)
# also remove SQC9*, which were run w/ a trapping column. The extra 10 or so minutes screws up
# some of the other prediction/alignment algos, so to be fair we are excluding those runs.
ev.f <- ev %>%
  mutate(Protein=sapply(strsplit(`Leading razor protein`, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })) %>%
  filter(!is.na(pep_new)) %>%
  filter(!grepl('SQC9', `Raw file`)) %>%
  dplyr::select(c('Modified sequence', 'Sequence', 'Leading razor protein', 'Protein', 'Raw file', 
                  'Retention time', 'Calibrated retention time',
                  'PEP', 'pep_new', 'pep_updated', 'residual', 'id', 'MS/MS scan number'),
                starts_with('Reporter intensity corrected'))
  #filter(PEP < 0.01)

## load RT estimation outputs  ------------------------------------------

## load ssrcalc output

hi_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrFA323_313_list1.txt')
hi_dat <- rbind(hi_dat, read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ssrFA7794_648_list2.txt'))

# append predicted HI to DART dataframe
ev.f$HI <- hi_dat$`HI (pred)`[match(ev.f$Sequence, hi_dat$Sequence)]

ssrcalc_coefs <- lm(ev.f$`Retention time` ~ ev.f$HI)$coefficients

ev.f$ssrcalc_error <- ev.f$`Retention time` - ((ev.f$HI * ssrcalc_coefs[2]) + ssrcalc_coefs[1])
ev.f$ssrcalc_RT_corrected <- ((ev.f$HI * ssrcalc_coefs[2]) + ssrcalc_coefs[1])


# load BioLCCC output

biolccc_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/biolccc_output.txt')
# i know acetylated and oxidized sequences will be thrown out, but don't need them for a quick comparison of the method
ev.f$biolccc_rt <- biolccc_dat$`Retention time, min`[match(ev.f$Sequence, biolccc_dat$Sequence)]

biolccc_coefs <- lm(ev.f$`Retention time` ~ ev.f$biolccc_rt)$coefficients

ev.f$biolccc_error <- ev.f$`Retention time` - ((ev.f$biolccc_rt * biolccc_coefs[2]) + biolccc_coefs[1])
ev.f$biolccc_RT_corrected <- ((ev.f$biolccc_rt * biolccc_coefs[2]) + biolccc_coefs[1])  

## load ELUDE output

elude_dat <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180613_5/elude_out.txt', skip=2)
elude_dat <- elude_dat %>% 
  arrange(Peptide) %>% group_by(Peptide) %>%
  summarise(RT=unique(Predicted_RT))

# match up peptide sequences between MaxQuant output and ELUDE output
ev.f$elude_RT <- elude_dat$RT[match(ev.f$Sequence, elude_dat$Peptide)]

# remove NA observations, and run linear regression between observed RTs and predicted RTs,
# to correct for any global systematic biases in the LC configuration (such as added dwell time, etc.)
elude_coefs <- lm(ev.f$`Retention time` ~ ev.f$elude_RT)$coefficients

# apply linear regression coefs and put back into ev.f data frame
ev.f$elude_error <- ev.f$`Retention time` - ((ev.f$elude_RT * elude_coefs[2]) + elude_coefs[1])
ev.f$elude_RT_corrected <- (ev.f$elude_RT * elude_coefs[2]) + elude_coefs[1]

# MaxQuant Match-Between-Runs (MBR) -- Calibrated Retention Time

exps <- sort(unique(ev.f$`Raw file`))
ev.mq <- data.frame()
for(exp in exps) {
  ev.a <- ev.f %>% 
    #filter(PEP < 0.01) %>% 
    filter(`Raw file`==exp) %>%
    filter(!is.na(`Calibrated retention time`))
  if(nrow(ev.a) == 0) return(c(0))
  coefs <- lm(ev.a$`Retention time` ~ ev.a$`Calibrated retention time`)$coefficients
  ev.a$RT_calibrated <- ((ev.a$`Calibrated retention time` * coefs[2]) + coefs[1])
  ev.mq <- rbind(ev.mq, ev.a)
}
rm(ev.a, coefs, exps)

# re-sort MBR dataframe by PSM id
ev.mq <- ev.mq %>% arrange(id)

ev.f$MBR_error <- ev.mq$`Retention time`-ev.mq$RT_calibrated
ev.f$MBR_RT_corrected <- ev.mq$RT_calibrated

# iRT

df_irt <- read_tsv('/gd/bayesian_RT/dat/irt.txt')

df_irt_f <- df_irt %>% 
  # remove IDs without predicted RTs
  filter(!is.na(PP.DeltaRT)) %>%
  # remove IDs that aren't assigned in MaxQuant
  filter(PEP.StrippedSequence %in% ev.f$Sequence) %>%
  # rename columns
  transmute(Sequence=PEP.StrippedSequence,
            iRT_error=PP.DeltaRT,
            `Raw file`=R.FileName,
            Proteins=PG.UniprotIds,
            iRT_RT=PP.EmpiricalRT,
            ScanNo=PSM.MS2ScanNumber,
            iRT_RT_corrected=PP.RTPredicted) %>%
  # remove ".raw" from the raw file name
  mutate(`Raw file`=substr(`Raw file`, 1, nchar(`Raw file`)-4)) %>%
  # create raw file-scan number-sequence identifier
  # don't want to compare different scans and definitely don't want
  # to compare scans ID'd to different sequences
  mutate(id=paste0(`Raw file`, '_', ScanNo, '_', Sequence)) %>%
  # remove duplicate identifiers
  distinct(id, .keep_all=T)

# create IDs from MaxQuant output
mq_id <- paste0(ev.f$`Raw file`, '_', ev.f$`MS/MS scan number`, '_', ev.f$Sequence)

ev.f <- cbind(ev.f, 
              df_irt_f[match(mq_id, df_irt_f$id), c('iRT_error', 'iRT_RT_corrected', 'iRT_RT')])

## refine PEPs with inferred RTs

# null distribution: normal distribution with center at mean(RT) and standard deviation of sd(RT)
# rt_minus: probability of observing peptide RT at random. 
# i.e., density of RT evaluated on the null distribution
rt_minus <- dnorm(ev.f$`Retention time`, 
                  mean=mean(ev.f$`Retention time`), 
                  sd=sd(ev.f$`Retention time`))

# rt_plus: probability of observing the peptide RT if the sequence is assigned correctly
# in this case, probability of observing RT with respect to ELUDE's predicted RT
# each unique peptide get's its own rt_plus distribution, where the mean is the ELUDE corrected RT,
# and the standard deviation is the standard deviation of the absolute ELUDE alignment error
rt_plus_ssrcalc <- dnorm(ev.f$`Retention time`, 
                         mean=ev.f$ssrcalc_RT_corrected, 
                         sd=sd(abs(ev.f$ssrcalc_error), na.rm=T))
rt_plus_biolccc <- dnorm(ev.f$`Retention time`, 
                         mean=ev.f$biolccc_RT_corrected, 
                         sd=sd(abs(ev.f$biolccc_error), na.rm=T))
rt_plus_elude <- dnorm(ev.f$`Retention time`, 
                       mean=ev.f$elude_RT_corrected, 
                       sd=sd(abs(ev.f$elude_error), na.rm=T))
rt_plus_irt <- dnorm(ev.f$iRT_RT, 
                       mean=ev.f$iRT_RT_corrected, 
                       sd=sd(abs(ev.f$iRT_error), na.rm=T))
rt_plus_MBR <- dnorm(ev.f$`Retention time`, 
                     mean=ev.f$MBR_RT_corrected, 
                     sd=sd(abs(ev.f$MBR_error), na.rm=T))

# MaxQuant PEP estimates are strange and artefactual sometimes. Make sure none exceed 1 or are at 0
pep <- ev.f$PEP
pep[pep > 1] <- 1
pep[pep == 0] <- .Machine$double.xmin

#                                         P(RT|delta=0)*P(delta=0)
# PEP.new = P(delta=0|RT) =   ___________________________________________________
#                             P(RT|delta=0)*P(delta=0) + P(RT|delta=1)*P(delta=1)

ssrcalc_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_ssrcalc * (1 - pep)))
ev.f$ssrcalc_pep <- ssrcalc_pep

biolccc_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_biolccc * (1 - pep)))
ev.f$biolccc_pep <- biolccc_pep

elude_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_elude * (1 - pep)))
ev.f$elude_pep <- elude_pep

irt_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_irt * (1 - pep)))
ev.f$irt_pep <- irt_pep

MBR_pep = (rt_minus * pep) / ((rt_minus * pep) + (rt_plus_MBR * (1 - pep)))
ev.f$MBR_pep <- MBR_pep


# save to file ------------------------------------------------------------

save(ev.f, file='/gd/bayesian_RT/dat/othermethods_20190422.rds')

## Define 2D PEP Scatter function, similar to Fig3A -----------------------------------------------

pep_scatter_plot <- function(pep, new_pep, filename, name, hi.color='black', title) {
  pep_log <- log10(pep)
  new_pep_log <- log10(new_pep)
  
  conf_limit <- 1e-8
  nbins <- 80
  x.bin <- seq(log10(conf_limit), 0, length=nbins)
  y.bin <- seq(log10(conf_limit), 0, length=nbins)
  
  freq <- as.data.frame(table(findInterval(pep_log, x.bin),
                              findInterval(new_pep_log, y.bin)))
  # remove any rows/columns over the 'nbins' limit
  freq <- freq %>%
    filter(as.numeric(Var1) <= nbins) %>%
    filter(as.numeric(Var2) <= nbins)
  
  freq[,1] <- as.numeric(freq[,1])
  freq[,2] <- as.numeric(freq[,2])
  
  freq2D <- diag(nbins)*0
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
  
  colfunc <- colorRampPalette(c('white', hi.color))
  
  pdf(file=paste0('manuscript/Figs/',filename), width=2.33, height=2.5)
  
  # layout(t(c(1, 2)), widths=c(7, 1))
  
  par(mar=c(2.25,2.25,2,0.5),
      pty='s', las=1,
      cex.axis=0.75, cex.lab=0.85, cex.main=1)
  
  cols <- colfunc(20)
  
  # Normal
  image(x.bin, y.bin, freq2D[-1,-1], col=cols,
        xlab=NA, ylab=NA,
        xaxs='i', yaxs='i',
        xaxt='n', yaxt='n', useRaster=F)
  
  abline(a=0, b=1, col='black')
  abline(h=-2, col='black', lty=2, lwd=1)
  abline(v=-2, col='black', lty=2, lwd=1)
  #segments(x0=-2, x1=-2, y0=log10(conf_limit), y1=-2, col='black', lty=2)
  #segments(x0=-2, x1=0, y0=-2, y1=-2, col='black', lty=2)
  
  rect(xleft=-2, xright=0, ybottom=log10(conf_limit), ytop=-2,
       border=NA, col=rgb(1,0,0,0.05))
  rect(xleft=log10(conf_limit), xright=-2, ybottom=-2, ytop=0,
       border=NA, col=rgb(0,0,1,0.05))
  
  text(-7.5, -1, 'Downgraded', cex=0.85, adj=c(0, 0.5))
  text(-1.2, -4.4, 'Upgraded', cex=0.85, adj=c(0, 0), srt=270)
  
  rng <- seq(-10, 0, 2)
  axis(1, tck=-0.02,  
       at=rng, labels=fancy_scientific(10^rng),
       mgp=c(0, 0.1, 0))
  axis(2, tck=-0.01, 
       at=rng, labels=fancy_scientific(10^rng),
       mgp=c(0, 0.2, 0), las=1)
  
  mtext('Spectra (MaxQuant)', 1, line=1, cex=0.85)
  mtext(name, 2, line=1.5, cex=0.85, las=3)
  mtext(paste0('Error Probability (PEP)\n', title), 3, line=0.06, cex=0.85, font=2)
  
  # par(mar=c(2.5, 0, 3, 1), pty='m')
  # image(matrix(seq(-1, 1, length.out=nbins), ncol=nbins), col=colfunc(nbins),
  #       xlab=NA, ylab=NA, xaxt='n', yaxt='n')
  # mtext('Density', side=3, line=0.25, cex=0.85)
  
  dev.off()
}

# Generate plots

# same colors as Fig 2
cols <- c(brewer.pal(9, 'Blues')[5], brewer.pal(9, 'BuGn')[5], brewer.pal(9, 'Blues')[8],
          brewer.pal(9, 'Reds')[4], brewer.pal(9, 'PuRd')[5], brewer.pal(9, 'Reds')[7])

pep_scatter_plot(pep, ev.f$ssrcalc_pep, 'pep_scatter_ssrcalc_v1.pdf', 'SSRCalc', cols[1], 'SSRCalc')
pep_scatter_plot(pep, ev.f$biolccc_pep, 'pep_scatter_biolccc_v1.pdf', 'BioLCCC', cols[2], 'BioLCCC')
pep_scatter_plot(pep, ev.f$elude_pep, 'pep_scatter_elude_v1.pdf', 'ELUDE', cols[3], 'ELUDE')
pep_scatter_plot(pep, ev.f$irt_pep, 'pep_scatter_irt_v1.pdf', 'iRT', cols[4], 'iRT')
pep_scatter_plot(pep, ev.f$MBR_pep, 'pep_scatter_MBR_v1.pdf', 'MaxQuant Calibrated RT', cols[5], 'MaxQuant')
pep_scatter_plot(pep, ev.f$pep_updated, 'pep_scatter_dart_v1.pdf', 'DART-ID', cols[6], 'DART-ID')

## Compare performance across these 4 methods -----------------------------------------------------

# first calculate q-values (FDR) from PEPs
pep_to_qval <- function(p) {
  (cumsum(p[order(p)]) / seq(1, length(p)))[order(order(p))]
}

ev.f$ssrcalc_qval <- pep_to_qval(ev.f$ssrcalc_pep)
ev.f$biolccc_qval <- pep_to_qval(ev.f$biolccc_pep)
ev.f$elude_qval <- pep_to_qval(ev.f$elude_pep)
ev.f$irt_qval <- pep_to_qval(ev.f$irt_pep)
ev.f$MBR_qval <- pep_to_qval(ev.f$MBR_pep)
ev.f$qval_updated <- pep_to_qval(ev.f$pep_updated)
ev.f$qval <- pep_to_qval(pep)

# fold-change of IDs as a function of the confidence threshold
x <- logseq(5e-4, 1, 100)

# frame to hold the results
df <- data.frame()
method.names <- c('Spectra', 'SSRCalc', 'BioLCCC', 'ELUDE', 'iRT', 'MaxQuant MBR', 'DART-ID')
counter <- 1
for(i in x) {
  cat('\r', counter, '/', length(x), '       ')
  flush.console()
  counter <- counter + 1
  
  ratios <- c(
    1, # Spectra
    sum(ev.f$ssrcalc_qval < i, na.rm=T) /    sum(ev.f$qval < i & !is.na(ev.f$ssrcalc_qval)),
    sum(ev.f$biolccc_qval < i, na.rm=T) /    sum(ev.f$qval < i & !is.na(ev.f$biolccc_qval)),
    sum(ev.f$elude_qval   < i, na.rm=T) /    sum(ev.f$qval < i & !is.na(ev.f$elude_qval)),
    sum(ev.f$irt_qval     < i, na.rm=T) /    sum(ev.f$qval < i & !is.na(ev.f$irt_qval)),
    sum(ev.f$MBR_qval     < i)          /    sum(ev.f$qval < i & !is.na(ev.f$MBR_qval)),
    sum(ev.f$qval_updated < i)          /    sum(ev.f$qval < i)
  )
  ident <- c(
    sum(ev.f$qval < i)                  /      nrow(ev.f),
    sum(ev.f$ssrcalc_qval < i, na.rm=T) /      sum(!is.na(ev.f$ssrcalc_qval)),
    sum(ev.f$biolccc_qval < i, na.rm=T) /      sum(!is.na(ev.f$biolccc_qval)),
    sum(ev.f$elude_qval   < i, na.rm=T) /      sum(!is.na(ev.f$elude_qval)),
    sum(ev.f$irt_qval     < i, na.rm=T) /      sum(!is.na(ev.f$irt_qval)),
    sum(ev.f$MBR_qval     < i)          /      sum(!is.na(ev.f$MBR_qval)),
    sum(ev.f$qval_updated < i)          /      nrow(ev.f)
  )
  
  df <- rbind(df, data.frame(
    x=as.numeric(i),
    ratio=as.numeric(ratios),
    ident=as.numeric(ident),
    Method=as.character(method.names)
  ))
}
df$Method <- factor(df$Method, levels=method.names)


## FDR fold change --------------------------------------------------------------------------------

pdf(file='manuscript/Figs/inference_fdr_increase.pdf', width=5, height=4)

par(mar=c(2.2,3,1.75,1),
    las=1, cex.axis=0.85, cex.lab=1, cex.main=1)

plot(0, 0, type='n', xlim=c(-3, -1), ylim=c(-15, 180),
     xlab=NA, ylab=NA, xaxs='i', yaxs='i', xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2, lwd=1)

# same colors as Fig 2, grey for Spectra
cols <- c(rgb(0.3,0.3,0.3), brewer.pal(9, 'Blues')[5], brewer.pal(9, 'BuGn')[5], brewer.pal(9, 'Blues')[8],
          brewer.pal(9, 'Reds')[4], brewer.pal(9, 'PuRd')[5], brewer.pal(9, 'Reds')[7])

for(i in 1:length(method.names)) {
  df_a <- df %>% filter(Method == method.names[i])
  lines(log10(df_a$x), (df_a$ratio-1)*100, col=cols[i], lty=1, lwd=2)
}

rng <- seq(-3, 0, by=0.5)
axis(1, tck=-0.02, at=rng, mgp=c(0, 0.1, 0),
     labels=c('0.1%', NA, '1%', NA, '10%', NA, '100%'))
axis(2, tck=-0.02, at=c(-25,seq(0, 175, 25)), mgp=c(0, 0.3, 0))


legend('topright', c('Spectra', 'SSRCalc', 'BioLCCC', 'ELUDE', 'iRT', 'MaxQuant MBR', 'DART-ID'),
       lwd=4, lty=1, col=cols, seg.len=1.2, ncol=1,
       bty='n', cex=0.8, x.intersp=0.6, y.intersp=1.2, inset=c(0.02, 0.0))

mtext('FDR Threshold', 1, line=1.2, cex=1)
mtext('% Increase', 2, line=1.7, cex=1, las=3)
mtext('Increase in PSMs', 3, line=0.2, cex=1, font=2)

dev.off()


# FDR Fold Change, DART + MBR only ----------------------------------------

pdf(file='manuscript/Figs/inference_fdr_increase_dart_MBR.pdf', width=2.33, height=2.5)

par(mar=c(2,2.5,1.5,1),
    las=1, cex.axis=0.85, cex.lab=1, cex.main=1)

plot(0, 0, type='n', xlim=c(-3, -1), ylim=c(-15, 180),
     xlab=NA, ylab=NA, xaxs='i', yaxs='i', xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2, lwd=1)

# same colors as Fig 2, grey for Spectra
cols <- c(rgb(0.3,0.3,0.3), brewer.pal(9, 'PuRd')[5], brewer.pal(9, 'Reds')[7])

method.names <- c('Spectra', 'MaxQuant MBR', 'DART-ID')

for(i in 1:length(method.names)) {
  df_a <- df %>% filter(Method == method.names[i])
  lines(log10(df_a$x), (df_a$ratio-1)*100, col=cols[i], lty=1, lwd=2)
}

rng <- seq(-3, 0, by=0.5)
axis(1, tck=-0.02, at=rng, mgp=c(0, 0.1, 0),
     labels=c('0.1%', NA, '1%', NA, '10%', NA, '100%'))
axis(2, tck=-0.02, at=c(-25,seq(0, 175, 25)), mgp=c(0, 0.3, 0))


legend('topright', c('Spectra', 'MaxQuant MBR', 'DART-ID'),
       lwd=3, lty=1, col=cols, seg.len=1, ncol=1,
       bg='white', bty='o', cex=0.7, x.intersp=0.6, y.intersp=1, inset=c(-0.0, -0.0))

mtext('FDR Threshold', 1, line=1, cex=0.85)
mtext('% Increase', 2, line=1.6, cex=0.85, las=3)
mtext('Increase in PSMs', 3, line=0.2, cex=0.85, font=2)

dev.off()

# DART vs. MQ updated PEPs ------------------------------------------------

dart_pep_log <- log10(ev.f$pep_updated)
mbr_pep_log <- log10(ev.f$MBR_pep)

conf_limit <- 1e-8
nbins <- 80
x.bin <- seq(log10(conf_limit), 0, length=nbins)
y.bin <- seq(log10(conf_limit), 0, length=nbins)

freq <- as.data.frame(table(findInterval(mbr_pep_log, x.bin),
                            findInterval(dart_pep_log, y.bin)))
# remove any rows/columns over the 'nbins' limit
freq <- freq %>%
  filter(as.numeric(Var1) <= nbins) %>%
  filter(as.numeric(Var2) <= nbins)

freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D <- diag(nbins)*0
freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

colfunc <- colorRampPalette(c('white', 'red'))

pdf(file=paste0('manuscript/Figs/pep_scatter_dart_v_mbr.pdf'), width=2.33, height=2.5)

par(mar=c(2.25,2.2,1.4,0.75),
    pty='s', las=1, cex.axis=0.7, cex.lab=0.7, cex.main=1)

cols <- colfunc(20)

# Normal
image(x.bin, y.bin, freq2D[-1,-1], col=cols, xlab=NA, ylab=NA,
      xaxs='i', yaxs='i', xaxt='n', yaxt='n', useRaster=F)

abline(a=0, b=1, col='black')
abline(h=-2, col='black', lty=2, lwd=1)
abline(v=-2, col='black', lty=2, lwd=1)

rng <- seq(-10, 0, 2)
axis(1, tck=-0.02, mgp=c(0, 0.1, 0),
     at=rng, labels=fancy_scientific(10^rng))
axis(2, tck=-0.01, mgp=c(0, 0.2, 0), las=1,
     at=rng, labels=fancy_scientific(10^rng))

mtext('MaxQuant Calibrated RT', 1, line=1, cex=0.85)
mtext('DART-ID', 2, line=1.3, cex=0.85, las=3)
mtext(paste0('Error Probability (PEP)'), 3, line=0.1, cex=0.85, font=2)

dev.off()


# Compare CVs b/n DART and MBR upgraded PSMs ------------------------------

ev.f <- ev.f %>%
  filter(!is.na(Protein)) %>%
  filter(!grepl("REV__|CON__", Protein)) %>%
  # only take SQC experiments
  filter(grepl("SQC", `Raw file`)) %>%
  # remove SQC9* experiments
  filter(!grepl('SQC9', `Raw file`))

# remove carrier, empty channels
data.cols <- grep('Reporter intensity corrected', colnames(ev.f))
ev.f <- ev.f[,-data.cols[c(1,2,3,4)]] 
  
# Normalize RI data within each experiment

ev_norm <- data.frame()
exps <- sort(unique(ev.f$`Raw file`))

for(i in 1:length(exps)) {
  e <- exps[i]
  print(e)
  
  # get subset of ev for this experiment
  .ev <- ev.f %>% filter(`Raw file` == e)
  
  # column indicecs of RI intensities
  ri_cols <- grep('Reporter intensity corrected', colnames(.ev))
  
  # normalize data
  .ev <- normalize_ri_data_table(.ev, ri_cols)
  
  # attach to output DF
  ev_norm <- rbind(ev_norm, .ev)
}
rm(.ev)


pep_thresh <- 1e-2
prot.psm.thresh <- 5

prots_orig <- ev_norm %>% 
  filter(PEP < pep_thresh) %>% 
  group_by(Protein) %>%
  summarise(n=n(), seqs=length(unique(`Modified sequence`))) %>%
  filter(n > prot.psm.thresh & seqs > 1) %>%
  arrange(desc(n)) %>%
  pull(Protein)

prots_new <- ev_norm %>%
  filter(PEP > pep_thresh & pep_new < pep_thresh) %>%
  group_by(Protein) %>%
  summarise(n=n(), seqs=length(unique(`Modified sequence`))) %>%
  filter(n > prot.psm.thresh & seqs > 1) %>%
  arrange(desc(n)) %>%
  pull(Protein)

# only take the intersection with original proteins
prots <- prots_new[prots_new %in% prots_orig]

# take high fold-change proteins (high signal-to-noise)
# run t-test for each protein
data.cols <- grep('Reporter intensity corrected', colnames(ev_norm))
prot_sigs <- sapply(1:length(prots), function(i) {
  e <- (ev_norm$Protein == prots[i])
  t.test(apply(ev_norm[e, data.cols[c(1, 3, 5)]], 1, mean), apply(ev_norm[e,data.cols[c(2, 4, 6)]], 1, mean),
         alternative='two.sided', var.equal=T)$p.value
})

# only take proteins with high fold-change
# i.e., significance of < 0.05 after applying bonferroni correction
prots <- prots[(prot_sigs*length(prots)) < 0.05]

prot_cvs_orig <- zeros(length(prots), length(data.cols))
prot_cvs_new  <- zeros(length(prots), length(data.cols))
prot_cvs_mbr  <- zeros(length(prots), length(data.cols))
prot_cvs_null <- zeros(length(prots), length(data.cols))

set.seed(1)

all_seqs <- unique(ev_norm %>% filter(PEP < pep_thresh) %>% pull(`Modified sequence`))

for(i in 1:length(prots)) {
  cat('\r', i, '/', length(prots), '-', prots[i], '                           ')
  flush.console()
  
  ev.a <- ev_norm %>% filter(Protein==prots[i])
  dmat <- data.matrix(ev.a %>% 
                        filter(PEP < pep_thresh) %>% 
                        group_by(`Modified sequence`) %>%
                        summarise_at(colnames(ev.a)[data.cols], mean, na.rm=T) %>%
                        select(-`Modified sequence`))
  prot_cvs_orig[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
  dmat <- data.matrix(ev.a %>% 
                        filter(PEP > pep_thresh & pep_updated < pep_thresh) %>% 
                        group_by(`Modified sequence`) %>%
                        summarise_at(colnames(ev.a)[data.cols], mean, na.rm=T) %>%
                        select(-`Modified sequence`))
  prot_cvs_new[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
  dmat <- data.matrix(ev.a %>% 
                        filter(PEP > pep_thresh & MBR_pep < pep_thresh) %>% 
                        group_by(`Modified sequence`) %>%
                        summarise_at(colnames(ev.a)[data.cols], mean, na.rm=T) %>%
                        select(-`Modified sequence`))
  prot_cvs_mbr[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
  
  # generate null distribution
  # get some random peptides
  random_seqs <- sample(all_seqs, size=10)
  dmat <- data.matrix(ev_norm %>% 
                        filter(`Modified sequence` %in% random_seqs) %>%
                        filter(PEP < pep_thresh) %>% 
                        group_by(`Modified sequence`) %>%
                        summarise_at(colnames(ev.a)[data.cols], mean, na.rm=T) %>%
                        select(-`Modified sequence`))
  prot_cvs_null[i,] <- apply(dmat, 2, sd) / apply(dmat, 2, mean)
}

cat("\nCollecing results...\n")
cvs_all <- data.frame()
cvs_all <- rbind(cvs_all, melt(prot_cvs_orig)  %>% mutate(Method="Spectra"))
cvs_all <- rbind(cvs_all, melt(prot_cvs_new)   %>% mutate(Method="DART-ID"))
cvs_all <- rbind(cvs_all, melt(prot_cvs_mbr)  %>% mutate(Method="MBR"))
cvs_all <- rbind(cvs_all, melt(prot_cvs_null)  %>% mutate(Method="Null"))

cvs_all$Method <- factor(cvs_all$Method, levels=c("Spectra", "DART-ID", "MBR", "Null"))

save(cvs_all, file='~/git/DART-ID_2018/dat/cvs_othermethods_20190117.rds')
# load data
load('dat/cvs_all_20180815.rds')

# plot CVs ----------------------------------------------------------------

methods <- as.character(unique(cvs_all$Method))
boxs <- list(DART=cvs_all$value[cvs_all$Method=='DART-ID'],
             MBR=cvs_all$value[cvs_all$Method=='MBR'],
             Spectra=cvs_all$value[cvs_all$Method=='Spectra'],
             Decoy=cvs_all$value[cvs_all$Method=='Null'])

p <- 
  ggplot(cvs_all) +
  geom_density_ridges(aes(x=value, y=rev(Method), group=rev(Method), fill=rev(Method)), 
                      rel_min_height=0.01, bandwidth=0.02) +
  scale_x_continuous(limits=c(0, 0.6)) +
  scale_y_discrete(limits=rev(levels(cvs_all$Method)), 
                   labels=rev(c('Decoy', 'Spectra', 'MBR', 'DART-ID')),
                   expand=c(0.01, 0)) +
  scale_fill_manual(values=c(cb[4], cb[1], cb[3], cb[2]), guide=F) +
  labs(y=NULL, x='CV of Relative Quantification  ',
       title='Consistency of Protein Quantification') +
  theme_ridges() + theme(
    plot.margin=margin(0.25, 0.3, 0.25, 0.1, 'cm'),
    axis.title.x=element_text(hjust=0.5, size=10),
    #axis.title.y=element_text(hjust=0.5, size=12),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    plot.title=element_text(size=10, hjust=0, vjust=1, lineheight=2, 
                            margin=margin(0,0,0.2,0,'cm'))
  )
ggsave('manuscript/Figs/cv_ridges_mbr.pdf', p, 'pdf', width=2.33, height=2.5, units='in')
