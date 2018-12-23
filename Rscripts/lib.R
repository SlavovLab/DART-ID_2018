#source('validate.lib.R')

load.ev <- function(ev, par.file='dat/params.Fit2.RData', include.REV=FALSE, include.CON=FALSE) {
  load(par.file)
  # remove abnormal LC experiments
  # load experiments from correlation testing in similar.lc.R
  exps.lc <- unlist(read_csv('dat/exps.corr.txt')[,2])
  names(exps.lc) <- NULL
  
  ev <- ev %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>% # Only use Elite experiments
    filter(`Raw file` %in% exps.lc) # Remove abnormal LC experiments
  
  ## Filter of PEP < .05
  ev.f <- ev %>% filter(PEP < 0.05) %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) # Only use Elite experiments

  # Remove Reverse matches
  if(!include.REV) {
    ev.f <- ev.f %>% filter(!grepl('REV*', `Leading razor protein`)) 
  }  
  # Remove Contaminants
  if(!include.CON) {
    ev.f <- ev.f %>% filter(!grepl('CON*',`Leading razor protein`))
  }
  
  ev.f <- ev.f %>% 
    filter(PEP < 0.05) %>%
    mutate(Protein=`Leading razor protein`) %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP", "Protein") %>%
    mutate(exp_id=`Raw file`) %>%  # new column - exp_id = numeric version of experiment file
    mutate_at("exp_id", funs(as.numeric(as.factor(.)))) %>%
    mutate(`Stan ID`=`Peptide ID`) %>%
    mutate_at("Stan ID", funs(as.numeric(as.factor(.))))
  ev.f <<- ev.f
  
  experiment_factors <<- as.factor(ev.f$`Raw file`)
  num_exps <<- length(unique(ev.f[["exp_id"]]))
  
  ## "true peptide id" matches peptide id in evidence file
  ## "peptide id" is index in 1:num_peptides for stan
  raw_peptide_id <<- ev.f[["Peptide ID"]]
  pep.id.list <<- unique(raw_peptide_id)
  num_peptides <<- length(unique(ev.f$`Stan ID`))
  
  # parse linear regression params from STAN output
  beta0 <<- pars[sprintf('beta_0[%i]', seq(1, num_exps))]
  beta1 <<- pars[sprintf('beta_1[%i]', seq(1, num_exps))]
  beta2 <<- pars[sprintf('beta_2[%i]', seq(1, num_exps))]
  
  # sigma params from Fit2
  split.point <<- pars[sprintf('split_point[%i]', seq(1, num_exps))]
  sigma.slope <<- pars[sprintf('sigma_slope[%i]', seq(1, num_exps))]
  sigma.intercept <<- pars[sprintf('sigma_intercept[%i]', seq(1, num_exps))]
  sigma.slope.global <<- pars['sigma_slope_global']
  
  # sigma params from Fit1
  sigmas <<- pars[grep('sigma\\[', names(pars))]
  
  mus <<- pars[grep('mu\\[', names(pars))]
}

clean.file.name <- function(name) {
  name <- gsub('#', '', name)
  name <- gsub('.*_NC_', '', name)
  name <- gsub('.raw', '', name)
  name <- gsub('set', '', name)
  name <- gsub('[-|+|=]', '_', name)
  return(name)
}

theme_bert <- function() {
  theme_bw(base_size=10, base_family="Helvetica") %+replace% 
  #theme_bw(base_size=16, base_family="Helvetica") %+replace% 
    theme(
      #panel.background  = element_rect(fill=NULL, color='black', size=0.5),
      #panel.background=element_blank(),
      #panel.border=element_blank(),
      #axis.line.x=element_line(color='black', size=0.25),
      #axis.line.y=element_line(color='black', size=0.25),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      #axis.line=element_line(size=0.5, color='#888888'),
      axis.line=element_blank(),
      axis.ticks=element_line(color='black', size=0.25),
      #axis.ticks=element_blank(),
      panel.grid=element_blank()
    )
}

fold.change.comp <- function(exps, begin=1e-5, end=1, num.steps=100, log=T, add.dummy=T) {
  library(pracma)
  # only use the same raw files between all evidence files
  common.exps <- Reduce(intersect, lapply(exps, function(exp) { exp$`Raw file` }))
  exps <- lapply(exps, function(exp) { exp[exp$`Raw file` %in% common.exps,]})
  
  # equally spaced steps in log space
  if(log) {
    x <- logseq(begin, end, n=num.steps)
  } else {
    x <- seq(begin, end, length.out=num.steps)
  }
  
  # frame to hold the results
  df <- data.frame()
  counter <- 0
  for(i in x) {
    ratios <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.new < i, na.rm=TRUE) / 
         sum(exp$PEP < i & !is.na(exp$PEP.new)))
    }))
    
    # add dummy maxquant experiment
    if(add.dummy) {
      ratios <- c(ratios, 1)
      exp.names <- c(names(exps), 'MaxQuant')
    } else {
      exp.names <- names(exps)
    }
    
    df <- rbind(df, data.frame(
      x=as.numeric(i),
      PEP=as.numeric(ratios),
      Method=as.character(exp.names)
    ))
  }
  return(df)
}

# fancy scientific scales
# from: https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  #l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  #l <- gsub("e", "%*%10^", l)
  # make sure +0 just turns into 0
  l <- gsub("\\+00", "00", l)
  # return this as an expression
  return(parse(text=l))
}

# from: https://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
emptyPlot <- function() {
  ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
}

bHeatmap <- c('#710e0d', '#8b42df', '#4c94dc', '#40df91', '#fef245')
vc <- c("#004358", "#1F8A70", "#BEDB39", "#FFE11A", "#FD7400")
av <- c("#20B2CF", "#85DB86", "#F29F05", "#F25C05", "#D92525")
#cb <- c('#20B2CF', '#FF6666', '#F68930', '#FFFFFF')
cb <- c('#20B2CF', '#FF6666', '#888888', '#FFFFFF')

#' Process Experiment Description Excel Sheet
process.desc <- function() {
  desc <- read.csv('dat/SingleCellExperiments_Description.csv')
  # only take the portion we need
  desc <- desc[1:308,]
  # remove rows that don't have experiment info
  desc <- desc[grepl('^\\d{2,3}[A-Z]{1}', desc$`Exp...`),]
  desc <- desc[,c(1, 3:12)]
  names(desc)[1] <- 'Exp'
  rownames(desc) <- NULL
  desc <- melt(desc, id.vars=c('Exp'), variable.name='Channel', value.name='Sample',
               factorsAsStrings=F)
  desc$Exps <- as.character(desc$Exp)
  # channel ID
  desc$ch <- as.numeric(desc$Channel)
  
  desc$Sample <- as.character(desc$Sample)
  desc$Sample[desc$Sample %in% c('', '0', 'Empty', 'N/A', 'PBS')] <- NA
  
  ## Quantity - how many cells were in each channel
  desc$Quantity <- str_extract(desc$Sample, '\\d+(e\\d)?(\\.\\dk?)?')
  # convert scientific notation
  desc$Quantity[grep('\\de\\d', desc$Quantity)] <- 
    as.numeric(desc$Quantity[grep('\\de\\d', desc$Quantity)])
  # convert "k" notation
  desc$Quantity[grep('\\d\\.\\dk', desc$Quantity)] <- 
    as.numeric(str_extract(desc$Quantity[grep('\\d+\\.\\dk', desc$Quantity)], '\\d+\\.\\d')) * 1000
  
  ## Type - what kind of cell it was
  # J = Jurkat
  # H = Hek293
  # U = U937
  # ES|EB = Mouse
  desc$Type <- str_extract(desc$Sample, '[JHUjhu]{1}|([Ee]{1}[SsBb]{1})')
  # Mixed type - more than one type of cell
  # most likely a carrier channel
  desc$Type[grepl('\\&|and', desc$Sample)] <- 'Mixed'
  
  # Uppercase = intact, single cell
  # Lowercase = diluted, lysed mixture
  desc$Diluted <- F
  desc$Diluted[grepl('[jhu]|es|eb', desc$Type)] <- T
  desc$Diluted[is.na(desc$Type)] <- NA
  
  # Re-normalize cell type for easier searching
  desc$Type <- toupper(desc$Type)
  
  return(desc)
}

fold.change.comp <- function(exps, range=c(1e-3, 1)) {
  library(pracma)
  # only use the same raw files between all evidence files
  #common.exps <- Reduce(intersect, lapply(exps, function(exp) { exp$`Raw file` }))
  #exps <- lapply(exps, function(exp) { exp[exp$`Raw file` %in% common.exps,]})
  
  num.steps <- 100
  # equally spaced steps in log space
  #x <- logseq(1e-10, 1, n=num.steps)
  x <- logseq(range[1], range[2], n=num.steps)
  # frame to hold the results
  df <- data.frame()
  counter <- 0
  for(i in x) {
    ratios <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.updated < i) / 
         sum(exp$PEP < i))
    }))
    identified_spectra <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP < i)) / nrow(exp)
    }))
    identified_spectra_new <- unlist(lapply(exps, function(exp) {
      (sum(exp$PEP.updated < i)) / nrow(exp)
    }))
    df <- rbind(df, data.frame(
      x=as.numeric(i),
      PEP=as.numeric(ratios),
      ident=as.numeric(identified_spectra),
      ident_new=as.numeric(identified_spectra_new),
      Method=as.character(names(exps))
    ))
  }
  return(df)
}

remove.channels <- function(ev) {
  # load excel description
  desc <- process.desc()
  desc.types <- dcast(desc, Exp~Type, value.var='Type', fun.aggregate=length)
  desc.types$Exp <- as.character(desc.types$Exp)
  
  ## remove carrier and empty channels
  # these channel indices vary between experiments, 
  # so we're gonna have to loop thru these and check for each one
  exps <- str_extract(clean.file.name(unique(ev$`Raw file`)), '[1-9][0-9][A-Z]')
  exps <- unique(exps[!is.na(exps)])
  
  for(i in exps) {
    #if(sum(grepl(i,ev$exp)) <= 0) next
    inds.remove <- desc[desc$Exp==i & (desc$Type=='MIXED' | desc$Type=='NA' | is.na(desc$Type)),'ch']
    
    cat('Channel(s)',inds.remove,'from experiment',i,'are blanks or carrier channels. Removing...\n')
    
    # make relative to the ev data frame
    inds.remove <- inds.remove + data.cols[1] - 1
    # set to NA
    ev[,inds.remove] <- NA
  }
  
  return(ev)
}

get.exp.ids <- function(raw.files) {
  exps <- str_extract(clean.file.name(raw.files), '[1-9][0-9][A-Z]')
  exps <- unique(exps[!is.na(exps)])
  return(exps)
}

filter.exps <- function(ev, pep.thresh=1e-2) {
  # ev <- read_tsv('dat/evidence_elite.txt')
  
  #removed.exps <- c()
  
  # confidence threshold - harsh to increase confidence
  # and decrease computation time
  # pep.thresh <- 1e-2
  
  # make sure we have enough observations in all experiments
  # lets say... need at least 100 to get a good idea of experiment quality
  #below.thresh <- ev %>% 
  #  group_by(`Raw file`) %>% 
  #  summarise(n=sum(PEP < pep.thresh)) %>%
  #  arrange(desc(n)) %>%
  #  filter(n < 100) %>%
  #  pull(`Raw file`)
  
  #ev <- ev %>% filter(!`Raw file` %in% below.thresh)
  
  raw.files <- unique(ev$`Raw file`)
  cor.mat <- matrix(nrow=length(raw.files), ncol=length(raw.files))
  
  # correlate high-confidence PSMs between experiments
  # only need the upper triangular of the correlation matrix
  for(i in 1:length(raw.files)) {
    ev.a <- ev %>% 
      filter(PEP < pep.thresh) %>%
      filter(`Raw file`==raw.files[i]) %>%
      select(c('Mod. peptide ID', 'Retention time'))
    
    
    for(j in i:length(raw.files)) {
      if(i == j) next
      
      cat('\r', i, '-', j, '       ')  
      flush.console()
      
      ev.b <- ev %>%
        filter(PEP < pep.thresh) %>%
        filter(`Raw file`==raw.files[j]) %>%
        filter(`Mod. peptide ID` %in% ev.a$`Mod. peptide ID`) %>%
        select(c('Mod. peptide ID', 'Retention time')) %>%
        group_by(`Mod. peptide ID`) %>%
        summarise(`Retention time`=mean(`Retention time`))
      
      # require at least 5 common points to do this analysis
      if(nrow(ev.b) < 5) {
        cor.mat[i,j] <- NA
        next
      }
      
      b <- ev.b$`Retention time`
      a <- ev.a %>%
        filter(`Mod. peptide ID` %in% ev.b$`Mod. peptide ID`) %>%
        group_by(`Mod. peptide ID`) %>%
        summarise(`Retention time`=mean(`Retention time`)) %>%
        pull(`Retention time`)
      
      cor.mat[i,j] <- cor(a, b)
    }
  }
  
  return(cor.mat)
}

plot.cor.mat <- function(cor.mat, show.text=F) {
  library(ggplot2)
  library(reshape2)
  
  # only take upper triangle. don't want to destroy the matrix structure,
  # so set the lower triangle to a dummy value so it can be removed later
  cor.mat[lower.tri(cor.mat)] <- 100
  
  # set diagnoal to 1
  diag(cor.mat) <- 1
  
  df <- melt(cor.mat)
  df <- df %>%
    filter(!value == 100)
  
  p <- ggplot(df, aes(y=Var1, x=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value='grey50',
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    theme(axis.text.y = element_text(size=8),
          axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 8, hjust = 1)
          #axis.text.x = element_blank()
    ) +
    labs(x=NULL, y=NULL) +
    coord_fixed()
  
  if(show.text) {
    p <- p + geom_text(aes(label=format(value, digits=2)))
  }
  
  p
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

normalize_ri_data <- function(dmat) {
  dmat[dmat==0] <- NA
  
  ## normalize data
  
  # first normalize by column, by median
  dmat <- t(t(dmat) / apply(dmat, 2, median, na.rm=T))
  # now normalize across rows, by mean
  dmat <- dmat / apply(dmat, 1, mean, na.rm=T)
  # remove rows without quant
  dmat <- dmat[!apply(apply(dmat, 1, is.na), 2, sum) > 0,]
  dmat
}
normalize_ri_data_table <- function(ev.f, dcols, remove.empty.rows=T, impute.data) {
  ev.f[,dcols][ev.f[,dcols]==0] <- NA
  
  ## normalize data
  
  # first normalize by column, by median
  ev.f[,dcols] <- t(t(ev.f[,dcols]) / apply(ev.f[,dcols], 2, median, na.rm=T))
  # now normalize across rows, by mean
  ev.f[,dcols] <- ev.f[,dcols] / apply(ev.f[,dcols], 1, mean, na.rm=T)
  # remove rows without quant
  if(remove.empty.rows) {
    ev.f <- ev.f[!apply(apply(ev.f[,dcols], 1, is.na), 2, sum) > 0,]
  }
  #if(impute.data) {
  #  library(VIM)
  #  # remove rows with more than half NA
  #  ev.f <- ev.f[!apply(apply(ev.f[,dcols], 1, is.na), 2, sum) > floor(length(dcols) / 2),]
  #  VIM::kNN(ev.f[,dcols], imp_var=F)
  #}
  return(ev.f)
}

extract_uniprot_id <- function(leading_protein) {
  sapply(strsplit(leading_protein, "\\|"), function(p) {
    if(length(unlist(p)) == 1) return(p[1])
    else if(length(unlist(p)) == 3) return(p[2])
    else return(p[1])
  })
}

# stolen from seqinr: https://cran.r-project.org/web/packages/seqinr/index.html
read.fasta <- function(file = system.file("sequences/ct.fasta.gz", package = "seqinr"), 
                       seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE,
                       set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE,
                       bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
                       endian = .Platform$endian, apply.mask = TRUE) {
  
  seqtype <- match.arg(seqtype) # default is DNA
  
  ##############################
  #
  # Regular flat FASTA text file
  #
  ##############################
  if(!bfa){ # Regular text file
    #
    # Read the fasta file as a vector of strings:
    #
    lines <- readLines(file)
    #
    # Remove comment lines starting with a semicolon ';'
    #
    if(legacy.mode){
      comments <- grep("^;", lines)
      if(length(comments) > 0) lines <- lines[-comments]
    }
    #
    # Get the line numbers where sequences names are:
    #
    ind <- which(substr(lines, 1L, 1L) == ">")
    #
    # Compute the total number of sequences:
    #
    nseq <- length(ind)
    if(nseq == 0){
      stop("no line starting with a > character found")
    }
    #
    # Localize sequence data:
    #
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    #
    # Read sequences:
    #
    sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], collapse = ""))
    if(seqonly) return(sequences)
    #
    # Read sequence names:
    #
    nomseq <- lapply(seq_len(nseq), function(i){
      firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
      substr(firstword, 2, nchar(firstword))
    })
    #
    # Turn DNA sequences in lower case letters if required:
    #
    if(seqtype == "DNA"){
      if(forceDNAtolower){
        sequences <- as.list(tolower(sequences))
      }
    }
    #
    # Turn it into a vector of single chars if required:
    #
    # if(as.string == FALSE) sequences <- lapply(sequences, s2c)
    #
    # Set sequence attributes when required:
    #
    if(set.attributes){
      for(i in seq_len(nseq)){
        Annot <- lines[ind[i]]
        if(strip.desc) Annot <- substr(Annot, 2L, nchar(Annot))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot,
                                           class = switch(seqtype, "AA" = "SeqFastaAA", "DNA" = "SeqFastadna"))
      }
    }
    #
    # Give the sequences names to the list elements:
    #
    names(sequences) <- nomseq
    return(sequences)
  }
  ##############################
  #
  # MAQ binary FASTA file
  #
  ##############################
  if(bfa){
    if(seqtype != "DNA") stop("binary fasta file available for DNA sequences only")
    # Open file in binary mode:
    mycon <- file(file, open = "rb")
    r2s <- words(4) # byte to tetranucleotide
    
    readOneBFARecord <- function(con, sizeof.longlong, endian, apply.mask){
      len <- readBin(con, n = 1, what = "int", endian = endian)
      if(length(len) == 0) return(NULL) # end of file 
      name <- readBin(con, n = 1, what = "character", endian = endian)
      ori_len <- readBin(con, n = 1, what = "int", endian = endian)
      len <- readBin(con, n = 1, what = "int", endian = endian)
      seq <- readBin(con, n = len*sizeof.longlong, what = "raw", size = 1, endian = endian)
      mask <- readBin(con, n = len*sizeof.longlong, what = "raw", size = 1, endian = endian)
      if(endian == "little"){
        neword <- sizeof.longlong:1 + 
          rep(seq(0, (len - 1)*sizeof.longlong, by = sizeof.longlong), 
              each = sizeof.longlong)
        # something like 8 7 6 5 4 3 2 1 16 15 14 13 12 11 10 9 ...
        seq <- seq[neword]
        mask <- mask[neword]
      }
      seq4 <- c2s(r2s[as.integer(seq) + 1])
      seq4 <- substr(seq4, 1, ori_len)
      if(apply.mask){
        mask4 <- c2s(r2s[as.integer(mask) + 1])
        mask4 <- substr(mask4, 1, ori_len)
        npos <- gregexpr("a", mask4, fixed = TRUE)[[1]]
        for(i in npos) substr(seq4, i, i + 1) <- "n"
      }
      return(list(seq = seq4, name = name))
    } # end readOneBFARecord
    
    sequences <- vector(mode = "list")
    nomseq <- vector(mode = "list")
    i <- 1
    repeat{
      res <- readOneBFARecord(mycon, sizeof.longlong, endian, apply.mask)
      if(is.null(res)) break
      sequences[[i]] <- res$seq
      nomseq[[i]] <- res$name
      i <- i + 1
    }
    close(mycon)
    nseq <- length(sequences)
    if(seqonly) return(sequences)
    # if(as.string == FALSE) sequences <- lapply(sequences, s2c)
    if(set.attributes){
      for(i in seq_len(nseq)){
        if(!strip.desc) Annot <- c2s(c(">", nomseq[[i]]))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = "SeqFastadna")
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
}




