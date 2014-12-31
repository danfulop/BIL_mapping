library(stringr)
library(plyr)
library(ggplot2)
library(circlize)
library(scam)

load(file="/Users/Dani/UCD/BILs/pen12t.Rdata")

# Use shape constrained monotonically increasing general additive models to convert physical to genetic distance, using scam package

fits <- vector('list', 12) # list for storing the scam model fits
for (i in 1:length(levels(pen12t$chromosome))) {
  chrom <- levels(pen12t$chromosome)[i]
  fit.name <- paste0(chrom,"fit")
  scamfit <- scam(formula=PEN12_Pos ~ 0 + s(position_bp, bs="mpi", k=30), data=pen12t[pen12t$chromosome==chrom,]) # use low k so that small pericentromeric bins show up more
  fits[[i]] <- scamfit
}
save(fits, file="/Users/Dani/UCD/BILs/scam_fits.Rdata") # save scam fits' list object
load("/Users/Dani/UCD/BILs/scam_fits.Rdata")

load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

# predict genetic distance on bin.stats' physical distance coordinates

chrom <- levels(bin.stats$chr)
bin.predict <- vector('list', 12) # list of tables for storing the scam model genetic distance predictions
for (h in 1:length(chrom) ) {
  chrom.tab <- bin.stats[bin.stats$chr==chrom[h], ]
  # predict genetic distance 3 times with: bin.mid, bin.start, and bin.end
  bin.mid <- data.frame(bin.mid=chrom.tab$bin.mid)
  newdata <- data.frame(position_bp=bin.mid$bin.mid)
  bin.mid$gen.bin.mid <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
  bin.start <- data.frame(bin.start=chrom.tab$bin.start)
  newdata <- data.frame(position_bp=bin.start$bin.start)
  bin.start$gen.bin.start <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
  bin.end <- data.frame(bin.end=chrom.tab$bin.end)
  newdata <- data.frame(position_bp=bin.end$bin.end)
  bin.end$gen.bin.end <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
  bin.predict[[h]] <- cbind(bin.mid, bin.start, bin.end) # cbind 3 predictions
  zero.cM <- colwise(function(x, chr.start) x - (chr.start * chr.start/x ) ) # function to proportionally zero the predicted gen. dists.
  bin.predict[[h]][, c(2,4,6)] <- zero.cM(bin.predict[[h]][, c(2,4,6)], chr.start = bin.start$gen.bin.start[1]) # zero the cM predictions
}
bin.predict <- do.call(rbind, bin.predict) # rbind bin.predict's 12 elements
gen.bin.stats <- cbind(bin.stats[, 1:2], bin.predict) # cbind bin.stats' 1st 2 cols (i.e. bin & chr) & bin.predict
save(gen.bin.stats, file="/Users/Dani/UCD/BILs/leaf_traits/gen.bin.stats.Rdata")
load("/Users/Dani/UCD/BILs/leaf_traits/gen.bin.stats.Rdata")

# Additive plots
#-----------
setwd("/Users/Dani/UCD/BILs/leaf_traits/genDist_additive_qtl_plots")

# Function to generate genetic distance plots of one whole dataset, with 1 plot per trait
plot.map <- function(map.dat, bin.stats, dat.name) {
  for(i in 1:length(map.dat) ) {
    if(map.dat[[i]]$n.coef==0) {
      next
    } else {
      trait.name <- names(map.dat)[i]
      nz.coef <- map.dat[[i]]$non.zero.coefs
      nz.coef$chr <- as.factor(substr(nz.coef$chr,7,10))
      ch.lev <- levels(nz.coef$chr)
      n.ch.lev <- length(ch.lev)
      max.coef <- max(nz.coef$coefs)
      min.coef <- min(nz.coef$coefs)
      nz.coef$color[nz.coef$coefs < 0] <- "magenta"
      nz.coef$color[nz.coef$coefs > 0] <- "green"
      tmp <- lapply(1:nrow(nz.coef), function(i) { # create 4 new columns with genetic distance information for the intervals
        gen.int0.95.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(nz.coef$int0.95[i], ":")[[1]][1] ]
        gen.int0.95.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(nz.coef$int0.95[i], ":")[[1]][2] ]
        gen.int0.90.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(nz.coef$int0.90[i], ":")[[1]][1] ]
        gen.int0.90.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(nz.coef$int0.90[i], ":")[[1]][2] ]
        c(gen.int0.95.start, gen.int0.95.end, gen.int0.90.start, gen.int0.90.end)
      })
      gen.dist.cols <- do.call(rbind, tmp)
      colnames(gen.dist.cols) <- c('gen.int0.95.start', 'gen.int0.95.end', 'gen.int0.90.start', 'gen.int0.90.end')
      nz.coef <- cbind(nz.coef, gen.dist.cols)
      dat <- merge(gen.bin.stats, nz.coef, all.x=T)
      dat$y.lo <- -5*(max(abs(dat$coefs), na.rm=T)/20)
      dat$y.hi <- -1*(max(abs(dat$coefs), na.rm=T)/20)
      qtl.plot.90 <- ggplot(dat) + geom_rect(aes(xmin=gen.bin.start, xmax=gen.bin.end, ymin=y.lo, ymax=y.hi), color="black", fill="white", size=0.1) +
        geom_rect(aes(xmin=gen.int0.90.start, xmax=gen.int0.90.end, ymin=y.hi, ymax=abs(coefs), fill=color), alpha=0.25) +
        geom_rect(aes(xmin=gen.bin.start, xmax=gen.bin.end, ymin=y.hi, ymax=abs(coefs), fill=color, color=color), size=0.05) + facet_grid(chr ~ .) +
        theme_bw(12) + labs(y="Absolute magnitude of coefficients", x="Genotypic bin structure in genetic distance (cM)") +
        scale_color_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") +
        scale_fill_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") + theme(legend.position="none")
      ggsave(filename = paste("genDist", dat.name, trait.name, "90corr","pdf", sep="."), qtl.plot.90, width=8, height=10.5)
      qtl.plot.95 <- ggplot(dat) + geom_rect(aes(xmin=gen.bin.start, xmax=gen.bin.end, ymin=y.lo, ymax=y.hi), color="black", fill="white", size=0.1) +
        geom_rect(aes(xmin=gen.int0.95.start, xmax=gen.int0.95.end, ymin=y.hi, ymax=abs(coefs), fill=color), alpha=0.25) +
        geom_rect(aes(xmin=gen.bin.start, xmax=gen.bin.end, ymin=y.hi, ymax=abs(coefs), fill=color, color=color), size=0.05) + facet_grid(chr ~ .) +
        theme_bw(12) + labs(y="Absolute magnitude of coefficients", x="Genotypic bin structure in genetic distance (cM)") +
        scale_color_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") +
        scale_fill_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") + theme(legend.position="none")
      ggsave(filename = paste("genDist", dat.name, trait.name, "95corr","pdf", sep="."), qtl.plot.95, width=8, height=10.5)
    }
  }
}

load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/comp.map.Rdata")
plot.map(comp.map, gen.bin.stats, dat.name="comp")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/circ.map.Rdata")
plot.map(circ.map, gen.bin.stats, dat.name="circ")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/sym.map.Rdata")
plot.map(sym.map, gen.bin.stats, dat.name="sym")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/asym.map.Rdata")
plot.map(asym.map, gen.bin.stats, dat.name="asym")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/FT.map.Rdata")
#FT.map <- list(FT=FT.map) # modify FT results so that they fit the structure of the plotting function, i.e. list with sublist(s)
plot.map(FT.map, gen.bin.stats, dat.name="FT")

#---------

# Function to generate physical distance plots of one whole dataset, with 1 plot per trait
map.dat=comp.map; dat.name="comp"; i=1
plot.map <- function(map.dat, bin.stats, dat.name) {
  for(i in 1:length(map.dat) ) {
    if(map.dat[[i]]$n.coef==0) {
      next
    } else {
      trait.name <- names(map.dat)[i]
      nz.coef <- map.dat[[i]]$non.zero.coefs
      nz.coef$chr <- as.factor(substr(nz.coef$chr,7,10))
      ch.lev <- levels(nz.coef$chr)
      n.ch.lev <- length(ch.lev)
      max.coef <- max(nz.coef$coefs)
      min.coef <- min(nz.coef$coefs)
      nz.coef$color <- NULL
      nz.coef$color[nz.coef$coefs < 0] <- "magenta"
      nz.coef$color[nz.coef$coefs > 0] <- "green"
      dat <- merge(bin.stats, nz.coef, all.x=T)
      dat$bin.mid <- as.numeric(as.character(dat$bin.mid)) / 1000000
      dat$bin.start <- as.numeric(as.character(dat$bin.start)) / 1000000
      dat$bin.end <- as.numeric(as.character(dat$bin.end)) / 1000000
      dat$int0.95.start <- as.numeric(as.character(dat$int0.95.start)) / 1000000
      dat$int0.95.end <- as.numeric(as.character(dat$int0.95.end)) / 1000000
      dat$int0.90.start <- as.numeric(as.character(dat$int0.90.start)) / 1000000
      dat$int0.90.end <- as.numeric(as.character(dat$int0.90.end)) / 1000000
      dat$y.lo <- -5*(max(abs(dat$coefs), na.rm=T)/20)
      dat$y.hi <- -1*(max(abs(dat$coefs), na.rm=T)/20)
      # Plot by physical distance
      qtl.plot.90 <- ggplot(dat) + geom_rect(aes(xmin=bin.start, xmax=bin.end, ymin=y.lo, ymax=y.hi), color="black", fill="white", size=0.1) +
        geom_rect(aes(xmin=int0.90.start, xmax=int0.90.end, ymin=y.hi, ymax=abs(coefs), fill=color), alpha=0.25) +
        geom_rect(aes(xmin=bin.start, xmax=bin.end, ymin=y.hi, ymax=abs(coefs), fill=color)) + facet_grid(chr ~ .) +
        theme_bw(12) + labs(y="Absolute magnitude of coefficients", x="Genotypic bin structure in physical distance (Mbp)") +
        scale_fill_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") + theme(legend.position="none")
      ggsave(filename = paste(dat.name, trait.name, "90corr","pdf", sep="."), qtl.plot.90, width=8, height=10.5)
      qtl.plot.95 <- ggplot(dat) + geom_rect(aes(xmin=bin.start, xmax=bin.end, ymin=y.lo, ymax=y.hi), color="black", fill="white", size=0.1) +
        geom_rect(aes(xmin=int0.95.start, xmax=int0.95.end, ymin=y.hi, ymax=abs(coefs), fill=color), alpha=0.25) +
        geom_rect(aes(xmin=bin.start, xmax=bin.end, ymin=y.hi, ymax=abs(coefs), fill=color)) + facet_grid(chr ~ .) +
        theme_bw(12) + labs(y="Absolute magnitude of coefficients", x="Genotypic bin structure in physical distance (Mbp)") +
        scale_fill_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") + theme(legend.position="none")
      ggsave(filename = paste(dat.name, trait.name, "95corr","pdf", sep="."), qtl.plot.95, width=8, height=10.5)
    }
  }
}

load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/comp.map.Rdata")
plot.map(comp.map, bin.stats, dat.name="comp")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/circ.map.Rdata")
plot.map(circ.map, bin.stats, dat.name="circ")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/sym.map.Rdata")
plot.map(sym.map, bin.stats, dat.name="sym")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/asym.map.Rdata")
plot.map(asym.map, bin.stats, dat.name="asym")
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/FT.map.Rdata")
#FT.map <- list(FT=FT.map) # modify FT results so that they fit the structure of the plotting function, i.e. list with sublist(s)
plot.map(FT.map, bin.stats, dat.name="FT")

#---------

# Epistatic plots
#-----------
# Modify bin.stats to have doubled information
load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )
# head(bin.stats)
# summary(bin.stats)
bin.stats$bin2 <- bin.stats$bin
bin.stats$chr2 <- bin.stats$chr
bin.stats$bin.mid2 <- bin.stats$bin.mid
bin.stats$bin.start2 <- bin.stats$bin.start
bin.stats$bin.end2 <- bin.stats$bin.end
names(bin.stats)[1:5] <- paste0(names(bin.stats)[1:5], "1")
head(bin.stats)
summary(bin.stats)
doubled.bin.stats <- bin.stats

load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) 
bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )
# create a data.frame with chromosome start and end info. to initialize the circos plot
circ.init <- lapply(1:nlevels(bin.stats$chr), function(i) {
  chr <- levels(bin.stats$chr)[i]
  chr.stats <- bin.stats[bin.stats$chr==chr, ]
  start <- chr.stats$bin.start[1]
  end <- chr.stats$bin.end[nrow(chr.stats)]
  c(chr, start, end)
} )
circ.init <- data.frame(do.call(rbind, circ.init), stringsAsFactors = FALSE)
colnames(circ.init) <- c("chr", "start", "end")
circ.init$chr <- as.factor(circ.init$chr)
circ.init$start <- as.numeric(circ.init$start)
circ.init$end <- as.numeric(circ.init$end)
#summary(circ.init)


load("/Users/Dani/UCD/BILs/final_epistatic_sparsenet_results/comp.epi.map.Rdata")
map.dat=comp.epi.map; dat.name="comp"; i=1
map.dat=asym.epi.map; dat.name="asym"; i=1

# Function to plot epistatic results
# Plot all QTL plot in a 2D plane as a heatmap. Plot by bin-number. Make additional tick marks to delimit chromosomes
# Make 2 versions for each trait, 0.95 and 0.90 correlation plots, to compare them
#---------
# TO DO
# 1) reduce space between tracks => DONE
# 2) reduce cell padding => DONE
# 3) add y-axis scale lines in white
# 4) implement circos.links for epistatic QTLs
# 5) adjust track size with track height in circos.genomicTrackPlotRegion() => DONE
# 6) increase size of the chromosome distance scale text and tick marks
# 7) change to genetic distance
plot.epi.map <- function(map.dat, doubled.bin.stats, dat.name, circ.init) {
  for(i in 1:length(map.dat) ) {
    if(map.dat[[i]]$n.coef==0) {
      next
    } else {
      trait.name <- names(map.dat)[i]
      nz.coef <- map.dat[[i]]$non.zero.coefs
      # separate into additive and epistatic QTL
      adt.nz.coef <- nz.coef[str_count(nz.coef$chr, "ch")==1, ]
      epi.nz.coef <- nz.coef[grep("//", nz.coef$chr), ]
      # for additive QTL duplicate the information
      adt.nz.coef[, c(3:5,8,9,11,12)] <- lapply(adt.nz.coef[, c(3:5,8,9,11,12)], function(f) as.numeric(as.character(str_trim(f))) ) # trim and coerce to numeric certain columns
      # for epistatic QTL deconvolute the data into separate chromosome and bin position columns
      epi.nz.coef$bin1 <- laply(epi.nz.coef$bin, function(x) str_split(x, "_x_")[[1]][1])
      epi.nz.coef$chr1 <- laply(epi.nz.coef$chr, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin.mid1 <- laply(epi.nz.coef$bin.mid, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin.start1 <- laply(epi.nz.coef$bin.start, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin.end1 <- laply(epi.nz.coef$bin.end, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$int0.951 <- laply(epi.nz.coef$int0.95, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$int0.95.start1 <- laply(epi.nz.coef$int0.95.start, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$int0.95.end1 <- laply(epi.nz.coef$int0.95.end, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$int0.901 <- laply(epi.nz.coef$int0.90, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$int0.90.start1 <- laply(epi.nz.coef$int0.90.start, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$int0.90.end1 <- laply(epi.nz.coef$int0.90.end, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin2 <- laply(epi.nz.coef$bin, function(x) str_split(x, "_x_")[[1]][2])
      epi.nz.coef$chr2 <- laply(epi.nz.coef$chr, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$bin.mid2 <- laply(epi.nz.coef$bin.mid, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$bin.start2 <- laply(epi.nz.coef$bin.start, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$bin.end2 <- laply(epi.nz.coef$bin.end, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$int0.952 <- laply(epi.nz.coef$int0.95, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$int0.95.start2 <- laply(epi.nz.coef$int0.95.start, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$int0.95.end2 <- laply(epi.nz.coef$int0.95.end, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$int0.902 <- laply(epi.nz.coef$int0.90, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$int0.90.start2 <- laply(epi.nz.coef$int0.90.start, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$int0.90.end2 <- laply(epi.nz.coef$int0.90.end, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef <- epi.nz.coef[c(13:17,6,18:ncol(epi.nz.coef))] # eliminate original columns except coefs, and place coefs in its proper order/position
      epi.nz.coef[, c(3:5,8,9,11,12,15:17,19,20,22,23)] <- lapply(epi.nz.coef[, c(3:5,8,9,11,12,15:17,19,20,22,23)], function(f) as.numeric(as.character(str_trim(f))) )
      epi.nz.coef$int951min <- laply(epi.nz.coef$int0.951, function(x) str_split(x, ":")[[1]][1])
      epi.nz.coef$int951max <- laply(epi.nz.coef$int0.951, function(x) str_split(x, ":")[[1]][2])
      epi.nz.coef$int901min <- laply(epi.nz.coef$int0.901, function(x) str_split(x, ":")[[1]][1])
      epi.nz.coef$int901max <- laply(epi.nz.coef$int0.901, function(x) str_split(x, ":")[[1]][2])
      epi.nz.coef$int952min <- laply(epi.nz.coef$int0.952, function(x) str_split(x, ":")[[1]][1])
      epi.nz.coef$int952max <- laply(epi.nz.coef$int0.952, function(x) str_split(x, ":")[[1]][2])
      epi.nz.coef$int902min <- laply(epi.nz.coef$int0.902, function(x) str_split(x, ":")[[1]][1])
      epi.nz.coef$int902max <- laply(epi.nz.coef$int0.902, function(x) str_split(x, ":")[[1]][2])
      # make 2 bed-like files for the links/chords for epistatic QTLs, ditto for each of the correlation intervals
      epi1.bed <- with(epi.nz.coef, data.frame(chr=chr1, start=bin.start1, end=bin.end1, coefs=coefs) )
      epi2.bed <- with(epi.nz.coef, data.frame(chr=chr2, start=bin.start2, end=bin.end2, coefs=coefs) )
      epi951.bed <- with(epi.nz.coef, data.frame(chr=chr1, start=int0.95.start1, end=int0.95.end1, coefs=coefs) )
      epi952.bed <- with(epi.nz.coef, data.frame(chr=chr2, start=int0.95.start2, end=int0.95.end2, coefs=coefs) )
      epi901.bed <- with(epi.nz.coef, data.frame(chr=chr1, start=int0.90.start1, end=int0.90.end1, coefs=coefs) )
      epi902.bed <- with(epi.nz.coef, data.frame(chr=chr2, start=int0.90.start2, end=int0.90.end2, coefs=coefs) )
      bin.bed <- with(bin.stats, data.frame(chr=chr, start=bin.start, end=bin.end) ) # BED-like data.frame specifying BIN structure
      adt.nz.coef$color[adt.nz.coef$coefs > 0] <- 1 # adding a color column may be "illegal"
      adt.nz.coef$color[adt.nz.coef$coefs < 0] <- -1
      adt.nz.coef$coefs <- abs(adt.nz.coef$coefs)
      adt.bed <- with(adt.nz.coef, data.frame(chr=chr, start=bin.start, end=bin.end, coefs=coefs, color=color) )
      adt.bed$transparency <- 0
      adt.int95.bed <- with(adt.nz.coef, data.frame(chr=chr, start=int0.95.start, end=int0.95.end, coefs=coefs, color=color) )
      adt.int95.bed$transparency <- 0.75
      adt.int90.bed <- with(adt.nz.coef, data.frame(chr=chr, start=int0.90.start, end=int0.90.end, coefs=coefs, color=color) )
      adt.int90.bed$transparency <- 0.75
      adt.dat.list95 <- list(adt.int95.bed, adt.bed)
      adt.dat.list90 <- list(adt.int90.bed, adt.bed)
      adt.ylim <- c(0, max(adt.nz.coef$coefs) )
      # start plotting
      par(mar = c(1, 1, 1, 1) ) #, lwd = 0.5, cex = 1) # investigate later whether all these params. are ideal
      circos.par("start.degree" = 90, "track.margin"=c(0.005,0)) # chromosome-01 is at top right, as opposed starting at 0 degrees i.e. below horizontal on the right
      circos.genomicInitialize(circ.init, sector.width=circ.init$end )
      # chromosome bounding box can be customized w.r.t. background and border, see circos.trackPlotRegion()      
      circos.genomicTrackPlotRegion(data = adt.dat.list90, track.height=0.28, cell.padding=c(0.01,0,0.01,0), ylim=adt.ylim, panel.fun = function(region, value, ...) {      
        i = getI(...) # assign the plotting iteration through the data list's elemetns to a counter
        if (i == 1) { # have 75% transparency if it's 1st plotting iteration, and 0 transp. otherwise
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("magenta", "black", "green"), transparency=0.75 )
        } else {
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("magenta", "black", "green"), transparency=0 )
        }
        circos.genomicRect(region, value, col=sign.col(value$color), ybottom=0, ytop.column=1, border=sign.col(value$color), ...)
      }, bg.col="gray87", bg.border=NA ) 
      circos.genomicTrackPlotRegion(data = bin.bed, ylim=c(0,1), track.height=0.07, cell.padding=c(0,0,0,0), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col="white", border="black", ybottom=0, ytop=1, lwd=0.5, ...)
      }, bg.border=NA )
      col.fun = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"), transparency=0.75 )
      col.vector <- col.fun(epi.nz.coef$coefs)
      circos.genomicLink(epi901.bed, epi902.bed, col=col.vector, border=col.vector, lwd=0.5)
      col.fun = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"), transparency=0 )
      col.vector <- col.fun(epi.nz.coef$coefs)
      circos.genomicLink(epi1.bed, epi2.bed, col=col.vector, border=col.vector, lwd=0.5)
      circos.clear()
    }
  }
}
#-----------
# Useful for adjusting circos plotting parameters
# circos.par("cell.padding") # [1] 0.02 0.00 0.02 0.00; bottom, left, top, right; 1st & 3rd are % radius, 2nd & 4th are degrees
# circos.par("track.height") # [1] 0.2 , i.e. 20% of radius
# circos.par("track.margin") # [1] 0.01 0.01 # bottom and top margins, left and right are controlled by gap.degree
# circos.par("gap.degree") # [1] 1
# circos.info()

# Useful for debugging inside trackPlotRegions
#         print(get.cell.meta.data("sector.index"))
#         print(get.cell.meta.data("xlim"))
#         print(get.cell.meta.data("ylim"))
#         print(region)
#         print(value)

























