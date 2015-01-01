library(stringr)
library(plyr)
library(ggplot2)
library(circlize)
library(scam)

load(file="/Users/Dani/UCD/BILs/pen12t.Rdata")

# Use shape constrained monotonically increasing general additive models to convert physical to genetic distance, using scam package

# fits <- vector('list', 12) # list for storing the scam model fits
# for (i in 1:length(levels(pen12t$chromosome))) {
#   chrom <- levels(pen12t$chromosome)[i]
#   fit.name <- paste0(chrom,"fit")
#   scamfit <- scam(formula=PEN12_Pos ~ 0 + s(position_bp, bs="mpi", k=30), data=pen12t[pen12t$chromosome==chrom,]) # use low k so that small pericentromeric bins show up more
#   fits[[i]] <- scamfit
# }
# save(fits, file="/Users/Dani/UCD/BILs/scam_fits.Rdata") # save scam fits' list object
load("/Users/Dani/UCD/BILs/scam_fits.Rdata")

load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

# predict genetic distance on bin.stats' physical distance coordinates
# chrom <- levels(bin.stats$chr)
# bin.predict <- vector('list', 12) # list of tables for storing the scam model genetic distance predictions
# for (h in 1:length(chrom) ) {
#   chrom.tab <- bin.stats[bin.stats$chr==chrom[h], ]
#   # predict genetic distance 3 times with: bin.mid, bin.start, and bin.end
#   bin.mid <- data.frame(bin.mid=chrom.tab$bin.mid)
#   newdata <- data.frame(position_bp=bin.mid$bin.mid)
#   bin.mid$gen.bin.mid <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
#   bin.start <- data.frame(bin.start=chrom.tab$bin.start)
#   newdata <- data.frame(position_bp=bin.start$bin.start)
#   bin.start$gen.bin.start <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
#   bin.end <- data.frame(bin.end=chrom.tab$bin.end)
#   newdata <- data.frame(position_bp=bin.end$bin.end)
#   bin.end$gen.bin.end <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
#   bin.predict[[h]] <- cbind(bin.mid, bin.start, bin.end) # cbind 3 predictions
#   zero.cM <- colwise(function(x, chr.start) x - (chr.start * chr.start/x ) ) # function to proportionally zero the predicted gen. dists.
#   bin.predict[[h]][, c(2,4,6)] <- zero.cM(bin.predict[[h]][, c(2,4,6)], chr.start = bin.start$gen.bin.start[1]) # zero the cM predictions
#   if(bin.predict[[h]][1, 4] != 0) { # re-check zeroing of genDist b/c of numerical issues that sometimes result in a tiny negative or positive value !=0 as the start of the chromosome
#     bin.predict[[h]][, c(2,4,6)] <- bin.predict[[h]][, c(2,4,6)] - bin.predict[[h]][1, 4]
#   } else next
# }
# bin.predict <- do.call(rbind, bin.predict) # rbind bin.predict's 12 elements
# gen.bin.stats <- cbind(bin.stats[, 1:2], bin.predict) # cbind bin.stats' 1st 2 cols (i.e. bin & chr) & bin.predict
# save(gen.bin.stats, file="/Users/Dani/UCD/BILs/leaf_traits/gen.bin.stats.Rdata")
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
#-----------
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
load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
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

# Create gen.circ.init to initialize genetic distance plots
gen.circ.init <- lapply(1:nlevels(gen.bin.stats$chr), function(i) {
  chr <- levels(gen.bin.stats$chr)[i]
  chr.stats <- gen.bin.stats[gen.bin.stats$chr==chr, ]
  start <- chr.stats$gen.bin.start[1]
  end <- chr.stats$gen.bin.end[nrow(chr.stats)]
  c(chr, start, end)
} )
gen.circ.init <- data.frame(do.call(rbind, gen.circ.init), stringsAsFactors = FALSE)
colnames(gen.circ.init) <- c("chr", "start", "end")
gen.circ.init$chr <- as.factor(gen.circ.init$chr)
gen.circ.init$start <- as.numeric(gen.circ.init$start)
gen.circ.init$end <- as.numeric(gen.circ.init$end)

load("/Users/Dani/UCD/BILs/final_epistatic_sparsenet_results/comp.epi.map.Rdata")
map.dat=comp.epi.map; dat.name="comp"; i=1

# Function to plot epistatic results
#---------
# ...to install and reload the hacked version to label chrom-axes in cM if length < 1000
# install.packages("/Users/Dani/UCD/R/circlize_0.2.0.tar.gz", repos=NULL, type="source")
# detach("package:circlize", unload=TRUE)
# library(circlize)
plot.epi.map <- function(map.dat, bin.stats, gen.bin.stats, dat.name, circ.init, gen.circ.init) {
  for(i in 1:length(map.dat) ) {
    if(map.dat[[i]]$n.coef==0) {
      next
    } else {
      trait.name <- names(map.dat)[i]
      nz.coef <- map.dat[[i]]$non.zero.coefs
      # separate into additive and epistatic QTL
      adt.nz.coef <- nz.coef[str_count(nz.coef$chr, "ch")==1, ]
      epi.nz.coef <- nz.coef[grep("//", nz.coef$chr), ]
      # setup additive data
      adt.nz.coef[, c(3:5,8,9,11,12)] <- lapply(adt.nz.coef[, c(3:5,8,9,11,12)], function(f) as.numeric(as.character(str_trim(f))) ) # trim and coerce to numeric certain columns
      adt.nz.coef$color[adt.nz.coef$coefs > 0] <- 1
      adt.nz.coef$color[adt.nz.coef$coefs < 0] <- -1
      adt.nz.coef$coefs <- abs(adt.nz.coef$coefs)
      # setup additive QTL BED objects for phys. dist. plots
      adt.bed <- with(adt.nz.coef, data.frame(chr=chr, start=bin.start, end=bin.end, coefs=coefs, color=color) )
      adt.int95.bed <- with(adt.nz.coef, data.frame(chr=chr, start=int0.95.start, end=int0.95.end, coefs=coefs, color=color) )
      adt.int90.bed <- with(adt.nz.coef, data.frame(chr=chr, start=int0.90.start, end=int0.90.end, coefs=coefs, color=color) )
      adt.dat.list95 <- list(adt.int95.bed, adt.bed)
      adt.dat.list90 <- list(adt.int90.bed, adt.bed)
      adt.ylim <- c(0, max(adt.nz.coef$coefs) )
      # add gen. dist. info. and setup additive QTL BED objects for gen. dist. plots
      tmp <- lapply(1:nrow(adt.nz.coef), function(i) { # create 6 new columns with genetic distance information for the intervals
        gen.bin.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==as.character(adt.nz.coef$bin[i]) ]
        gen.bin.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==as.character(adt.nz.coef$bin[i]) ]
        gen.int0.95.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(adt.nz.coef$int0.95[i], ":")[[1]][1] ]
        gen.int0.95.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(adt.nz.coef$int0.95[i], ":")[[1]][2] ]
        gen.int0.90.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(adt.nz.coef$int0.90[i], ":")[[1]][1] ]
        gen.int0.90.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(adt.nz.coef$int0.90[i], ":")[[1]][2] ]
        c(gen.bin.start, gen.bin.end, gen.int0.95.start, gen.int0.95.end, gen.int0.90.start, gen.int0.90.end)
      })
      gen.dist.cols <- do.call(rbind, tmp)
      colnames(gen.dist.cols) <- c('gen.bin.start', 'gen.bin.end','gen.int0.95.start', 'gen.int0.95.end', 'gen.int0.90.start', 'gen.int0.90.end')
      adt.nz.coef <- cbind(adt.nz.coef, gen.dist.cols)
      # setup additive QTL BED objects for gen. dist. plots
      gen.adt.bed <- with(adt.nz.coef, data.frame(chr=chr, start=gen.bin.start, end=gen.bin.end, coefs=coefs, color=color) )
      gen.adt.int95.bed <- with(adt.nz.coef, data.frame(chr=chr, start=gen.int0.95.start, end=gen.int0.95.end, coefs=coefs, color=color) )
      gen.adt.int90.bed <- with(adt.nz.coef, data.frame(chr=chr, start=gen.int0.90.start, end=gen.int0.90.end, coefs=coefs, color=color) )
      gen.adt.dat.list95 <- list(gen.adt.int95.bed, gen.adt.bed)
      gen.adt.dat.list90 <- list(gen.adt.int90.bed, gen.adt.bed)
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
      # add gen. dist. info. for epistatic QTL
      tmp <- lapply(1:nrow(epi.nz.coef), function(i) { # create new columns with genetic distance information for the intervals
        gen.bin.start1 <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==as.character(epi.nz.coef$bin1[i]) ]
        gen.bin.end1 <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==as.character(epi.nz.coef$bin1[i]) ]
        gen.int0.95.start1 <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(epi.nz.coef$int0.951[i], ":")[[1]][1] ]
        gen.int0.95.end1 <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(epi.nz.coef$int0.951[i], ":")[[1]][2] ]
        gen.int0.90.start1 <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(epi.nz.coef$int0.901[i], ":")[[1]][1] ]
        gen.int0.90.end1 <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(epi.nz.coef$int0.901[i], ":")[[1]][2] ]
        gen.bin.start2 <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==as.character(epi.nz.coef$bin2[i]) ]
        gen.bin.end2 <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==as.character(epi.nz.coef$bin2[i]) ]
        gen.int0.95.start2 <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(epi.nz.coef$int0.952[i], ":")[[1]][1] ]
        gen.int0.95.end2 <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(epi.nz.coef$int0.952[i], ":")[[1]][2] ]
        gen.int0.90.start2 <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(epi.nz.coef$int0.902[i], ":")[[1]][1] ]
        gen.int0.90.end2 <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(epi.nz.coef$int0.902[i], ":")[[1]][2] ]
        c(gen.bin.start1, gen.bin.end1, gen.int0.95.start1, gen.int0.95.end1, gen.int0.90.start1, gen.int0.90.end1, gen.bin.start2,
          gen.bin.end2, gen.int0.95.start2, gen.int0.95.end2, gen.int0.90.start2, gen.int0.90.end2)
      })
      gen.dist.cols <- do.call(rbind, tmp)
      colnames(gen.dist.cols) <- c('gen.bin.start1', 'gen.bin.end1','gen.int0.95.start1', 'gen.int0.95.end1', 'gen.int0.90.start1', 'gen.int0.90.end1',
                                   'gen.bin.start2', 'gen.bin.end2','gen.int0.95.start2', 'gen.int0.95.end2', 'gen.int0.90.start2', 'gen.int0.90.end2')
      epi.nz.coef <- cbind(epi.nz.coef, gen.dist.cols)
      # make 2 bed-like files for the links/chords for epistatic QTLs, ditto for each of the correlation intervals
      gen.epi1.bed <- with(epi.nz.coef, data.frame(chr=chr1, start=gen.bin.start1, end=gen.bin.end1, coefs=coefs) )
      gen.epi2.bed <- with(epi.nz.coef, data.frame(chr=chr2, start=gen.bin.start2, end=gen.bin.end2, coefs=coefs) )
      gen.epi951.bed <- with(epi.nz.coef, data.frame(chr=chr1, start=gen.int0.95.start1, end=gen.int0.95.end1, coefs=coefs) )
      gen.epi952.bed <- with(epi.nz.coef, data.frame(chr=chr2, start=gen.int0.95.start2, end=gen.int0.95.end2, coefs=coefs) )
      gen.epi901.bed <- with(epi.nz.coef, data.frame(chr=chr1, start=gen.int0.90.start1, end=gen.int0.90.end1, coefs=coefs) )
      gen.epi902.bed <- with(epi.nz.coef, data.frame(chr=chr2, start=gen.int0.90.start2, end=gen.int0.90.end2, coefs=coefs) )
      # bin structure BED objects
      bin.bed <- with(bin.stats, data.frame(chr=chr, start=bin.start, end=bin.end) ) # BED-like data.frame specifying BIN structure
      gen.bin.bed <- with(gen.bin.stats, data.frame(chr=chr, start=gen.bin.start, end=gen.bin.end) ) # BED-like data.frame specifying gen.bin structure
      # set up plotting parameters to reuse, such as track.margin, etc.
      init.tkHt = 0.05
      gap.degree = 2
      bottom.margin = 0.005
      adt.tkHt = 0.20 # additive track height
      cpad = 0.01
      bin.tkHt = 0.07
      transp = 0.75
      transp2 = 0.3
      bkgd.col = "gray90" # gray95 works well too, but it's not very visible in my monitor
      init.cex = 1.2
      major.by=20
      # Physical distance plots
      # 0.90 correlation
      par(mar = c(1, 1, 1, 1), cex = init.cex ) #, lwd = 0.5, cex = 1) # investigate later whether all these params. are ideal
      circos.par("start.degree" = 90, "track.margin"=c(bottom.margin,0), "gap.degree"=gap.degree) # chromosome-01 is at top right, as opposed starting at 0 degrees i.e. below horizontal on the right
      circos.genomicInitialize(circ.init, sector.width=circ.init$end, track.height=init.tkHt )
      par(cex = 1 )
      circos.genomicTrackPlotRegion(data = adt.dat.list90, track.height=adt.tkHt, cell.padding=c(cpad,0,cpad,0), ylim=adt.ylim, panel.fun = function(region, value, ...) {      
        
        i = getI(...) # assign the plotting iteration through the data list's elements to a counter
        if (i == 1) { # have 75% transparency if it's 1st plotting iteration, and 0 transp. otherwise
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("magenta", "black", "green"), transparency=transp )
        } else {
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("magenta", "black", "green"), transparency=0 )
        }
        circos.genomicRect(region, value, col=sign.col(value$color), ybottom=0, ytop.column=1, border=sign.col(value$color), ...)
      }, bg.col=bkgd.col, bg.border=NA ) 
      circos.genomicTrackPlotRegion(data = bin.bed, ylim=c(0,1), track.height=bin.tkHt, cell.padding=c(0,0,0,0), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col="white", border="black", ybottom=0, ytop=1, lwd=0.5, ...)
      }, bg.border=NA )
      col.fun = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"), transparency=transp )
      col.vector <- col.fun(epi.nz.coef$coefs)
      circos.genomicLink(epi901.bed, epi902.bed, col=col.vector, border=col.vector, lwd=0.5)
      col.fun = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"), transparency=0 )
      col.vector <- col.fun(epi.nz.coef$coefs)
      circos.genomicLink(epi1.bed, epi2.bed, col=col.vector, border=col.vector, lwd=0.5)
      circos.clear()
      # Genetic distance plots
      # 0.90 correlation
      par(mar = c(1, 1, 1, 1), cex = init.cex ) #, lwd = 0.5, cex = 1) # investigate later whether all these params. are ideal
      circos.par("start.degree" = 90, "track.margin"=c(bottom.margin,0), "gap.degree"=gap.degree) # chromosome-01 is at top right, as opposed starting at 0 degrees i.e. below horizontal on the right
      circos.genomicInitialize(gen.circ.init, sector.width=gen.circ.init$end, track.height=init.tkHt, major.by=major.by )
      par(cex = 1 )
      circos.genomicTrackPlotRegion(data = gen.adt.dat.list90, track.height=adt.tkHt, cell.padding=c(cpad,0,cpad,0), ylim=adt.ylim, panel.fun = function(region, value, ...) {      
        i = getI(...) # assign the plotting iteration through the data list's elemetns to a counter
        if (i == 1) { # have 75% transparency if it's 1st plotting iteration, and 0 transp. otherwise
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("magenta", "black", "green"), transparency=transp )
        } else {
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("magenta", "black", "green"), transparency=0 )
        }
        circos.genomicRect(region, value, col=sign.col(value$color), ybottom=0, ytop.column=1, border=sign.col(value$color), ...)
      }, bg.col=bkgd.col, bg.border=NA ) 
      circos.genomicTrackPlotRegion(data = gen.bin.bed, ylim=c(0,1), track.height=bin.tkHt, cell.padding=c(0,0,0,0), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col="white", border="black", ybottom=0, ytop=1, lwd=0.5, ...)
      }, bg.border=NA )
      col.fun = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"), transparency=transp )
      col.vector <- col.fun(epi.nz.coef$coefs)
      circos.genomicLink(gen.epi901.bed, gen.epi902.bed, col=col.vector, border=col.vector, lwd=0.5)
      col.fun = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"), transparency=0 )
      col.vector <- col.fun(epi.nz.coef$coefs)
      circos.genomicLink(gen.epi1.bed, gen.epi2.bed, col=col.vector, border=col.vector, lwd=0.5)
      # add legends
      legend.col = colorRamp2(breaks=c(min(epi.nz.coef$coefs),0,max(epi.nz.coef$coefs)), colors=c("magenta", "black", "green"))
      num.col.vector <- c( seq(max(epi.nz.coef$coefs), 0, length.out=7), seq(0, min(epi.nz.coef$coefs), length.out=7)[2:7])
      legend.col.vector <- legend.col(num.col.vector)
      text.col.vector <- as.character(round(num.col.vector, 2) )
      text.col.vector[c(2,3,5,6,8,9,11,12)] <- ""
      legend(x=0.77, y=0.98, legend=text.col.vector,  fill=legend.col.vector, border=NA, bty="n", bg="white", y.intersp=0.5, x.intersp=0.2, cex=0.6)
      legend(x=0.755, y=0.98, legend=rep("", 13),  fill=legend.col.vector, border=NA, bty="n", bg="white", y.intersp=0.5, x.intersp=0.2, cex=0.6)
      legend(x=0.74, y=0.98, legend=rep("", 13),  fill=legend.col.vector, border=NA, bty="n", bg="white", y.intersp=0.5, x.intersp=0.2, cex=0.6)
      text(x=0.893, y=1, labels="Epistatic QTL", cex=0.7)
      text.col.vector <- c("green", "magenta")
      legend(x=-0.98, y=0.98, legend=c("", ""), fill=text.col.vector, border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=-0.98, y=0.962, legend=c("", ""), fill=text.col.vector, border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=-1, y=0.962, legend=c("", ""), fill=text.col.vector, border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=-1, y=0.98, legend=c("positive", "negative"), fill=text.col.vector, border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      text(x=-0.844, y=1, labels="Additive QTL", cex=0.7)
      circos.clear()
    }
  }
}
#-----------
# TO DO
# 1) reduce space between tracks => DONE
# 2) reduce cell padding => DONE
# 3) add y-axis scale lines in white for additive QTLs; prob. requires maing space in the top btw. chr. 01 & 12 to fit the text of the y-axis
# 4) implement circos.links for epistatic QTLs => DONE
# 5) adjust track size with track height in circos.genomicTrackPlotRegion() => DONE
# 6) increase size of the chromosome distance scale text and tick marks => DONE
# 7) change to genetic distance => DONE
# 8) change tick labels units to cM => DONE
# 9) check gen.dist borders => DONE, but I don't know why there are points outside the plotting regions
#10) add legends


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

























