library(stringr)
library(plyr)
library(ggplot2)

# use bin.stats and epi.bin.stats to setup the framework of the plots

# Additive plots
#-----------
load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names

setwd("/Users/Dani/UCD/BILs/leaf_traits/additive_qtl_plots")

# DEPRECATED Function to generate plots of one whole dataset, with 1 plot per chromosome
#------
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
      for(j in 1:n.ch.lev) {
        chr.stats <- bin.stats[bin.stats$chr==ch.lev[j],]
        chr.coef <- nz.coef[nz.coef$chr==ch.lev[j],]
        chr.dat <- merge(chr.stats, chr.coef, all.x=T)
        chr.dat$bin.mid <- as.numeric(as.character(chr.dat$bin.mid))
        chr.dat$bin.start <- as.numeric(as.character(chr.dat$bin.start))
        chr.dat$bin.end <- as.numeric(as.character(chr.dat$bin.end))
        chr.dat$y.lo <- -1.5*(max(abs(chr.dat$coefs), na.rm=T)/20)
        chr.dat$y.hi <- -0.5*(max(abs(chr.dat$coefs), na.rm=T)/20)
        # Plot by physical distance
        qtl.plot <- ggplot(chr.dat) + geom_segment(aes(y=y.lo, yend=y.hi, x=bin.start, xend=bin.start), size=0.1) +
          geom_segment(aes(y=y.lo, yend=y.hi, x=bin.end, xend=bin.end), size=0.1) + 
          geom_segment(aes(y=y.lo, yend=y.lo, x=min(bin.start), xend=max(bin.end) ) , size=0.1) +
          geom_segment(aes(y=y.hi, yend=y.hi, x=min(bin.start), xend=max(bin.end) ) , size=0.1) +
          geom_point(aes(x=bin.mid, y=abs(coefs), color=color), size=3) + 
          theme_bw(16) + labs(y="Absolute magnitude of coefficients", x="Physical distance (bp)")
        if ( all(chr.dat$coefs > 0, na.rm=TRUE) ) {
          qtl.plot <- qtl.plot + scale_color_identity("Sign of coefficients", labels=c("positive"), guide="legend")
        } else if ( all(chr.dat$coefs < 0, na.rm=TRUE) ) {
          qtl.plot <- qtl.plot + scale_color_identity("Sign of coefficients", labels=c("negative"), guide="legend")
        } else {
          qtl.plot <- qtl.plot + scale_color_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend")
        }
        ggsave(filename = paste(dat.name, trait.name, ch.lev[j],"pdf", sep="."), qtl.plot, width=15, height=7.5)
      }
    }
  }
}
#------
# **ADD bin-number plotting** to this and also epistatic plotting below?
# Function to generate plots of one whole dataset, with 1 plot per trait
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
      dat$y.lo <- -3*(max(abs(dat$coefs), na.rm=T)/20)
      dat$y.hi <- -1*(max(abs(dat$coefs), na.rm=T)/20)
      # Plot by physical distance
      qtl.plot <- ggplot(dat) + geom_segment(aes(y=y.lo, yend=y.hi, x=bin.start, xend=bin.start), size=0.05, alpha=0.7) +
        geom_segment(aes(y=y.lo, yend=y.hi, x=bin.end, xend=bin.end), size=0.05, alpha=0.7) + 
        geom_segment(aes(y=y.lo, yend=y.lo, x=bin.start, xend=bin.end ), size=0.05) +
        geom_segment(aes(y=y.hi, yend=y.hi, x=bin.start, xend=bin.end ), size=0.05) +
        geom_point(aes(x=bin.mid, y=abs(coefs), color=color), size=3) + facet_wrap(~ chr, ncol=2) +
        theme_grey(16) + labs(y="Absolute magnitude of coefficients", x="Genotypic bin structure in physical distance (Mbp)") +
        scale_color_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend")
      ggsave(filename = paste(dat.name, trait.name,"pdf", sep="."), qtl.plot, width=30, height=15)
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
FT.map <- list(FT=FT.map) # modify FT results so that they fit the structure of the plotting function, i.e. list with sublist(s)
plot.map(FT.map, bin.stats, dat.name="FT")

#---------

# Epistatic plots
#-----------
# Modify bin.stats to have doubled information
head(bin.stats)
summary(bin.stats)
bin.stats$bin.mid <- as.numeric(as.character(bin.stats$bin.mid))
bin.stats$bin.start <- as.numeric(as.character(bin.stats$bin.start))
bin.stats$bin.end <- as.numeric(as.character(bin.stats$bin.end))
bin.stats$bin2 <- bin.stats$bin
bin.stats$chr2 <- bin.stats$chr
bin.stats$bin.mid2 <- bin.stats$bin.mid
bin.stats$bin.start2 <- bin.stats$bin.start
bin.stats$bin.end2 <- bin.stats$bin.end
names(bin.stats)[1:5] <- paste0(names(bin.stats)[1:5], "1")

# Function to plot epistatic results
plot.epi.map <- function(map.dat, bin.stats, dat.name) {
  for(i in 1:length(map.dat) ) {
    if(map.dat[[i]]$n.coef==0) {
      next
    } else {
      trait.name <- names(map.dat)[i]
      nz.coef <- map.dat[[i]]$non.zero.coefs
      # separate into additive and epistatic QTL. ID add. by SL2.40, and epi. by //
      adt.nz.coef <- nz.coef[grep("SL2.40", nz.coef$chr), ]
      epi.nz.coef <- nz.coef[grep("//", nz.coef$chr), ]
      # for additive QTL duplicate the information
      adt.nz.coef$chr <- as.factor(substr(adt.nz.coef$chr,7,10))
      
      # for epistatic QTL deconvolute the data into separate chromosome and bin position columns
      epi.nz.coef$bin1 <- laply(epi.nz.coef$bin, function(x) str_split(x, "_x_")[[1]][1])
      epi.nz.coef$bin2 <- laply(epi.nz.coef$bin, function(x) str_split(x, "_x_")[[1]][2])
      epi.nz.coef$chr1 <- laply(epi.nz.coef$chr, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$chr2 <- laply(epi.nz.coef$chr, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$bin.mid1 <- laply(epi.nz.coef$bin.mid, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin.mid2 <- laply(epi.nz.coef$bin.mid, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$bin.start1 <- laply(epi.nz.coef$bin.start, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin.start2 <- laply(epi.nz.coef$bin.start, function(x) str_split(x, "//")[[1]][2])
      epi.nz.coef$bin.end1 <- laply(epi.nz.coef$bin.end, function(x) str_split(x, "//")[[1]][1])
      epi.nz.coef$bin.end2 <- laply(epi.nz.coef$bin.end, function(x) str_split(x, "//")[[1]][2])
      # I can plot by either physical distance or bin-number; bin-number prob. works better for a number of reasons
      # Make additional tick marks to delimit chromosomes
      
      # plot all QTL plot in a 2D plane as a heatmap
    }
  }
}

#-----------
load("/Users/Dani/UCD/BILs/leaf_traits/final_epistatic_sparsenet_results/comp.epi.map.Rdata")
























