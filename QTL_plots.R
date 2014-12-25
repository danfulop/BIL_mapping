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
      qtl.plot <- ggplot(dat) + geom_rect(aes(xmin=bin.start, xmax=bin.end, ymin=y.lo, ymax=y.hi), color="black", fill="white", size=0.05) +
        geom_rect(aes(xmin=bin.start, xmax=bin.end, ymin=y.hi, ymax=abs(coefs), fill=color)) + facet_grid(chr ~ .) +
        theme_bw(8) + labs(y="Absolute magnitude of coefficients", x="Genotypic bin structure in physical distance (Mbp)") +
        scale_fill_identity("Sign of coefficients", labels=c("positive", "negative"), guide="legend") + theme(legend.position="none")
      ggsave(filename = paste(dat.name, trait.name,"pdf", sep="."), qtl.plot, width=8, height=10.5)
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
head(bin.stats)
summary(bin.stats)

map.dat=comp.epi.map; dat.name="comp"; i=1
map.dat=asym.epi.map; dat.name="asym"; i=1

# Function to plot epistatic results
# Plot all QTL plot in a 2D plane as a heatmap. Plot by bin-number. Make additional tick marks to delimit chromosomes
# Make 2 versions for each trait, 0.95 and 0.90 correlation plots, to compare them
#---------
plot.epi.map <- function(map.dat, bin.stats, dat.name) {
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
      adt.nz.coef[, c(3:5,8,9,11,12)] <- lapply(adt.nz.coef[, c(3:5,8,9,11,12)], function(f) as.numeric(as.character(str_trim(f))) )
      adt.nz.coef <- cbind(adt.nz.coef, adt.nz.coef[-6]) # exclude coefs column from 2nd copy
      names(adt.nz.coef)[c(1:5,7:12)] <- paste0(names(adt.nz.coef)[c(1:5,7:12)], "1")
      names(adt.nz.coef)[13:ncol(adt.nz.coef)] <- paste0(names(adt.nz.coef)[13:ncol(adt.nz.coef)], "2")
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
      # I may need to split the int0.95 and int0.90 labels into first and last bins
      # What I actually need is 1 dataset for each set of data to be plotted, as many as 3: selected bin, 0.95 interval, and 0.90 interval
      # ...and where all the pertinent bins have the correct coefficient value
      # at least that's the way to make it work with geom_raster()
      # Perhaps geom_tile() is more flexible...
      # I could use geom_tile() with a free scale ...and then use the coordinates instead of bin-numbers
      # OR, I could used geom_rect instead w/ free scale and w/ coordinates or bins
      # ** I may need a matrix of points!! **
      # Use BIN-number, and not full BIN label for axes
      nz.coef <- rbind(adt.nz.coef, epi.nz.coef)
      head(nz.coef)
      nz.coef$int951min <- laply(nz.coef$int0.951, function(x) str_split(x, ":")[[1]][1])
      nz.coef$int951max <- laply(nz.coef$int0.951, function(x) str_split(x, ":")[[1]][2])
      nz.coef$int901min <- laply(nz.coef$int0.901, function(x) str_split(x, ":")[[1]][1])
      nz.coef$int901max <- laply(nz.coef$int0.901, function(x) str_split(x, ":")[[1]][2])
      nz.coef$int952min <- laply(nz.coef$int0.952, function(x) str_split(x, ":")[[1]][1])
      nz.coef$int952max <- laply(nz.coef$int0.952, function(x) str_split(x, ":")[[1]][2])
      nz.coef$int902min <- laply(nz.coef$int0.902, function(x) str_split(x, ":")[[1]][1])
      nz.coef$int902max <- laply(nz.coef$int0.902, function(x) str_split(x, ":")[[1]][2])
      head(nz.coef)
      # save nz.coef, just b/c it's a convenient view of the data
      
      plot.dat <- merge(nz.coef, bin.stats) #, all=TRUE) # merge nz.coef with doubled bin.stats
      plot.dat$coefs[is.na(plot.dat$coefs)] <- 0 # make all NA coefficients == 0
      bin.stats.p <- bin.stats
      bin.stats.p$coefs <- NA
      head(plot.dat, 20)
      plot.dat
      epi.plot <- ggplot(bin.stats, aes(x=bin1, y=bin2)) + geom_raster()
      epi.plot <- epi.plot + geom_rect(data=nz.coef, aes(xmin=int901min, xmax=int901max, ymin=int902min, ymax=int902max, fill=coefs), alpha=0.4)
      epi.plot <- epi.plot + geom_rect(data=nz.coef, aes(xmin=int951min, xmax=int951max, ymin=int952min, ymax=int952max, fill=coefs))
      epi.plot <- epi.plot + geom_rect(data=nz.coef, aes(xmin=int902min, xmax=int902max, ymin=int901min, ymax=int901max, fill=coefs), alpha=0.4)
      epi.plot <- epi.plot + geom_rect(data=nz.coef, aes(xmin=int952min, xmax=int952max, ymin=int951min, ymax=int951max, fill=coefs))
      epi.plot <- epi.plot + scale_fill_gradient2(low="magenta", mid="black", high="green")
      epi.plot
      
      str_split(nz.coef[1, 'int0.951'], ":")[[1]][1]
      str_split(nz.coef[1, 'int0.951'], ":")[[1]][2]
      # 1) lay down the "structure" w/ geom_raster or geom_tile
      # 2) then use geom_rect to draw 0.90 and 0.95 intervals in 2D
    }
  }
}
#-----------

load("/Users/Dani/UCD/BILs/final_epistatic_sparsenet_results/comp.epi.map.Rdata")
























