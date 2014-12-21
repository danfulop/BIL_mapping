# use bin.stats and epi.bin.stats to setup the framework of the plots

# Additive plots
#-----------
load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names

setwd("/Users/Dani/UCD/BILs/leaf_traits/additive_qtl_plots")

# Function to generate plots of one whole dataset
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
      for(j in 1:n.ch.lev) {
        chr.stats <- bin.stats[bin.stats$chr==ch.lev[j],]
        chr.coef <- nz.coef[nz.coef$chr==ch.lev[j],]
        chr.dat <- merge(chr.stats, chr.coef, all.x=T)
        chr.dat$bin.mid <- as.numeric(as.character(chr.dat$bin.mid))
        chr.dat$bin.start <- as.numeric(as.character(chr.dat$bin.start))
        chr.dat$bin.end <- as.numeric(as.character(chr.dat$bin.end))
        chr.dat$y.lo <- -(max(abs(chr.dat$coefs), na.rm=T)/20)
        chr.dat$y.hi <- 0
        # Plot by physical distance
        qtl.plot <- ggplot(chr.dat) + geom_segment(aes(y=y.lo, yend=y.hi, x=bin.start, xend=bin.start), size=0.1) +
          geom_segment(aes(y=y.lo, yend=y.hi, x=bin.end, xend=bin.end), size=0.1) + 
          geom_segment(aes(y=y.lo, yend=y.lo, x=min(bin.start), xend=max(bin.end) ) , size=0.1) +
          geom_segment(aes(y=y.hi, yend=y.hi, x=min(bin.start), xend=max(bin.end) ) , size=0.1) +
          geom_point(aes(x=bin.mid, y=abs(coefs)), shape=17, size=2, color="red") + theme_bw(16) +
          labs(y="Absolute magnitude of coefficients", x="Physical distance (bp)")
        ggsave(filename = paste(dat.name, trait.name, ch.lev[j],"pdf", sep="."), qtl.plot, width=15, height=7.5)
      }
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

#-----------