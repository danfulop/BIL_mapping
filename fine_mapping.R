# Script to fine map IL-QTL using BIL-QTL

library(stringr)
library(IRanges)
library(foreach)

setwd("~/UCD/BILs/fine_mapping_n_eqtl_enrichment/")
files <- list.files(paste0(getwd(), "/05.bin.pvalues"))
common.traits <- c(1:3,5,14,24:27) # traits in common btw. ILs and BILs
files <- files[common.traits] # subset by needed traits
trait.names <- sapply(files, function(x) str_match(x, ".+\\.(.+)\\.txt")[1,2])

# Iterate through the wanted traits and cbind them into a single data.frame and/or list
setwd(paste0(getwd(), "/05.bin.pvalues"))
il.qtl.df <- foreach(i=1:length(files), .combine = cbind) %do% {
  dat <- read.delim(files[i])
  colnames(dat)[2:6] <- paste0(trait.names[i], "_", colnames(dat)[2:6])
  if (i > 1) dat <- dat[2:6]
  dat
}
dim(il.qtl.df)
head(il.qtl.df)

# write.csv(il.qtl.df[,1], "bin.csv", row.names=F, quote=F)

il.qtl.list <- lapply(1:length(files), function(i) {
  dat <- read.delim(files[i])
  dat
} )
names(il.qtl.list) <- trait.names
str(il.qtl.list)
head(il.qtl.list[[1]], 20)

# IL-BIN coordinates
il.bin <- read.csv("../bin2_coord.csv")
head(il.bin)
il.bin$chr <- paste0("SL2.40ch", str_pad(str_replace(as.character(il.bin$bin), "d\\.([0-9]+)[A-Z]", "\\1"), 2, pad="0"))
head(il.bin)
# which((il.bin$end-il.bin$start < 0) ) # 70
# merge with each traits' QTL data frame
il.qtl.list <- lapply(il.qtl.list, function(x) merge( il.bin, x, by="bin", sort=FALSE) )
str(il.qtl.list)
head(il.qtl.list[[1]], 20)
# keep only significant QTL per trait
il.sig.qtl <- lapply(il.qtl.list, function(x) subset(x, adj.pval < 0.05) )
str(il.sig.qtl)
head(il.sig.qtl[[1]])

# put info. directly into RangedData so it can have the associated value columns
il.sig.rd <- lapply(il.sig.qtl, function(x) RangedData(IRanges(start=x$start, end=x$end, names=x$bin2), space=x$chr, values=x[c(1,6:10)]) )
il.sig.rd[[1]]

# DO likewise with BIL QTL data to put into a RangedData list
load("~/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
dim(bin.stats)
head(bin.stats)
# bin.statsin.statsbin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
# num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
# bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

# Load BIL QTL results
load("~/UCD/BILs/final_additive_sparsenet_results/comp.map.Rdata")
# make a list of non.zero.coefs data.frames with 
comp.bil.qtl <- lapply(comp.map, function(x) x$non.zero.coefs )
names(comp.bil.qtl) <- paste0("comp.", names(comp.bil.qtl)) 
rm(comp.map)

load("~/UCD/BILs/final_additive_sparsenet_results/FT.map.Rdata")
ft.bil.qtl <- lapply(FT.map, function(x) x$non.zero.coefs )
rm(FT.map)

load("~/UCD/BILs/final_additive_sparsenet_results/circ.map.Rdata")
circ.bil.qtl <- lapply(circ.map[2:5], function(x) x$non.zero.coefs ) # exclude Area as it has no QTL
names(circ.bil.qtl) <- paste0("circ.", names(circ.bil.qtl)) 
rm(circ.map)

# reorder BIL data
comp.bil.qtl <- comp.bil.qtl[c(4,2,1,3)]
circ.bil.qtl <- circ.bil.qtl[c(2,1,3:4)]

# check that the trait order matches between IL and BIL data, and reorder BIL data
names(il.sig.rd)
names(comp.bil.qtl)
names(ft.bil.qtl)
names(circ.bil.qtl)

# join BIL datasets
bil.qtl <- c(comp.bil.qtl, ft.bil.qtl, circ.bil.qtl)
head(bil.qtl[[1]])

# put into a list of RangedData objects
bil.rd <- lapply(bil.qtl, function(x) RangedData(IRanges(start=x$int0.90.start, end=x$int0.90.end, names=x$bin), space=x$chr, values=x[c(3:9)]) )
bil.rd[[1]]
length(bil.rd)

# find overlaps
ovs <- lapply(1:9, function(i) {
  findOverlaps(query = bil.rd[[i]], subject = il.sig.rd[[i]])
})
names(ovs) <- names(bil.rd)

overlappingBins <- lapply(1:9, function(i) {
  cbind(il.sig.qtl[[i]][as.matrix(ovs[[i]])[,2], ], bil.qtl[[i]][as.matrix(ovs[[i]])[,1], ])
})
names(overlappingBins) <- names(bil.rd)

# APPEND "il." and "bil." to appropriate columns to distinguish
for (i in 1:9) {
  names(overlappingBins[[i]])[1:2] <- paste0("il.", names(overlappingBins[[i]])[1:2]) 
  names(overlappingBins[[i]])[11] <- paste0("bil.", names(overlappingBins[[i]])[11])
}

setwd("../")
save(overlappingBins, file="overlappingBins.Rdata")

# Save overlappingBins' tables
for (i in 1:9) {
  filename <- paste0(names(overlappingBins)[i], "_il.bil.overlappingBins.csv")
  write.csv(overlappingBins[[i]], file=filename, quote=F, row.names=F)
}

#--------
#Script to generate fine mapping circular plots showing IL-QTl and BIL-QTl together
# ...to install and reload the hacked version to label chromosomes' axes in cM if length < 1000
# install.packages("/Users/Dani/UCD/R/circlize_0.2.0.tar.gz", repos=NULL, type="source") 
library(circlize)
library(plyr)
library(scam)
library(scales)
library(stringr)

setwd("~/UCD/BILs/fine_mapping_n_eqtl_enrichment/")
load("~/UCD/BILs/leaf_traits/gen.bin.stats.Rdata") # load gen.bin info
head(gen.bin.stats)

load("~/UCD/BILs/scam_fits.Rdata") # load monotonic spline smoothed regression conversion of bp to cM

load("~/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

il.bin <- read.csv("bin2_coord.csv")
head(il.bin)
# make a chromosome column in the same style as in BIL's gen.bin.stats
il.bin$chr <- paste0("ch", str_pad(str_replace(as.character(il.bin$bin), "d\\.([0-9]+)[A-Z]", "\\1"), 2, pad="0"))
head(il.bin)
class(il.bin$chr)

# convert IL-bin phys. coordinates to genetic distance
chrom <- unique(il.bin$chr)
bin.predict <- vector('list', 12) # list of tables for storing the scam model genetic distance predictions
for (h in 1:length(chrom) ) {
  chrom.tab <- il.bin[il.bin$chr==chrom[h], ]
  # predict genetic distance 3 times with: bin.mid, bin.start, and bin.end
  bin.start <- data.frame(bin.start=chrom.tab$start)
  newdata <- data.frame(position_bp=bin.start$bin.start)
  bin.start$gen.bin.start <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
  bin.end <- data.frame(bin.end=chrom.tab$end)
  newdata <- data.frame(position_bp=bin.end$bin.end)
  bin.end$gen.bin.end <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
  bin.predict[[h]] <- cbind(bin.start, bin.end) # cbind 3 predictions
  bil.tab <- bin.stats[bin.stats$chr==chrom[h], ]
  bil.start <- data.frame(bin.start=bil.tab$bin.start)
  newdata <- data.frame(position_bp=bil.start$bin.start)
  bil.start$gen.bin.start <- predict(fits[[h]], newdata=newdata, type="response", newdata.guaranteed=TRUE)
  # *** zero it differently!!  =>  do so w.r.t. the BIL data so that they match!!
  # I need the unaltered cM prediction (based on BIL data) for the start of the chromosome
  zero.cM <- colwise(function(x, chr.start) x - (chr.start * chr.start/x ) ) # function to proportionally zero the predicted gen. dists.
  bin.predict[[h]][, c(2,4)] <- zero.cM(bin.predict[[h]][, c(2,4)], chr.start = bil.start$gen.bin.start[1]) # zero the cM predictions
  bil.zero <- bil.start$gen.bin.start[1] - (bil.start$gen.bin.start[1] * bil.start$gen.bin.start[1] / bil.start$gen.bin.start[1])
  if(bil.zero != 0) { # re-check zeroing of genDist b/c of numerical issues that sometimes result in a tiny negative or positive value !=0 as the start of the chromosome
    bin.predict[[h]][, c(2,4)] <- bin.predict[[h]][, c(2,4)] - bil.zero
  } else next
}
bin.predict <- do.call(rbind, bin.predict) # rbind bin.predict's 12 elements

il.bin.gen <- cbind(il.bin[, c(1:2,5)], bin.predict) # cbind il.bin' 1st 2 cols (i.e. bin & chr) & bin.predict
#save(il.bin.gen, file="il.bin.gen.Rdata")
names(il.bin.gen)[c(4,6)]
head(il.bin.gen)
# check that the last BIL-bin boundary for each chromosome occurs downstream of the IL one
sapply(1:12, function(i) {
  tmp = chrom[i]
  bil.chr.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$chr==tmp][length(gen.bin.stats$gen.bin.end[gen.bin.stats$chr==tmp])]
  il.chr.end <- il.bin.gen$gen.bin.end[il.bin.gen$chr==tmp][length(il.bin.gen$gen.bin.end[il.bin.gen$chr==tmp])]
  bil.chr.end > il.chr.end
}) # Yes, they do.

# cbind genetic distance predictions to IL QTL data
il.qtl <- lapply(il.sig.qtl, function(x) merge(il.bin.gen[c(1:3,5,7)], x[c(2:4,6:10)], by="bin2", sort=F))
setwd("../")
save(il.qtl, file="il.qtl.Rdata")
bil.qtl
gen.bin.stats
save(bil.qtl, file="bil.qtl.Rdata")

#-----------
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

# Function to plot fine mapping results, iterating through a list of mapping results of phenotypic traits for both ILs and BILs
#---------
load("il.bin.gen.Rdata")
load("il.qtl.Rdata")
load("bil.qtl.Rdata")

plot.fine.map <- function(il.qtl, bil.qtl, gen.bin.stats, il.bin.gen, gen.circ.init, n.dig, start.degree, gap.degree) {  
  for(j in 1:length(bil.qtl) ) {
    dat.name <- names(bil.qtl)[[j]]
    il.dat <- il.qtl[[j]]
    bil.dat <- bil.qtl[[j]]
    il.dat$chr <- as.factor(il.dat$chr)
    bil.dat$chr <- as.factor(substr(bil.dat$chr,7,10))
    il.dat$color[il.dat$slope < 0] <- -1
    il.dat$color[il.dat$slope > 0] <- 1
    bil.dat$color[bil.dat$coefs < 0] <- -1
    bil.dat$color[bil.dat$coefs > 0] <- 1
    bil.dat$coefs <- abs(bil.dat$coefs) # take magnitude of coefficients so that they're all plotted "upwards"
    il.dat$slope <- abs(il.dat$slope)
    tmp <- lapply(1:nrow(bil.dat), function(i) { # create 4 new columns with genetic distance information for the intervals
      gen.bin.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==bil.dat$bin[i] ]
      gen.bin.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==bil.dat$bin[i] ]
      gen.int0.90.start <- gen.bin.stats$gen.bin.start[gen.bin.stats$bin==str_split(bil.dat$int0.90[i], ":")[[1]][1] ]
      gen.int0.90.end <- gen.bin.stats$gen.bin.end[gen.bin.stats$bin==str_split(bil.dat$int0.90[i], ":")[[1]][2] ]
      c(gen.bin.start, gen.bin.end, gen.int0.90.start, gen.int0.90.end)
    })
    gen.dist.cols <- do.call(rbind, tmp)
    colnames(gen.dist.cols) <- c('gen.bin.start', 'gen.bin.end', 'gen.int0.90.start', 'gen.int0.90.end')
    bil.dat <- cbind(bil.dat, gen.dist.cols)
    il.lines <- cbind(gen.circ.init, l1=rep(0, 12), l2=rep(max(il.dat$slope)/2, 12), l3=rep(max(il.dat$slope), 12) ) # y-axis dotted lines for monogenic results ring
    bil.lines <- cbind(gen.circ.init, l1=rep(0, 12), l2=rep(max(bil.dat$coefs)/2, 12), l3=rep(max(bil.dat$coefs), 12) ) # y-axis dotted lines for monogenic results ring
    # setup QTL BED objects for gen. dist. plots
    il.gen.bed <- with(il.dat, data.frame(chr=chr, start=gen.bin.start, end=gen.bin.end, coefs=slope, color=color, row.names=NULL) )
    bil.gen.bed <- with(bil.dat, data.frame(chr=chr, start=gen.bin.start, end=gen.bin.end, coefs=coefs, color=color, row.names=NULL) )
    bil.gen.int90.bed <- with(bil.dat, data.frame(chr=chr, start=gen.int0.90.start, end=gen.int0.90.end, coefs=coefs, color=color, row.names=NULL) )
    bil.list <- list(bil.gen.int90.bed, bil.gen.bed)
    il.ylim <- c(0, max(il.dat$slope) )
    bil.ylim <- c(0, max(bil.dat$coefs) )
    # Setup some plotting parameters
    col.top = "green4"
    col.bottom = "magenta4"
    transp = 0.5 # transparency value used for the interval around the QTL peak; note: transparency = 1 - alpha. 0.75 was the old value, when the mid color was black
    text.col.vector2 <- c(col.top, col.bottom) # label vector for additive QTL legend
    # bin structure BED objects
    il.gen.bin.bed <- with(il.bin.gen, data.frame(chr=chr, start=gen.bin.start, end=gen.bin.end) ) # BED-like data.frame specifying IL gen.bin structure
    gen.bin.bed <- with(gen.bin.stats, data.frame(chr=chr, start=gen.bin.start, end=gen.bin.end) ) # BED-like data.frame specifying BIL gen.bin structure
    # set up plotting parameters to reuse, such as track.margin, etc.
    init.tkHt = 0.05 # track height (in proportion of the circle's radius) of the chromosomal distance track used to initialize the circos plot
    bottom.margin = 0.005 # bottom margin between tracks (in proportion of the circle's radius)
    qtl.tkHt = 0.2 # track height for IL & BIL QTL (in proportion of the circle's radius)
    cpad = 0.01 # cell padding on left and right
    bin.tkHt = 0.05 # BIN structure track height (in proportion of the circle's radius)
    bkgd.col = "gray91" # gray95 works well too, but it's not very visible in some monitors for which gray90 is better
    init.cex = 1.35 # font size scalar for initializing the circos plot
    bin.col=alpha("black", 0.5) # color used for the BINs' rectangles, i.e. black w/ 50% trans
    transp = 0.5 # transparency value used for the interval around the QTL peak; note: transparency = 1 - alpha. 0.75 was the old value, when the mid color was black
    # Individual circular plot function
    circ.plot <- function (plot.path, init.cex, start.degree, bottom.margin, gap.degree, gen.circ.init, init.tkHt, major.by=20, il.lines, bil.lines, qtl.tkHt, cpad, il.ylim, bil.ylim, bkgd.col, 
                           bil.list, gen.bin.bed, il.gen.bed, bin.tkHt, bin.col, text.col.vector2, bil.dat, il.dat) {
      #pdf(file=plot.path, width=8, height=8)
      svg(file=plot.path, width=8, height=8)
      par(mar = c(1, 1, 1, 1), cex = init.cex ) # initialize plotting window and some plotting params; note: if margins aren't all equal then plot won't be circular
      circos.par("start.degree" = start.degree, "track.margin"=c(bottom.margin,0), "gap.degree"=gap.degree) # chromosome-01 is at top right, as opposed starting at 0 degrees i.e. below horizontal on the right
      circos.genomicInitialize(gen.circ.init, sector.width=gen.circ.init$end, track.height=init.tkHt, major.by=major.by ) # initialize circos plot w/ chromosomal axes
      par(cex = 1 ) # change the font scalar back to 1
      # First plot the light grey background and Y-axis dotted lines
      circos.genomicTrackPlotRegion(data = bil.lines, track.height=qtl.tkHt, cell.padding=c(cpad,0,cpad,0), ylim=bil.ylim, panel.fun = function(region, value, ...) {
        cell.xlim = get.cell.meta.data("cell.xlim")
        circos.lines(x=cell.xlim, y=c(value$l1, value$l1), lty=2, lwd=0.3, col="gray50")
        circos.lines(x=cell.xlim, y=c(value$l2, value$l2), lty=2, lwd=0.3, col="gray50" )
        circos.lines(x=cell.xlim, y=c(value$l3, value$l3), lty=2, lwd=0.3, col="gray50" )
      }, bg.col=bkgd.col, bg.border=NA )
      # Then plot the QTL rectangles or bar plots
      circos.genomicTrackPlotRegion(data = bil.list, track.height=qtl.tkHt, track.index=2, cell.padding=c(cpad,0,cpad,0), ylim=bil.ylim, panel.fun = function(region, value, ...) {      
        i = getI(...) # assign the plotting iteration through the data list's elements to a counter
        if (i == 1) { # have 75% transparency if it's 1st plotting iteration, and 0 transp. otherwise
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c(col.bottom, "white", col.top), transparency=transp ) # QTL interval
        } else {
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c(col.bottom, "white", col.top), transparency=0 ) # QTL peak
        }
        circos.genomicRect(region, value, col=sign.col(value$color), ybottom=0, ytop.column=1, border=sign.col(value$color), ...)
      }, bg.border=NA )
      # Plot BIL BIN structure
      circos.genomicTrackPlotRegion(data = gen.bin.bed, ylim=c(0,1), track.height=bin.tkHt, cell.padding=c(0,0,0,0), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col="white", border=bin.col, ybottom=0, ytop=1, lwd=0.5, ...)
      }, bg.border=NA )
      # First plot the light grey background and Y-axis dotted lines
      circos.genomicTrackPlotRegion(data = il.lines, track.height=qtl.tkHt, cell.padding=c(cpad,0,cpad,0), ylim=il.ylim, panel.fun = function(region, value, ...) {
        cell.xlim = get.cell.meta.data("cell.xlim")
        circos.lines(x=cell.xlim, y=c(value$l1, value$l1), lty=2, lwd=0.3, col="gray50")
        circos.lines(x=cell.xlim, y=c(value$l2, value$l2), lty=2, lwd=0.3, col="gray50" )
        circos.lines(x=cell.xlim, y=c(value$l3, value$l3), lty=2, lwd=0.3, col="gray50" )
      }, bg.col=bkgd.col, bg.border=NA )
      # Then plot the QTL rectangles or bar plots
      circos.genomicTrackPlotRegion(data = il.gen.bed, track.height=adt.tkHt, track.index=4, cell.padding=c(cpad,0,cpad,0), ylim=il.ylim, panel.fun = function(region, value, ...) {      
        i = getI(...) # assign the plotting iteration through the data list's elements to a counter
        if (i == 1) { # have 75% transparency if it's 1st plotting iteration, and 0 transp. otherwise
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("blue3", "white", "yellow3"), transparency=transp ) # QTL interval
        } else {
          sign.col = colorRamp2(breaks=c(-1,0,1), colors=c("blue3", "white", "yellow3"), transparency=0 ) # QTL peak
        }
        circos.genomicRect(region, value, col=sign.col(value$color), ybottom=0, ytop.column=1, border=sign.col(value$color), ...)
      }, bg.border=NA )
      # Plot IL BIN structure
      circos.genomicTrackPlotRegion(data = il.gen.bin.bed, ylim=c(0,1), track.height=bin.tkHt, cell.padding=c(0,0,0,0), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col="white", border=bin.col, ybottom=0, ytop=1, lwd=0.5, ...)
      }, bg.border=NA )
      # add legends. for all of these I use several calls to legend() appropriately spaced to acheive the desired color box size; only 1 of these has the legend labels.
      text(x=-1.07, y=1, labels="A. BIL QTL", cex=0.9, pos=4)
      text(x=-1.07, y=0.83, labels="B. BIL bins", cex=0.9, pos=4)
      x.il.labs = 0.7
      text(x=x.il.labs, y=1, labels="C. IL QTL", cex=0.9, pos=4)
      text(x=x.il.labs, y=0.83, labels="D. IL bins", cex=0.9, pos=4)
      # Track labels A = IL QTL, B = IL bins, C = BIL QTL, and D = BIL bins
      text(x=0, y=0.894, offset=0, labels="A", cex=1.1)
      text(x=0, y=0.7146, offset=0, labels="B", cex=1.1)
      text(x=0, y=0.634, offset=0, labels="C", cex=1.1)
      text(x=0, y=0.455, offset=0, labels="D", cex=1.1)
      # BIL color legend
      legend(x=-1.04, y=1, legend=c("", ""), fill=c(col.top, col.bottom), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=-1.04, y=0.982, legend=c("", ""), fill=c(col.top, col.bottom), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=-1.06, y=0.982, legend=c("", ""), fill=c(col.top, col.bottom), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=-1.06, y=1, legend=c("positive", "negative"), fill=c(col.top, col.bottom), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      # IL color legend
      legend(x=x.il.labs+0.03, y=1, legend=c("", ""), fill=c("yellow3", "blue3"), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=x.il.labs+0.03, y=0.982, legend=c("", ""), fill=c("yellow3", "blue3"), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=x.il.labs+0.01, y=0.982, legend=c("", ""), fill=c("yellow3", "blue3"), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      legend(x=x.il.labs+0.01, y=1, legend=c("positive", "negative"), fill=c("yellow3", "blue3"), border=NA, bty="n", bg="white", cex=0.7, x.intersp=1.2, y.intersp=1.2, xjust=0)
      # Y-axis labels for BIL QTL tracks
      if (max(bil.dat$coefs) < 2e-2 ) { # use scientific notation when appropriate
        text(x=0, y=0.936, offset=0, labels=as.character(scientific(max(bil.dat$coefs), 2)), cex=0.6)
        text(x=0, y=0.8478, offset=0, labels=as.character(scientific(max(bil.dat$coefs)/2, 2)), cex=0.6)
      } else {
        text(x=0, y=0.936, offset=0, labels=as.character(round(max(bil.dat$coefs), n.dig)), cex=0.6)
        text(x=0, y=0.8478, offset=0, labels=as.character(round(max(bil.dat$coefs)/2, n.dig)), cex=0.6)
      }
      text(x=0, y=0.756, offset=0, labels="0", cex=0.6)
      # Y-axis labels for IL QTL track
      if (max(il.dat$slope) < 2e-2 ) { # use scientific notation when appropriate
        text(x=0, y=0.6776, offset=0, labels=as.character(scientific(max(il.dat$slope), 2)), cex=0.6)
        text(x=0, y=0.5872, offset=0, labels=as.character(scientific(max(il.dat$slope)/2, 2)), cex=0.6)
      } else {
        text(x=0, y=0.6776, offset=0, labels=as.character(round(max(il.dat$slope), n.dig)), cex=0.6)
        text(x=0, y=0.5872, offset=0, labels=as.character(round(max(il.dat$slope)/2, n.dig)), cex=0.6)
      }
      text(x=0, y=0.495, offset=0, labels="0", cex=0.6)
      dev.off()
      circos.clear()
    }
    genDist.plot.path="~/UCD/BILs/fine_mapping_n_eqtl_enrichment/fine_mapping_plots/"
    #plot.path <- paste0(genDist.plot.path, "/", dat.name, "_fine_mapping", ".pdf")
    plot.path <- paste0(genDist.plot.path, "/", dat.name, "_fine_mapping", ".svg")
    circ.plot(plot.path, init.cex, start.degree, bottom.margin, gap.degree, gen.circ.init, init.tkHt, major.by=20, il.lines, bil.lines, qtl.tkHt, cpad, il.ylim, bil.ylim, bkgd.col,
              bil.list, gen.bin.bed, il.gen.bed, bin.tkHt, bin.col, text.col.vector2, bil.dat, il.dat)
  }
}

plot.fine.map(il.qtl, bil.qtl, gen.bin.stats, il.bin.gen, gen.circ.init, n.dig=2, start.degree=86, gap.degree=c(rep(1.5, 11), 8) )

#---------
# ADD Flowering Time fine-mapping results




