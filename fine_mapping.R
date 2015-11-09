# Script to fine map IL-QTL using BIL-QTL

library(stringr)
library(IRanges)
library(foreach)
library(reshape2)

setwd("~/UCD/BILs/fine_mapping_n_eqtl_enrichment/")
files <- list.files(paste0(getwd(), "/05.bin.pvalues"))
common.traits <- c(1:3,5,24:27) # traits in common btw. ILs and BILs
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
load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
dim(bin.stats)
head(bin.stats)
# bin.statsin.statsbin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
# num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
# bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

# Load BIL QTL results
load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/comp.map.Rdata")
# make a list of non.zero.coefs data.frames with 
comp.bil.qtl <- lapply(comp.map, function(x) x$non.zero.coefs )
names(comp.bil.qtl) <- paste0("comp.", names(comp.bil.qtl)) 
rm(comp.map)

load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/circ.map.Rdata")
circ.bil.qtl <- lapply(circ.map[2:5], function(x) x$non.zero.coefs ) # exclude Area as it has not QTL
names(circ.bil.qtl) <- paste0("circ.", names(circ.bil.qtl)) 
rm(circ.map)

# check that the trait order matches between IL and BIL data, and reorder BIL data
names(il.sig.rd)
names(comp.bil.qtl)
names(circ.bil.qtl)

# reorder BIL data
comp.bil.qtl <- comp.bil.qtl[c(4,2,1,3)]
circ.bil.qtl <- circ.bil.qtl[c(2,1,3:4)]

# join BIL datasets
bil.qtl <- c(comp.bil.qtl, circ.bil.qtl)
head(bil.qtl[[1]])

# put into a list of RangedData objects
bil.rd <- lapply(bil.qtl, function(x) RangedData(IRanges(start=x$int0.90.start, end=x$int0.90.end, names=x$bin), space=x$chr, values=x[c(3:9)]) )
bil.rd[[1]]
length(bil.rd)

# find overlaps
ovs <- lapply(1:8, function(i) {
  findOverlaps(query = bil.rd[[i]], subject = il.sig.rd[[i]])
})
names(ovs) <- names(bil.rd)

overlappingBins <- lapply(1:8, function(i) {
  cbind(il.sig.qtl[[i]][as.matrix(ovs[[i]])[,2], ], bil.qtl[[i]][as.matrix(ovs[[i]])[,1], ])
})
names(overlappingBins) <- names(bil.rd)

# APPEND "il." and "bil." to appropriate columns to distinguish
for (i in 1:8) {
  names(overlappingBins[[i]])[1:2] <- paste0("il.", names(overlappingBins[[i]])[1:2]) 
  names(overlappingBins[[i]])[11] <- paste0("bil.", names(overlappingBins[[i]])[11])
}

setwd("../")
save(overlappingBins, file="overlappingBins.Rdata")

# Save overlappingBins' tables
for (i in 1:8) {
  filename <- paste0(names(overlappingBins)[i], "_il.bil.overlappingBins.csv")
  write.csv(overlappingBins[[i]], file=filename, quote=F, row.names=F)
}

# Make circular plots with IL-QTL + BIL-QTL for fine-mapping viz!
