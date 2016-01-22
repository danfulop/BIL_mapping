load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata")
head(bin.stats)
nrow(bin.stats)

load("/Users/Dani/UCD/BILs/leaf_traits/gen.bin.stats.Rdata")
head(gen.bin.stats)
nrow(gen.bin.stats)

bins.info <- read.table("/Users/Dani/UCD/BILs/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like", header = T, sep = "\t")
dim(bins.info)

bins.info[1:10, 1:6]
summary(bins.info[, 1:6])

#remove SL2.40 from chromosome label
bins.info$chr <- substr(bins.info$chr,7,10)


# convert the bp coordinates in bins.file to cM genetic distance
library(scam)
library(plyr)

# load(file="/Users/Dani/UCD/BILs/pen12t.Rdata")

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

# load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
# bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
# num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
# bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

# predict genetic distance on bins.file's physical distance coordinates
chrom <- unique(bins.info$chr)
bin.predict <- vector('list', 12) # list of tables for storing the scam model genetic distance predictions
for (h in 1:length(chrom) ) {
  chrom.tab <- bins.info[bins.info$chr==chrom[h], ]
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
  if(bin.predict[[h]][1, 4] != 0) { # re-check zeroing of genDist b/c of numerical issues that sometimes result in a tiny negative or positive value !=0 as the start of the chromosome
    bin.predict[[h]][, c(2,4,6)] <- bin.predict[[h]][, c(2,4,6)] - bin.predict[[h]][1, 4]
  } else next
}
bin.predict <- do.call(rbind, bin.predict) # rbind bin.predict's 12 elements
head(bin.predict)
tail(bin.predict)

gen.bins.info <- bins.info
gen.bins.info[2:4] <- bin.predict[c(2,4,6)] # replace bp with cM
gen.bins.info[1:10, 1:6]
tail(gen.bins.info[1:6])

save(gen.bins.info, file="gen.bins.info.Rdata")
