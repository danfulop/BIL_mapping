# Script to fine map IL-QTL using BIL-QTL

library(stringr)
library(foreach)

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

il.qtl.list <- lapply(1:length(files), function(i) {
  dat <- read.delim(files[i])
  dat
} )
names(il.qtl.list) <- trait.names
str(il.qtl.list)

# treat QTL at adjacent BINs as a single QTL