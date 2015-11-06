# Script to fine map IL-QTL using BIL-QTL

library(stringr)


setwd("~/UCD/BILs/fine_mapping_n_eqtl_enrichment/")
files <- list.files(paste0(getwd(), "/05.bin.pvalues"))
common.traits <- c(1:3,5,24:27) # traits in common btw. ILs and BILs
files <- files[common.traits] # subset by needed traits
trait.names <- sapply(files, function(x) str_match(x, ".+\\.(.+)\\.txt")[1,2])


