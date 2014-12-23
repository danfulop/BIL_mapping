# Calculate correlation among BINs
# Add 0.95 correlation interval to mapping results
library(hopach) # genomics clustering package w/ convenient correlation/distance matrix function

load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information

load("/Users/Dani/UCD/BILs/leaf_traits/genotab.recoded.Rdata")
bin.code <- genotab.recoded
bin.code <- bin.code[1:(nrow(bin.code)-1), 1:(ncol(bin.code)-2)] # remove BIL & FinBIL cols, and M82 row
bin.code <- t(bin.code) # transpose
dim(bin.code)
head(bin.code)
tail(bin.code)
any(is.na(bin.code)) # FALSE

bin.cor <- distancematrix(bin.code, d="cor") # Calculate correlation distance matrix, this is = 1 - correlation!
bin.cor <- as.matrix(bin.cor)
colnames(bin.cor) <- rownames(bin.code)
rownames(bin.cor) <- rownames(bin.code)
bin.cor[420:432, 420:432]
