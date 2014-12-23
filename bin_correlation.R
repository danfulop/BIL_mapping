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
length(which(bin.cor == 0)) # 1117
bin.cor[which(bin.cor == 0)] # these are perfect correlations
length(which(bin.cor < 0)) # 2096
bin.cor[which(bin.cor < 0)] # ALL these cells are -99999, which is odd and seems like a numerical issue.
bin.cor[which(bin.cor < 0)] <- 1 # making these dist = 1, so that the correlation is 0
bin.cor <- 1 - bin.cor # convert from distnce to correlation matrix
save(bin.cor, file="/Users/Dani/UCD/BILs/leaf_traits/bin.cor.Rdata")
load("/Users/Dani/UCD/BILs/leaf_traits/bin.cor.Rdata")

load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names

# Function to add correlation intervals to mapping results
#--------
add.cor <- function(map.dat, bin.cor, bin.stats) {
  for(i in 1:length(map.dat) ) {
    if(map.dat[[i]]$n.coef==0) {
      next
    } else {
      for (j in 1:nrow(map.dat[[i]]$non.zero.coefs) ) {
        tmp95 <- which(bin.cor[rownames(bin.cor) == map.dat[[i]]$non.zero.coefs$bin[j], ] > 0.95)
        first.bin.tmp95 <- tmp95[1]
        name.first.bin.tmp95 <- names(tmp95)[1]
        last.bin.tmp95 <- tmp95[length(tmp95)]
        name.last.bin.tmp95 <- names(tmp95)[length(tmp95)]
        # ADD columns int0.95 & int0.90 (names -- concat e.g. 421:430), int0.95.start, int0.90.start, int0.95.end, int0.90.end
        map.dat[[i]]$non.zero.coefs$int0.95[j] <- paste(name.first.bin.tmp95, name.last.bin.tmp95, sep=":")
        map.dat[[i]]$non.zero.coefs$int0.95.start[j] <- as.numeric(as.character(bin.stats$bin.start[first.bin.tmp95]))
        map.dat[[i]]$non.zero.coefs$int0.95.end[j] <- as.numeric(as.character(bin.stats$bin.end[last.bin.tmp95]))
        tmp90 <- which(bin.cor[rownames(bin.cor) == map.dat[[i]]$non.zero.coefs$bin[j], ] > 0.9)
        first.bin.tmp90 <- tmp90[1]
        name.first.bin.tmp90 <- names(tmp90)[1]
        last.bin.tmp90 <- tmp90[length(tmp90)]
        name.last.bin.tmp90 <- names(tmp90)[length(tmp90)]
        map.dat[[i]]$non.zero.coefs$int0.90[j] <- paste(name.first.bin.tmp90, name.last.bin.tmp90, sep=":")
        map.dat[[i]]$non.zero.coefs$int0.90.start[j] <- as.numeric(as.character(bin.stats$bin.start[first.bin.tmp90]))
        map.dat[[i]]$non.zero.coefs$int0.90.end[j] <- as.numeric(as.character(bin.stats$bin.end[last.bin.tmp90]))
      }
    }
  }
  map.dat
}
#----------

# Run the add.cor function on all the additive results
#-------
setwd("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/")
load("comp.map.Rdata")
comp.map <- add.cor(map.dat=comp.map, bin.cor, bin.stats)
save(comp.map, file="comp.map.Rdata")
load("circ.map.Rdata")
circ.map <- add.cor(map.dat=circ.map, bin.cor, bin.stats)
save(circ.map, file="circ.map.Rdata")
load("sym.map.Rdata")
sym.map <- add.cor(map.dat=sym.map, bin.cor, bin.stats)
save(sym.map, file="sym.map.Rdata")
load("asym.map.Rdata")
asym.map <- add.cor(map.dat=asym.map, bin.cor, bin.stats)
save(asym.map, file="asym.map.Rdata")
load("FT.map.Rdata")
FT.map <- list(FT=FT.map) # modify FT results so that they fit the structure of the plotting function, i.e. list with sublist(s)
FT.map <- add.cor(map.dat=FT.map, bin.cor, bin.stats)
save(FT.map, file="FT.map.Rdata")
#--------

