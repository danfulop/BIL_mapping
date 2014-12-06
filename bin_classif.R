setwd("/Users/Dani/UCD/BILs/")
bg <- read.delim("bin-genotype.BILs.2014-07-25.txt")
summary(bg)
head(bg)
names(bg)

bg$binN <- rownames(bg) # use row numbers as BIN number and assign them to a column
summary(bg[469:471])
head(bg[c(1:5,469:471)])
names(bg)
summary(bg$chr)
table(bg[c(1,471)]) # this table shows the chromosome that each BIN is in. This info could be collated more succintly, as there are no BINs that "span" 2 or more chromosomes

bgl <- as.list(bg[,5:470]) # convert only the BIN genotype columns into a list
length(unique(bgl)) # 434, total number of **unique** BILs
length(which(duplicated(bgl)==TRUE)) # 32 BILs with 1 or more identical BIN genotypes

bgm <- as.matrix(bg)
#bgm <- as.matrix(bg[,5:470])
head(bgm)
#bgm.uniq <- unique(bgm) # subset unique BILs
dim(bgm)
#dim(bgm.uniq) # same dimension, so why did I do this?
#names(bgm.uniq)
#head(bgm.uniq)
bgmT <- t(bgm)
duplicated(bgmT)
which(duplicated(bgmT))
length(which(duplicated(bgmT))) #32 **for this set I could one-by-one find out which BIL is identical ...iterate through other BILs, not in the duplicated set.
head(bgmT)
#TRANSPOSE and then do unique(bgm)
dupl <- NULL
dupv <- which(duplicated(bgmT)) # duplicate vector
#dupvuq <- which(duplicated(bgmT))
#kl = 1

# For each BIL that's duplicated, this loop finds identical BILs, puts all their names into a list, and collates these in a list of lists
for (i in 1:length(dupv)) {
  dup <- bgmT[dupv[i],]
  #bgmTsub <- bgmT[-dupv[i],]
  dltmp <- NULL
  for (j in 1:nrow(bgmT)) {
    if (all(bgmT[j,] == dup)) {
      dltmp <- c(dltmp, rownames(bgmT)[j])
    } else {
      next
    }
    dupl[i] <- list(dltmp)
    #dupl[i] <- list(rownames(bgmT)[dupv[i]]=dltmp)
  }
  #add addition to the list of lists
  #dupl <- c(dupl, dltmp)
}
names(dupl) <- rownames(bgmT)[which(duplicated(bgmT))]
dupl
length(dupl) # duplicates list has 32 items, as expected
dupl.uq <- unique(dupl)
length(dupl.uq) # 28 ...i.e. only a few BILs w/ more than 2 identical lines

dupl.uq  
#save(dupl.uq,file="BIL_replicated_sets.Rdata")
summary(bg[,c(1:4,471)])
head(bg[,c(1:4,471)])

head(bgmT)
head(rownames(bgmT))

#bgmT_rn <- bgmT
#bgmT_rn$FinBIL <- rownames(bgmT_rn)

bgmT <- as.data.frame(bgmT)
bgmT$BIL <- rownames(bgmT) 


for (i in 1:length(dupl.uq)) {
  rep <- unlist(dupl.uq[[i]])
  for (j in 2:length(rep)) {
    rep2rename <- unlist(dupl.uq[[i]][j])
    rownames(bgmT)[rownames(bgmT) == rep2rename] <- paste0(unlist(dupl.uq[[i]][1]),".rep.",dupl.uq[[i]][j]) # rename BIL a concat. name with the name of it's lower-numbered duplicate 1st
  }
}
bgmT$FinBIL <- substr(rownames(bgmT), 1, 7) # FinBIL = final BIL names. i.e. keep just 1 BIL ID for each set of duplicated BILs
summary(as.factor(bgmT$FinBIL))

# Figure out if/which BILs are entirely M82?
v=NULL
for (i in 1:nrow(bgmT)) {
  v[i] <- all(bgmT[i,1:1048] == "M82")
}

which(v) # 181 214 226 255 column indices of BILs that are "mostly" M82, i.e.
# no bins w/ PEN b/c of bin constraints at the start and end of chromosomes
# 207 (181) & 244 (214) may have a tiny bit of PEN at the ends
# BIL_261 (226 index) is fully M82? ...check w/ Mike
# BIL_302 (255 index) is fully M82? ...check w/ Mike 

#bgmT[c(177, 210, 222, 251),1049:1050]
bgmT$FinBIL[c(177, 210, 222, 251)] <- "M82"  # Why are these indices different than the ones above? ~> they were determined to be M82 beforehand

dim(bgmT)
bil.reassign <- bgmT[5:470,1049:1050]
#save(bil.reassign, file="~/UCD/BILs/bil.reassign.Rdata")

names(bgmT)
head(bgmT)
bgmT <- as.data.frame(bgmT)
colnames(bgmT)[1:1048] <- paste0("BIN_", substr(colnames(bgmT)[1:1048], 2, 5))
summary(bgmT)
tail(bgmT)
head(bgmT)
dim(bgmT)
# Clean up some tidbits
bgmT[1:4,1049:1050] <- NA # the 1st 4 rows on on the BIL (1049) and FinBIL (1050) columns should be NA
bgmT <- bgmT[1:470,] # the last row, i.e. 471 should be removed because it's redundant w/ the BIN-# columns, which make up almost the whole data frame
save(bgmT, file="/Users/Dani/UCD/BILs/bgmT.Rdata")
