setwd("/Users/Dani/UCD/BILs/")
bg <- read.delim("bin-genotype.BILs.2014-07-25.txt")
summary(bg)
head(bg)
names(bg)

bg$binN <- rownames(bg)
summary(bg[469:471])
head(bg[c(1:5,469:471)])
names(bg)
summary(bg$chr)
table(bg[c(1,471)])


bgl <- as.list(bg[,5:470])
length(unique(bgl)) # 434
length(duplicated(bgl)) #466

bgm <- as.matrix(bg)
#bgm <- as.matrix(bg[,5:470])
head(bgm)
bgm.uniq <- unique(bgm)
dim(bgm)
dim(bgm.uniq)
names(bgm.uniq)
head(bgm.uniq)
bgmT <- t(bgm)
duplicated(bgmT)
which(duplicated(bgmT))
length(which(duplicated(bgmT))) #32 **for this set I could one-by-one find out which BIL is identical ...iterate through other BILs, not in the duplicated set.
head(bgmT)
#TRANSPOSE and then do unique(bgm)
dupl <- NULL
dupv <- which(duplicated(bgmT))
dupvuq <- which(duplicated(bgmT))
#kl = 1
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
length(dupl)
dupl.uq <- unique(dupl)
length(dupl.uq) #28

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
    rownames(bgmT)[rownames(bgmT) == rep2rename] <- paste0(unlist(dupl.uq[[i]][1]),".rep.",dupl.uq[[i]][j])
  }
}
bgmT$FinBIL <- substr(rownames(bgmT), 1, 7)
summary(as.factor(bgmT$FinBIL))
v=NULL
for (i in 1:nrow(bgmT)) {
  v[i] <- all(bgmT[i,1:1048] == "M82")
}

which(v) # 181 214 226 255 column indices of BILs that are "mostly" M82, i.e.
# no bins w/ PEN b/c of bin constraints at the start and end of chromosomes
# 207 & 244 are 
# BIL_261 (226 index) is fully M82
# BIL_302 (255 index) is fully M82

#bgmT[c(177, 210, 222, 251),1049:1050]
bgmT$FinBIL[c(177, 210, 222, 251)] <- "M82"


bil.reassign <- bgmT[,1049:1050]
save(bil.reassign, file="bil.reassign.Rdata")

names(bgmT)
head(bgmT)
bgmT <- as.data.frame(bgmT)
colnames(bgmT) <- paste0("BIN_", substr(colnames(bgmT), 2, 5))
summary(bgmT)
tail(bgmT)
head(bgmT)
save(bgmT, file="/Users/Dani/UCD/BILs/bgmT.Rdata")
