setwd("/Users/Dani/UCD/BILs/")
#bg <- read.delim("bin-genotype.BILs.2014-07-25.txt")
bg <- read.delim("~/UCD/BILs/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like.txt")
summary(bg)
head(bg)
names(bg)

dim(bg) # 1049  443
bg$binN <- rownames(bg) # use row numbers as BIN number and assign them to a column
dim(bg) # 1049  444
summary(bg[439:444])
head(bg[c(1:5,439:444)])
# names(bg)
# summary(bg$chr)
# table(bg[c(5,443)]) # this table shows the chromosome that each BIN is in. This info could be collated more succintly, as there are no BINs that "span" 2 or more chromosomes
# 
# bgl <- as.list(bg[,5:443]) # convert only the BIN genotype columns into a list
# length(bgl) # 439
# length(unique(bgl)) # 412, total number of **unique** BILs
# length(which(duplicated(bgl)==TRUE)) # 27 BILs with 1 or more identical BIN genotypes

bgm <- as.matrix(bg)
#bgm <- as.matrix(bg[,5:470])
#head(bgm)
#bgm.uniq <- unique(bgm) # subset unique BILs
#dim(bgm) # 1049  444
#dim(bgm.uniq) # same dimension, so why did I do this?
#names(bgm.uniq)
#head(bgm.uniq)
bgmT <- t(bgm)

# duplicated(bgmT)
# which(duplicated(bgmT))
# length(which(duplicated(bgmT))) # 27 **for this set I could one-by-one find out which BIL is identical ...iterate through other BILs, not in the duplicated set.
# head(bgmT)

#TRANSPOSE and then do unique(bgm)
# dupl <- NULL
# dupv <- which(duplicated(bgmT)) # duplicate vector
#dupvuq <- which(duplicated(bgmT))
#kl = 1

# ** Skip portions of script with duplicate removal **

# For each BIL that's duplicated, this loop finds identical BILs, puts all their names into a list, and collates these in a list of lists
# for (i in 1:length(dupv)) {
#   dup <- bgmT[dupv[i],]
#   #bgmTsub <- bgmT[-dupv[i],]
#   dltmp <- NULL
#   for (j in 1:nrow(bgmT)) {
#     if (all(bgmT[j,] == dup)) {
#       dltmp <- c(dltmp, rownames(bgmT)[j])
#     } else {
#       next
#     }
#     dupl[i] <- list(dltmp)
#     #dupl[i] <- list(rownames(bgmT)[dupv[i]]=dltmp)
#   }
#   #add addition to the list of lists
#   #dupl <- c(dupl, dltmp)
# }
# names(dupl) <- rownames(bgmT)[which(duplicated(bgmT))]
# dupl
# length(dupl) # duplicates list has 32 items, as expected
# dupl.uq <- unique(dupl)
# length(dupl.uq) # 28 ...i.e. only a few BILs w/ more than 2 identical lines

# dupl.uq  
# #save(dupl.uq,file="BIL_replicated_sets.Rdata")
# summary(bg[,c(1:4,471)])
# head(bg[,c(1:4,471)])
# 
# head(bgmT)
# head(rownames(bgmT))

#bgmT_rn <- bgmT
#bgmT_rn$FinBIL <- rownames(bgmT_rn)

bgmT <- as.data.frame(bgmT)
bgmT$BIL <- rownames(bgmT) 

# Skip this too
# for (i in 1:length(dupl.uq)) {
#   rep <- unlist(dupl.uq[[i]])
#   for (j in 2:length(rep)) {
#     rep2rename <- unlist(dupl.uq[[i]][j])
#     rownames(bgmT)[rownames(bgmT) == rep2rename] <- paste0(unlist(dupl.uq[[i]][1]),".rep.",dupl.uq[[i]][j]) # rename BIL a concat. name with the name of it's lower-numbered duplicate 1st
#   }
# }
bgmT$FinBIL <- substr(rownames(bgmT), 1, 7) # FinBIL = final BIL names. i.e. keep just 1 BIL ID for each set of duplicated BILs
summary(as.factor(bgmT$FinBIL))

# Figure out if/which BILs are entirely M82?
# v=NULL
# for (i in 1:nrow(bgmT)) {
#   v[i] <- all(bgmT[i,1:1048] == "M82")
# }
# 
# which(v) # 181 214 226 255 column indices of BILs that are "mostly" M82, i.e.
# # no bins w/ PEN b/c of bin constraints at the start and end of chromosomes
# # 207 (181) & 244 (214) may have a tiny bit of PEN at the ends
# # BIL_261 (226 index) is fully M82? ...check w/ Mike
# # BIL_302 (255 index) is fully M82? ...check w/ Mike 
# 
# #bgmT[c(177, 210, 222, 251),1049:1050]
# bgmT$FinBIL[c(177, 210, 222, 251)] <- "M82"  # Why are these indices different than the ones above? ~> they were determined to be M82 beforehand

dim(bgmT) # 444 1051
bgmT[1:4,1:6]
bgmT[441:444,c(1:6,1048:1051)]

names(bgmT)
head(bgmT)
#bgmT <- as.data.frame(bgmT)
colnames(bgmT)[1:1049] <- paste0("BIN_", substr(colnames(bgmT)[1:1049], 2, 5))
summary(bgmT)
tail(bgmT)
head(bgmT)
#dim(bgmT)
# Clean up some tidbits
bgmT[1:4,1050:1051] <- NA # the 1st 4 rows on on the BIL (1050) and FinBIL (1051) columns should be NA
bgmT <- bgmT[1:443,] # the last row, i.e. 471 should be removed because it's redundant w/ the BIN-# columns, which make up almost the whole data frame
bgmT[1:4,1:6]
bgmT[1:6,1048:1051]
bgmT[439:443,c(1:6,1048:1051)]

# ADD M82 to this data frame
m82.geno <- rep("M82", 1051)
bgmT <- rbind(bgmT, m82.geno)
rownames(bgmT)[444] <- "M82"
pen.geno <- rep("PEN", 1051)
bgmT <- rbind(bgmT, pen.geno)
rownames(bgmT)[445] <- "PEN"

bgmT[439:445,c(1:6,1048:1051)]
bgmT[1:6,c(1:6,1048:1051)]

save(bgmT, file="/Users/Dani/UCD/BILs/bgmT.Rdata")
load("/Users/Dani/UCD/BILs/bgmT.Rdata")
bil.reassign <- bgmT[5:445,1050:1051]
save(bil.reassign, file="~/UCD/BILs/bil.reassign.Rdata")
