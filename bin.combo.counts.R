library(parallel)
library(ggplot2)

setwd("/Users/Dani/UCD/BILs/leaf_traits/")
load("epitab.Rdata")
load("epi.bin.stats.Rdata")

# epitab[1:10, 1:10]
# dim(epitab)
# epitab[1:10, 1:10]
# epitab[1:10, (ncol(epitab)-10):ncol(epitab)]
epiCols <- epitab[1:(nrow(epitab)-1), 1050:(ncol(epitab)-2)]

length(epiCols[,1]) # 439
length(which(epiCols[,1]==0)) # 427
length(which(epiCols[,1]!=0)) # 12

epi.bin.combo.count <- mclapply(1:ncol(epiCols), function(x) length(which(epiCols[x]!=0)), mc.cores=4)
epi.bin.combo.count <- do.call(c, epi.bin.combo.count)
names(epi.bin.combo.count) <- colnames(epiCols)
epi.bin.combo.count <- data.frame(cbind(bin.combo=colnames(epiCols), count=epi.bin.combo.count))
head(epi.bin.combo.count)
epi.bin.combo.count$count <- as.numeric(as.character(epi.bin.combo.count$count))
setwd("/Users/Dani/UCD/BILs/")
save(epi.bin.combo.count, file="epi.bin.combo.count.Rdata")
load("epi.bin.combo.count.Rdata")
write.csv(epi.bin.combo.count, file="epi.bin.combo.count.csv")
hist((epi.bin.combo.count$count))
qplot(x=count, data=epi.bin.combo.count, geom="histogram", main="BIN combination counts in BILs", xlab="Number of BILs with a given BIN combination")

count.tab <- table(epi.bin.combo.count$count)
count.tab [1]
count.dat <- cbind(count=as.numeric(names(count.tab)), n.times=as.numeric(count.tab))
count.dat <- as.data.frame(count.dat)
save(count.dat, file="count.dat.Rdata")
load("count.dat.Rdata")
write.csv(count.dat, file="count.dat.csv")

count.dat <- rbind(c(count=0, n.times=243385), count.dat)

n1 <- nls(n.times ~ n0 * exp(l*count), data=count.dat, start=list(n0=200000, l=-1) )
n1
summary(n1)
fitted(n1)
coef(n1)[2]
avg <- (1/coef(n1)[2])*-1
avg

# p1 <- glm(n.times ~ count, family=poisson(), data=count.dat)
# p1
# summary(p1)
# fitted.values(p1)
# effects(p1)
# exp(effects(p1))
# residuals(p1)
# exp(residuals(p1))
# 
# exp(11.5617205) # 105000.5
# exp(-0.2947666) # ~0.75



