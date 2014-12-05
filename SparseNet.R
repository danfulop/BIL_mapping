# SparseNet
#install.packages("sparsenet")
library(reshape2)
library(sparsenet)

load("/Users/Dani/UCD/BILs/leaf_traits/ctall.Rdata") # ** replace this w/ a table of predicted values excluding random effects **
load("/Users/Dani/UCD/BILs/bgmT.Rdata")

head(bgmT,10)[c(1047:1050)] # the 1st 4 rows on on the BIL (1049) and FinBIL (1050) columns should be NA
tail(bgmT,10)[c(1047:1048)]
dim(bgmT) # 471 So, the genotype data is in rows --> 5:470
rownames(bgmT)[5:470]
bgmAlt<- bgmT[5:470,]
#bgmAlt$binN <- paste0("BIN_", bgmAlt$binN)
head(bgmAlt)
tail(bgmAlt)
rownames(bgmAlt)
bgmAlt$geno <- rownames(bgmAlt) # this is *not* a good idea!


# mgeno <- melt(bgmAlt, id.vars="binN")
# head(mgeno)
# names(mgeno)[2] <- "BIL"
# head(mgeno)

head(ctall)
all.comp.dat <- merge(ctall, bgmAlt, by="geno")
head(all.comp.dat)[,1057:1058]
names(all.comp.dat)[2] <- "compAll"
rownames(all.comp.dat) <- all.comp.dat$geno

acd <- all.comp.dat[c(1:2,11:ncol(all.comp.dat))]
head(acd)[1:15]
summary(acd)
acd <- droplevels(acd)
head(acd)[1:15]
summary(acd)
levels(acd[3])


recode <- function(x) {
  if (x == "M82") {
    x <- 0
  } else if (x == "PEN") {
    x <- 2
  } else if (x == "HET") {
    x <- 1
  }
}

# levels(acd$BIN_1)
# str <- names(acd)[3]
# str
# levels(acd$assign(str))
# 
# 
# recode <- function(x) {
#   x <- as.factor(x)
#   x <- relevel(x, ref="M82")
# }

acd.recoded <- acd
acd.recoded[3:ncol(acd.recoded)] <- apply(acd.recoded[3:ncol(acd.recoded)], c(1,2), recode)
acd.recoded <- droplevels(acd.recoded)
head(acd.recoded)[1:6,1:15]
summary(acd.recoded)



# for (i in 3:ncol(acd)) {
#   for (j in 1:nrow(acd)) {
#     if (acd[i,j] == "M82") {
#        acd[i,j] <- 0
#     } else if (acd[i,j] == "PEN") {
#        acd[i,j] <- 2
#     } else if (acd[i,j] == "HET") {
#        acd[i,j] <- 1
#     }
#   }
# }


acd <- acd.recoded 
acd <- droplevels(acd)
summary(acd[1:10])
dim(acd)

acdm <- as.matrix(acd[3:ncol(acd)], rownames.force=F)
colnames(acdm) <- NULL
compAll <- as.vector(as.matrix(acd[2], rownames.force=F))
names(compAll)
acdat <- list(x=as.matrix(acdm), y=compAll)

names(acdat)

library(sparsenet)
?sparsenet
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, pmax=358, lambda.min.ratio=0.1, nlambda=20, warm="both", nfolds=287, trace.it = TRUE)

# Try running lambda0() with the above x and y to test...
x=acdat$x
y=acdat$y
lambda0 = 0

cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, pmax=358, nlambda=50, ngamma=9, max.gamma=150, min.gamma=1.000001, lambda.min.ratio=0.01, nlambda=20, warm="both", nfolds=287, trace.it = TRUE)

#cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, trace.it = TRUE)
# lambda.manual = exp(seq(from = log(max.lambda), to = log(max.lambda * lambda.min.ratio), length = nlambda))
lambda.manual = exp(seq(from = log(1), to = log(1 * 0.00001), length = 150))
#lambda0=0
#lambda0=NULL

#cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, pmax=1000, ngamma=9, max.gamma=10000, min.gamma=1.000001, lambda=lambda.manual, warm="both", nfolds=357, trace.it = TRUE)
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, ngamma=9, max.gamma=10000, min.gamma=1.000001, lambda=lambda.manual, warm="both", nfolds=357, trace.it = TRUE)
par(mfrow=c(1,1))
plot(cv.sp)
cv.sp$parms.1se; cv.sp$parms.min
cv.sp$parms.min[[1]] # gamma=10000
cv.sp$parms.min[[2]] # lambda=0.03606177

lambda.manual = exp(seq(from = log(0.15), to = log(0.015), length = 150))
cv.sp_narrow <- cv.sparsenet(x=acdat$x, y=acdat$y, ngamma=20, max.gamma=40, min.gamma=1.000001, lambda=lambda.manual, warm="both", nfolds=100, trace.it = TRUE)
par(mfrow=c(1,1))
plot(cv.sp_narrow)
cv.sp_narrow$parms.1se; cv.sp_narrow$parms.min
str(cv.sp_narrow)
cv.sp_narrow$which.min; cv.sp_narrow$which.1se
cv.sp_narrow$nzero

lambda.manual = exp(seq(from = log(0.1353353), to = log(0.01831564), length = 150))
spn <- sparsenet(x=acdat$x, y=acdat$y, pmax=130, ngamma=9, lambda=lambda.manual, warm="both")
spn
par(mfrow=c(3,3))
plot(spn)
plot(spn, xvar="lambda")

lambda.manual = exp(seq(from = log(0.15), to = log(0.015), length = 150))
spn <- sparsenet(x=acdat$x, y=acdat$y, pmax=150, max.gamma=10, ngamma=9, lambda=lambda.manual, warm="both")
spn
par(mfrow=c(3,3))
plot(spn)
plot(spn, xvar="lambda")

lambda.manual = exp(seq(from = log(0.15), to = log(0.015), length = 150))
#lambda.manual = exp(seq(from = log(0.07368285), to = log(3.619442e-02), length = 50))
#gamma.manual = seq(from=9.900000e+35, to=5.15260563, length = 9)
spn.min <- sparsenet(x=acdat$x, y=acdat$y, pmax=150, ngamma=9, max.gamma=10, min.gamma=5.15260563, lambda=lambda.manual, warm="both")
#spn.min <- sparsenet(x=acdat$x, y=acdat$y, pmax=150, ngamma=9, max.gamma=10, lambda=lambda.manual, warm="both")
spn.min
str(spn.min)
names(spn.min)
spn.min$parms
par(mfrow=c(3,3))
plot(spn.min)
plot(spn.min, xvar="lambda")
plot(spn.min, xvar="norm")

spn.min$coefficients$g9  # preferred set of parameters
dim(spn.min$coefficients$g9$beta) # [1] 1048   74
spn.min$coefficients$g9$beta[,74]
class(spn.min$coefficients$g9$beta[,74])
allComp.results.additive <-  cbind(bg[c(1:4,471)], spn.min$coefficients$g9$beta[,74])
head(allComp.results.additive)
names(allComp.results.additive)[6] <- "coef"
allComp.results.additive$chr <- substr(allComp.results.additive$chr,7,11)
allComp.results.additive$bin.len <- allComp.results.additive$bin.end - allComp.results.additive$bin.start
head(allComp.results.additive)
summary(allComp.results.additive)
hist(allComp.results.additive$coef[allComp.results.additive$coef!=0])
sort(allComp.results.additive$coef[allComp.results.additive$coef!=0])

res1 <- allComp.results.additive[allComp.results.additive$chr=="ch01",]
summary(res1)
head(res1)
dim(res1)
hist(res1$coef[res1$coef!=0])
res1$coef[res1$coef!=0]
res1$bin.start.arb <- (as.numeric(res1$binN) -1) * 10
res1$bin.end.arb <- (as.numeric(res1$binN)) * 10
res1$bin.mid.arb <- res1$bin.start.arb + 5
res1.nz <- res1[res1$coef!=0,]
class(res1.nz)
res1.nz
library(ggplot2)
# code for physical dist. plot
ggplot(res1) + geom_segment(aes(y=-0.2, yend=0.2, x=bin.start, xend=bin.start)) + geom_segment(aes(y=-0.2, yend=0.2, x=bin.end, xend=bin.end)) + 
  geom_segment(aes(y=-0.2, yend=-0.2, x=min(bin.start), xend=max(bin.end) ) ) + geom_segment(aes(y=0.2, yend=0.2, x=min(bin.start), xend=max(bin.end) ) ) +
  geom_point(data=res1.nz, aes(x=bin.mid, y=coef), shape=17, size=3 , color="red") + theme_bw() + labs(y="Coefficients", x="Bins")

# REMOVE x-axis ticks from this plot, as they're meaningless
ggplot(res1) + geom_segment(aes(y=-0.1, yend=0.1, x=bin.start.arb, xend=bin.start.arb), size=0.1) + geom_segment(aes(y=-0.1, yend=0.1, x=bin.end.arb, xend=bin.end.arb), size=0.1) + 
  geom_segment(aes(y=-0.1, yend=-0.1, x=min(bin.start.arb), xend=max(bin.end.arb) ), size=0.1) + geom_segment(aes(y=0.1, yend=0.1, x=min(bin.start.arb), xend=max(bin.end.arb) ), size=0.1) +
  geom_text(aes(label=binN, x=bin.mid.arb, y=0), size=3.5) + 
  geom_point(data=res1.nz, aes(x=bin.mid.arb, y=coef), shape=17, size=3 , color="red") + theme_bw() + labs(y="Coefficients", x="Bins")


res2 <- allComp.results.additive[allComp.results.additive$chr=="ch02",]
summary(res2)
head(res2)
dim(res2)
hist(res2$coef[res2$coef!=0])
res2$coef[res2$coef!=0]
res2$bin.start.arb <- (as.numeric(res2$binN) -1) * 10
res2$bin.end.arb <- (as.numeric(res2$binN)) * 10
res2$bin.mid.arb <- res2$bin.start.arb + 5
res2.nz <- res2[res2$coef!=0,]
class(res2.nz)
res2.nz
library(ggplot2)
# code for physical dist. plot
ggplot(res2) + geom_segment(aes(y=-0.2, yend=0.2, x=bin.start, xend=bin.start)) + geom_segment(aes(y=-0.2, yend=0.2, x=bin.end, xend=bin.end)) + 
  geom_segment(aes(y=-0.2, yend=-0.2, x=min(bin.start), xend=max(bin.end) ) ) + geom_segment(aes(y=0.2, yend=0.2, x=min(bin.start), xend=max(bin.end) ) ) +
  geom_point(data=res2.nz, aes(x=bin.mid, y=coef), shape=17, size=3 , color="red") + theme_bw() + labs(y="Coefficients", x="Bins")

ggplot(res2) + geom_segment(aes(y=-0.1, yend=0.1, x=bin.start.arb, xend=bin.start.arb), size=0.1) + geom_segment(aes(y=-0.1, yend=0.1, x=bin.end.arb, xend=bin.end.arb), size=0.1) + 
  geom_segment(aes(y=-0.1, yend=-0.1, x=min(bin.start.arb), xend=max(bin.end.arb) ), size=0.1) + geom_segment(aes(y=0.1, yend=0.1, x=min(bin.start.arb), xend=max(bin.end.arb) ), size=0.1) +
  geom_text(aes(label=binN, x=bin.mid.arb, y=0), size=3.5) + 
  geom_point(data=res2.nz, aes(x=bin.mid.arb, y=coef), shape=17, size=3 , color="red") + theme_bw() + labs(y="Coefficients", x="Bins")

res7 <- allComp.results.additive[allComp.results.additive$chr=="ch07",]
summary(res7)
head(res7)
dim(res7)
hist(res7$coef[res7$coef!=0])
res7$coef[res7$coef!=0]
res7$bin.start.arb <- (as.numeric(res7$binN) -1) * 10
res7$bin.end.arb <- (as.numeric(res7$binN)) * 10
res7$bin.mid.arb <- res7$bin.start.arb + 5
res7.nz <- res7[res7$coef!=0,]
class(res7.nz)
res7.nz
library(ggplot2)
# code for physical dist. plot
ggplot(res7) + geom_segment(aes(y=-0.2, yend=0.2, x=bin.start, xend=bin.start)) + geom_segment(aes(y=-0.2, yend=0.2, x=bin.end, xend=bin.end)) + 
  geom_segment(aes(y=-0.2, yend=-0.2, x=min(bin.start), xend=max(bin.end) ) ) + geom_segment(aes(y=0.2, yend=0.2, x=min(bin.start), xend=max(bin.end) ) ) +
  geom_point(data=res7.nz, aes(x=bin.mid, y=coef), shape=17, size=3 , color="red") + theme_bw() + labs(y="Coefficients", x="Bins")

# REMOVE x-axis ticks from this plot, as they're meaningless
ggplot(res7) + geom_segment(aes(y=-0.1, yend=0.1, x=bin.start.arb, xend=bin.start.arb), size=0.1) + geom_segment(aes(y=-0.1, yend=0.1, x=bin.end.arb, xend=bin.end.arb), size=0.1) + 
  geom_segment(aes(y=-0.1, yend=-0.1, x=min(bin.start.arb), xend=max(bin.end.arb) ), size=0.1) + geom_segment(aes(y=0.1, yend=0.1, x=min(bin.start.arb), xend=max(bin.end.arb) ), size=0.1) +
  geom_text(aes(label=binN, x=bin.mid.arb, y=0), size=3.5) + 
  geom_point(data=res7.nz, aes(x=bin.mid.arb, y=coef), shape=17, size=3 , color="red") + theme_bw() + labs(y="Coefficients", x="Bins")


# Make chromosome-specific bin-numbers?

## RETRY with the predicted values after fitting ONLY random effects.  ...and do like ~200-fold CV  ...that way all of the genotypes are present at all times.

##--
train.data=gendata(100,1000,nonzero=30,rho=0.3,snr=3)
fit=sparsenet(train.data$x,train.data$y)
par(mfrow=c(3,3))
plot(fit)
plot(fit, xvar="lambda")
par(mfrow=c(1,1))
fitcv=cv.sparsenet(train.data$x,train.data$y,trace.it=TRUE)
plot(fitcv)
