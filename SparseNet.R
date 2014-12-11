# SparseNet
#install.packages("sparsenet")
library(reshape2)
library(sparsenet)

setwd("/Users/Dani/UCD/BILs/leaf_traits/")

load("/Users/Dani/UCD/BILs/bgmT.Rdata")
dim(bgmT)
head(bgmT)[,1:6]
head(bgmT)[,1046:1050]
head(bgmT)[,1:6]
dim(bgmT)
genotab <- bgmT[5:nrow(bgmT),] # just BIN genotypes, 1st 4 rows are BIN stats

# Recode BIN genotypes before merging and save data frame 
# Function to recode genotype information into 0s, 1s, and 2s
recode <- function(x) {
  if (x == "M82") {
    x <- 0
  } else if (x == "PEN") {
    x <- 2
  } else if (x == "HET") {
    x <- 1
  }
}

genotab.recoded <- genotab
genotab.recoded[1:(ncol(genotab)-2)] <- apply(genotab.recoded[1:(ncol(genotab)-2)], c(1,2), recode)
genotab.recoded <- droplevels(genotab.recoded)
head(genotab.recoded)[1:6,1:15]
head(genotab.recoded)[1:6,(ncol(genotab.recoded)-5):ncol(genotab.recoded)]
summary(genotab.recoded) # the mean value for the genotype/BIN columns is an indication of what proportion of the population has non-M82 "allele"

lambda.manual = exp(seq(from = -6, to = 0, length = 50))

# Function to fit SparseNet regression -- mapping function
map.fx <- function(df, colN, gt.tab, lambda.manual) {
  trait.name <- names(df)[colN] # (list element) name for focal trait
  trait.dat <- df[[get("trait.name")]] # Assign focal trait data.frame
  trait.dat <- merge(gt.tab, trait.dat, by.x="BIL", by.y="genotype") # merge genotype and trait data
  trait.dat <- trait.dat[-which(names(trait.dat)=="FinBIL.y")] # remove redundant "FinBIL.y"
  names(trait.dat)[which(names(trait.dat)=="FinBIL.x")] <- "FinBIL" # rename FinBIL.x into FinBIL
  trait.dat <- droplevels(trait.dat)
  geno.mat <- as.matrix(trait.dat[2:(ncol(trait.dat)-5)], rownames.force=F)
  response <- as.vector(as.matrix(trait.dat['predPlusResid'], rownames.force=F)) # response vector for sparsenet
  cv.sp <- cv.sparsenet(x=geno.mat, y=response, lambda=lambda.manual)
  plot(cv.sp)
  #cv.sp$parms.min
  #cv.sp$parms.1se 
  
  # SAVE PLOTS, run through all traits, examine plots
  
}

# Complexity data
load("comp.pred.Rdata")
x <- comp.pred
y=1
comp.pred[[1]]
head(comp.pred[[1]], 50)

# Sym PCs
load("sym.pred.Rdata")
#sym.pred[[1]]
df <- sym.pred
colN=1
gt.tab <- genotab.recoded

# Asym PCs
load("asym.pred.Rdata")
asym.pred[[1]]
df <- asym.pred
colN=2
gt.tab <- genotab.recoded


names(comp.pred)
ctall <- comp.pred[[4]]
head(ctall)

all.comp.dat <- merge(bgmAlt, ctall, by.y="genotype", by.x="BIL")
dim(all.comp.dat)
head(all.comp.dat)[,1045:1055]
tail(colnames(all.comp.dat), 10)
colnames(all.comp.dat)[1052]
all.comp.dat <- all.comp.dat[-1052]
tail(colnames(all.comp.dat), 10)
names(all.comp.dat)[2] <- "compAll"
rownames(all.comp.dat) <- all.comp.dat$geno
colnames(all.comp.dat)[1050] <- "FinBIL"

acd <- all.comp.dat[c(1:1050,ncol(all.comp.dat))]
head(acd)[1:15]
summary(acd)
acd <- droplevels(acd)
head(acd)[1:15]
summary(acd)
dim(acd)
levels(acd[1050]) # NULL, not a factor
acd <- acd[c(1,1050:1051,2:1049)]
colnames(acd)


acd.recoded <- acd
acd.recoded[4:ncol(acd.recoded)] <- apply(acd.recoded[4:ncol(acd.recoded)], c(1,2), recode) # it would be best to do this recoding on the genotype data frame before merging it with the trait data
acd.recoded <- droplevels(acd.recoded)
head(acd.recoded)[1:6,1:15]
summary(acd.recoded)

acd <- acd.recoded 
acd <- droplevels(acd)
summary(acd[1:10])
dim(acd)

genodat <- acd[-c(3:4)]
summary(genodat[1:10]) # The first BIN is BIN 2!  Is that okay?  Did I lose a BIN somewhere along the way? YEP!!

ctall <- comp.pred[[4]]

acdm <- as.matrix(acd[4:ncol(acd)], rownames.force=F)
colnames(acdm) <- NULL
compAll <- as.vector(as.matrix(acd[3], rownames.force=F)) # response vector for sparsenet
names(compAll)
acdat <- list(x=as.matrix(acdm), y=compAll)

names(acdat)
dim(acdat$x)

library(sparsenet)
?sparsenet

# Hack lambda0() so that it returns a value
# lambda0 <- function (x, y, weights = rep(1, N), exclude = NULL) 
# {
#   if (length(exclude)) 
#     x = x[, -exclude]
#   N = length(y)
#   ybar = weighted.mean(y, weights)
#   yvar = weighted.mean((y - ybar)^2, weights)
#   y = (y - ybar)/sqrt(yvar)
#   weights = weights/N
#   xbar = t(weights) %*% x
#   xvar = t(weights) %*% (x^2) - xbar^2
#   grad = abs(t(y * weights) %*% x)/sqrt(xvar)
#   max(grad, na.rm=TRUE)
# }
# x=acdat$x
# y=acdat$y
# lambda0 <- lambda0(x=x, y=y)
# min.lambda <- log(lambda0) / 0.0001
# 
# lambda.manual = exp(seq(from = log(lambda0), to = min.lambda, length = 50))
# cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual)
# plot(cv.sp)

lambda.manual = exp(seq(from = -200, to = 0, length = 50))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual)
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se

lambda.manual = exp(seq(from = -100, to = 0, length = 50))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual)
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se

cv.sp$parms.min[[2]] - cv.sp$parms.1se[[2]] # if this difference is = zero, then iterate


lambda.manual = exp(seq(from = -10, to = -2, length = 50))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual)
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se 
# If I had to stick to such a wide lambda range, I should make the length longer so as to
# allow a better estimate of gamma.1se

lambda.manual = exp(seq(from = -6, to = -2, length = 100))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual)
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se

lambda.manual = exp(seq(from = -6, to = -2, length = 100))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual, nfolds=30, trace.it=TRUE)
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se
# increase the nfolds until the paths are smooth
# 

lambda.manual = exp(seq(from = -6, to = -2, length = 200))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual, ngamma=100, nfolds=10, trace.it=TRUE, warm="both")
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se
# The paths don't get any smoother with higher nfolds nor longer lengths.
# Probably what would make it smoother is more gamma values in the range
# 10-fold may be best
# YEP. more values of both gamma and lambda w/ 10-fold CV is the way to go
# ...although thinking about what portion on the BILs I need to get good validation makes
# me think that lower fold CV is better b/c each set has a larger proportion of the BILs,
# e.g. 6- or 5-fold CV seems good.
# maybe 10-fold is a bit better, b/c ther CV error is lower  ...could do that w/ n=200 for
# lambda and gamma

lambda.manual = exp(seq(from = -6, to = -2, length = 100))
cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda=lambda.manual, ngamma=18, nfolds=6, trace.it=TRUE, warm="both")
plot(cv.sp)
cv.sp$parms.min
cv.sp$parms.1se

cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, lambda.min.ratio=0.0001, nlambda=50, nfolds=10, trace.it = TRUE, pmax=1048, warm="both")

sp <- sparsenet(x=acdat$x, y=acdat$y)

cv.sp <- cv.sparsenet(x=acdat$x, y=acdat$y, pmax=358, lambda.min.ratio=0.1, nlambda=20, warm="both", nfolds=287, trace.it = TRUE)


# ** 10 gamma values and 50 lambda may be reasonable ...it's the number used in the CRAN announcement of the package
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
save(cv.sp_narrow, file="cv.sp_narrow.Rdata")
par(mfrow=c(1,1))
plot(cv.sp_narrow) # now that I'm including the residuals *and* the plant replicates 
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
