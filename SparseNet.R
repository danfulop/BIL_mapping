#install.packages("sparsenet")
library(reshape2)
library(sparsenet)
library(ggplot2)
library(foreach)
library(doParallel)

setwd("/Users/Dani/UCD/BILs/leaf_traits/")

load("/Users/Dani/UCD/BILs/bgmT.Rdata")
genotab <- bgmT[5:nrow(bgmT),] # just BIN genotypes, 1st 4 rows are BIN stats
genotab <- droplevels(genotab)
genotab <- genotab[genotab$BIL!="PEN",] # remove PEN
# save(genotab, file="genotab.Rdata")
bin.stats <- bgmT[1:4,1:(ncol(bgmT)-2)]
bin.stats <- droplevels(bin.stats)
bin.stats <- data.frame(t(bin.stats))
bin.stats$bin <- rownames(bin.stats)
bin.stats$bin <- factor(bin.stats$bin, levels=bin.stats$bin)
bin.stats <- bin.stats[,c(5,1:4)]
#save(bin.stats, file="bin.stats.Rdata")
load("bin.stats.Rdata")

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
genotab.recoded[1:(ncol(genotab.recoded)-2)] <- apply(genotab.recoded[1:(ncol(genotab.recoded)-2)], c(1,2), recode)
genotab.recoded <- droplevels(genotab.recoded)
#save(genotab.recoded, file="genotab.recoded.Rdata")

##-----------
load("genotab.recoded.Rdata")

# Function to fit SparseNet regression -- mapping function
map.fx <- function(datl, ln, gt.tab, bin.stats, lambda.manual) {
  trait.name <- names(datl)[ln] # (list element) name for focal trait
  trait.dat <- datl[[get("trait.name")]] # Assign focal trait data.frame
  trait.dat <- merge(gt.tab, trait.dat, by.x="BIL", by.y="genotype") # merge genotype and trait data
  trait.dat <- trait.dat[-which(names(trait.dat)=="FinBIL.y")] # remove redundant "FinBIL.y"
  names(trait.dat)[which(names(trait.dat)=="FinBIL.x")] <- "FinBIL" # rename FinBIL.x into FinBIL
  trait.dat <- droplevels(trait.dat)
  geno.mat <- as.matrix(trait.dat[2:(ncol(trait.dat)-5)], rownames.force=F) # genotype matrix, i.e. Xs
  response <- as.vector(as.matrix(trait.dat['predPlusResid'], rownames.force=F)) # response vector for sparsenet, i.e. Ys
  tmp <- vector("list", length=10) # temp list for storing 1se parameters
  for(j in 1:10) { # redo cross-validation 10 times to get more "stable" tuning parameters
    cv.sp <- cv.sparsenet(x=geno.mat, y=response, lambda=lambda.manual, ngamma=18, nfolds=6, warm="both")
    tmp[[j]] <- cv.sp$parms.1se 
  }
  mean.gamma <- mean(unlist(tmp)[seq(1,19,2)]) # mean gamma
  mean.lambda <- mean(unlist(tmp)[seq(2,20,2)]) # mean lambda
  # Fit sparsenet one more time with the mean tuning parameters
  lambda.seq <- seq(from=exp(-5.5), to=mean.lambda, length=10) # lambda sequence ending in mean.lambda
  gamma.seq <- seq(from=150 , to=mean.gamma, length=9) # gamma sequence ending in mean.gamma
  sp.fit <- sparsenet(x=geno.mat, y=response, lambda=lambda.seq, gamma=gamma.seq, warm="both") # full dataset sparsenet fit
  coefs <- sp.fit$coefficients$g9$beta[,1] # save preferred set of coefficients
  coefs <- cbind(bin.stats, coefs) # combine with bin information
  non.zero.coefs <- coefs[coefs!=0] # non-zero coefficients
  n.coef <- length(non.zero.coefs) # number of non-zero coefficients
  results <- list(coefs=coefs, non.zero.coefs=non.zero.coefs, n.coef=n.coef, gamma=mean.gamma, lambda=mean.lambda, sp.fit=sp.fit)
  results
}

# SAVE PLOTS, run through all traits, examine plots => if I reuse this code I'll have to add back dat.name to the function's variables
#   plot.path <- paste0("/Users/Dani/UCD/BILs/CVplots/", dat.name, ".", trait.name, ".CVplot", ".pdf")
#   pdf(file=plot.path)
#   plot(cv.sp)
#   dev.off()

## Fit SparseNet regressions
lambda.manual = exp(seq(from = -5.5, to = -1.5, length = 50))
# Complexity data 
load("comp.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
# Initialize lists in which to save the ggplot objects and the data.frames with model-fitted means 
comp.map <- vector("list", length=length(comp.pred))
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
system.time(comp.map <- foreach(i=1:length(comp.pred), .options.multicore=mcoptions) %dopar% { # run loop
  comp.map[i] <- map.fx(datl=comp.pred, ln=i, gt.tab=genotab.recoded, bin.stats=bin.stats, lambda.manual=lambda.manual)
})
names(comp.map) <- names(comp.pred)
save(comp.map, file="comp.map.Rdata")

# Circ. data
load("circ.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
# Initialize lists in which to save the ggplot objects and the data.frames with model-fitted means 
circ.map <- vector("list", length=length(circ.pred))
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
system.time(circ.map <- foreach(i=1:length(circ.pred), .options.multicore=mcoptions) %dopar% { # run loop
  circ.map[i] <- map.fx(datl=circ.pred, ln=i, gt.tab=genotab.recoded, bin.stats=bin.stats, lambda.manual=lambda.manual)
})
names(circ.map) <- names(circ.pred)
save(circ.map, file="circ.map.Rdata")

# Sym PCs
load("sym.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
# Initialize lists in which to save the ggplot objects and the data.frames with model-fitted means 
sym.map <- vector("list", length=length(sym.pred))
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
system.time(sym.map <- foreach(i=1:length(sym.pred), .options.multicore=mcoptions) %dopar% { # run loop
  sym.map[i] <- map.fx(datl=sym.pred, ln=i, gt.tab=genotab.recoded, bin.stats=bin.stats, lambda.manual=lambda.manual)
})
names(sym.map) <- names(sym.pred)
save(sym.map, file="sym.map.Rdata")

# Asym PCs
load("asym.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
# Initialize lists in which to save the ggplot objects and the data.frames with model-fitted means 
asym.map <- vector("list", length=length(asym.pred))
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
system.time(asym.map <- foreach(i=1:length(asym.pred), .options.multicore=mcoptions) %dopar% { # run loop
  asym.map[i] <- map.fx(datl=asym.pred, ln=i, gt.tab=genotab.recoded, bin.stats=bin.stats, lambda.manual=lambda.manual)
})
names(asym.map) <- names(asym.pred)
save(asym.map, file="asym.map.Rdata")


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

## RETRY with the predicted values after fitting ONLY random effects.

##--
train.data=gendata(100,1000,nonzero=30,rho=0.3,snr=3)
fit=sparsenet(train.data$x,train.data$y)
par(mfrow=c(3,3))
plot(fit)
plot(fit, xvar="lambda")
par(mfrow=c(1,1))
fitcv=cv.sparsenet(train.data$x,train.data$y,trace.it=TRUE)
plot(fitcv)
