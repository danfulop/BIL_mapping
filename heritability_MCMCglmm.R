# Script to calculate broad-sense heritabilities of the BIL traits

# Pseudo-code:
# Use a random effects model on the BILs w/o M82:
# trait ~ 1 + (1 | BIL) + (1 | block) + (1 | expt)
    # do I need to use sum to zero contrast?!
        # maybe not since there are no fixed effects
    # I could do away with the block effect because there are so few blocks that it's hard to properly estimate a block effect

# does the mixed crossing scheme present a problem?

library(lme4)
library(ggplot2)
library(MCMCglmm)
library(coda)
library(parallel)
library(doParallel)
library(foreach)

#----
# FUNCTIONS
# Function to calculate broad-sense heritability, repeatability, and "stochastic index" from lme4 models
paramMLEs <- function(x, y, randform) {  # Include as input column index of response values, i.e. y
  resp <- names(x)[y] # assign the column name of the 'y' index
  formula <- as.formula( paste0( get("resp"), " ~ 1 + (1|FinBIL) + (1|block) ", get("randform") ) )
  fit <- lmer(formula=formula, data=x, REML=TRUE)
  tab <- as.data.frame(VarCorr(fit)) # save random effects table
  gen.idx <- which(tab$grp == "FinBIL")
  if (gen.idx == 3) {
    h2 <- with(tab, vcov[3] / sum(vcov))
    r2 <- with(tab, sum(vcov[1:3]) / sum(vcov))
    stoch <- with(tab, sum(vcov[4:5]) / sum(vcov[1:2]))
  } else {
    h2 <- with(tab, vcov[2] / sum(vcov))
    r2 <- with(tab, sum(vcov[1:2]) / sum(vcov))
    stoch <- with(tab, sum(vcov[3:4]) / vcov[1])
  }
  total.var <- with(tab, sum(vcov))
  res <- c(h2=h2, r2=r2, stoch=stoch, total.var=total.var)
  return(res)
}

# Function to calculate random effects MLEs to use for building "empirical bayes" weakly informative priors
mle.fx <- function(x, y, randform) {  # Include as input column index of response values, i.e. y
  resp <- names(x)[y] # assign the column name of the 'y' index
  formula <- as.formula( paste0( get("resp"), " ~ 1 + (1|FinBIL) + (1|block) ", get("randform") ) )
  fit <- lmer(formula=formula, data=x, REML=TRUE)
  tab <- as.data.frame(VarCorr(fit)) # save random effects table
  tab
}

# Function to run MCMC models
# dat=dataset, leaflet=logical var. for leaflet.type, mle.tab=MLEs to use for tayloring the priors
# fit1 <- MCMCglmm(fixed = fixed.formula, random = random.formula, data=x, prior=prior, nitt=26000, thin=20, burnin=6000, pr=T)
# fit2 <- MCMCglmm(fixed = fixed.formula, random = random.formula, data=x, prior=prior, nitt=26000, thin=20, burnin=6000, pr=T)
mc.fx <- function(x, y, leaflet, mle.tab) {
  resp <- names(x)[y] # assign the column name of the 'y' index
  fixed.formula <- as.formula( paste0( get("resp"), " ~ 1 ") )
  mles <- mle.tab$vcov
  names(mles) <- mle.tab$grp
  if (leaflet) {
    random.formula <- as.formula(" ~ FinBIL + block + plant + leaflet.type ")
    prior <- list(G=list(G1=list(V=1, nu=1, alpha.mu=mles["FinBIL"], alpha.V=(25*sqrt(mles["FinBIL"]))^2), G2=list(V=1, nu=1, alpha.mu=mles["block"], alpha.V=(25*sqrt(mles["block"]))^2), G3=list(V=1, nu=1, alpha.mu=mles["plant"], alpha.V=(25*sqrt(mles["plant"]))^2), G4=list(V=1, nu=1, alpha.mu=mles["leaflet.type:plant"], alpha.V=(25*sqrt(mles["leaflet.type:plant"]))^2)), R=list(V=1, nu=0.002) )  
  } else {
    random.formula <- as.formula(" ~ FinBIL + block + plant ")
    prior <- list(G=list(G1=list(V=1, nu=1, alpha.mu=mles["FinBIL"], alpha.V=(25*sqrt(mles["FinBIL"]))^2), G2=list(V=1, nu=1, alpha.mu=mles["block"], alpha.V=(25*sqrt(mles["block"]))^2), G3=list(V=1, nu=1, alpha.mu=mles["plant"], alpha.V=(25*sqrt(mles["plant"]))^2)), R=list(V=1, nu=0.002) )  
  }
  fit1 <- MCMCglmm(fixed = fixed.formula, random = random.formula, data=x, prior=prior, nitt=120000, thin=100, burnin=20000, pr=T)
  fit2 <- MCMCglmm(fixed = fixed.formula, random = random.formula, data=x, prior=prior, nitt=120000, thin=100, burnin=20000, pr=T)
  diag <- gelman.diag(as.mcmc.list(list(fit1vcv=fit1$VCV, fit2vcv=fit2$VCV)), autoburnin=FALSE, multivariate=TRUE )
  mpsrf <- diag$mpsrf
  mc.list <- list(fit1=fit1, fit2=fit2, diag=diag, mpsrf=mpsrf)
  mc.list
}

# Function to take a pair of MCMC runs and calculate H^2, R^2, and stoch. idx. estimates
mc.resTab <- function(mc1, mc2, leaflet, mle.tab){
  ranefTab <- as.data.frame(rbind(mc1, mc2))
  ranefTab$h2 <- apply(ranefTab, 1, function(x) x["FinBIL"] / sum(x))
  if (leaflet) {
    ranefTab$r2 <- apply(ranefTab, 1, function(x) (x["FinBIL"] + x["plant"] + x["leaflet.type"]) / sum(x) )
    ranefTab$stoch <- apply(ranefTab, 1, function(x) (x["block"] + x["units"]) / (x["plant"] + x["leaflet.type"]) )
  } else {
    ranefTab$r2 <- apply(ranefTab, 1, function(x) (x["FinBIL"] + x["plant"]) / sum(x) )
    ranefTab$stoch <- apply(ranefTab, 1, function(x) (x["block"] + x["units"]) / x["plant"] )
  }
  res.fx <- function (x) {
    mean <- mean(x)
    median <- median(x)
    mode <- as.numeric(posterior.mode(as.mcmc(x)))
    CI <- HPDinterval(as.mcmc(x))
    res <- c(mean=mean, median=median, mode=mode, lowerCI=CI[1], upperCI=CI[2])
    res
  }
  h2.res <- res.fx(ranefTab$h2)
  r2.res <- res.fx(ranefTab$r2)
  stoch.res <- res.fx(ranefTab$stoch)
  h2.mle <- mle.tab[1]
  r2.mle <- mle.tab[2]
  stoch.mle <- mle.tab[3]
  results <- c(h2.mle, h2=h2.res, r2.mle, r2=r2.res, stoch.mle, stoch=stoch.res)
  names(results)[c(1,7,13)] <- paste0(names(results)[c(1,7,13)], ".mle")
  results
}
#----
# Load datasets and run analyses

# Load Leaf complexity dataset
setwd("/Users/Dani/UCD/BILs/leaf_traits/")
# load data tables
load("comp.rn.Rdata")
head(comp.rn)
levels(comp.rn$FinBIL)
# REMOVE M82 & PEN
comp.noparents <- comp.rn[comp.rn$FinBIL!="M82" & comp.rn$FinBIL!="PEN",]
comp.noparents <- droplevels(comp.noparents)
levels(comp.noparents$FinBIL)
head(comp.noparents)

# run loops
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
comp.resTab <- foreach(i=1:4, .combine="rbind", .options.multicore=mcoptions, .errorhandling="pass") %dopar% { # run loop
  y = i + 4
  paramMLEs(x=comp.noparents, y=y, randform="+ (1 | plant)" )
}
rownames(comp.resTab) <- names(comp.noparents)[5:8] # name the list elements by their original trait names
comp.resTab
write.csv(comp.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MLEs/comp.resTab.csv", quote=F)

comp.mles <- foreach(i=1:4, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 4
  mle.fx(x=comp.noparents, y=y, randform="+ (1 | plant)" )
}
names(comp.mles) <- names(comp.noparents)[5:8] # name the list elements by their original trait names

comp.mc <- foreach(i=1:4, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 4
  mle.tab=comp.mles[[i]]
  mc.fx(x=comp.noparents, y=y, leaflet=FALSE, mle.tab=mle.tab)
}
names(comp.mc) <- names(comp.noparents)[5:8]
lapply(comp.mc, function(x) x[[3]][[2]] ) # print out MPSRF's to check convergence
save(comp.mc, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/comp.mc.Rdata")

comp.mc.resTab <- foreach(i=1:4, .combine="rbind", .options.multicore=mcoptions) %dopar% { # run loop
  mc=comp.mc[[i]]
  mle.tab=comp.resTab[i,]
  mc.resTab(mc1 = mc[[1]][[3]], mc2 = mc[[2]][[3]], leaflet=FALSE, mle.tab=mle.tab)
}
rownames(comp.mc.resTab) <- names(comp.noparents)[5:8] # name the list elements by their original trait names
comp.mc.resTab
write.csv(comp.mc.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/comp.mc.resTab.csv", quote=F)

rm(comp.mc) # remove to free up memory
#-----
load("labels.Rdata")
circ <- read.delim("BIL.circAR.all.2011.txt")
circ <- merge(circ, labels, by="plant")
circ <- droplevels(circ)
names(circ)[2] <- "leaflet.type"
circ <- circ[circ$FinBIL!="M82" & circ$FinBIL!="PEN",]
circ <- droplevels(circ)

# ***CONSIDER*** standardizing (scale and mean centering) circ traits
head(circ)
circ[,3:7] <- scale(circ[,3:7])

# run loop
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
circ.resTab <- foreach(i=1:5, .combine="rbind", .options.multicore=mcoptions, .errorhandling="pass") %dopar% { # run loop
  y = i + 2
  paramMLEs(x=circ, y=y, randform="+ (1 | plant / leaflet.type)" )
}
rownames(circ.resTab) <- names(circ)[3:7] # name the list elements by their original trait names
circ.resTab
write.csv(circ.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MLEs/circ.resTab.csv", quote=F)

circ.mles <- foreach(i=1:5, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 2
  mle.fx(x=circ, y=y, randform="+ (1 | plant / leaflet.type)" )
}
names(circ.mles) <- names(circ)[3:7]

circ$leaflet.type <- paste0(circ$plant, circ$leaflet.type)
circ.mc <- foreach(i=1:5, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 2
  mle.tab=circ.mles[[i]]
  mc.fx(x=circ, y=y, leaflet=TRUE, mle.tab=mle.tab)
}
names(circ.mc) <- names(circ)[3:7]
lapply(circ.mc, function(x) x[[3]][[2]] ) # print out MPSRF's to check convergence
save(circ.mc, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/circ.mc.Rdata")

circ.mc.resTab <- foreach(i=1:5, .combine="rbind", .options.multicore=mcoptions) %dopar% { # run loop
  mc=circ.mc[[i]]
  mle.tab=circ.resTab[i,]
  mc.resTab(mc1 = mc[[1]][[3]], mc2 = mc[[2]][[3]], leaflet=TRUE, mle.tab=mle.tab)
}
rownames(circ.mc.resTab) <- names(circ)[3:7] # name the list elements by their original trait names
circ.mc.resTab
write.csv(circ.mc.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/circ.mc.resTab.csv", quote=F)

rm(circ.mc)
#----
sym <- read.delim("BIL_Spcascores.txt")
sym <- merge(sym, labels, by="plant")
sym <- droplevels(sym)
names(sym)[2] <- "leaflet.type"
sym <- sym[with(sym, FinBIL!="M82" & FinBIL!="PEN"),]
sym <- droplevels(sym)

sym[,3:11] <- scale(sym[,3:11]) # scale & center

# run lme4 loop
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
sym.resTab <- foreach(i=1:9, .combine="rbind", .options.multicore=mcoptions, .errorhandling="pass") %dopar% { # run loop
  y = i + 2
  paramMLEs(x=sym, y=y, randform="+ (1 | plant / leaflet.type)" )
}
rownames(sym.resTab) <- names(sym)[3:11] # name the list elements by their original trait names
sym.resTab
write.csv(sym.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MLEs/sym.resTab.csv", quote=F)

sym.mles <- foreach(i=1:9, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 2
  mle.fx(x=sym, y=y, randform="+ (1 | plant / leaflet.type)" )
}
names(sym.mles) <- names(sym)[3:11]

sym$leaflet.type <- paste0(sym$plant, sym$leaflet.type) # relabel leaf.type (by prepending the plantID) to enable nested ranefs in MCMCglmm
sym.mc <- foreach(i=1:9, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 2
  mle.tab=sym.mles[[i]]
  mc.fx(x=sym, y=y, leaflet=TRUE, mle.tab=mle.tab)
}
names(sym.mc) <- names(sym)[3:11] # name the list elements by their original trait names
lapply(sym.mc, function(x) x[[3]][[2]] ) # print out MPSRF's to check convergence
save(sym.mc, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC//Users/Dani/UCD/BILs/h2_et_al_MCMC/sym.mc.Rdata")

sym.mc.resTab <- foreach(i=1:9, .combine="rbind", .options.multicore=mcoptions) %dopar% { # run loop
  mc=sym.mc[[i]]
  mle.tab=sym.resTab[i,]
  mc.resTab(mc1 = mc[[1]][[3]], mc2 = mc[[2]][[3]], leaflet=TRUE, mle.tab=mle.tab)
}
rownames(sym.mc.resTab) <- names(sym)[3:11] # name the list elements by their original trait names
sym.mc.resTab
write.csv(sym.mc.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/sym.mc.resTab.csv", quote=F)

rm(sym.mc)
#----
asym <- read.delim("BIL_ASpcascores.txt")
asym <- merge(asym, labels, by="plant")
asym <- droplevels(asym)
names(asym)[2] <- "leaflet.type"
asym <- asym[asym$leaflet.type != "t",] # remove terminal leaflets, because they shouldn't be included in asymmetric shape analysis
asym <- droplevels(asym) # drop levels again to remove terminal leaflet type level
asym[3:9] <- abs(asym[3:9]) # take the absolute value of the PCs so that the asymmetry traits are enantiomer-independent (i.e. left and right lateral leaflets are treated equally)
asym <- asym[with(asym, FinBIL!="M82" & FinBIL!="PEN"),]
asym <- droplevels(asym)

asym[,3:9] <- scale(asym[,3:9]) # scale & center

# run loop
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
asym.resTab <- foreach(i=1:7, .combine="rbind", .options.multicore=mcoptions, .errorhandling="pass") %dopar% { # run loop
  y = i + 2
  paramMLEs(x=asym, y=y, randform="+ (1 | plant)" )
}
rownames(asym.resTab) <- names(asym)[3:9] # name the list elements by their original trait names
asym.resTab
write.csv(asym.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/asym.resTab.csv", quote=F)

asym.mles <- foreach(i=1:7, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 2
  mle.fx(x=asym, y=y, randform="+ (1 | plant)" )
}
names(asym.mles) <- names(asym)[3:9] # name the list elements by their original trait names

asym.mc <- foreach(i=1:7, .options.multicore=mcoptions) %dopar% { # run loop
  y = i + 2
  mle.tab=asym.mles[[i]]
  mc.fx(x=asym, y=y, leaflet=FALSE, mle.tab=mle.tab)
}
names(asym.mc) <- names(asym)[3:9] # name the list elements by their original trait names
save(asym.mc, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/asym.mc.Rdata")
lapply(asym.mc, function(x) x[[3]][[2]] ) # print out MPSRF's to check convergence

asym.mc.resTab <- foreach(i=1:7, .combine="rbind", .options.multicore=mcoptions) %dopar% { # run loop
  mc=asym.mc[[i]]
  mle.tab=asym.resTab[i,]
  mc.resTab(mc1 = mc[[1]][[3]], mc2 = mc[[2]][[3]], leaflet=FALSE, mle.tab=mle.tab)
}
rownames(asym.mc.resTab) <- names(asym)[3:9] # name the list elements by their original trait names
asym.mc.resTab
write.csv(asym.mc.resTab, file="/Users/Dani/UCD/BILs/h2_et_al_MCMC/asym.mc.resTab.csv", quote=F)

rm(asym.mc)
#-----
# deprecated code
# head(ranefTab)
# 
# c.prior <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=25^2), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=25^2)), R=list(V=1, nu=0.002) )
# mc.pe.aPC1c <- MCMCglmm(fixed = PC1 ~ 1, random = ~ FinBIL + plant, data=asym, prior=c.prior, nitt=39000, thin=30, burnin=9000, pr=T)
# summary(mc.pe.aPC1c)
# plot(mc.pe.aPC1c)
# 
# d.prior <- list(G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=25^2), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=25^2)), R=list(V=1, nu=0.002) )
# mc.pe.aPC1d <- MCMCglmm(fixed = PC1 ~ 1, random = ~ FinBIL + block + plant, data=asym, prior=d.prior, nitt=65000, thin=50, burnin=15000, pr=T)
# summary(mc.pe.aPC1d)
# plot(mc.pe.aPC1d)
# 
# mc.pe.aPC1f <- MCMCglmm(fixed = PC1 ~ 1, random = ~ FinBIL + block + plant, data=asym, prior=e.prior, nitt=65000, thin=50, burnin=15000, pr=T)

# empBay.resTab <- paramTab(mc.pe.aPC1e$VCV, mc.pe.aPC1e$VCV)
# head(empBay.resTab)
# 
# h2.res <- res.fx(empBay.resTab$h2)
# r2.res <- res.fx(empBay.resTab$r2)
# stoch.res <- res.fx(empBay.resTab$stoch)
# hist(log10(empBay.resTab$stoch))

# Pseudo-code:
# 1) Write function to store lme4 ranef tables in a list
# 2) Write function to run the MCMC analysis on a whole dataset, save the 2 model fits, and save the gelman.diag() result for each trait
#       Thus, use a list of lists with 3 elements to the trait-specific sublists
#       If all gelman.diag are okay, then proceed to step 2
# 3) Write another function that uses paramTab and res.fx internally to iterate through the above lists and generate
#    tables of results
# ** REMEMBER to relabel leaflet.type w/ unique labels in order to do nesting in a simple manner **

# aPC1 <- lmer(PC1 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
# mles <- as.data.frame(VarCorr(aPC1))$vcov; names(mles) <- as.data.frame(VarCorr(aPC1))$grp
# e.prior <- list(G=list(G1=list(V=1, nu=2, alpha.mu=mles["FinBIL"], alpha.V=25^2), G2=list(V=1, nu=2, alpha.mu=mles["block"], alpha.V=25^2), G3=list(V=1, nu=2, alpha.mu=mles["plant"], alpha.V=25^2)), R=list(V=1, nu=0.002) )
# mc.pe.aPC1e <- MCMCglmm(fixed = PC1 ~ 1, random = ~ FinBIL + block + plant, data=asym, prior=e.prior, nitt=65000, thin=50, burnin=15000, pr=T)

# maybe make it a "round" 120K w/ 20K burn-in and thin=50

# x = circ
# y = 3
# leaflet = TRUE
# mle.tab = circ.mles[[1]]
# rm(list=c("x", "y", "leaflet", "mle.tab", "resp", "fixed.formula", "mles", "random.formula", "prior", "fit1"))
# lm.area <- lmer(Area ~ 1+(1|FinBIL)+(1|plant/leaflet.type)+(1|block), data=x, REML=TRUE)
# summary(lm.area)