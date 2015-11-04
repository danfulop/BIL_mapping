library(reshape2)
library(sparsenet)
library(ggplot2)
library(foreach)
library(doParallel)
library(plyr)

setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/danielfulop/leaf_traits/")

# Load data, recode genotypes, etc
#------
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
#----
load("bin.stats.Rdata")

# Function to recode genotype information into 0s, 1s, and 2s
#-----
recode <- function(x) {
  if (x == "M82") {
    x <- 0
  } else if (x == "PEN") {
    x <- 2
  } else if (x == "HET") {
    x <- 1
  }
}
#-----
genotab.recoded <- genotab
genotab.recoded[1:(ncol(genotab.recoded)-2)] <- apply(genotab.recoded[1:(ncol(genotab.recoded)-2)], c(1,2), recode)
genotab.recoded <- droplevels(genotab.recoded)
#save(genotab.recoded, file="genotab.recoded.Rdata")
#-------
load("genotab.recoded.Rdata")

# Function to fit SparseNet regression -- mapping function
#----
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
  tmp2 <- vector("list", length=10) # temp list for storing min parameters
  for(j in 1:10) { # redo cross-validation 10 times to get more "stable" tuning parameters
    cv.sp <- cv.sparsenet(x=geno.mat, y=response, lambda=lambda.manual, ngamma=36, nfolds=6, warm="both") # cross-validation
    tmp[[j]] <- cv.sp$parms.1se # CV params. for 1 std. dev. away from min.
    tmp2[[j]] <- cv.sp$parms.min # CV params.
  }
  mean.gamma <- mean(unlist(tmp)[seq(1,19,2)]) # mean gamma
  mean.lambda <- mean(unlist(tmp)[seq(2,20,2)]) # mean lambda
  min.cv.lambda <- mean(unlist(tmp2)[seq(2,20,2)]) # min cv lambda  
  cv.lambda.seq <- cv.sp$sparsenet.fit$lambda
  start.lambda.seq <- cv.lambda.seq[which(cv.lambda.seq[] > min.cv.lambda)[length(which(cv.lambda.seq[] > min.cv.lambda))] ] # store lambda value just beyond cv.lambda.min
  # Fit sparsenet one more time with the mean tuning parameters
  lambda.seq <- exp(seq(from=log(start.lambda.seq), to=log(mean.lambda), length=10)) # lambda sequence ending in mean.lambda
  sp.fit <- sparsenet(x=geno.mat, y=response, lambda=lambda.seq, min.gamma=mean.gamma, ngamma=9, warm="both") # full dataset sparsenet fit
  coefs <- sp.fit$coefficients$g9$beta[,1] # save preferred set of coefficients
  coefs <- cbind(bin.stats, coefs) # combine with bin information
  non.zero.coefs <- coefs[coefs$coefs!=0,] # non-zero coefficients
  n.coef <- nrow(non.zero.coefs) # number of non-zero coefficients **this prob. fails b/c it should be nrow, b/c now it's a data.frame
  results <- list(coefs=coefs, non.zero.coefs=non.zero.coefs, n.coef=n.coef, gamma=mean.gamma, lambda=mean.lambda, start.lambda.seq=start.lambda.seq, sp.fit=sp.fit)
  results
}
#----

# Fit SparseNet regressions
#-------
lambda.manual = exp(seq(from = -5.5, to = -1.5, length = 100))
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
registerDoParallel(cores=5) # register parallel backend
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
registerDoParallel(cores=9) # register parallel backend
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
registerDoParallel(cores=7) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
system.time(asym.map <- foreach(i=1:length(asym.pred), .options.multicore=mcoptions) %dopar% { # run loop
  asym.map[i] <- map.fx(datl=asym.pred, ln=i, gt.tab=genotab.recoded, bin.stats=bin.stats, lambda.manual=lambda.manual)
})
names(asym.map) <- names(asym.pred)
save(asym.map, file="asym.map.Rdata")
#------

# Construct genotypic matrix for epistatic searches
#----
geno.mat <- as.matrix(genotab.recoded[,1:(ncol(genotab.recoded)-2)])
tmp <- mclapply(combn(ncol(geno.mat), 2, simplify=FALSE), function(i) (geno.mat[, i[1]] * geno.mat[, i[2]])/2, mc.cores=4) # multiply BIN vectors
geno.epi.mat <- do.call(cbind, tmp) # cbind to make them into a matrix
tmp.names <- mclapply(combn(length(colnames(geno.mat)), 2, simplify=FALSE), function(i) paste0(colnames(geno.mat)[i][1], "_x_", colnames(geno.mat)[i][2]), mc.cores=4) # multiply BIN vectors
epi.names <- do.call(c, tmp.names)
colnames(geno.epi.mat) <- epi.names
ncol(geno.epi.mat) # 549676
geno.epi.mat <- geno.epi.mat[, colSums(geno.epi.mat) != 0] # Remove epistatic "bins" with all zeroes
ncol(geno.epi.mat) # 306291 => 243385 all zero columns
geno.mat <- cbind(geno.mat, geno.epi.mat) # append epistatic "bin" matrix to original geno.mat
epitab <- as.data.frame(geno.mat) # convert to data.frame and append BIL + FinBIL columns
epitab <- cbind(epitab,  genotab.recoded[,c(ncol(genotab.recoded)-1,ncol(genotab.recoded))])
# save(epitab, file="epitab.Rdata")
#-----
load("epitab.Rdata")

# Construct "epi.bin.stats"
#-----
# make a bin-ID column
tmp.bin <- mclapply(combn(length(bin.stats$bin), 2, simplify=FALSE), function(i) paste0(bin.stats$bin[i][1], "_x_", bin.stats$bin[i][2]), mc.cores=4)
epi.bin <- do.call(c, tmp.bin)
# make a chr column, similar syntax to epi.names, but trimming SL2.40. Use "//" delimiter
tmp.chr <- mclapply(combn(length(bin.stats$chr), 2, simplify=FALSE), function(i) paste0(substr(bin.stats$chr[i][1], 7, 10), "//", substr(bin.stats$chr[i][2], 7, 10)), mc.cores=4)
epi.chr <- do.call(c, tmp.chr)
# join the bin info columns with "//" delimiter
tmp.bin.mid <- mclapply(combn(length(bin.stats$bin.mid), 2, simplify=FALSE), function(i) paste0( str_trim(bin.stats$bin.mid[i][1]), "//", str_trim(bin.stats$bin.mid[i][2])), mc.cores=4)
epi.bin.mid <- do.call(c, tmp.bin.mid)
tmp.bin.start <- mclapply(combn(length(bin.stats$bin.start), 2, simplify=FALSE), function(i) paste0( str_trim(bin.stats$bin.start[i][1]), "//", str_trim(bin.stats$bin.start[i][2])), mc.cores=4)
epi.bin.start <- do.call(c, tmp.bin.start)
tmp.bin.end <- mclapply(combn(length(bin.stats$bin.end), 2, simplify=FALSE), function(i) paste0( str_trim(bin.stats$bin.end[i][1]), "//", str_trim(bin.stats$bin.end[i][2])), mc.cores=4)
epi.bin.end <- do.call(c, tmp.bin.end)
epi.bin.mat <- cbind(epi.bin, epi.chr, epi.bin.mid, epi.bin.start, epi.bin.end) # cbind above 5 columns
epi.bin.mat <- epi.bin.mat[epi.bin.mat[,1] %in% intersect(epi.bin.mat[,1], colnames(geno.epi.mat)),] # trim rows to match set of colnames(geno.epi.mat)
colnames(epi.bin.mat) <- substr(colnames(epi.bin.mat), 5, 13)
epi.bin.stats <- rbind(bin.stats, epi.bin.mat, deparse.level = 0) # rbind these joined columns to bin.stats
save(epi.bin.stats, file="epi.bin.stats.Rdata")
#-----
load("epi.bin.stats.Rdata")

# Function to check for the lambda range
#------
plot.fx <- function(datl, ln, gt.tab, bin.stats, lambda.manual, dat.name) {
  trait.name <- names(datl)[ln] # (list element) name for focal trait
  trait.dat <- datl[[get("trait.name")]] # Assign focal trait data.frame
  trait.dat <- merge(gt.tab, trait.dat, by.x="BIL", by.y="genotype") # merge genotype and trait data
  trait.dat <- trait.dat[-which(names(trait.dat)=="FinBIL.y")] # remove redundant "FinBIL.y"
  names(trait.dat)[which(names(trait.dat)=="FinBIL.x")] <- "FinBIL" # rename FinBIL.x into FinBIL
  trait.dat <- droplevels(trait.dat)
  geno.mat <- as.matrix(trait.dat[2:(ncol(trait.dat)-5)], rownames.force=F) # genotype matrix, i.e. Xs
  response <- as.vector(as.matrix(trait.dat['predPlusResid'], rownames.force=F)) # response vector for sparsenet, i.e. Ys
  cv.sp <- cv.sparsenet(x=geno.mat, y=response, ngamma=9, lambda=lambda.manual, nfolds=6, warm="both")
  plot.path <- paste0(getwd(), "/epiCVplots/", dat.name, ".", trait.name, ".CVplot", ".pdf")
  pdf(file=plot.path)
  plot(cv.sp)
  dev.off()
}
#------

# Check lambda ranges
#-------
lambda.manual = exp(seq(from = -5.5, to = -1.5, length = 20))
# Complexity data 
load("comp.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
nofun <- function(a,b) NULL
system.time(foreach(i=1:length(comp.pred), .options.multicore=mcoptions, .combine='nofun') %dopar% { # run loop
  plot.fx(datl=comp.pred, ln=i, gt.tab=epitab, bin.stats=bin.stats, lambda.manual=lambda.manual, dat.name="comp")
})

# Circ. data
load("circ.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
registerDoParallel(cores=5) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
nofun <- function(a,b) NULL
system.time(foreach(i=1:length(circ.pred), .options.multicore=mcoptions, .combine='nofun') %dopar% { # run loop
  plot.fx(datl=circ.pred, ln=i, gt.tab=epitab, bin.stats=bin.stats, lambda.manual=lambda.manual, dat.name="circ")
})

# Sym PCs
load("sym.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
registerDoParallel(cores=9) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
nofun <- function(a,b) NULL
system.time(foreach(i=1:length(sym.pred), .options.multicore=mcoptions, .combine='nofun') %dopar% { # run loop
  plot.fx(datl=sym.pred, ln=i, gt.tab=epitab, bin.stats=bin.stats, lambda.manual=lambda.manual, dat.name="sym")
})

# Asym PCs
load("asym.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
registerDoParallel(cores=7) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
nofun <- function(a,b) NULL
system.time(foreach(i=1:length(asym.pred), .options.multicore=mcoptions, .combine='nofun') %dopar% { # run loop
  plot.fx(datl=asym.pred, ln=i, gt.tab=epitab, bin.stats=bin.stats, lambda.manual=lambda.manual, dat.name="asym")
})
#------

# Function to fit SparseNet regression **with epistasis** -- mapping function
#----
epi.map.fx <- function(datl, ln, gt.tab, bin.stats, lambda.manual) {
  trait.name <- names(datl)[ln] # (list element) name for focal trait
  trait.dat <- datl[[get("trait.name")]] # Assign focal trait data.frame
  trait.dat <- merge(gt.tab, trait.dat, by.x="BIL", by.y="genotype") # merge genotype and trait data
  trait.dat <- trait.dat[-which(names(trait.dat)=="FinBIL.y")] # remove redundant "FinBIL.y"
  names(trait.dat)[which(names(trait.dat)=="FinBIL.x")] <- "FinBIL" # rename FinBIL.x into FinBIL
  trait.dat <- droplevels(trait.dat)
  geno.mat <- as.matrix(trait.dat[2:(ncol(trait.dat)-5)], rownames.force=F) # genotype matrix, i.e. Xs
  response <- as.vector(as.matrix(trait.dat['predPlusResid'], rownames.force=F)) # response vector for sparsenet, i.e. Ys
  tmp <- vector("list", length=10) # temp list for storing 1se parameters
  tmp2 <- vector("list", length=10) # temp list for storing min parameters
  registerDoParallel(cores=10) # register parallel backend
  mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options
  nofun <- function(a,b) NULL
  foreach(j=1:10, .options.multicore=mcoptions, .combine='nofun') %dopar% { # redo cross-validation 10 times to get more "stable" tuning parameters
    cv.sp <- cv.sparsenet(x=geno.mat, y=response, lambda=lambda.manual, ngamma=18, nfolds=6, warm="both")
    tmp[[j]] <- cv.sp$parms.1se # 1 std. dev. error params.
    tmp2[[j]] <- cv.sp$parms.min # min CV error params.
  }
  mean.gamma <- mean(unlist(tmp)[seq(1,19,2)]) # mean gamma
  mean.lambda <- mean(unlist(tmp)[seq(2,20,2)]) # mean lambda
  min.cv.lambda <- mean(unlist(tmp2)[seq(2,20,2)]) # min cv lambda  
  cv.lambda.seq <- cv.sp$sparsenet.fit$lambda
  start.lambda.seq <- cv.lambda.seq[which(cv.lambda.seq[] > min.cv.lambda)[length(which(cv.lambda.seq[] > min.cv.lambda))] ] # store lambda value just beyond cv.lambda.min
  # Fit sparsenet one more time with the mean tuning parameters
  lambda.seq <- exp(seq(from=log(start.lambda.seq), to=log(mean.lambda), length=10)) # lambda sequence ending in mean.lambda
  sp.fit <- sparsenet(x=geno.mat, y=response, lambda=lambda.seq, min.gamma=mean.gamma, ngamma=9, warm="both") # full dataset sparsenet fit
  coefs <- sp.fit$coefficients$g9$beta[,1] # save preferred set of coefficients
  coefs <- cbind(bin.stats, coefs) # combine with bin information
  non.zero.coefs <- coefs[coefs$coefs!=0,] # non-zero coefficients
  n.coef <- nrow(non.zero.coefs) # number of non-zero coefficients **this prob. fails b/c it should be nrow, b/c now it's a data.frame
  results <- list(coefs=coefs, non.zero.coefs=non.zero.coefs, n.coef=n.coef, gamma=mean.gamma, lambda=mean.lambda, start.lambda.seq=start.lambda.seq, sp.fit=sp.fit)
  results
}
# if I were to use the results from a single CV instead, use code below
# which.parms <- cv.sp$which.1se
# coefs <- cv.sp$sparsenet.fit$coefficients[[ which.parms[2] ]]$beta[, which.parms[1] ] # access to the 1.se coefficients
#----

# Lambda ranges
#----
comp.lambda <- list(comp.pri = c(-4, -2), comp.int = c(-4.2, -2), comp.sec = c(-4, -1.5), comp.all = c(-4.2, -2) )
circ.lambda <- list(circ.Area = c(-3.5, -1.5), circ.Circ = c(-4, -2), circ.AR = c(-4.5, -2), circ.Round = c(-4.8, -2.3), circ.Solidity = c(-4.5, -2.5) )
sym.lambda <- list(sym.PC1 = c(-4, -2), sym.PC2 = c(-4.5, -2.5), sym.PC3 = c(-4.5, -2.5), sym.PC4 = c(-4, -2), sym.PC5 = c(-4, -2.5), 
                   sym.PC6 = c(-3.7, -2.2), sym.PC7 = c(-3.8, -1.8), sym.PC8 = c(-3.7, -2), sym.PC9 = c(-3.5, -1.5) )
asym.lambda <- list(asym.PC1 = c(-3.5, -1.5), asym.PC2 = c(-3.5, -1.5), asym.PC3 = c(-3.5, -1.5), asym.PC4 = c(-3.5, -1.5), asym.PC5 = c(-3.5, -1.5),
                    asym.PC6 = c(-3.5, -1.5), asym.PC7 = c(-3.5, -1.5) )
#----

# Fit epistatic SparseNet regressions
#-------
registerDoParallel(cores=5) # register parallel backend
mcoptions <- list(preschedule=FALSE, set.seed=FALSE) # multi-core options
# Complexity data 
load("comp.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
comp.epi.map <- vector("list", length=length(comp.pred)) # Initialize results list
comp.epi.map <- foreach(i=1:length(comp.pred), .options.multicore=mcoptions) %dopar% { # run loop
  lambda.manual <- exp(seq(from=comp.lambda[[i]][1], to=comp.lambda[[i]][2], length=(abs(comp.lambda[[i]][1] - comp.lambda[[i]][2]))*10 ))
  comp.epi.map[i] <- epi.map.fx(datl=comp.pred, ln=i, gt.tab=epitab, bin.stats=epi.bin.stats, lambda.manual=lambda.manual)
}
names(comp.epi.map) <- names(comp.pred)
save(comp.epi.map, file="comp.epi.map.Rdata")

# Circ. data
load("circ.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
circ.epi.map <- vector("list", length=length(circ.pred)) # Initialize results list
circ.epi.map <- foreach(i=1:length(circ.pred), .options.multicore=mcoptions) %dopar% { # run loop
  lambda.manual <- exp(seq(from=circ.lambda[[i]][1], to=circ.lambda[[i]][2], length=(abs(circ.lambda[[i]][1] - circ.lambda[[i]][2]))*10 ))
  circ.epi.map[i] <- epi.map.fx(datl=circ.pred, ln=i, gt.tab=epitab, bin.stats=epi.bin.stats, lambda.manual=lambda.manual)
}
names(circ.epi.map) <- names(circ.pred)
save(circ.epi.map, file="circ.epi.map.Rdata")

# Sym PCs
load("sym.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
sym.epi.map <- vector("list", length=length(sym.pred)) # Initialize results list
sym.epi.map <- foreach(i=1:length(sym.pred), .options.multicore=mcoptions) %dopar% { # run loop
  lambda.manual <- exp(seq(from=sym.lambda[[i]][1], to=sym.lambda[[i]][2], length=(abs(sym.lambda[[i]][1] - sym.lambda[[i]][2]))*10 ))
  sym.epi.map[i] <- epi.map.fx(datl=sym.pred, ln=i, gt.tab=epitab, bin.stats=epi.bin.stats, lambda.manual=lambda.manual)
}
names(sym.epi.map) <- names(sym.pred)
save(sym.epi.map, file="sym.epi.map.Rdata")

# Asym PCs
load("asym.pred.Rdata") # load LMM predicted response values + residuals (1 value per plant)
asym.epi.map <- vector("list", length=length(asym.pred)) # Initialize results list
asym.epi.map <- foreach(i=1:length(asym.pred), .options.multicore=mcoptions) %dopar% { # run loop
  lambda.manual <- exp(seq(from=asym.lambda[[i]][1], to=asym.lambda[[i]][2], length=(abs(asym.lambda[[i]][1] - asym.lambda[[i]][2]))*10 ))
  asym.epi.map[i] <- epi.map.fx(datl=asym.pred, ln=i, gt.tab=epitab, bin.stats=epi.bin.stats, lambda.manual=lambda.manual)
}
names(asym.epi.map) <- names(asym.pred)
save(asym.epi.map, file="asym.epi.map.Rdata")
#------
