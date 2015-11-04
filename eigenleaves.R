# Script to generate EFD-PC's QTL images for BIL paper
# devtools::install_github("vbonhomme/Momocs")

library(Momocs)
library(ggplot2)
library(ade4)
library(stringr)

setwd("~/UCD/BILs/EFD_figures/")
outlines <- nef2Coe("neffile") # matrix of harmonics ...still need Coe object
dim(outlines)
outlines[1:6,1:10]
# colnames(outlines)
# class(outlines) # matrix

fac <- rownames(outlines)
fac <- data.frame(leaflet=fac, plant=substr(fac,1,4), block=substr(fac,1,1), type=substr(fac,5,5), plant.n.type=substr(fac,1,5))
head(fac)

bil.coe <- OutCoe(coe=outlines, fac=fac, method="efourier", norm=T) # convert from matrix to OutCoe object
# take mean of leaflet types per leaf, and then take mean of leaflets per plant
mshape.by.type <- mshapes(bil.coe, fac="plant.n.type")
names(mshape.by.type) # [1] "Coe" "shp"
dim(mshape.by.type$Coe) # [1] 4571   80
mshape.by.type$Coe$coe[1:6,1:6]
mshape.by.plant <- mshapes(mshape.by.type$Coe, fac="plant")
dim(mshape.by.plant$Coe) # [1] 1523   80
mshape.by.plant$Coe$coe[1:6,1:6]
str(mshape.by.plant$Coe)

mean.coe.mat <- as.matrix(mshape.by.plant$Coe$coe)
dim(mean.coe.mat)
head(mean.coe.mat)
dim(mean.coe.mat)
rownames(mean.coe.mat)

load("/Users/Dani/UCD/BILs/leaf_traits/labels.Rdata") # BIL line IDs + field IDs
dim(labels)
head(labels)
length(unique(labels$plant))
# reorder labels by plant ID
labels <- with(labels, labels[order(plant),])
head(labels)

setdiff(rownames(mean.coe.mat), labels$plant)
length(setdiff(rownames(mean.coe.mat), labels$plant)) # 270
setdiff(labels$plant, rownames(mean.coe.mat))
length(setdiff(labels$plant, rownames(mean.coe.mat))) # 94

# length(labels$plant %in% rownames(mean.coe.mat))
# length(rownames(mean.coe.mat) %in% labels$plant)
sub.mcoe <- subset(mean.coe.mat, subset = rownames(mean.coe.mat) %in% labels$plant) # subset mean Coe data my plants in labels
dim(sub.mcoe) # [1] 1253   80
head(sub.mcoe)

sub.labs <- subset(labels, subset = labels$plant %in% rownames(sub.mcoe))
dim(sub.labs)
head(sub.labs)

# Get means by BIL lines
sub.mcoe <- as.data.frame(sub.mcoe)
sub.mcoe$plant <- rownames(sub.mcoe)
mcoe.n.labs <- merge(sub.mcoe, sub.labs[c(1,3)], by="plant")
dim(mcoe.n.labs) # 1] 1253   82
head(mcoe.n.labs)
mshape.by.plant2 <- OutCoe(coe=as.matrix(mcoe.n.labs[2:81]), fac=mcoe.n.labs[c(1,82)], method="efourier", norm=T) # convert from matrix to OutCoe object, this time with BIL IDs
mshape.by.bil <- mshapes(mshape.by.plant2, fac="genotype") # get mean shapes by BIL

load("/Users/Dani/UCD/BILs/leaf_traits/genotab.Rdata")
dim(genotab)
genotab[1:6, 1:6]

setdiff(sub.labs$genotype, rownames(genotab)) # "PEN"
# sub.labs$genotype %in% rownames(genotab)
sub.labs <- sub.labs[sub.labs$genotype!="PEN",] # remove PEN
head(sub.labs)

load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/sym.map.Rdata")
names(sym.map$PC1)
sym.map$PC1$non.zero.coefs

pc1gen <- genotab[,sym.map$PC1$non.zero.coefs$bin]
dim(pc1gen) # [1] 440  18
pc1gen[1:6,1:6]
pc1gen$genotype <- rownames(pc1gen)
# pc1gen <- merge(pc1gen, sub.labs[c(1,3)], by="genotype")
head(pc1gen)
bil.mat <- mshape.by.bil$Coe$coe
dim(bil.mat)
head(bil.mat)
bil.mat <- as.data.frame(bil.mat)
bil.mat$genotype <- rownames(bil.mat)
bil.mat.n.pc1qtl <- merge(bil.mat, pc1gen, by="genotype")
head(bil.mat.n.pc1qtl)
dim(bil.mat.n.pc1qtl)
pc1.coe <- OutCoe(coe=as.matrix(bil.mat.n.pc1qtl[2:81]), fac=bil.mat.n.pc1qtl[82:99], method="efourier", norm=T)

# Function to plot EFD outlines for a given QTL
ef.qtl.plot <- function(pc.coe) {
  for(i in 1:ncol(pc.coe$fac)) {
    bin.name <- colnames(pc.coe$fac)[i]
    m.shp.coe <- mshapes(pc.coe, fac=bin.name)
    m82.shp <- as.data.frame(m.shp.coe$shp$M82)
    m82.shp <- rbind(m82.shp, m82.shp[1,])
    pen.shp <- as.data.frame(m.shp.coe$shp$PEN)
    pen.shp <- rbind(pen.shp, pen.shp[1,])
    m.shp.plot <- ggplot(aes(x=x, y=y), data=m82.shp) + geom_path(color="magenta", size=2) + geom_path(aes(x=x, y=y), color="green", size=2, data=pen.shp) + 
      coord_fixed() + labs(title=bin.name, x="length", y="width") + xlim(-1,1) + ylim(-0.8, 0.8)
    ggsave(filename = paste0(bin.name, ".pdf"), plot=m.shp.plot, width=7.5, height=6)
  }
}

setwd("~/UCD/BILs/EFD_figures/PC1")
ef.qtl.plot(pc1.coe)

setwd("~/UCD/BILs/EFD_figures/PC2")
pc2gen <- genotab[,sym.map$PC2$non.zero.coefs$bin]
pc2gen$genotype <- rownames(pc2gen)
bil.mat.n.pc2qtl <- merge(bil.mat, pc2gen, by="genotype")
pc2.coe <- OutCoe(coe=as.matrix(bil.mat.n.pc2qtl[2:81]), fac=bil.mat.n.pc2qtl[82:ncol(bil.mat.n.pc2qtl)], method="efourier", norm=T)
ef.qtl.plot(pc2.coe)

setwd("~/UCD/BILs/EFD_figures/PC3")
pc3gen <- genotab[,sym.map$PC3$non.zero.coefs$bin]
pc3gen$genotype <- rownames(pc3gen)
bil.mat.n.pc3qtl <- merge(bil.mat, pc3gen, by="genotype")
pc3.coe <- OutCoe(coe=as.matrix(bil.mat.n.pc3qtl[2:81]), fac=bil.mat.n.pc3qtl[82:ncol(bil.mat.n.pc3qtl)], method="efourier", norm=T)
ef.qtl.plot(pc3.coe)

setwd("~/UCD/BILs/EFD_figures/PC4")
pc4gen <- genotab[,sym.map$PC4$non.zero.coefs$bin]
pc4gen$genotype <- rownames(pc4gen)
bil.mat.n.pc4qtl <- merge(bil.mat, pc4gen, by="genotype")
pc4.coe <- OutCoe(coe=as.matrix(bil.mat.n.pc4qtl[2:81]), fac=bil.mat.n.pc4qtl[82:ncol(bil.mat.n.pc4qtl)], method="efourier", norm=T)
ef.qtl.plot(pc4.coe)

load("/Users/Dani/UCD/BILs/final_additive_sparsenet_results/circ.map.Rdata")
names(circ.map) # "Area"     "Circ."    "AR"       "Round"    "Solidity"
# SKIP area because it has not QTL
setwd("~/UCD/BILs/EFD_figures/circ")
circ.gen <- genotab[,circ.map$Circ.$non.zero.coefs$bin]
circ.gen$genotype <- rownames(circ.gen)
circ.qtl <- merge(bil.mat, circ.gen, by="genotype")
circ.coe <- OutCoe(coe=as.matrix(circ.qtl[2:81]), fac=circ.qtl[82:ncol(circ.qtl)], method="efourier", norm=T)
ef.qtl.plot(circ.coe)

setwd("~/UCD/BILs/EFD_figures/AR")
AR.gen <- genotab[,circ.map$AR$non.zero.coefs$bin]
AR.gen$genotype <- rownames(AR.gen)
AR.qtl <- merge(bil.mat, AR.gen, by="genotype")
AR.coe <- OutCoe(coe=as.matrix(AR.qtl[2:81]), fac=AR.qtl[82:ncol(AR.qtl)], method="efourier", norm=T)
ef.qtl.plot(AR.coe)

setwd("~/UCD/BILs/EFD_figures/roundness")
Round.gen <- genotab[,circ.map$Round$non.zero.coefs$bin]
Round.gen$genotype <- rownames(Round.gen)
Round.qtl <- merge(bil.mat, Round.gen, by="genotype")
Round.coe <- OutCoe(coe=as.matrix(Round.qtl[2:81]), fac=Round.qtl[82:ncol(Round.qtl)], method="efourier", norm=T)
ef.qtl.plot(Round.coe)

setwd("~/UCD/BILs/EFD_figures/solidity")
Solidity.gen <- genotab[,circ.map$Solidity$non.zero.coefs$bin]
Solidity.gen$genotype <- rownames(Solidity.gen)
Solidity.qtl <- merge(bil.mat, Solidity.gen, by="genotype")
Solidity.coe <- OutCoe(coe=as.matrix(Solidity.qtl[2:81]), fac=Solidity.qtl[82:ncol(Solidity.qtl)], method="efourier", norm=T)
ef.qtl.plot(Solidity.coe)

# Plot scatter plots of: PC3 and Aspect Ratio, and PC3 and Solidity
# 491 & 533
setwd("~/UCD/BILs/EFD_figures")
mshp491 <- as.data.frame(mshape.by.bil$shp$BIL_491)
mshp491 <- rbind(mshp491, mshp491[1,])
mshp491.plot <- ggplot(aes(x=x, y=y), data=mshp491) + geom_path(color="black", size=2) +
  coord_fixed() + labs(x="", y="") + xlim(-1,1) + ylim(-0.8, 0.8) + theme_classic() + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(filename = paste0("mshp491", ".pdf"), mshp491.plot, width=7.5, height=6)

mshp533 <- as.data.frame(mshape.by.bil$shp$BIL_533)
mshp533 <- rbind(mshp533, mshp533[1,])
mshp533.plot <- ggplot(aes(x=x, y=y), data=mshp533) + geom_path(color="black", size=2) +
  coord_fixed() + labs(x="", y="") + xlim(-1,1) + ylim(-0.8, 0.8) + theme_classic() + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(filename = paste0("mshp533", ".pdf"), mshp533.plot, width=7.5, height=6)

mshp461 <- as.data.frame(mshape.by.bil$shp$BIL_461)
mshp461 <- rbind(mshp461, mshp461[1,])
mshp461.plot <- ggplot(aes(x=x, y=y), data=mshp461) + geom_path(color="black", size=2) +
  coord_fixed() + labs(x="", y="") + xlim(-1,1) + ylim(-0.8, 0.8) + theme_classic() + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(filename = paste0("mshp461", ".pdf"), mshp461.plot, width=7.5, height=6)

mshp484 <- as.data.frame(mshape.by.bil$shp$BIL_484)
mshp484 <- rbind(mshp484, mshp484[1,])
mshp484.plot <- ggplot(aes(x=x, y=y), data=mshp484) + geom_path(color="black", size=2) +
  coord_fixed() + labs(x="", y="") + xlim(-1,1) + ylim(-0.8, 0.8) + theme_classic() + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(filename = paste0("mshp484", ".pdf"), mshp484.plot, width=7.5, height=6)
