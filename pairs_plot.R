library(GGally)
library(ggplot2)
library(gpairs)
library(scales)
library(stringr)

setwd("/Users/Dani/UCD/BILs/model-fitted-means/")
load("comp.list.Rdata")
load("circ.list.Rdata")
load("sym.list.Rdata")
load("asym.list.Rdata")
load("FT.means.Rdata")
FT.list <- list(FT = FT.means)

join.means <- function(mlist, dat.name) {
  for(i in 1:length(mlist)) {
    trait.name <- names(mlist)[i]
    if(i == 1) {
      df <- mlist[[i]]$fittedMeans[, c('geno', 'Estimate')] 
    } else {
      df <- merge(df, mlist[[i]]$fittedMeans[, c('geno', 'Estimate')], by="geno")
    }
    names(df)[ncol(df)] <- paste0(dat.name, ".", trait.name)
  }
  df
} 

comp.means <- join.means(comp.list, "comp")
circ.means <- join.means(circ.list, "circ")
sym.means <- join.means(sym.list, "sym")
asym.means <- join.means(asym.list, "asym")
FT.means <- join.means(FT.list, "FT")

jdat <- cbind(comp.means, circ.means[2:ncol(circ.means)], sym.means[2:5], asym.means[2:5])
jdat <- merge(jdat, FT.means, by="geno", all=T)
head(jdat)
summary(jdat)
dim(jdat)
jdat$circ.Area <- jdat$circ.Area / 10000

jdat2 <- cbind(comp.means, circ.means[2:ncol(circ.means)], sym.means[2:5])
jdat2$circ.Area <- jdat2$circ.Area / 10000

colnames(circ.means)
jdat3 <- cbind(circ.means, sym.means[2:5])

ggpairs(data = na.omit(jdat), columns = 2:19, upper = list(continuous = "cor"), lower = list(continuous = "smooth", params = c(alpha = 0.25, color = "blue")), diag = list(continuous = "density"))
ggpairs(data = na.omit(jdat), columns = 2:19, upper = list(continuous = "cor"), lower = list(continuous = "smooth", params = c(alpha = 0.25, color = "blue")), diag = list(continuous = "bar"))

pdf(file="gpairs_all_16x9.pdf", width=16, height=9)
gpairs(na.omit(jdat)[2:19], upper.pars=list(scatter='corrgram'), lower.pars=list(scatter='loess'), outer.labels="all", scatter.pars=list(pch=16, col=alpha("black", 0.25)), diag.pars=list(fontsize=12) )
dev.off()

pdf(file="gpairs_all_10x7.5.pdf", width=10, height=7.5)
gpairs(na.omit(jdat)[2:19], upper.pars=list(scatter='corrgram'), lower.pars=list(scatter='loess'), outer.labels="all", scatter.pars=list(pch=16, col=alpha("black", 0.25)), diag.pars=list(fontsize=7), axis.pars=list(n.ticks=3, fontsize=6) )
dev.off()

pdf(file="gpairs_comp.circ_7.5x7.5.pdf", width=7.5, height=7.5)
gpairs(jdat2[2:10], upper.pars=list(scatter='corrgram'), lower.pars=list(scatter='loess'), outer.labels="all", scatter.pars=list(pch=16, col=alpha("black", 0.25)), diag.pars=list(fontsize=10) )
dev.off()

pdf(file="gpairs_circ.sym_7.5x7.5.pdf", width=7.5, height=7.5)
gpairs(jdat2[6:14], upper.pars=list(scatter='corrgram'), lower.pars=list(scatter='loess'), outer.labels="all", scatter.pars=list(pch=16, col=alpha("black", 0.25)), diag.pars=list(fontsize=10) )
dev.off()

pdf(file="gpairs_circ.sym_7.5x7.5.pdf", width=7.5, height=7.5)
gpairs(jdat2[6:14], upper.pars=list(scatter='corrgram'), lower.pars=list(scatter='loess'), outer.labels="all", scatter.pars=list(pch=16, col=alpha("black", 0.25)), diag.pars=list(fontsize=10) )
dev.off()

# Plot scatter plots of: PC3 and Aspect Ratio, and PC3 and Solidity
colnames(jdat3)
head(jdat3)
AR_vs_symPC3 <- ggplot(data=jdat3, aes(x=circ.AR, y=sym.PC3)) + geom_point(alpha=0.5) + labs(x="Aspect Ratio", y="symmetric PC3")
ggsave("AR_vs_symPC3.pdf", AR_vs_symPC3, width=7.5, height=7.5)
AR_vs_symPC3_text <- ggplot(data=jdat3, aes(x=circ.AR, y=sym.PC3)) + geom_text(aes(label=substr(geno,5,7)), size=2.8) + labs(x="Aspect Ratio", y="symmetric PC3")
ggsave("AR_vs_symPC3_text.pdf", AR_vs_symPC3_text, width=7.5, height=7.5)
Solid_vs_symPC3 <- ggplot(data=jdat3, aes(x=circ.Solidity, y=sym.PC3)) + geom_point(alpha=0.5) + labs(x="Solidity", y="symmetric PC3")
ggsave("Solidity_vs_symPC3.pdf", Solid_vs_symPC3, width=7.5, height=7.5)
Solid_vs_symPC3_text <- ggplot(data=jdat3, aes(x=circ.Solidity, y=sym.PC3)) + geom_text(aes(label=substr(geno,5,7)), size=2.8) + labs(x="Solidity", y="symmetric PC3")
ggsave("Solidity_vs_symPC3_text.pdf", Solid_vs_symPC3_text, width=7.5, height=7.5)


#--------------
# jdat$geno[which(is.na(jdat$FT.FT))] # PEN     BIL_157 BIL_498 BIL_321 BIL_525 BIL_416 BIL_294 BIL_327
# # which of the above have introgressions in chromosome 6? => NONE!
# lateFT.vector <- as.character(jdat$geno[which(is.na(jdat$FT.FT))]) # vector of very late-flowering BILs that were still in the dataset after final genotyping
# lateFT.vector <- lateFT.vector[lateFT.vector!='PEN']
# lateFT.vector
# 
# names(FT.list)
# FT.list[[1]][[1]]
# head(FT.list[[1]][[1]])
# mean.sd <- mean(FT.list[[1]][[1]]$StdErr) # [1] 6.145809 ~ 6.15  ...to be used to "impute" data for very late-flowering BILs
