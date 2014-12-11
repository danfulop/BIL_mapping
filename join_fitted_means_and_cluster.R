setwd("~/UCD/BILs/")
# join all fitted means' data.frames => may be able to do it w/ ldply
collate.results <- function(l) {
  for(i in 1:length(l) ) {
    if (i==1) {
      tmp <- l[[i]][[1]][c('Estimate', 'geno')]
      names(tmp)[1] <- paste0(names(l)[i], '.Estimate')
      collated.df <- tmp 
    } else {
      tmp <- l[[i]][[1]][c('Estimate', 'geno')] # dump table w/ estimate + geno into a tmp data.frame
      names(tmp)[1] <- paste0(names(l)[i], '.Estimate') # tag trait name onto "Estimate" => for sym and asym PCs, add "sym" or "asym"
      collated.df <- merge (collated.df, tmp, by='geno') # add/merge to output data.frame => create it if trait #1 within list, otherwise merge
    }
  }
  collated.df
}

# complexity data
comp.means <- collate.results(comp.list)
head(comp.means)
# circ. data
circ.means <- collate.results(circ.list)
head(circ.means)
# Sym. EFD-PCs
sym.means <- collate.results(sym.list)
head(sym.means)
names(sym.means)[2:ncol(sym.means)] <- paste0("sym.", names(sym.means)[2:ncol(sym.means)])
# Asym. EFD-PCs
asym.means <- collate.results(asym.list)
head(asym.means)
names(asym.means)[2:ncol(asym.means)] <- paste0("asym.", names(asym.means)[2:ncol(asym.means)])

# merge 4 resulting data.frames
all.traits <- merge(comp.means, circ.means, by='geno')
all.traits <- merge(all.traits, sym.means, by='geno')
all.traits <- merge(all.traits, asym.means, by='geno')
head(all.traits)
row.names(all.traits) <- all.traits$geno
all.traits.nolabs <- all.traits[2:ncol(all.traits)]
head(all.traits.nolabs)

# do simple hierarchical clustering
hclust.all <- hclust(dist(all.traits.nolabs))
save(hclust.all, file="hclust.all.Rdata")
png(filename="hclust.all.png", width=5000, height=1500)
plot(hclust.all)
dev.off()
pdf(file="hclust.all.pdf", width=65, height=15)
plot(hclust.all)
dev.off()

# Try pvclust
library(pvclust)
pvclust.all <- pvclust(all.traits.nolabs)
save(pvclust.all, file="pvclust.all.Rdata")
plot(pvclust.all)

# transpose matrix/DF so that pvclust clusters BILs and not traits
t.all <- t(all.traits.nolabs)
head(t.all)
# distance: correlation, cluster method: average
pvclust.all.bils <- pvclust(t.all)
save(pvclust.all.bils, file="pvclust.all.bils.Rdata")
pdf(file="pvclust.all.bils.pdf", width=65, height=15)
plot(pvclust.all.bils)
dev.off()
# distance: correlation, cluster method: ward
pvclust.ward.all.bils <- pvclust(t.all, method.hclust="ward.D2")
save(pvclust.ward.all.bils, file="pvclust.ward.all.bils.Rdata")
pdf(file="pvclust.ward.all.bils.pdf", width=65, height=15)
plot(pvclust.ward.all.bils)
dev.off()
# choose clusters
trait.clusters <- pvpick(pvclust.ward.all.bils)
trait.clusters
