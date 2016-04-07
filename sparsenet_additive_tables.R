# output sparsenet additive QTL search tables
setwd("~/UCD/BILs/final_additive_sparsenet_results/")

output.tables <- function(dat, trait.class) {
  for (i in 1:length(dat)) {
    trait.name <- names(dat)[i]
    tab.name <- paste0(trait.class, ".", trait.name, ".csv")
    tab <- dat[[i]]$non.zero.coefs
    write.csv(tab, file=tab.name, quote=F, row.names=F)
  }
}

load("comp.map.Rdata")
output.tables(comp.map, "comp")

load("circ.map.Rdata")
output.tables(circ.map, "shape")

load("FT.map.Rdata")
output.tables(FT.map, "FT")

load("sym.map.Rdata")
output.tables(sym.map, "sym")

load("asym.map.Rdata")
output.tables(asym.map, "asym")
