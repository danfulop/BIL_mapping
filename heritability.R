# Script to calculate broad-sense heritabilities of the BIL traits

# Pseudo-code:
# Use a random effects model on the BILs w/o M82:
# trait ~ 1 + (1 | BIL) + (1 | block) + (1 | expt)
    # do I need to use sum to zero contrast?!
        # maybe not since there are no fixed effects
    # I could do away with the block effect because there are so few blocks that it's hard to properly estimate a block effect

# does the mixed crossing scheme present a problem?

library(lme4)
library(blme)
library(ggplot2)
library(MCMCglmm)

# library(plyr)
# library(coefplot2)
# library(fdrtool)
# library(parallel)
# library(doParallel)
# library(foreach)

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

comp.all <- lmer(all ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=comp.noparents)
summary(comp.all)
VarCorr(comp.all)
est <- as.data.frame(VarCorr(comp.all))
h2.comp.all <- est[2,4] / sum(est$vcov) # 0.3954368
stoch.comp.all <- with(as.data.frame(VarCorr(comp.all)), (vcov[3] + vcov[4]) / vcov[1] ) # 0.8723195
comp.all.profCI <- confint(comp.all)
comp.all.profCI
comp.all.bootCI <- confint(comp.all, method="boot")
comp.all.bootCI^2

comp.sec <- lmer(sec ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=comp.noparents)
h2.comp.sec <- with(as.data.frame(VarCorr(comp.sec)), vcov[2] / sum(vcov) ) # 0.3204147
stoch.comp.sec <- with(as.data.frame(VarCorr(comp.sec)), (vcov[3] + vcov[4]) / vcov[1] ) # 0.9417766
comp.sec.profCI <- confint(comp.sec)

comp.int <- lmer(int ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=comp.noparents)
h2.comp.int <- with(as.data.frame(VarCorr(comp.int)), vcov[2] / sum(vcov) ) # 0.4142432
stoch.comp.int <- with(as.data.frame(VarCorr(comp.int)), (vcov[3] + vcov[4]) / vcov[1] ) # 1.205608

comp.pri <- lmer(pri ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=comp.noparents)
h2.comp.pri <- with(as.data.frame(VarCorr(comp.pri)), vcov[2] / sum(vcov) ) # 0.2276865
stoch.comp.pri <- with(as.data.frame(VarCorr(comp.pri)), (vcov[3] + vcov[4]) / vcov[1] ) # 2.180645

#-----
load("labels.Rdata")
circ <- read.delim("BIL.circAR.all.2011.txt")
circ <- merge(circ, labels, by="plant")
circ <- droplevels(circ)
names(circ)[2] <- "leaflet.type"
circ <- circ[circ$FinBIL!="M82" & circ$FinBIL!="PEN",]
circ <- droplevels(circ)

area <- lmer(Area ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant / leaflet.type), data=circ)
h2.area <- with(as.data.frame(VarCorr(area)), vcov[3] / sum(vcov) ) # 0.06015235
stoch.area <- with(as.data.frame(VarCorr(area)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 0.5866483

round <- lmer(Round ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant / leaflet.type), data=circ)
h2.round <- with(as.data.frame(VarCorr(round)), vcov[3] / sum(vcov) ) # 0.2796237
stoch.round <- with(as.data.frame(VarCorr(round)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 3.798839

AR <- lmer(AR ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant / leaflet.type), data=circ)
h2.AR <- with(as.data.frame(VarCorr(AR)), vcov[3] / sum(vcov) ) # 0.2707402
stoch.AR <- with(as.data.frame(VarCorr(AR)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 3.989171

circularity <- lmer(Circ. ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant / leaflet.type), data=circ)
h2.circ <- with(as.data.frame(VarCorr(circularity)), vcov[3] / sum(vcov) ) # 0.2697808
stoch.circ <- with(as.data.frame(VarCorr(circularity)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 3.127448

solid <- lmer(Solidity ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant / leaflet.type), data=circ)
h2.solid <- with(as.data.frame(VarCorr(solid)), vcov[3] / sum(vcov) ) # 0.2557488
stoch.solid <- with(as.data.frame(VarCorr(solid)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 4.359837

#----
sym <- read.delim("BIL_Spcascores.txt")
sym <- merge(sym, labels, by="plant")
sym <- droplevels(sym)
names(sym)[2] <- "leaflet.type"
sym <- sym[with(sym, FinBIL!="M82" & FinBIL!="PEN"),]
sym <- droplevels(sym)

# relabel leaf.type (by prepending the plantID) to enable nested ranefs in MCMCglmm
mc.sPC1 <- MCMCglmm(fixed = PC1 ~ 1, random = ~    )

sPC1 <- lmer(PC1 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
summary(sPC1)
h2.sPC1 <- with(as.data.frame(VarCorr(sPC1)), vcov[3] / sum(vcov) ) # 0.05943493
with(as.data.frame(VarCorr(sPC1)), (vcov[1] + vcov[2]) /  sum(vcov) ) # 0.1471853
with(as.data.frame(VarCorr(sPC1)), sum(vcov) ) # 0.001141159
rt.sPC1 <- with(as.data.frame(VarCorr(sPC1)), (vcov[1] + vcov[2]) / vcov[3] ) # 2.476411
stoch.sPC1 <- with(as.data.frame(VarCorr(sPC1)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 5.390345

#----
x <- sym
y <- 3
randform <- "+ (1| plant/leaflet.type)"

x <- comp.noparents
y <- 5
randform <- "+ (1| plant)"

# Function to calculate broad-sense heritability, repeatability, and "stochastic index" from lme4 models
paramMLEs <- function(x, y, i, randform) {  # Include as input column index of response values, i.e. y
  resp <- names(x)[y] # assign the column name of the 'y' index
  formula <- as.formula( paste0( get("resp"), " ~ 1 + (1|FinBIL) + (1|block) ", get("randform") ) )
  fit <- lmer(formula, data=x, REML=TRUE)
  tab <- as.data.frame(VarCorr(fit)) # save random effects table
  gen.idx <- which(tab$grp == "FinBIL")
  if (gen.idx == 3) {
    h2 <- with(tab, vcov[3] / sum(vcov))
    r2 <- with(tab, sum(vcov[1:3]) / sum(vcov))
    stoch <- with(tab, sum(vcov[4:5]) / sum(vcov[1:2]))
    
  } else if (gen.idx == 2) {
    h2 <- with(tab, vcov[2] / sum(vcov))
    r2 <- with(tab, sum(vcov[1:2]) / sum(vcov))
    stoch <- with(tab, sum(vcov[3:4]) / vcov[1])
  }
  res <- c(h2=h2, r2=r2, stoch=stoch)
}

#----

sPC2 <- lmer(PC2 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC2 <- with(as.data.frame(VarCorr(sPC2)), vcov[3] / sum(vcov) ) # 0.1365813
with(as.data.frame(VarCorr(sPC2)), (vcov[1] + vcov[2]) /  sum(vcov) ) # 0.1215611
with(as.data.frame(VarCorr(sPC2)), sum(vcov) ) # 0.0006895134
rt.sPC2 <- with(as.data.frame(VarCorr(sPC2)), (vcov[1] + vcov[2]) / vcov[3] ) # 0.8900276
stoch.sPC2 <- with(as.data.frame(VarCorr(sPC2)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 6.102754

sPC3 <- lmer(PC3 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC3 <- with(as.data.frame(VarCorr(sPC3)), vcov[3] / sum(vcov) ) # 0.1527682
with(as.data.frame(VarCorr(sPC3)), (vcov[1] + vcov[2]) / sum(vcov) ) # 0.1408828
with(as.data.frame(VarCorr(sPC3)), sum(vcov) ) # 0.0005015638
rt.sPC3 <- with(as.data.frame(VarCorr(sPC3)), (vcov[1] + vcov[2]) / vcov[3] ) # 0.9221996
stoch.sPC3 <- with(as.data.frame(VarCorr(sPC3)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 5.013737

alt.sPC3 <- lmer(PC3 ~ 1 + (1 | FinBIL / plant / leaflet.type) + (1 | block), data=sym)
summary(sPC3)
summary(alt.sPC3)

noint.sPC3 <- lmer(PC3 ~ 0 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
as.data.frame(VarCorr(noint.sPC3))
as.data.frame(VarCorr(sPC3))

sPC4 <- lmer(PC4 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC4 <- with(as.data.frame(VarCorr(sPC4)), vcov[3] / sum(vcov) ) # 0.03440149
with(as.data.frame(VarCorr(sPC4)), (vcov[1] + vcov[2]) / sum(vcov) ) # 0.1226908
with(as.data.frame(VarCorr(sPC4)), sum(vcov) ) # 0.0004484479
rt.sPC4 <- with(as.data.frame(VarCorr(sPC4)), (vcov[1] + vcov[2]) / vcov[3] ) # 3.566438
stoch.sPC4 <- with(as.data.frame(VarCorr(sPC4)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 6.870179

sPC5 <- lmer(PC5 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC5 <- with(as.data.frame(VarCorr(sPC5)), vcov[3] / sum(vcov) ) # 0.04993123
rt.sPC5 <- with(as.data.frame(VarCorr(sPC5)), (vcov[1] + vcov[2]) / vcov[3] ) # 1.418052
stoch.sPC5 <- with(as.data.frame(VarCorr(sPC5)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 12.41809

sPC6 <- lmer(PC6 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC6 <- with(as.data.frame(VarCorr(sPC6)), vcov[3] / sum(vcov) ) # 0.03320082
rt.sPC6 <- with(as.data.frame(VarCorr(sPC6)), (vcov[1] + vcov[2]) / vcov[3] ) # 0.7797706
stoch.sPC6 <- with(as.data.frame(VarCorr(sPC6)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 36.34399

sPC7 <- lmer(PC7 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC7 <- with(as.data.frame(VarCorr(sPC7)), vcov[3] / sum(vcov) ) # 0.03098909
rt.sPC7 <- with(as.data.frame(VarCorr(sPC7)), (vcov[1] + vcov[2]) / vcov[3] ) # 1.132871
stoch.sPC7 <- with(as.data.frame(VarCorr(sPC7)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 26.60194

sPC8 <- lmer(PC8 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC8 <- with(as.data.frame(VarCorr(sPC8)), vcov[3] / sum(vcov) ) # 0.03651483
rt.sPC8 <- with(as.data.frame(VarCorr(sPC8)), (vcov[1] + vcov[2]) / vcov[3] ) # 1.365528
stoch.sPC8 <- with(as.data.frame(VarCorr(sPC8)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 18.32304

sPC9 <- lmer(PC9 ~ 1 + (1 | FinBIL) + (1 | block) + (1 |  plant / leaflet.type), data=sym)
h2.sPC9 <- with(as.data.frame(VarCorr(sPC9)), vcov[3] / sum(vcov) ) # 0.02540988
rt.sPC9 <- with(as.data.frame(VarCorr(sPC9)), (vcov[1] + vcov[2]) / vcov[3] ) # 2.050967
stoch.sPC9 <- with(as.data.frame(VarCorr(sPC9)), (vcov[4] + vcov[5]) / (vcov[1] + vcov[2]) ) # 17.70082


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

# Maybe keep terminal leaflet data for H^2
# OR analyse H^2 sep. for term. and lats.
# INCLUDE the total variance estimate as a means of assessing total info. amount when comparing sym and asym PCs

aPC1 <- lmer(PC1 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC1 <- with(as.data.frame(VarCorr(aPC1)), vcov[2] / sum(vcov) ) # 0.01965591
with(as.data.frame(VarCorr(aPC1)), sum(vcov) ) # 0.0005644242
wgv.aPC1 <- with(as.data.frame(VarCorr(aPC1)), vcov[1] / sum(vcov) ) # 0.007763238
rt.aPC1 <- with(as.data.frame(VarCorr(aPC1)), vcov[1] / vcov[2] ) # 0.3949569
stoch.aPC1 <- with(as.data.frame(VarCorr(aPC1)), (vcov[3] + vcov[4]) / vcov[1] ) # 125.2803

aPC2 <- lmer(PC2 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC2 <- with(as.data.frame(VarCorr(aPC2)), vcov[2] / sum(vcov) ) # 0.03478533
with(as.data.frame(VarCorr(aPC2)), sum(vcov) ) # 0.000351139
wgv.aPC2 <- with(as.data.frame(VarCorr(aPC2)), vcov[1] / sum(vcov) ) # 0.03466213
rt.aPC2 <- with(as.data.frame(VarCorr(aPC2)), vcov[1] / vcov[2] ) # 0.9964583
stoch.aPC2 <- with(as.data.frame(VarCorr(aPC2)), (vcov[3] + vcov[4]) / vcov[1] ) # 26.84638

aPC3 <- lmer(PC3 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC3 <- with(as.data.frame(VarCorr(aPC3)), vcov[2] / sum(vcov) ) # 0.02453142
with(as.data.frame(VarCorr(aPC3)), sum(vcov) ) # 0.0001977388
wgv.aPC3 <- with(as.data.frame(VarCorr(aPC3)), vcov[1] / sum(vcov) ) # 0.02969915
rt.aPC3 <- with(as.data.frame(VarCorr(aPC3)), vcov[1] / vcov[2] ) # 1.210658
stoch.aPC3 <- with(as.data.frame(VarCorr(aPC3)), (vcov[3] + vcov[4]) / vcov[1] ) # 31.845

aPC4 <- lmer(PC4 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC4 <- with(as.data.frame(VarCorr(aPC4)), vcov[2] / sum(vcov) ) # 0.03109231
with(as.data.frame(VarCorr(aPC4)), sum(vcov) ) # 0.0001583888
rt.aPC4 <- with(as.data.frame(VarCorr(aPC4)), vcov[1] / vcov[2] ) # 0.7685719
stoch.aPC4 <- with(as.data.frame(VarCorr(aPC4)), (vcov[3] + vcov[4]) / vcov[1] ) # 39.54572

aPC5 <- lmer(PC5 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC5 <- with(as.data.frame(VarCorr(aPC5)), vcov[2] / sum(vcov) ) # 0.0461871
stoch.aPC5 <- with(as.data.frame(VarCorr(aPC5)), (vcov[3] + vcov[4]) / vcov[1] ) # 65.90329

aPC6 <- lmer(PC6 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC6 <- with(as.data.frame(VarCorr(aPC6)), vcov[2] / sum(vcov) ) # 0.03155996
stoch.aPC6 <- with(as.data.frame(VarCorr(aPC6)), (vcov[3] + vcov[4]) / vcov[1] ) # 42.57669

aPC7 <- lmer(PC7 ~ 1 + (1 | FinBIL) + (1 | block) + (1 | plant), data=asym)
h2.aPC7 <- with(as.data.frame(VarCorr(aPC7)), vcov[2] / sum(vcov) ) # 0.02135341
stoch.aPC7 <- with(as.data.frame(VarCorr(aPC7)), (vcov[3] + vcov[4]) / vcov[1] ) # 172.7999

