setwd("~/UCD/BILs/IL_BIL_flowering/Field_flowering_time/")

data <- read.csv("Flowering_Time_6.csv")

head(data)

summary(data)

any(duplicated(data$stake))

data$block <- substr(data$stake,1,1)

data.BIL <- data[!data$block %in% c("A","B","C","D","E","F","G","H","I","J"),]

head(data.BIL)

summary(data)

data.BIL$plant <- grep(pattern="[0-9]*",data.BIL$stake)

data.BIL$row <- ceiling(data.BIL$plant/37)

head(data.BIL,n=100)

data.BIL$col <- ((data.BIL$plant-1) %% 37) +1

head(data.BIL,n=100)

summary(data.BIL)

#convert columns and rows to factors
#note that it might be best to use ordered factors...needs to be investigated.

data.BIL$col <- factor(data.BIL$col)

data.BIL$row <- factor(data.BIL$row)

data.BIL$block <- factor(data.BIL$block)

names(data.BIL)[2] <- "BIL"

data.BIL$BIL <- droplevels(data.BIL$BIL)

sort(table(data.BIL$BIL)) #looks good!

library(lme4) # load in package with the mixed-effects modeling function

head(data.BIL)

data.BIL$BIL <- relevel(data.BIL$BIL,ref="M82")

## restructure data.BIL to
FT.BIL <- data.BIL
FT.BIL <- FT.BIL[, colnames(FT.BIL) != 'plant']
summary(FT.BIL)
colnames(FT.BIL)[1] <- 'plant'
FT.BIL$plant <- as.character(FT.BIL$plant)
# ADD zeroes to too short BIL IDs
len3idx <- which(nchar(FT.BIL$plant) == 3) # length == 3 BIL names indices
len2idx <- which(nchar(FT.BIL$plant) == 2) # length == 2 BIL names indices
FT.BIL$plant[len3idx] <- paste0("BIL_0", substr(FT.BIL$plant[len3idx], 2, 3))
FT.BIL$plant[len2idx] <- paste0("BIL_00", substr(FT.BIL$plant[len2idx], 2, 2))
FT.BIL$plant <- as.factor(FT.BIL$plant)
summary(FT.BIL)
# Fix BIL labels, i.e. prepend "BIL_"
FT.BIL$BIL <- as.character(FT.BIL$BIL)
FT.BIL$BIL[which(FT.BIL$BIL!="M82" & FT.BIL$BIL!="PEN")] <- paste0("BIL_", FT.BIL$BIL[which(FT.BIL$BIL!="M82" & FT.BIL$BIL!="PEN")])
FT.BIL$BIL <- as.factor(FT.BIL$BIL)
FT.BIL$BIL <- relevel(FT.BIL$BIL, ref="M82")
summary(FT.BIL)
dim(FT.BIL) # 1665    6
sort(unique(FT.BIL$BIL[which(is.na(FT.BIL$days))]))
FT.BIL_all <- FT.BIL

FT.BIL <- na.omit(FT.BIL)
dim(FT.BIL) # 1539    6
FT.BIL <- droplevels(FT.BIL)
save(FT.BIL, file="FT.BIL.Rdata")

# Deprecated code
#--------
# ## ADD simulated data for very late-flowering BILs
# # Sample from rnorm(n, mean = 128 days, sd = 6.15 -- the average within line st.dev.)
# # 128 days is 3*stdev + 1 day beyond the last measured day (i.e. 109 days).
# lateFT.vector # vector of very late-flowering BILs that were still in the dataset after final genotyping
# # [1] "BIL_157" "BIL_498" "BIL_321" "BIL_525" "BIL_416" "BIL_294" "BIL_327"
# # the lateFT.vector vector comes from the pairs_plot.R script
# for(i in 1:length(lateFT.vector) ) {
#   plants <- as.character(FT.BIL_all$plant[FT.BIL_all$BIL==lateFT.vector[i] ])
#   for(j in 1:length(plants) ) {
#     FT.BIL_all$days[FT.BIL_all$plant == plants[j] ] <- rnorm(1, mean = 128, sd = 6.15)
#   }
# }
# 
# FT.BIL <- FT.BIL_all
# FT.BIL <- na.omit(FT.BIL)
# dim(FT.BIL) # 1549    6  , i.e. not all 7 BILs "imputed" above have plants that were alive it seems, as the total only increased by 10 and not by 14
# FT.BIL <- droplevels(FT.BIL)
# save(FT.BIL, file="FT.BIL.Rdata")
#---------

## Modeling column-wise to get a sense for the compaction effect

mcCols <- lmer(days ~ col + (1|block) + (1|BIL), data = data.BIL)
mcCols
cols2 <- as.data.frame(fixef(mcCols))
rownames(cols2)[1] <- "col1"
names(cols2)[1] <- "coef"
cols2$val[1] <- cols2$coef[1]
cols2$val[2:37] <- cols2$coef[1] + cols2$coef[2:37]
cols2$col <- as.factor(rownames(cols2))

summary(cols2)

#library(ggplot2)

ggplot(cols2, aes(x=col, y=val) ) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="column", y="Flowering Time (days)", title="inc. BIL-gt as random effect")

-----------

lme1 <- lmer(days ~ BIL + (1|block) + (1|row) + (1|col),data=data.BIL)
summary(lme1)

lme2 <- update(lme1,.~. - (1|row))
summary(lme2)

anova(lme1,lme2)

lme3 <- update(lme2,.~. - (1|block))

anova(lme2,lme3)

lme4 <- update(lme2,.~. - (1|col)) #note : must have 1 random effect to use lmer
anova(lme2,lme4) #so leave column in the model.  Final model is lme3

#p-values for fixed effects
#need another package
library(languageR)

lme3.pvals <- pvals.fnc(lme3,addPlot=F) #takes 5 or 10 minutes

pvals <- lme3.pvals$fixed

head(pvals)
summary(pvals)
pvals$Estimate <- as.numeric(pvals$Estimate)
pvals$pMCMC <- as.numeric(pvals$pMCMC)
pvals$pMCMC[1] <- 1 #the M82 pval is meaningless...it is whether M82 is different from 0.  Set to 1 so that it is not flagged

#add the intercept to all of the ILs:

pvals$Estimate[-1] <- pvals$Estimate[-1] + pvals$Estimate[1]

pvals$gt <- rownames(pvals)
pvals$gt[1] <- "M82"
pvals$gt.plot <- sub("IL[0-9][[:print:]]*","IL",pvals$gt)
head(pvals)

#convert p-values to text
pvals$pMCMC.txt <- ""
pvals$pMCMC.txt[pvals$pMCMC<.05] <- "*"
pvals$pMCMC.txt[pvals$pMCMC<.01] <- "**"
pvals$pMCMC.txt[pvals$pMCMC<.001] <- "***"

summary(pvals)

library(ggplot2)

p <- ggplot(pvals,aes(x=gt,y=Estimate,fill=gt.plot)) #set up the basic plot, put it in object "p"
p <- p + geom_bar(stat="identity") #define that want a bar plot and do not need any stats done
p <- p + scale_fill_brewer(palette="Set2") # nicer colors
p <- p + opts(title="Flowering time in ILs")

    #the next lines set up the x-axis
p <- p + scale_x_discrete(name="genotype", #title for the x axis
                          limits=pvals$gt[order(pvals$Estimate)]) #use "limits" to define the order

p <- p + opts(axis.text.x = theme_text(angle=-90)) #orient the text vertically

p <- p + geom_text(mapping=aes(label=pMCMC.txt,angle=90,vjust=.8))

p #plot it!





