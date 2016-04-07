# plot of punctate panel for finest submapping

library(ggplot2)
library(stringr)
library(plyr)
library(reshape2)

setwd("~/UCD/BILs")

#-------
punc <- read.csv("punctate_finemap.csv", header = F)
dim(punc) #20 14

colnames(punc) <- paste0("bin", str_pad(0:23, 2, pad="0"))

punc.coords <- as.data.frame(t(punc[1:2,]))
punc.coords$bin <- row.names(punc.coords)
colnames(punc.coords) <- c("start", "end", "bin")
punc.coords <- punc.coords[-1,]
punc.coords$start <- as.numeric(as.character(punc.coords$start))
punc.coords$end <- as.numeric(as.character(punc.coords$end))
summary(punc.coords)
punc.coords$start.n <- 0:22
punc.coords$end.n <- 1:23

punc.geno <- as.data.frame(punc[3:20,])
colnames(punc.geno)[1] <- "bil"
punc.geno <- droplevels(punc.geno)
summary(punc.geno)
dim(punc.geno)
punc.geno$ymax <- 18:1
punc.geno$ymin <- 17:0

m.punc.geno <- melt(punc.geno, id.vars=c(1,25:26), variable.name="bin")

j.punc <- merge(punc.coords, m.punc.geno, by="bin")
colnames(j.punc)[9] <- "geno"
head(j.punc)
j.punc$bil.lab <- ""
j.punc$bil.lab[1:18] <- as.character(j.punc$bil[1:18])
j.punc$bil.lab <- as.factor(j.punc$bil.lab)
summary(j.punc)

ggplot(j.punc, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=geno)) + geom_rect() + theme_classic(20) + xlim(61900000,64834305) +
  scale_fill_manual(values=c("magenta4", "green4")) + geom_text(data=j.punc[1:18,], aes(label=bil, y=ymin+0.5), x=61898100, hjust="left", size=7)

ggplot(j.punc, aes(xmin=start.n, xmax=end.n, ymin=ymin, ymax=ymax, fill=geno)) + geom_rect(color="black") + theme_classic(20) + xlim(-2,23) +
  scale_fill_manual(values=c("magenta", "green")) + geom_text(data=j.punc[1:18,], aes(label=bil, y=ymin+0.5), x=-0.25, hjust="right", size=7)

punc.geno$bil

#------
# Use the actual boundary files instead

punc.reord <- read.csv("punctate_finemap_reordered_trimmed.csv", head=F)

bil.order <- as.character(punc.reord[,1])
bil.order
files <- paste0(bil.order, ".boundaries")
files

setwd("~/UCD/BILs/punctate_boundaries/extracted_chr10_boundaries/")

punc.boundaries <- data.frame(start=numeric(), end=numeric(), geno=character(), bil=character(), ymax=numeric(), ymin=numeric(), stringsAsFactors=FALSE)
# punc.boundaries <- data.frame(matrix(,ncol=6))
for (i in 1:length(files)) {
  bil.id <- substr(files[i], 1, 7)
  bound <- read.delim(files[i], header=F, col.names=c("chr", "start", "end", "geno"), stringsAsFactors=FALSE, colClasses=c("character", "numeric", "numeric", "character"))[2:4]
  bound$bil <- as.character(bil.id)
  bound$ymax <- 26.5 - (i - 1) - ((i - 1) * 0.5)
  bound$ymin <- 26.5 - i - ((i -1) * 0.5)
  punc.boundaries <- rbind(punc.boundaries, bound)
}

# trim the starts so they don't get excluded by xlim
punc.bound.trim <- punc.boundaries
punc.bound.trim$lab <- ""

# gotta also remove lines prior to 64 Mbp?
for (i in 1:nrow(punc.bound.trim)) {
  if (punc.bound.trim$start[i] < 64000000) punc.bound.trim$start[i] = 64000000
  idx <- grep(punc.bound.trim$bil[i], punc.bound.trim$bil)[1]
  punc.bound.trim$lab[idx] <- punc.bound.trim$bil[i]
}

punc.bound.trim$y <- punc.bound.trim$ymin + 0.5
punc.bound.trim$col <- "black"
punc.bound.trim$col[1:31] <- "purple"

punc.bound.trim[punc.bound.trim$bil=="BIL_009",]
punc.bound.trim[punc.bound.trim$bil=="BIL_384",]
punc.bound.trim[punc.bound.trim$bil=="BIL_188",]

punc.plot <-  ggplot(punc.bound.trim, aes(xmin=start, xmax=end, ymax=ymax, ymin=ymin, y=y, fill=geno)) + 
  geom_rect(color="black") + xlim(63953000, 64829650) + 
  geom_text(aes(label=lab, color=col), x=63997000, hjust="right", size=9, fontface = "bold") +
  scale_fill_manual(values=c("magenta", "green")) + theme_classic(28) +
  scale_color_manual(values=c("black", "mediumpurple1")) +
  theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(), axis.title=element_blank()) +
  geom_segment(x=64395176, xend=64395176, y=-0.5, yend=27, linetype=2) +
  geom_segment(x=64509358, xend=64509358, y=-0.5, yend=27, linetype=2)

setwd("~/UCD/BILs/")
ggsave("punctate_fineMap_subpanel.pdf", punc.plot, width=556/20, height=170/15)
ggsave("punctate_fineMap_subpanel.svg", punc.plot, width=556/20, height=170/15)

# 64443068 - 64407814
# 64509358 - 64395176
# 64509358 - 64407814
# 64509358 - 64448843

