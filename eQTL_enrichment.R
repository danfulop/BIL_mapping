# Daniel Fulop 11/29/2015
# Script for finding overlaps between IL eQTL of leaf development genes and BIL epi-QTL from leaf complexity and gross leaflet morphology traits (i.e. excluding EFD-PC traits)
# the leaf development genes are the gene network expanded list from Yasu's PNAS paper plus the eQTL cluster module that's enriched for leaf development GO terms

library(stringr)
library(plyr)
library(IRanges)
library(foreach)
library(doParallel)

setwd("~/UCD/BILs/fine_mapping_n_eqtl_enrichment/eQTL_enrichment/")

# input master eQTL table
eqtl.tab <- read.delim("~/UCD/BILs/eQTL_data/No6.29.added.cis.trans.NOpseudoRENAME.GeneAnnot.tsv")
head(eqtl.tab)
names(eqtl.tab)

sub.eqtl <- eqtl.tab[c(1:2,4, 10:12, 16:17, 24)] # keep only needed columns
summary(sub.eqtl)

trans.eqtl <- sub.eqtl[sub.eqtl$cis.trans=="true.trans", 1:8] # keep only trans-eQTL
trans.eqtl <- droplevels(trans.eqtl)
names(trans.eqtl)[2] <- "chr.gene"
trans.eqtl$chr.cor.bin <- paste0("ch", str_pad(str_replace(as.character(trans.eqtl$cor.bin), "d\\.([0-9]+)[A-Z]", "\\1"), 2, pad="0"))

# EXCLUDE eQTL that have both ends in the same chromosome
nrow(trans.eqtl) # 3630
# trans.eqtl <- trans.eqtl[trans.eqtl$chr.gene != trans.eqtl$chr.cor.bin, ]
nrow(trans.eqtl) # 3162 
head(trans.eqtl)

rm(sub.eqtl)

# rows of 1st occurrence of genes, to extract the relevant info.
itag.idx <- which(!duplicated(trans.eqtl$itag))

gene.tab <- trans.eqtl[itag.idx, c(1:2,7:8)] # table of genes to use for overlap
nrow(gene.tab) # 1801 genes, down from 3630 because of multiple eQTL per gene (and down from 2115 when eQTL w/ both ends in same chr. hadn't been excluded)
gene.tab <- droplevels(gene.tab)
summary(gene.tab) 
head(gene.tab)
tail(gene.tab)

# Trim gene IDs of the gene model specification, i.e. last 4 characters
gene.tab$itag <- substr(gene.tab$itag, 1, 14)

# input Yasu's gene network expanded leaf development table
yasu.df <- read.csv("../../eQTL_data/literature_curated_leaf_larger.csv")
dim(yasu.df) # 697  5
head(yasu.df)

# keep only the geneID/itag vector and filter by having an eQTL
yasu.list <- intersect(yasu.df$Gene.ID, gene.tab$itag)
length(yasu.list) # 112 (down from 122 before chr1==chr2 filter)

# Input Steven's leaf module
leaf.mod <- read.csv("../../eQTL_data/annotation_cluster_16_list_SL.csv")
dim(leaf.mod) # 108  3
head(leaf.mod)

# Keep only genes with trans-eQTL
leaf.mod.list <- intersect(leaf.mod$ITAG, gene.tab$itag)
length(leaf.mod.list) # 56 ( down from 66 before chr1==chr2 filter)

# Join the Yasu's and Steven's lists
jt.list <- union(yasu.list, leaf.mod.list)
length(jt.list) # 153 (down from 173)

# Repeat w/ photosynthesis module
photo.mod <- read.csv("../../eQTL_data/annotation_cluster_11_list_photosynthesis.csv")
dim(photo.mod) # 715  6
photo.mod$ITAG <- substr(photo.mod$ITAG, 1, 14)
photo.mod.list <- intersect(photo.mod$ITAG, gene.tab$itag)
length(photo.mod.list) # 295 (down from 329 before chr1==chr2 filter)
(329 - 295) / 329 * 100 # ~10% filtered out

# Join the photo & joint lists
jt.photo.list <- union(jt.list, photo.mod.list)
length(jt.photo.list) # 373 (down from 421 before chr1==chr2 filter)
(421 - 373) / 421 * 100 # ~11% filtered out

# I need to filter the *both* gene.tab and cor.bin.tab by the joint list (jt.list) of leaf development genes
# *and* also filter trans.eqtl by jt.list (but keeping duplicates) so as to be able to form the list of gene:cor.bin labels

# filter gene.tab by jt.list
gene.tab.leaf <- gene.tab[which(gene.tab$itag %in% jt.list), ]
dim(gene.tab.leaf) # 153  4
head(gene.tab.leaf)

# filter gene.tab by jt.photo.list
gene.tab.photoleaf <- gene.tab[which(gene.tab$itag %in% jt.photo.list), ]
dim(gene.tab.photoleaf) # 373  4
head(gene.tab.photoleaf)

# get correlated bins for both joint lists
trans.eqtl$itag <- substr(trans.eqtl$itag, 1, 14)
cor.bins.leaf <- unique(trans.eqtl$cor.bin[trans.eqtl$itag %in% jt.list])
cor.bins.photoleaf <- unique(trans.eqtl$cor.bin[trans.eqtl$itag %in% jt.photo.list])

# attach IL bin coordinate information
il.bin <- read.csv("../bin2_coord.csv")
head(il.bin)
nrow(il.bin) # 126 sub-bins
summary(il.bin)

cor.bin.leaf.tab <- merge(il.bin, as.data.frame(cor.bins.leaf), by.x="bin", by.y="cor.bins.leaf")
nrow(cor.bin.leaf.tab) # 32 (down from 39 before chr1==chr2 filter)
head(cor.bin.leaf.tab)
cor.bin.leaf.tab <- droplevels(cor.bin.leaf.tab)
summary(cor.bin.leaf.tab) # 24 bins, 26 sub-bins (down from 31 bins, 33 sub-bins)
# generate chromosome column from bin names
cor.bin.leaf.tab$chr <- paste0("ch", str_pad(str_replace(as.character(cor.bin.leaf.tab$bin), "d\\.([0-9]+)[A-Z]", "\\1"), 2, pad="0"))
# sort by chromosome & sub-bin
cor.bin.leaf.tab <- cor.bin.leaf.tab[with(cor.bin.leaf.tab, order(chr, bin2)), ]
head(cor.bin.leaf.tab)

cor.bin.photoleaf.tab <- merge(il.bin, as.data.frame(cor.bins.photoleaf), by.x="bin", by.y="cor.bins.photoleaf")
nrow(cor.bin.photoleaf.tab) # 48 (down from 55)
head(cor.bin.photoleaf.tab)
cor.bin.photoleaf.tab <- droplevels(cor.bin.photoleaf.tab)
summary(cor.bin.photoleaf.tab) # 38 bins, 42 sub-bins (down from 45 bins, 49 sub-bins)
# generate chromosome column from bin names
cor.bin.photoleaf.tab$chr <- paste0("ch", str_pad(str_replace(as.character(cor.bin.photoleaf.tab$bin), "d\\.([0-9]+)[A-Z]", "\\1"), 2, pad="0"))
# sort by chromosome &  sub-bin
cor.bin.photoleaf.tab <- cor.bin.photoleaf.tab[with(cor.bin.photoleaf.tab, order(chr, bin2)), ]
head(cor.bin.photoleaf.tab)

# Generate gene:cor.bin pairs for each eQTL in the joint lists
leaf.trans.eqtl <- trans.eqtl[trans.eqtl$itag %in% jt.list, ]
nrow(leaf.trans.eqtl) # 238 (down from 274)
leaf.gene.cor.bin <- with(leaf.trans.eqtl, paste(itag, cor.bin, sep=":"))

photoleaf.trans.eqtl <- trans.eqtl[trans.eqtl$itag %in% jt.photo.list, ]
nrow(photoleaf.trans.eqtl) # 546 (down from 624)
photoleaf.gene.cor.bin <- with(photoleaf.trans.eqtl, paste(itag, cor.bin, sep=":"))

# Input BIL's bin info
load("/Users/Dani/UCD/BILs/leaf_traits/bin.stats.Rdata") # load bin information
bin.stats$chr <- as.factor(substr(bin.stats$chr,7,10)) # Trim "SL2.40" from chromosome names
num.trim.fx <- colwise(function(x) as.numeric(str_trim(x) ) ) # make a column-wise function to trim factor columns of white space and convert them to numeric. Coercion to char isn't needed since trimming does that already
bin.stats[3:5] <- num.trim.fx(bin.stats[3:5] )

# Load BIL bin correlation data to use during permutations
load("/Users/Dani/UCD/BILs/leaf_traits/bin.cor.Rdata")
dim(bin.cor) # 1049 1049

# load epistatic mapping results and concatenate them
load("/Users/Dani/UCD/BILs/final_epistatic_sparsenet_results/comp.epi.map.Rdata")
load("/Users/Dani/UCD/BILs/final_epistatic_sparsenet_results/circ.epi.map.Rdata")

# skip single locus QTL from epistatic model search. lapply to keep only non.zero.coefs' epistatic rows & add a trait column
extract_epi.qtl <- function(epi.map) {
  epi.qtl.list <- lapply(1:length(epi.map), function(i) {
    trait.name <- names(epi.map)[i]
    epi.qtl <- epi.map[[i]][[2]][grep("_x_", epi.map[[i]][[2]]$bin), c(1:2,4:6,10:12)]
    epi.qtl$trait <- trait.name
    epi.qtl
  })
  do.call(rbind, epi.qtl.list)
}

comp.epi.qtl <- extract_epi.qtl(comp.epi.map)
circ.epi.qtl <- extract_epi.qtl(circ.epi.map[2:5]) # Skip traits 1, area, because it has not QTL

rm(comp.epi.map, circ.epi.map)

jt.epi.qtl <- rbind(comp.epi.qtl, circ.epi.qtl) # join epi.qtl

# Deconvolute / divide table into single ends of the epi.qtl => Keep only 0.9 correlation info
epi.split <- jt.epi.qtl[c(1:2,6:9)]
epi.split$bin1 <- laply(epi.split$bin, function(x) str_split(x, "_x_")[[1]][1])
epi.split$bin2 <- laply(epi.split$bin, function(x) str_split(x, "_x_")[[1]][2])
epi.split$chr1 <- laply(epi.split$chr, function(x) str_split(x, "//")[[1]][1])
epi.split$chr2 <- laply(epi.split$chr, function(x) str_split(x, "//")[[1]][2])
epi.split$start1 <- laply(epi.split$int0.90.start, function(x) str_split(x, "//")[[1]][1])
epi.split$start2 <- laply(epi.split$int0.90.start, function(x) str_split(x, "//")[[1]][2])
epi.split$end1 <- laply(epi.split$int0.90.end, function(x) str_split(x, "//")[[1]][1])
epi.split$end2 <- laply(epi.split$int0.90.end, function(x) str_split(x, "//")[[1]][2])
head(epi.split)
nrow(epi.split) # 285

# Filter epi.qtl to exclude QTL w/ both ends in same chromosome
# epi.split <- epi.split[epi.split$chr1 != epi.split$chr2, ]
head(epi.split)
nrow(epi.split) # 186

(285 - 186) / 285 * 100 # 34.7% filtered out => YES, it was worth filtering out the QTL w/ both ends on the same chromosome,
# because as expected the proportion filtered out is much larger in the epiQTL

# Put eQTL and epiQTL data into RangedData objects
leaf.gene.rd <- with(gene.tab.leaf, RangedData(IRanges(start=begin, end=end, names=itag), space=chr.gene, itag))
photoleaf.gene.rd <- with(gene.tab.photoleaf, RangedData(IRanges(start=begin, end=end, names=itag), space=chr.gene, itag))

leaf.cor.bin.rd <- with(cor.bin.leaf.tab, RangedData(IRanges(start=start, end=end, names=bin2), space=chr, bin, bin2 ))
photoleaf.cor.bin.rd <- with(cor.bin.photoleaf.tab, RangedData(IRanges(start=start, end=end, names=bin2), space=chr, bin, bin2 ))

epi1 <- with(epi.split, RangedData(IRanges(start=as.numeric(start1), end=as.numeric(end1), names=paste0(bin,":",trait)), space=chr1, epi.qtl=bin, trait))
epi2 <- with(epi.split, RangedData(IRanges(start=as.numeric(start2), end=as.numeric(end2), names=paste0(bin,":",trait)), space=chr2, epi.qtl=bin, trait))

# find overlaps
# overlaps with gene lists
epi1.lfg.ov <- findOverlaps(query=leaf.gene.rd, subject=epi1)
epi1.lfg.ov.df <- cbind(as.data.frame(leaf.gene.rd[as.matrix(epi1.lfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi1[as.matrix(epi1.lfg.ov)[,2], ])[c(1:3,6:7)])
names(epi1.lfg.ov.df)[1] <- "chr"
names(epi1.lfg.ov.df)[1:3] <- paste0("il.gene.", names(epi1.lfg.ov.df)[1:3])
names(epi1.lfg.ov.df)[5] <- "bil.chr1"
names(epi1.lfg.ov.df)[6:7] <- paste0("bil.", names(epi1.lfg.ov.df)[6:7], "1")
head(epi1.lfg.ov.df, 10)

epi2.lfg.ov <- findOverlaps(query=leaf.gene.rd, subject=epi2)
epi2.lfg.ov.df <- cbind(as.data.frame(leaf.gene.rd[as.matrix(epi2.lfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi2[as.matrix(epi2.lfg.ov)[,2], ])[c(1:3,6:7)])
names(epi2.lfg.ov.df)[1] <- "chr"
names(epi2.lfg.ov.df)[1:3] <- paste0("il.gene.", names(epi2.lfg.ov.df)[1:3])
names(epi2.lfg.ov.df)[5] <- "bil.chr2"
names(epi2.lfg.ov.df)[6:7] <- paste0("bil.", names(epi2.lfg.ov.df)[6:7], "2")
head(epi2.lfg.ov.df, 10)

# overlaps with cor.bins
epi1.lfcb.ov <- findOverlaps(query=leaf.cor.bin.rd, subject=epi1)
epi1.lfcb.ov.df <- cbind(as.data.frame(leaf.cor.bin.rd[as.matrix(epi1.lfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi1[as.matrix(epi1.lfcb.ov)[,2], ])[c(1:3,6:7)] )
names(epi1.lfcb.ov.df)[1] <- "chr"
names(epi1.lfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi1.lfcb.ov.df)[c(1:3)])
names(epi1.lfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi1.lfcb.ov.df)[c(4:5)])
names(epi1.lfcb.ov.df)[6] <- "bil.chr1"
names(epi1.lfcb.ov.df)[7:8] <- paste0("bil.", names(epi1.lfcb.ov.df)[7:8], "1")
head(epi1.lfcb.ov.df)

epi2.lfcb.ov <- findOverlaps(query=leaf.cor.bin.rd, subject=epi2)
epi2.lfcb.ov.df <- cbind(as.data.frame(leaf.cor.bin.rd[as.matrix(epi2.lfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi2[as.matrix(epi2.lfcb.ov)[,2], ])[c(1:3,6:7)] )
names(epi2.lfcb.ov.df)[1] <- "chr"
names(epi2.lfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi2.lfcb.ov.df)[c(1:3)])
names(epi2.lfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi2.lfcb.ov.df)[c(4:5)])
names(epi2.lfcb.ov.df)[6] <- "bil.chr2"
names(epi2.lfcb.ov.df)[7:8] <- paste0("bil.", names(epi2.lfcb.ov.df)[7:8], "2")
head(epi2.lfcb.ov.df)

# find epi.qtl with overlaps on both ends
leaf.gene1.corbin2 <- intersect(epi1.lfg.ov.df$epi.qtl, epi2.lfcb.ov.df$epi.qtl)
leaf.gene2.corbin1 <- intersect(epi2.lfg.ov.df$epi.qtl, epi1.lfcb.ov.df$epi.qtl)
length(leaf.gene1.corbin2) # 17
length(leaf.gene2.corbin1) # 16

# merge into single data.frames
leaf.gene1.corbin2.df <- merge(epi1.lfg.ov.df[epi1.lfg.ov.df$epi.qtl %in% leaf.gene1.corbin2, ], epi2.lfcb.ov.df[epi2.lfcb.ov.df$epi.qtl %in% leaf.gene1.corbin2, ], by.x=8:9, by.y=9:10 )
leaf.gene1.corbin2.df$gene_corbin <- paste(leaf.gene1.corbin2.df$itag, leaf.gene1.corbin2.df$il.bin, sep=":")
head(leaf.gene1.corbin2.df)

leaf.gene2.corbin1.df <- merge(epi2.lfg.ov.df[epi2.lfg.ov.df$epi.qtl %in% leaf.gene2.corbin1, ], epi1.lfcb.ov.df[epi1.lfcb.ov.df$epi.qtl %in% leaf.gene2.corbin1, ], by.x=8:9, by.y=9:10 )
leaf.gene2.corbin1.df$gene_corbin <- paste(leaf.gene2.corbin1.df$itag, leaf.gene2.corbin1.df$il.bin, sep=":")
head(leaf.gene2.corbin1.df)

# combine both sets of results
leaf.gene.corbin.df <- rbind(leaf.gene1.corbin2.df, leaf.gene2.corbin1.df)
nrow(leaf.gene.corbin.df) # 63

# check against gene:cor.bin vector of eQTL in the leaf development gene set
intersect(leaf.gene.corbin.df$gene_corbin, leaf.gene.cor.bin) # Solyc06g007350:d.11B  ...and this one should be suspect because it's in chromosome 6!!
# [1] "Solyc06g007350:d.6B"  "Solyc06g007350:d.11B" if I *do not* exclude QTL with both ends in the same chromosome

leaf.eqtl.epiqtl.overlap <- leaf.gene.corbin.df[leaf.gene.corbin.df$gene_corbin %in% intersect(leaf.gene.corbin.df$gene_corbin, leaf.gene.cor.bin), ]

sub.eqtl.tab <- eqtl.tab
sub.eqtl.tab$itag <- substr(sub.eqtl.tab$itag, 1, 14) # trim gene names to remove gene model number
sub.eqtl.tab <- unique(sub.eqtl.tab[, c(1, 25:33)]) # keep only the information that's useful for annotating the genes

leaf.eqtl.epiqtl.overlap <- merge(leaf.eqtl.epiqtl.overlap, sub.eqtl.tab, by="itag") # annotate the w/ gene function
nrow(leaf.eqtl.epiqtl.overlap) # 4, 2 w/ chromosome filtering

# save versions including and excluding QTL with both ends on the same chromosome
# save(leaf.eqtl.epiqtl.overlap, file="leaf.eqtl.epiqtl.overlap_no.filter.Rdata" )
# write.table(leaf.eqtl.epiqtl.overlap, row.names = F, sep = "\t", quote = F, file="leaf.eqtl.epiqtl.overlap_no.filter.tsv")
# 
# save(leaf.eqtl.epiqtl.overlap, file="leaf.eqtl.epiqtl.overlap_same.chr.filter.Rdata" )
# write.table(leaf.eqtl.epiqtl.overlap, row.names = F, sep = "\t", quote = F, file="leaf.eqtl.epiqtl.overlap_same.chr.filter.tsv")

#------
# Repeat with photosynthesis + leaf modules

# find overlaps
# overlaps with gene lists
epi1.phlfg.ov <- findOverlaps(query=photoleaf.gene.rd, subject=epi1)
epi1.phlfg.ov.df <- cbind(as.data.frame(photoleaf.gene.rd[as.matrix(epi1.phlfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi1[as.matrix(epi1.phlfg.ov)[,2], ])[c(1:3,6:7)])
names(epi1.phlfg.ov.df)[1] <- "chr"
names(epi1.phlfg.ov.df)[1:3] <- paste0("il.gene.", names(epi1.phlfg.ov.df)[1:3])
names(epi1.phlfg.ov.df)[5] <- "bil.chr1"
names(epi1.phlfg.ov.df)[6:7] <- paste0("bil.", names(epi1.phlfg.ov.df)[6:7], "1")
head(epi1.phlfg.ov.df, 10)

epi2.phlfg.ov <- findOverlaps(query=photoleaf.gene.rd, subject=epi2)
epi2.phlfg.ov.df <- cbind(as.data.frame(photoleaf.gene.rd[as.matrix(epi2.phlfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi2[as.matrix(epi2.phlfg.ov)[,2], ])[c(1:3,6:7)])
names(epi2.phlfg.ov.df)[1] <- "chr"
names(epi2.phlfg.ov.df)[1:3] <- paste0("il.gene.", names(epi2.phlfg.ov.df)[1:3])
names(epi2.phlfg.ov.df)[5] <- "bil.chr2"
names(epi2.phlfg.ov.df)[6:7] <- paste0("bil.", names(epi2.phlfg.ov.df)[6:7], "2")
head(epi2.phlfg.ov.df, 10)

# overlaps with cor.bins
epi1.phlfcb.ov <- findOverlaps(query=photoleaf.cor.bin.rd, subject=epi1)
epi1.phlfcb.ov.df <- cbind(as.data.frame(photoleaf.cor.bin.rd[as.matrix(epi1.phlfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi1[as.matrix(epi1.phlfcb.ov)[,2], ])[c(1:3,6:7)] )
names(epi1.phlfcb.ov.df)[1] <- "chr"
names(epi1.phlfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi1.phlfcb.ov.df)[c(1:3)])
names(epi1.phlfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi1.phlfcb.ov.df)[c(4:5)])
names(epi1.phlfcb.ov.df)[6] <- "bil.chr1"
names(epi1.phlfcb.ov.df)[7:8] <- paste0("bil.", names(epi1.phlfcb.ov.df)[7:8], "1")
head(epi1.phlfcb.ov.df)

epi2.phlfcb.ov <- findOverlaps(query=photoleaf.cor.bin.rd, subject=epi2)
epi2.phlfcb.ov.df <- cbind(as.data.frame(photoleaf.cor.bin.rd[as.matrix(epi2.phlfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi2[as.matrix(epi2.phlfcb.ov)[,2], ])[c(1:3,6:7)] )
names(epi2.phlfcb.ov.df)[1] <- "chr"
names(epi2.phlfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi2.phlfcb.ov.df)[c(1:3)])
names(epi2.phlfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi2.phlfcb.ov.df)[c(4:5)])
names(epi2.phlfcb.ov.df)[6] <- "bil.chr2"
names(epi2.phlfcb.ov.df)[7:8] <- paste0("bil.", names(epi2.phlfcb.ov.df)[7:8], "2")
head(epi2.phlfcb.ov.df)

# find epi.qtl with overlaps on both ends
photoleaf.gene1.corbin2 <- intersect(epi1.phlfg.ov.df$epi.qtl, epi2.phlfcb.ov.df$epi.qtl)
photoleaf.gene2.corbin1 <- intersect(epi2.phlfg.ov.df$epi.qtl, epi1.phlfcb.ov.df$epi.qtl)
length(photoleaf.gene1.corbin2) # 50
length(photoleaf.gene2.corbin1) # 49

# merge into single data.frames
photoleaf.gene1.corbin2.df <- merge(epi1.phlfg.ov.df[epi1.phlfg.ov.df$epi.qtl %in% photoleaf.gene1.corbin2, ], epi2.phlfcb.ov.df[epi2.phlfcb.ov.df$epi.qtl %in% photoleaf.gene1.corbin2, ], by.x=8:9, by.y=9:10 )
photoleaf.gene1.corbin2.df$gene_corbin <- paste(photoleaf.gene1.corbin2.df$itag, photoleaf.gene1.corbin2.df$il.bin, sep=":")
head(photoleaf.gene1.corbin2.df)

photoleaf.gene2.corbin1.df <- merge(epi2.phlfg.ov.df[epi2.phlfg.ov.df$epi.qtl %in% photoleaf.gene2.corbin1, ], epi1.phlfcb.ov.df[epi1.phlfcb.ov.df$epi.qtl %in% photoleaf.gene2.corbin1, ], by.x=8:9, by.y=9:10 )
photoleaf.gene2.corbin1.df$gene_corbin <- paste(photoleaf.gene2.corbin1.df$itag, photoleaf.gene2.corbin1.df$il.bin, sep=":")
head(photoleaf.gene2.corbin1.df)

# combine both sets of results
photoleaf.gene.corbin.df <- rbind(photoleaf.gene1.corbin2.df, photoleaf.gene2.corbin1.df)
nrow(photoleaf.gene.corbin.df) # 308
intersect(photoleaf.gene.corbin.df$gene_corbin, photoleaf.gene.cor.bin) # "Solyc06g007350:d.11B" "Solyc07g065860:d.10F" "Solyc09g059240:d.8B"
# [1] "Solyc06g007350:d.6B"  "Solyc06g008500:d.6B"  "Solyc06g007590:d.6B"  "Solyc06g008160:d.6B"  "Solyc06g007350:d.11B" "Solyc07g065860:d.10F" "Solyc02g062610:d.2C" 
# [8] "Solyc09g059240:d.8B"  IF I don't exclude QTL with both ends on the same chromosome

photoleaf.eqtl.epiqtl.overlap <- photoleaf.gene.corbin.df[photoleaf.gene.corbin.df$gene_corbin %in% intersect(photoleaf.gene.corbin.df$gene_corbin, photoleaf.gene.cor.bin), ]

photoleaf.eqtl.epiqtl.overlap <- merge(photoleaf.eqtl.epiqtl.overlap, sub.eqtl.tab, by="itag") # annotate the w/ gene function
nrow(photoleaf.eqtl.epiqtl.overlap) # 15, 6 w/ chromosome filter

# save(photoleaf.eqtl.epiqtl.overlap, file="photoleaf.eqtl.epiqtl.overlap_no.filter.Rdata" )
# write.table(photoleaf.eqtl.epiqtl.overlap, row.names = F, sep = "\t", quote = F, file="photoleaf.eqtl.epiqtl.overlap_no.filter.tsv")
# 
# save(photoleaf.eqtl.epiqtl.overlap, file="photoleaf.eqtl.epiqtl.overlap_same.chr.filter.Rdata" )
# write.table(photoleaf.eqtl.epiqtl.overlap, row.names = F, sep = "\t", quote = F, file="photoleaf.eqtl.epiqtl.overlap_same.chr.filter.tsv")

#-----
# Permutations

n.perm=10000; n.epiqtl=nrow(epi.split)

registerDoParallel(cores=4) # register parallel backend
mcoptions <- list(preschedule=TRUE, set.seed=FALSE) # multi-core options

system.time(count_ovs <- foreach(i=1:n.perm, .options.multicore=mcoptions, .combine="rbind") %dopar% {
  perm.epi1 <- sample(colnames(bin.cor), size=n.epiqtl)
  perm.epi2 <- sample(colnames(bin.cor), size=n.epiqtl)
  perm.epiqtl <- paste0(perm.epi1, "_x_", perm.epi2)
  int1 <- adply(perm.epi1, 1, .id=NULL, function(x) {
    chr <- as.character(bin.stats$chr[as.character(bin.stats$bin) == x])
    tmp <- which(bin.cor[rownames(bin.cor) == x, ] > 0.90)
    start.bin <- names(tmp)[1]
    end.bin <-  names(tmp)[length(tmp)]
    bins <- paste0(start.bin, ":", end.bin)
    start <- bin.stats$bin.start[bin.stats$bin == start.bin]
    end <- bin.stats$bin.end[bin.stats$bin == end.bin]
    c(chr=chr, bins=bins, start=start, end=end)
  })
  int2 <- adply(perm.epi2, 1, .id=NULL, function(x) {
    chr <- as.character(bin.stats$chr[as.character(bin.stats$bin) == x])
    tmp <- which(bin.cor[rownames(bin.cor) == x, ] > 0.90)
    start.bin <- names(tmp)[1]
    end.bin <-  names(tmp)[length(tmp)]
    bins <- paste0(start.bin, ":", end.bin)
    start <- bin.stats$bin.start[bin.stats$bin == start.bin]
    end <- bin.stats$bin.end[bin.stats$bin == end.bin]
    c(chr=chr, bins=bins, start=start, end=end)
  })
  epi1.df <- as.data.frame(cbind(epi.qtl=perm.epiqtl, bin=perm.epi1, int1))
  epi2.df <- as.data.frame(cbind(epi.qtl=perm.epiqtl, bin=perm.epi2, int2))
  epi1.rd <- with(epi1.df, RangedData(IRanges(start=as.numeric(start), end=as.numeric(end), names=epi.qtl), space=chr, epi.qtl=epi.qtl, bin=bin, int=bins, chr=chr))
  epi2.rd <- with(epi2.df, RangedData(IRanges(start=as.numeric(start), end=as.numeric(end), names=epi.qtl), space=chr, epi.qtl=epi.qtl, bin=bin, int=bins, chr=chr))
  
  # find overlaps
  # overlaps with gene lists
  epi1.lfg.ov <- findOverlaps(query=leaf.gene.rd, subject=epi1.rd)
  epi1.lfg.ov.df <- cbind(as.data.frame(leaf.gene.rd[as.matrix(epi1.lfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi1.rd[as.matrix(epi1.lfg.ov)[,2], ])[c(1:3,6:7)])
  names(epi1.lfg.ov.df)[1] <- "chr"
  names(epi1.lfg.ov.df)[1:3] <- paste0("il.gene.", names(epi1.lfg.ov.df)[1:3])
  names(epi1.lfg.ov.df)[5] <- "bil.chr1"
  names(epi1.lfg.ov.df)[6:7] <- paste0("bil.", names(epi1.lfg.ov.df)[6:7], "1")
  
  epi2.lfg.ov <- findOverlaps(query=leaf.gene.rd, subject=epi2.rd)
  epi2.lfg.ov.df <- cbind(as.data.frame(leaf.gene.rd[as.matrix(epi2.lfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi2.rd[as.matrix(epi2.lfg.ov)[,2], ])[c(1:3,6:7)])
  names(epi2.lfg.ov.df)[1] <- "chr"
  names(epi2.lfg.ov.df)[1:3] <- paste0("il.gene.", names(epi2.lfg.ov.df)[1:3])
  names(epi2.lfg.ov.df)[5] <- "bil.chr2"
  names(epi2.lfg.ov.df)[6:7] <- paste0("bil.", names(epi2.lfg.ov.df)[6:7], "2")
  
  # overlaps with cor.bins
  epi1.lfcb.ov <- findOverlaps(query=leaf.cor.bin.rd, subject=epi1.rd)
  epi1.lfcb.ov.df <- cbind(as.data.frame(leaf.cor.bin.rd[as.matrix(epi1.lfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi1.rd[as.matrix(epi1.lfcb.ov)[,2], ])[c(1:3,6:7)] )
  names(epi1.lfcb.ov.df)[1] <- "chr"
  names(epi1.lfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi1.lfcb.ov.df)[c(1:3)])
  names(epi1.lfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi1.lfcb.ov.df)[c(4:5)])
  names(epi1.lfcb.ov.df)[6] <- "bil.chr1"
  names(epi1.lfcb.ov.df)[7:8] <- paste0("bil.", names(epi1.lfcb.ov.df)[7:8], "1")
  
  epi2.lfcb.ov <- findOverlaps(query=leaf.cor.bin.rd, subject=epi2.rd)
  epi2.lfcb.ov.df <- cbind(as.data.frame(leaf.cor.bin.rd[as.matrix(epi2.lfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi2.rd[as.matrix(epi2.lfcb.ov)[,2], ])[c(1:3,6:7)] )
  names(epi2.lfcb.ov.df)[1] <- "chr"
  names(epi2.lfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi2.lfcb.ov.df)[c(1:3)])
  names(epi2.lfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi2.lfcb.ov.df)[c(4:5)])
  names(epi2.lfcb.ov.df)[6] <- "bil.chr2"
  names(epi2.lfcb.ov.df)[7:8] <- paste0("bil.", names(epi2.lfcb.ov.df)[7:8], "2")
  
  # find epi.qtl with overlaps on both ends
  leaf.gene1.corbin2 <- intersect(epi1.lfg.ov.df$epi.qtl, epi2.lfcb.ov.df$epi.qtl)
  leaf.gene2.corbin1 <- intersect(epi2.lfg.ov.df$epi.qtl, epi1.lfcb.ov.df$epi.qtl)
  
  # merge into single data.frames
  leaf.gene1.corbin2.df <- merge(epi1.lfg.ov.df[epi1.lfg.ov.df$epi.qtl %in% leaf.gene1.corbin2, ], epi2.lfcb.ov.df[epi2.lfcb.ov.df$epi.qtl %in% leaf.gene1.corbin2, ], by.x=8:9, by.y=9:10 )
  leaf.gene1.corbin2.df$gene_corbin <- paste(leaf.gene1.corbin2.df$itag, leaf.gene1.corbin2.df$il.bin, sep=":")
  
  leaf.gene2.corbin1.df <- merge(epi2.lfg.ov.df[epi2.lfg.ov.df$epi.qtl %in% leaf.gene2.corbin1, ], epi1.lfcb.ov.df[epi1.lfcb.ov.df$epi.qtl %in% leaf.gene2.corbin1, ], by.x=8:9, by.y=9:10 )
  leaf.gene2.corbin1.df$gene_corbin <- paste(leaf.gene2.corbin1.df$itag, leaf.gene2.corbin1.df$il.bin, sep=":")
  
  # combine both sets of results
  leaf.gene.corbin.df <- rbind(leaf.gene1.corbin2.df, leaf.gene2.corbin1.df)
  
  # check against gene:cor.bin vector of eQTL in the leaf development gene set
  leaf.eqtl.epiqtl.overlap <- leaf.gene.corbin.df[leaf.gene.corbin.df$gene_corbin %in% intersect(leaf.gene.corbin.df$gene_corbin, leaf.gene.cor.bin), ]
  n.leaf <- nrow(leaf.eqtl.epiqtl.overlap) 
  
  # Repeat with photosynthesis + leaf modules
  
  # find overlaps
  # overlaps with gene lists
  epi1.phlfg.ov <- findOverlaps(query=photoleaf.gene.rd, subject=epi1.rd)
  epi1.phlfg.ov.df <- cbind(as.data.frame(photoleaf.gene.rd[as.matrix(epi1.phlfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi1.rd[as.matrix(epi1.phlfg.ov)[,2], ])[c(1:3,6:7)])
  names(epi1.phlfg.ov.df)[1] <- "chr"
  names(epi1.phlfg.ov.df)[1:3] <- paste0("il.gene.", names(epi1.phlfg.ov.df)[1:3])
  names(epi1.phlfg.ov.df)[5] <- "bil.chr1"
  names(epi1.phlfg.ov.df)[6:7] <- paste0("bil.", names(epi1.phlfg.ov.df)[6:7], "1")
  
  epi2.phlfg.ov <- findOverlaps(query=photoleaf.gene.rd, subject=epi2.rd)
  epi2.phlfg.ov.df <- cbind(as.data.frame(photoleaf.gene.rd[as.matrix(epi2.phlfg.ov)[,1], ])[c(1:3,6)], as.data.frame(epi2.rd[as.matrix(epi2.phlfg.ov)[,2], ])[c(1:3,6:7)])
  names(epi2.phlfg.ov.df)[1] <- "chr"
  names(epi2.phlfg.ov.df)[1:3] <- paste0("il.gene.", names(epi2.phlfg.ov.df)[1:3])
  names(epi2.phlfg.ov.df)[5] <- "bil.chr2"
  names(epi2.phlfg.ov.df)[6:7] <- paste0("bil.", names(epi2.phlfg.ov.df)[6:7], "2")
  
  # overlaps with cor.bins
  epi1.phlfcb.ov <- findOverlaps(query=photoleaf.cor.bin.rd, subject=epi1.rd)
  epi1.phlfcb.ov.df <- cbind(as.data.frame(photoleaf.cor.bin.rd[as.matrix(epi1.phlfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi1.rd[as.matrix(epi1.phlfcb.ov)[,2], ])[c(1:3,6:7)] )
  names(epi1.phlfcb.ov.df)[1] <- "chr"
  names(epi1.phlfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi1.phlfcb.ov.df)[c(1:3)])
  names(epi1.phlfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi1.phlfcb.ov.df)[c(4:5)])
  names(epi1.phlfcb.ov.df)[6] <- "bil.chr1"
  names(epi1.phlfcb.ov.df)[7:8] <- paste0("bil.", names(epi1.phlfcb.ov.df)[7:8], "1")
  
  epi2.phlfcb.ov <- findOverlaps(query=photoleaf.cor.bin.rd, subject=epi2.rd)
  epi2.phlfcb.ov.df <- cbind(as.data.frame(photoleaf.cor.bin.rd[as.matrix(epi2.phlfcb.ov)[,1], ])[c(1:3,6:7)], as.data.frame(epi2.rd[as.matrix(epi2.phlfcb.ov)[,2], ])[c(1:3,6:7)] )
  names(epi2.phlfcb.ov.df)[1] <- "chr"
  names(epi2.phlfcb.ov.df)[c(1:3)] <- paste0("il.corbin.", names(epi2.phlfcb.ov.df)[c(1:3)])
  names(epi2.phlfcb.ov.df)[c(4:5)] <- paste0("il.", names(epi2.phlfcb.ov.df)[c(4:5)])
  names(epi2.phlfcb.ov.df)[6] <- "bil.chr2"
  names(epi2.phlfcb.ov.df)[7:8] <- paste0("bil.", names(epi2.phlfcb.ov.df)[7:8], "2")
  
  # find epi.qtl with overlaps on both ends
  photoleaf.gene1.corbin2 <- intersect(epi1.phlfg.ov.df$epi.qtl, epi2.phlfcb.ov.df$epi.qtl)
  photoleaf.gene2.corbin1 <- intersect(epi2.phlfg.ov.df$epi.qtl, epi1.phlfcb.ov.df$epi.qtl)
  
  # merge into single data.frames
  photoleaf.gene1.corbin2.df <- merge(epi1.phlfg.ov.df[epi1.phlfg.ov.df$epi.qtl %in% photoleaf.gene1.corbin2, ], epi2.phlfcb.ov.df[epi2.phlfcb.ov.df$epi.qtl %in% photoleaf.gene1.corbin2, ], by.x=8:9, by.y=9:10 )
  photoleaf.gene1.corbin2.df$gene_corbin <- paste(photoleaf.gene1.corbin2.df$itag, photoleaf.gene1.corbin2.df$il.bin, sep=":")
  
  photoleaf.gene2.corbin1.df <- merge(epi2.phlfg.ov.df[epi2.phlfg.ov.df$epi.qtl %in% photoleaf.gene2.corbin1, ], epi1.phlfcb.ov.df[epi1.phlfcb.ov.df$epi.qtl %in% photoleaf.gene2.corbin1, ], by.x=8:9, by.y=9:10 )
  photoleaf.gene2.corbin1.df$gene_corbin <- paste(photoleaf.gene2.corbin1.df$itag, photoleaf.gene2.corbin1.df$il.bin, sep=":")
  
  # combine both sets of results
  photoleaf.gene.corbin.df <- rbind(photoleaf.gene1.corbin2.df, photoleaf.gene2.corbin1.df)
  photoleaf.eqtl.epiqtl.overlap <- photoleaf.gene.corbin.df[photoleaf.gene.corbin.df$gene_corbin %in% intersect(photoleaf.gene.corbin.df$gene_corbin, photoleaf.gene.cor.bin), ]
  n.photoleaf <- nrow(photoleaf.eqtl.epiqtl.overlap) 
  c(n.leaf=n.leaf, n.photoleaf=n.photoleaf)
})

dim(count_ovs)
head(count_ovs, 50)

save(count_ovs, file="count_ovs.Rdata")
load("count_ovs.Rdata")
length(which(count_ovs[,1] >= 2)) / 10000 # 0.0084
length(which(count_ovs[,1] >= 4)) / 10000 # 0.0028
length(which(count_ovs[,2] >= 6)) / 10000 # 0.0044
length(which(count_ovs[,2] >= 15)) / 10000 # 0

#----------
# Check if distribution of number of genes among epiQTL bins differs significantly from that of all the bins

bil.gene.bin <- read.csv("../../fine_mapping_n_eqtl_enrichment/gene_BIN_table.csv")
bil.gene.bin <- bil.gene.bin[-1]
head(bil.gene.bin)

gene.num <- ddply(bil.gene.bin, .variables="BIN", function(x) nrow(x))
names(gene.num)[2] <- "num.genes"
head(gene.num)
 
hist(gene.num$num.genes)

head(epi.split)
epi.bins <- c(epi.split$bin1, epi.split$bin2)
length(epi.bins) # 570
gene.num.epiqtl <- gene.num[gene.num$BIN %in% epi.bins, ]
nrow(gene.num.epiqtl) # 280, so a few BINs are part of more than 1 epiQTL

hist(gene.num.epiqtl$num.genes)

ks.test(gene.num$num.genes, gene.num.epiqtl$num.genes) # D = 0.068636, p-value = 0.2521 => no need to redo permutations w/ adjusted sampling probs.


