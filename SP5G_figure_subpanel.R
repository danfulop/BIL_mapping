# Subpanel plot for SP5G figure showing genes within the QTL interval

library(ggplot2)

setwd("~/UCD/BILs/")
region <- read.csv("SP5G_region.csv")
head(region)

# which gene is SELF PRUNING 5G?
grep("SELF", region$Human.readable.description) #17
region$color <- "white"
region$color[17] <- "black"
region$begin[17] + ((region$end[17] - region$begin[17]) / 2 ) # 63037316
region$begin.Mbp <- region$begin / 1000000
region$end.Mbp <- region$end / 1000000
region$begin.Mbp[17] + ((region$end.Mbp[17] - region$begin.Mbp[17]) / 2 ) # 63.03732

# ggplot(region, aes(xmin=begin, xmax=end, fill=color)) + geom_rect(ymin=0, ymax=1, color="black") + scale_fill_manual(values=c("black", "white")) + 
#   theme_classic(16) + ylim(0,1.022) + geom_text(aes(x=63037316, y=1.022), label="SP5G", size=7, parse=T) +
#   theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(), axis.title=element_blank()) +
#   scale_x_reverse(labels = scales::comma)

region$lab <- ""
region$lab[17] <- "SP5G"
region$mid.Mbp <- region$begin.Mbp + ( (region$end.Mbp - region$begin.Mbp) / 2 )
region$y <- 1.05
head(region)

sp5g.plot <- ggplot(region, aes(xmin=begin.Mbp, xmax=end.Mbp, x=mid.Mbp, y=y, fill=color)) + geom_rect(ymin=0, ymax=1, color="black") + 
  scale_fill_manual(values=c("black", "white")) + theme_classic(54) + ylim(0,2) + 
  theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(), axis.title=element_blank()) +
  scale_x_reverse()
# geom_text(aes(label=lab), size=24, parse=T, vjust="bottom")
# + geom_text(x=63.03732, y=1.05, label="SP5G", size=24, parse=T, vjust="bottom")

ggsave("sp5g.plot.svg", sp5g.plot, "svg", width=250/10, height=60/10, units="in", limitsize=F)
ggsave("sp5g.plot.pdf", sp5g.plot, "pdf", width=250/10, height=60/10, units="in", limitsize=F)

sp5g.plot <- ggplot(region, aes(xmin=begin.Mbp, xmax=end.Mbp, fill=color)) + geom_rect(ymin=0, ymax=1, color="black", size=0.5) + scale_fill_manual(values=c("black", "white")) + 
  theme_classic(6) + ylim(0,1.6) + geom_text(aes(x=63.03732, y=1.2, label=lab), size=3, parse=T) +
  theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(), axis.title=element_blank()) +
  scale_x_reverse(labels = scales::comma)
ggsave("sp5g.plot_small.svg", sp5g.plot, "svg", width=250/80, height=40/80, units="in", limitsize=F)