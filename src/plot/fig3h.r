#####################
## repeat region associate accessory chromosome and 2-speed regions

cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
library(regioneR)
library(GenomicFeatures)

repeat.gr <- toGRanges(read.table("d:/data/XN21/xn1genome/ltr/LTR.bed", header = F))


df.core <- read.table("d:/data/XN21/xn1genome/panGenome/graphgenome/gaf/coreGene/core_graph_len.bed")
df.dis <- read.table("d:/data/XN21/xn1genome/panGenome/graphgenome/gaf/coreGene/Dispensable_graph_len.bed")

core.gr <- toGRanges(df.core)
dis.gr <- toGRanges(df.dis)

# calculate overlapped base of repeat && fast/slow regions
ovBase.core <- overlapRegions(repeat.gr, core.gr, get.bases = T)$ov.bases
# (sum(ovBase.fast)*10000)/sum(width(fast.gr))
m.core <- mean(ovBase.core)
ovBase.dis <- overlapRegions(repeat.gr, dis.gr, get.bases = T)$ov.bases
# (sum(ovBase.slow)*10000)/sum(width(slow.gr))
m.dis <- mean(ovBase.dis)

# chisq.test
p <- c(1/2, 1/2)
x <- c(m.core, m.dis)
chisq.test(x = x, p = p,correct =F)
## plot bar
cc <- c(rgb(100,186,159,maxColorValue = 255), rgb(148,57,111, maxColorValue = 255))
par(mar=c(1.5,2,5,2))
barplot(c(m.core/1000, m.dis/1000), horiz = T, space = 0.5,col = cc, las=1, xaxt="n")
legend(x = 0.4, y = 0.5, legend = c("Core", "Dispensable"), 
       col = cc, xpd = T, pch=15, bty = "n", ncol = 2)
axis(side = 3, cex=0.5)
text(x = 0.4, y=1.7, "*** X-squared p<2.2e-16")
title("Number of LTRs / 10 kb")