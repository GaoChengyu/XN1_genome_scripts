### permutation test of the location correlation secretome to closest repeat
library(GenometriCorr)
library(regioneR)
library(tidyverse)

repeat.gr <- toGRanges(read.table("d:/data/XN21/xn1genome/ltr/LTR.bed", header = F))

fast.repeat.gr<-overlapRegions(repeat.gr, fast.gr, get.bases = T)%>%
  select(1,2,3)%>%distinct()%>%toGRanges()

slow.repeat.gr<-overlapRegions(repeat.gr, slow.gr, get.bases = T)%>%
  select(1,2,3)%>%distinct()%>%toGRanges()

gene1.gr <-toGRanges(read.table("d:/data/XN21/xn1genome/panGenome/twospeed/result/fast.speed.gene.bed", header = F))

gs <- read.table("d:/data/XN21/xn1genome/centro/chrom.sizes", header = F)
gs.len <- gs$V2

### test

geCor <- GenometriCorr::GenometriCorrelation(query = gene1.gr, reference = fast.repeat.gr,
                                             permut.number = 100, keep.distributions = TRUE,
                                             chromosomes.to.proceed = gs$V1, chromosomes.length = gs.len)
GenometriCorr::graphical.report(geCor, pdffile = "d:/data/XN21/xn1genome/panGenome/twospeed/result/repeatPermutationTest_fast.SM.pdf")
