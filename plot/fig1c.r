library(circlize)
library(Rsamtools)
library(GenomicFeatures)
## prepare data
xn1.fai <- read.table("d:/data/xn1/xn1genome/ChrInfo/XN1.Chr.fasta.fai")
canu.gs =  xn1.fai$V2





chro <- xn1.fai$V1
starts = rep(0, length(canu.gs))
ends = canu.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])

#### start plot######################################### 
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree=c(rep(1,14), 5))
circos.genomicInitialize(data = genoCir[1:15,],
                         sector.names = chro,
                         labels.cex = 0.5, track.height = 0.05, plotType = "labels")
### a: ideagram of 15 Chrs
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(chro[1:15] == sector.index)
  circos.axis(labels.cex = 0.5,direction = "outside", labels.niceFacing = T, labels = "", minor.ticks = 5, lwd = 0.8, 
              major.at = c(0, canu.gs[idx]), major.tick.length = 0.4)
}, track.height = 0.05, bg.border = NA)

### b: plot gap and telo 
gap.bed<-read.table( "d:/data/xn1/xn1genome/ChrInfo/genomeGap.txt",header = F)
colnames(gap.bed)<-c("chr","start","end")
gap.bed$value<-gap.bed$end-gap.bed$start
#awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4}' chr.telo.txt|awk 'BEGIN{FS=OFS="\t"}{for (f=1; f <= NF; f+=1) {if (f !=1 && $f != "nd") {print $1,$f}}}' > telo.bed
telo.bed<-read.table( "d:/data/xn1/xn1genome/ChrInfo/telo.txt",header = F)
telo.bed1<-telo.bed[c(1,3)]
colnames(telo.bed1)<-c("chr","start")
telo.bed2<-telo.bed[c(1,4)]
colnames(telo.bed2)<-c("chr","start")
telo.bed<-rbind(telo.bed1,telo.bed2)
telo.bed<-drop_na(telo.bed)

telo.bed$end<-telo.bed$start
telo.bed$value<-1

all.bed.list<-list(gap.bed,telo.bed)


circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.03,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                circos.genomicRect(region, value, col = "black", border = "black", lwd=1)}
                                else{
                                  circos.genomicRect(region, value, col = "red", border = "red", lwd=1)
                                }
                                })

### c: reads mean coverage barplot 
read.df <- read.table("d:/data/xn1/pacbio/assembleQuality/GAP_fill/combine/minCov.txt.txt",header = F)
read.df<-read.df[,1:4]
colnames(read.df) <- c("chr", "start", "end", "value")
read.df$value <- log2(read.df$value + 1)
circos.genomicTrackPlotRegion(read.df,track.height=0.08, bg.border=NA, ylim = c(min(read.df$value), max(read.df$value)),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ybottom = 0, ytop.column = 1, 
                                                    col = "SlateBlue1", border = NA, lwd = 0.5, area = T)
                              })


####plot gene heatmap
txdb <- makeTxDbFromGFF("d:/data/xn1/xn1genome/ChrInfo/XN1.genome.gff", format = "gff")
ge <- genes(txdb)
ge.df <- data.frame(chr = seqnames(ge), start= start(ge), end=end(ge))%>%filter(chr!="ChrUnanchored")
## calculate gene density in slide windows
bed.gene <- genomicDensity(region = ge.df, window.size = 10000, overlap = F)
cc1 = colorRamp2(breaks = c(0, 1), colors = c("white","#E79397"))
circos.genomicTrackPlotRegion(bed.gene, track.height=0.08, bg.border=NA,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = cc1(value),border = NA, lwd=0.1)
                              })

####plot LTR
ltr<-read.delim("d:/data/xn1/xn1genome/ltr/XN1.Chr.fasta.out.gff3",header = F)
Copia_LTR_retrotransposon<-filter(ltr,V3=="Copia_LTR_retrotransposon")[c(1,4,5)]
Copia.df <- data.frame(chr=Copia_LTR_retrotransposon$V1, start=Copia_LTR_retrotransposon$V4, end=Copia_LTR_retrotransposon$V5)
bed.Copia <- genomicDensity(region = Copia.df, window.size = 10000, overlap = F)###circos.clear()

circos.genomicTrackPlotRegion(bed.Copia, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){

                              circos.genomicLines(region, value, col = "#f09640", border = NA,lwd = 0.02,
                                                     ybottom = 0, ytop.column = 1)
                              })

Gypsy_LTR_retrotransposon<-filter(ltr,V3=="Gypsy_LTR_retrotransposon")[c(1,4,5)]
Gypsy.df <- data.frame(chr=Gypsy_LTR_retrotransposon$V1, start=Gypsy_LTR_retrotransposon$V4, end=Gypsy_LTR_retrotransposon$V5)
bed.Gypsy <- genomicDensity(region = Gypsy.df, window.size = 10000, overlap = F)
vvs <- bed.Gypsy$value

circos.genomicTrackPlotRegion(bed.Gypsy, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                
                                circos.genomicLines(region, value, col = "#23749f", border = NA,lwd = 0.02,
                                                    ybottom = 0, ytop.column = 1)
                              })


### e: plot GC content
## GC content
GC = 0.4389398933180977

gc_content <- read.delim('d:/data/xn1/xn1genome/ChrInfo/binGC.txt',header = F)

bed.gc <- data.frame(chr = gc_content$V1, start = gc_content$V2, 
                     end = gc_content$V3,value=gc_content$V4-GC)


vv <- bed.gc$value
vv[which(vv %in% boxplot.stats(vv)$out)] <- mean(vv) 
bed.gc$value <- vv
bed.gc.list <- list(bed.gc[bed.gc$value>0, ], 
                    bed.gc[bed.gc$value<0, ])
circos.genomicTrackPlotRegion(bed.gc.list, track.height = 0.08, bg.border = NA, ylim = c(min(vv),max(vv)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  # circos.genomicLines(region, value, col = "forestgreen", border = NA,
                                  #                     baseline = 0,area = T)
                                  circos.genomicRect(region, value, col = "forestgreen", border = NA,
                                                     ybottom = 0, ytop.column = 1)
                                }else{
                                  # circos.genomicLines(region, value, col = "firebrick1", border = NA,
                                  #                     baseline = 0,area = T)
                                  circos.genomicRect(region, value, col = "firebrick1", border = NA,
                                                     ybottom = 0, ytop.column = 1)
                                }})
