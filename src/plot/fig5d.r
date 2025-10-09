library(ggtree)
library(ggnewscale)
library(tidyverse)
library(ggtreeExtra)
library(jjPlot)



tree<-read.tree("all_class_name.txt.tree")
tree <- groupClade(tree, .node = c("Ascomycota","Basidiomycota"))


annotation<-read.csv("num.CSV",header = T)

mat<-read.delim("utral_stat.tsv",header = T,row.names = 1)
mat<-as.data.frame(t(mat))
colnames(mat)<-c(seq(1:(ncol(mat)-1)),"species_sum")
mat$name<-rownames(mat)

matdf<-pivot_longer(mat,cols = c(-ncol(mat),-(ncol(mat)-1)),names_to = "type")
matdf$Presence<-matdf$value/matdf$species_sum
matdf$Absence<-1-matdf$value/matdf$species_sum
matdf<-matdf[order(matdf$type),]
matdf$seq<-seq(1:nrow(matdf))

matdf<-pivot_longer(matdf,cols = c(5,6),names_to = "ngroup",values_to = "nvalue")
matdf$type<-factor(matdf$type,levels = seq(1:23))

p<-ggtree(tree, aes(color=group), branch.length = 'none', 
          layout = "roundrect",
          ladderize=F,size=0.1) + 
  geom_tiplab(fontface = "bold.italic",size=3.5) + 
  #geom_nodelab(aes(subset=!isTip,label=node))+
  scale_color_manual(values=c("black","orange","darkgreen"))

p<-ggtree::rotate(p,14)


p1<-p+
  geom_fruit(
    data=annotation,
    geom=geom_text,
    mapping=aes(y=name, x=num,label=paste(num,sep = "")), 
    angle=0.1,
    size=3,  
    pwidth=0.2,
    offset = 4
    
  )

p2 <- p1 + 
  new_scale_fill() +  
  geom_fruit(
    data=matdf,
    geom=geom_jjPointPie,  
    mapping=aes(y=name, x=type, group=seq,fill=ngroup, pievar=nvalue),  
    color="white",
    line.size=0.01,
    width=0.3,
    pwidth=12,
    offset = 12
    
  ) + 
  scale_fill_manual(
    values=c(Presence="#3e6da5", Absence="#bfbebe"), 
    limits=c("Presence","Absence"),
    name="" 
  )+
  guides(color="none")+ 
  theme(
    legend.justification = c("right", "top") 
  )