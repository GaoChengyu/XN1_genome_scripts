library(ggtree)
library(ggtreeExtra)
library(tidytree)
library(tidyverse)
setwd("d:/data/XN1/xn1genome/CG/ortho2/")

tree<-read.tree("tree/spp_r8s_ultrametric.txt")
cafenum<-read.table("cafe/k2p/Gamma_clade_results.txt",sep = "\t",header = F)
cafenum$label<-gsub("<.+>","",cafenum$V1)
cafenum$info<-paste("+",cafenum$V2,"/","-",cafenum$V3,sep = "")%>%as.character()
cafe<-select(cafenum,c(4,5))%>%filter(label!="")

tibbletree<-tidytree::as_tibble(tree)%>%left_join(cafe)
tibbletree[is.na(tibbletree$info),][,5]<-""



tree<-as.treedata(tibbletree)

geneInfo<-read.delim("OrthoFinder/statMatrix.txt",header = T,row.names = 1)
geneInfo<-t(geneInfo)%>%as.data.frame()%>% select(c(2,3,9))
geneInfo$Non<-geneInfo[,1]-geneInfo[,3]
geneInfo$species<-row.names(geneInfo)
geneInfo<-geneInfo%>%select(c(2,3,4,5))%>%
  rename_with(~c("Not assigned","Species-species orthogroups","Non species-species orthogroups"),c(1,2,3))%>%
  pivot_longer(cols = c(1,2,3),names_to = "Type")


p<-ggtree(tree,layout = "roundrect")+
  geom_tiplab(align = T, fontface="italic")+
  geom_text(aes(label=tibbletree$info),vjust=-.3,hjust=1)+
  theme_tree2()+
  geom_treescale()

p1<-p+geom_fruit(
    data = geneInfo,
    geom = geom_col,
    mapping = aes(y=species,x=value,fill=Type),
    orientation="y",
    pwidth = 1,
    width=0.3,
    offset = 0.9,
    axis.params = list(
      axis="x",
      line.size=0.5,
      line.color="black",
      text.size=1.8,
      title="No. of Genes"))+ theme(legend.position = "top")+
  scale_x_continuous(breaks = c(0,50,100,150),expand = c(0,0),labels = c("150","100","50","0"))+
  xlab("MYA")



p2<-p1+geom_fruit(
  data =Pergenesinorthogroups,
  geom = geom_bar,
  stat="identity",
  mapping = aes(y=V1,x= as.numeric(V2)),
  fill="#F8AC8C",
  pwidth = 0.4,
  width=0.2,
  offset = 0.4,
  axis.params = list(
    axis="x",
    line.size=0.7,
    line.color="black",
    text.size=1.8,
    title="Genes in orthogroups (%)"
  )
  
)

p3<-p2+geom_fruit(
  data =peruniqgenes,
  geom = geom_bar,
  stat="identity",
  mapping = aes(y=V1,x= as.numeric(V2),fill=V3),
  pwidth = 0.4,
  width=0.2,
  offset = 0.7,
  axis.params = list(
    axis="x",
    line.size=0.7,
    line.color="black",
    text.size=1.8,
    title="Uniqe Genes (%)"
  )
  
)