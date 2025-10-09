library(ComplexHeatmap)
mat_R<-read.csv('d:/data/XN1/xn1genome/DEG/TPM/geneTPMmean2.CSV',header=T,row.names = 1)
col1 <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("Biotroph",3),rep("Necrotroph",2)),
                       col = list(Stage = c("CK"="white","Biotroph"="#9467bd","Necrotroph"="#e7ba52")))

right_annotation=rowAnnotation(cluster=c("DcHET-1","DcHET-2","DcHET-1","DcHET-2","DcHET-1","DcHET-1"),
                               col = list(cluster = c("DcHET-1"="#f7e6a6","DcHET-2"="#9CD1C8")))

listde<-read.table("d:/data/XN1/xn1genome/CG/ortho2/analysis/最大扩张家族/exp1.list",header = F)
listde<-c(listde[1])
listde<-unlist(listde)
listde<-substr(listde,1,9)
mat<-mat_R[rownames(mat_R) %in% listde,]
row_indices <- match(listde, rownames(mat))  
mat<-mat[row_indices,]
mat<-as.matrix(mat)
mat=t(scale(t(mat)))
mat<-na.omit(mat)

colnames(mat)<-c("0h","12h","36h","72h","192h","264h")


Heatmap(mat, name = "z-score",
        #row_split = factor(c(rep("Biotrophic module", 2), rep("Necrotrophic module", 10),rep("Transition module", 2))),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        #cluster_row = F,
        show_row_names =F,
        col = col1,
        top_annotation = ha,
        right_annotation = right_annotation
)
