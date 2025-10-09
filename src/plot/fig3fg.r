##########
###### SV pie plot 
del = 1631
dup = 561
ins = 1021
pie.data <- data.frame(group = c("DEL", "DUP", "INS"), 
                       ratio=c(del, dup, ins))
aa <- pie.data$ratio/sum(pie.data$ratio)
labs <- paste0(round(aa,3)*100,"%")
p <- ggpubr::ggpie(pie.data, x = "ratio", label = labs[c(1,2,3)], fill = "group", lab.pos = "in",
                   color="white", palette = RColorBrewer::brewer.pal(6, "Set2"), legend.title="")
pp <- ggpubr::ggpar(p, legend = "right", tickslab = F)


################################
# count SV in different regions
sv_stats <- data.frame(Region=c("Intergenic", "Intron", "Exon", "CDS", "Promoter"),
                       Count=c(24,105,40, 21,90))


region_colors <- c(
  "Intergenic" = "#4da0a0",
  "Intron" = "#6d8a8a",
  "Exon" = "#9b3a74",
  "CDS" = "#d14a8e",
  "Promoter" = "#ff7f00"
)

sv_stats$Region <- factor(sv_stats$Region, levels = c("Intergenic", "Intron", "Exon", "CDS", "Promoter"))


ggplot(sv_stats, aes(x = Region, y = Count, fill = Region)) +
  geom_bar(stat = "identity",width = 0.6) +
  scale_fill_manual(values = region_colors) +
  labs(title = "SV Distribution Across Genomic Regions",
       x = "Genomic Region",
       y = "Number of SVs / Mb") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
