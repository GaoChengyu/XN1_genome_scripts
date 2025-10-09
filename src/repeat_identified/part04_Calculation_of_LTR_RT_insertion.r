####Calculation of LTR-RT insertion
library(ape)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)
library(systemPipeR)

ltr <- read.table("d:/data/XN21/xn1genome/ltr/XN1.Chr.fa.pass.list")
genome <- "d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fa"

LTR_insert <- function(ltr, genome){
  ltr.df <- ltr
  a1 <- strsplit(ltr.df$V1, split="[.:]", perl=T)
  chrom <- unlist(lapply(a1, function(x) x[1]))
  s_lLTR <- as.numeric(unlist(lapply(a1, function(x) x[2])))
  e_rLTR <- as.numeric(unlist(lapply(a1, function(x) x[4])))
  
  a2 <- strsplit(ltr.df$V7, split="[.:]", perl=T)
  e_lLTR <- as.numeric(unlist(lapply(a2, function(x) x[2]))) - 1
  s_rLTR <- as.numeric(unlist(lapply(a2, function(x) x[4]))) + 1
  
  genome.fasta <- Rsamtools::FaFile(genome)
  l.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_lLTR, end = e_lLTR))
  r.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_rLTR, end = e_rLTR))
  
  l.ltr.seq <- getSeq(genome.fasta, l.ltr)
  r.ltr.seq <- getSeq(genome.fasta, r.ltr)
  
  pair_align <- pairwiseAlignment(l.ltr.seq, r.ltr.seq)
  seq1 <- as.character(pattern(pair_align))
  seq2 <- as.character(subject(pair_align))
  res <- list()
  # res2 <- list()
  for (i in 1:length(seq1)) {
    x <- list(strsplit(seq1, "")[[i]], strsplit(seq2, "")[[i]])
    xx <- as.DNAbin(x)
    d <- dist.dna(xx)[1]
    # res2[[i]] <- d
    # fungal substitution rate  nucleotides per site per year
    r <- 1.05*10e-9 
    divergent_time <- d/(2*r)
    res[[i]] <- divergent_time/1000000
  }
  rr1 <- unlist(res)
  return(rr1)
}
ydj.ltr.res <- density(LTR_insert(ltr, genome))
ydj.ltr.res$x[which.max(ydj.ltr.res$y)]
options(scipen = 200)
p3 <- ~plot(ydj.ltr.res, main = "Time of insertion of LTR",col="purple", 
            lwd=1.5, xlab="Age (MYA)",ylim=c(0,0.5), xlim=c(0,3)
)
abline(v=ydj.ltr.res$x[which.max(ydj.ltr.res$y)], col = "red")
cowplot::plot_grid(plotlist = list(p3))