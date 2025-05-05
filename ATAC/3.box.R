# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(DESeq2)
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(VennDiagram)
library(circlize)
library(enrichplot)
library(ReactomePA)
library(ComplexHeatmap)

# Get the TxDb object for mouse mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

files = list(cmp_rep1 = "ENCFF832UUS.bed",
             cmp_rep2 = "ENCFF343PTQ.bed",
             cfue_rep1 = "ENCFF796ZSB.bed",
             cfue_rep2 = "ENCFF599ZDJ.bed")

#plotAnnoBar and plotDistToTSS
peakAnnoList <- lapply(files, annotatePeak, 
                       TxDb=txdb,
                       tssRegion=c(-3000, 3000))
pdf("feature_distribution.pdf", width=10, height=8)
plotAnnoBar(peakAnnoList,title = " Feature Distribution")
dev.off()

pdf("tss_distribution.pdf", width=10, height=8)
plotDistToTSS(peakAnnoList, title = "Feature Distribution relative to TSS")
dev.off()
