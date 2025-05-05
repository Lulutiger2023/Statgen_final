library("DESeq2")
library("apeglm")
library("pheatmap")
library("RColorBrewer")
library("pheatmap")
library("vsn")
library(org.Mm.eg.db)
library(ggplot2)
library(EnhancedVolcano)
library(ggVennDiagram)
library(grid)
library(devtools)
library(enrichR)
library(viridis)

#HSC files
hsc_1 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/1_HSC/rna/ENCFF064MKY.tsv", header = TRUE, sep = "\t")
hsc_2 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/1_HSC/rna/ENCFF247FEJ.tsv", header = TRUE, sep = "\t")

#CMP files
cmp_1 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/2_CMP/rna/ENCFF623OLU.tsv", header = TRUE, sep = "\t")
cmp_2 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/2_CMP/rna/ENCFF691MHW.tsv", header = TRUE, sep = "\t")

#CFUE files
cfue_1 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/3_CFUE/rna/ENCFF655LMK.tsv", header = TRUE, sep = "\t") 
cfue_2 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/3_CFUE/rna/ENCFF667IDY.tsv", header = TRUE, sep = "\t")  

#ERY files
ery_1 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/4_ERY/rna/ENCFF342WUL.tsv", header = TRUE, sep = "\t") 
ery_2 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/4_ERY/rna/ENCFF858JHF.tsv", header = TRUE, sep = "\t") 

#HSC counts
hsc_1_counts <- hsc_1$expected_count
hsc_2_counts <- hsc_2$expected_count

#CMP counts
cmp_1_counts <- cmp_1$expected_count
cmp_2_counts <- cmp_2$expected_count

#CFUE counts
cfue_1_counts <- cfue_1$expected_count
cfue_2_counts <- cfue_2$expected_count  

#ERY counts
ery_1_counts <- ery_1$expected_count
ery_2_counts <- ery_2$expected_count

#Making DDS object
#countData
counts <- data.frame(
  hsc_1 = hsc_1_counts,
  hsc_2 = hsc_2_counts,
  cmp_1 = cmp_1_counts,
  cmp_2 = cmp_2_counts,
  cfue_1 = cfue_1_counts,
  cfue_2 = cfue_2_counts,
  ery_1 = ery_1_counts,
  ery_2 = ery_2_counts,
  row.names = cmp_1$gene_id
)

#colData
samples <- c("hsc_1", "hsc_2", "cmp_1", "cmp_2", "cfue_1", "cfue_2", "ery_1", "ery_2")
broad_category_samples <-rep(c("HSC","CMP","CFUE", "ERY"),each=2)
coldata <- data.frame(row.names = samples, cell_stage = broad_category_samples)

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData = coldata,
                                  design = ~ cell_stage)

#Pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## Performing differential gene expression analysis
dds <- DESeq(dds)
res <- results(dds)

##Transformation
vsd <- vst(dds, blind=FALSE)


##Making sample-to-sample distance matrix
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

colors <- inferno(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Hierarchial Clustering",
         labels_row = rep("", nrow(sampleDistMatrix))) 


##Making heatmap of top 500 most expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[, "cell_stage", drop = FALSE])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,  main = "Top 500 Most Expressed Genes")

##PCA
plotPCA(vsd, intgroup = "cell_stage") +
  ggtitle("PCA Plot of RNA-seq Samples by Cell Stage")

