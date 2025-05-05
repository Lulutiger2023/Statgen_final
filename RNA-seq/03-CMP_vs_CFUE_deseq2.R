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
library(limma)

#CMP files
cmp_1 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/2_CMP/rna/ENCFF623OLU.tsv", header = TRUE, sep = "\t")
cmp_2 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/2_CMP/rna/ENCFF691MHW.tsv", header = TRUE, sep = "\t")

#CFUE files
cfue_1 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/3_CFUE/rna/ENCFF655LMK.tsv", header = TRUE, sep = "\t") 
cfue_2 <- read.table("/Users/ishikaverma/Documents/STAT555-Project/data/3_CFUE/rna/ENCFF667IDY.tsv", header = TRUE, sep = "\t")  

#CMP files
cmp_1_counts <- cmp_1$expected_count
cmp_2_counts <- cmp_2$expected_count

#CFUE files
cfue_1_counts <- cfue_1$expected_count
cfue_2_counts <- cfue_2$expected_count  

#countData
cmp_cfue_count_data <- data.frame(
  cmp_1 = cmp_1_counts,
  cmp_2 = cmp_2_counts,
  cfue_1 = cfue_1_counts,
  cfue_2 = cfue_2_counts,
  row.names = cmp_1$gene_id
)

#colData
samples <- c("cmp_1", "cmp_2", "cfue_1", "cfue_2")
broad_category_samples <-rep(c("cmp","cfue"),each=2)
cmp_cfue_coldata <- data.frame(row.names = samples, cell_stage = broad_category_samples)

#dds
dds_cmp_cfue <- DESeqDataSetFromMatrix(countData = round(cmp_cfue_count_data),
                                       colData = cmp_cfue_coldata, design = ~ cell_stage )
levels(dds_cmp_cfue$cell_stage)
## Pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds_cmp_cfue) >= 10) >= smallestGroupSize
dds_cmp_cfue <- dds_cmp_cfue[keep,]

#deseq
dds_cmp_cfue <- DESeq(dds_cmp_cfue)
res_cmp_cfue <- results(dds_cmp_cfue)
mcols(res_cmp_cfue)$description
summary(res_cmp_cfue)

#lfc_shrinkage
resultsNames(dds_cmp_cfue)
resLFC_cmp_cfue <- lfcShrink(dds_cmp_cfue, coef="cell_stage_cmp_vs_cfue", type="apeglm")
resLFC_cmp_cfue


#Top DE genes
sum(resLFC_cmp_cfue$padj < 0.05, na.rm=TRUE)
sum(resLFC_cmp_cfue$padj < 0.05 & resLFC_cmp_cfue$log2FoldChange > 1, na.rm = TRUE) #upregulated
sum(resLFC_cmp_cfue$padj < 0.05 & resLFC_cmp_cfue$log2FoldChange < -1, na.rm = TRUE) #Downreguated
up_deseq2 <- subset(resLFC_cmp_cfue, padj < 0.05 & log2FoldChange >= 1)
up_deseq2_genes <- rownames(up_deseq2)
length(up_deseq2_genes)
down_deseq2 <- subset(resLFC_cmp_cfue, padj < 0.05 & log2FoldChange <= -1)
down_deseq2_genes <- rownames(down_deseq2)
length(down_deseq2_genes)

#Mapping ENSEMBL Ids to gene symbols

clean_ids <- sub("\\..*", "", rownames(resLFC_cmp_cfue))
clean_ids_dds <- sub("\\..*", "", rownames(res_cmp_cfue))

gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = clean_ids,
                       column = c('SYMBOL'),
                      keytype = "ENSEMBL")
gene_symbols_dds <- mapIds(org.Mm.eg.db,
                       keys = clean_ids_dds,
                       column = c('SYMBOL'),
                       keytype = "ENSEMBL")



gene_symbols <- gene_symbols[match(clean_ids, names(gene_symbols))]
gene_symbols_dds <- gene_symbols_dds[match(clean_ids_dds, names(gene_symbols_dds))]

resLFC_cmp_cfue$symbol <- gene_symbols
res_cmp_cfue$symbol <- gene_symbols_dds


# Volcano Plot

png("CMP_vs_CFUE_volcano_deseq2_before_LFCskrinkage.png", width = 2000, height = 1600, res = 300)
EnhancedVolcano(res_cmp_cfue,
                lab = res_cmp_cfue$symbol, #not shrunken
                x = 'log2FoldChange',
                xlim = c(-5, 5),
                pCutoff = 0.05,
                FCcutoff = 2,
                title = 'Volcano plot: CMP vs CFUE',
                subtitle = 'DESeq2',
                caption = 'Thresholds: padj < 0.05, |log2FC| > 1',
                pointSize = 2.0,
                labSize = 3.0,
                colAlpha = 0.8,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                y = 'padj',
                ylab = "-Log10(padj)",
                legendLabels = c("NS", "Log2 FC", "padj", "padj and log2 FC"))

dev.off()

## MA plot
plotMA(resLFC_cmp_cfue, alpha = 0.05)
title("MA Plot: HSC vs CMP")


#Transformation
vsd_cmp_cfue <- vst(dds_cmp_cfue, blind=FALSE)
#Heatmap of count matrix
select <- order(rowMeans(counts(dds_cmp_cfue,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds_cmp_cfue)[, "cell_stage", drop = FALSE])
pheatmap(assay(vsd_cmp_cfue)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,  main = "Top 50 Expressed Genes in CMP & CFUE")


#Plot venn diagram to compare genes in DESeq2 and Limma

#Upregulated genes
up_limma_genes <- readLines("/Users/ishikaverma/Documents/STAT555-Project/up_genes_limma_lfc1.txt")
down_limma_genes <- readLines("/Users/ishikaverma/Documents/STAT555-Project/down_genes_limma_lfc1.txt")

genes_list_venn_diagram <- list(A = up_deseq2_genes, B= up_limma_genes, C = down_deseq2_genes, D = down_limma_genes)

up_plot <- ggVennDiagram(genes_list_venn_diagram[1:2]) + scale_fill_gradient(low="grey90",high = "red")+
  ggtitle("Upregulated Genes- DESeq2 vs Limma") 
up_plot


down_plot <- ggVennDiagram(genes_list_venn_diagram[3:4]) + scale_fill_gradient(low="grey90",high = "green")+
  ggtitle("Downregulated Genes- DESeq2 vs Limma")
down_plot

# Make DESeq2 vs Limma log2Foldchnage plot
load("efit.Rdata")
common_deseq2_limma_genes <- intersect(rownames(res_cmp_cfue), rownames(efit$coefficients))
logfc_df <- data.frame(
  gene = common_deseq2_limma_genes,
  log2FC_deseq2 = res_cmp_cfue[common_deseq2_limma_genes, "log2FoldChange"],
  log2FC_limma = efit$coefficients[common_deseq2_limma_genes, "CMPvsCFUE"]
)

ggplot(logfc_df, aes(x = log2FC_deseq2, y = log2FC_limma)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Comparison of log2FoldChange: DESeq2 vs Limma",
    x = "log2FC (DESeq2)",
    y = "log2FC (Limma-voom)"
  )

# Make DESeq2 vs Limma -log10(padj)

res_limma_df <- topTable(efit, coef = "CMPvsCFUE", number = Inf)

original_ids_in_deseq2 <- rownames(res_cmp_cfue)

clean_ids <- sub("\\..*", "", original_ids_in_deseq2)

res_limma_df$ensembl_clean <- rownames(res_limma_df)
res_limma_df$ensembl_full <- original_ids_in_deseq2[match(res_limma_df$ensembl_clean, clean_ids)]

padj_df <- data.frame(
  gene = common_deseq2_limma_genes,
  padj_deseq2 = res_cmp_cfue[common_deseq2_limma_genes, "padj"],
  padj_limma = res_limma_df[common_deseq2_limma_genes, "adj.P.Val"]
)

padj_df$logpadj_deseq2 <- -log10(padj_df$padj_deseq2)
padj_df$logpadj_limma <- -log10(padj_df$padj_limma)

ggplot(padj_df, aes(x = logpadj_deseq2, y = logpadj_limma )) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Comparison of -log10 of p.adjusted values: DESeq2 vs Limma",
    x = "-log10(padj) (DESeq2)",
    y = "-log10(padj) (Limma-voom)"
  )

# Savinf upregulated and downregulated genes for functional analysis

top_genes_logfc_1 <- subset(res_cmp_cfue, padj < 0.05 & abs(log2FoldChange) >= 1)

up_genes <- subset(top_genes_logfc_1, log2FoldChange >= 1)
down_genes <- subset(top_genes_logfc_1, log2FoldChange <= -1)

up_ensembl_ids <- sub("\\..*", "", rownames(up_genes))
down_ensembl_ids <- sub("\\..*", "", rownames(down_genes))

write.table(up_ensembl_ids,
            file = "up_genes_logfc1_cmp_vs_cfue_deseq2_ensembl_ids.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(down_ensembl_ids,
            file = "down_genes_logfc1_cmp_vs_cfue_deseq2_ensembl_ids.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


