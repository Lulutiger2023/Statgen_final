library(rtracklayer)
library(GenomicRanges)
library(limma)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(EnrichedHeatmap)

# 1. Create a Venn diagram showing peak overlaps
# =============================================

# First, we need to load the peak files
cmp_rep1 <- import.bed("ENCFF832UUS.bed")
cmp_rep2 <- import.bed("ENCFF343PTQ.bed")
cfue_rep1 <- import.bed("ENCFF796ZSB.bed")
cfue_rep2 <- import.bed("ENCFF599ZDJ.bed")

# Create a named list of peak sets
peak_sets <- list(
  "CMP_Rep1" = cmp_rep1,
  "CMP_Rep2" = cmp_rep2,
  "CFUE_Rep1" = cfue_rep1,
  "CFUE_Rep2" = cfue_rep2
)

# Function to create a logical overlap matrix
create_overlap_matrix <- function(peak_list) {
  # Create a union of all peaks
  all_peaks <- Reduce(c, peak_list)
  all_peaks <- reduce(all_peaks)  # Merge overlapping peaks
  
  # Create a logical matrix showing which sets each peak in the union belongs to
  overlap_matrix <- matrix(FALSE, nrow = length(all_peaks), ncol = length(peak_list))
  colnames(overlap_matrix) <- names(peak_list)
  
  for (i in seq_along(peak_list)) {
    # Find overlaps between all_peaks and the current peak set
    overlaps <- findOverlaps(all_peaks, peak_list[[i]])
    # Mark the corresponding positions in the matrix as TRUE
    overlap_matrix[queryHits(overlaps), i] <- TRUE
  }
  
  return(overlap_matrix)
}

# Create the overlap matrix
overlap_matrix <- create_overlap_matrix(peak_sets)

# Use limma's vennDiagram function to create a Venn diagram
pdf("ATAC_peak_overlaps_venn.pdf", width=10, height=10)
vennDiagram(overlap_matrix, circle.col = brewer.pal(4, "Set1"))
dev.off()


# Additional ATAC-seq Visualizations
# This script provides additional visualizations for ATAC-seq analysis of CMP vs CFUE

# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(limma)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(EnrichedHeatmap)

# 1. Create a Venn diagram showing peak overlaps
# =============================================

# First, we need to load the peak files
cmp_rep1 <- import.bed("ENCFF832UUS.bed")
cmp_rep2 <- import.bed("ENCFF343PTQ.bed")
cfue_rep1 <- import.bed("ENCFF796ZSB.bed")
cfue_rep2 <- import.bed("ENCFF599ZDJ.bed")

# Create a named list of peak sets
peak_sets <- list(
  "CMP_Rep1" = cmp_rep1,
  "CMP_Rep2" = cmp_rep2,
  "CFUE_Rep1" = cfue_rep1,
  "CFUE_Rep2" = cfue_rep2
)

# Function to create a logical overlap matrix
create_overlap_matrix <- function(peak_list) {
  # Create a union of all peaks
  all_peaks <- Reduce(c, peak_list)
  all_peaks <- reduce(all_peaks)  # Merge overlapping peaks
  
  # Create a logical matrix showing which sets each peak in the union belongs to
  overlap_matrix <- matrix(FALSE, nrow = length(all_peaks), ncol = length(peak_list))
  colnames(overlap_matrix) <- names(peak_list)
  
  for (i in seq_along(peak_list)) {
    # Find overlaps between all_peaks and the current peak set
    overlaps <- findOverlaps(all_peaks, peak_list[[i]])
    # Mark the corresponding positions in the matrix as TRUE
    overlap_matrix[queryHits(overlaps), i] <- TRUE
  }
  
  return(overlap_matrix)
}

# Create the overlap matrix
overlap_matrix <- create_overlap_matrix(peak_sets)

# Use limma's vennDiagram function to create a Venn diagram
pdf("ATAC_peak_overlaps_venn.pdf", width=10, height=10)
vennDiagram(overlap_matrix, circle.col = brewer.pal(4, "Set1"))
dev.off()