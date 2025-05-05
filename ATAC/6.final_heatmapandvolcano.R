# ---- ATAC-seq Heatmap with Gene Symbols ----

# Load required packages
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Select significant differential peaks from ATAC-seq results
sig_atac <- atac_results_df %>%
  filter(padj < 0.05 & !is.na(padj)) %>%
  filter(!is.na(gene_symbol)) %>%  # Only keep peaks with gene annotation
  arrange(padj)

# If there are too many significant peaks, select the top 100 most significant ones
if(nrow(sig_atac) > 100) {
  sig_atac <- head(sig_atac, 100)
}

# Get normalized count data
atac_norm_counts <- counts(dds_atac, normalized = TRUE)

# Extract data needed for the heatmap
atac_heatmap_data <- atac_norm_counts[sig_atac$peak_id, ]

# Log2 transform the count data
atac_heatmap_data <- log2(atac_heatmap_data + 1)

# Create row annotation - using gene symbols
row_anno <- data.frame(
  Gene = sig_atac$gene_symbol,
  LogFC = sig_atac$log2FoldChange,
  row.names = sig_atac$peak_id
)

# Create column annotation
col_anno <- data.frame(
  Condition = c("CMP", "CMP", "CFUE", "CFUE"),
  row.names = colnames(atac_heatmap_data)
)

# Define colors
ann_colors <- list(
  Condition = c(CMP = "#377EB8", CFUE = "#E41A1C"),
  LogFC = colorRampPalette(c("blue", "white", "red"))(100)
)

# Create heatmap with specific customizations
# 1. Sort rows by LogFC to group red blocks together
# 2. Enable clustering for samples
# 3. Show gene symbols

# First, sort the data by log2FoldChange to group similar patterns
sig_atac <- sig_atac[order(sig_atac$log2FoldChange, decreasing = TRUE), ]

# Reorder the heatmap data to match the sorted sig_atac
atac_heatmap_data <- atac_norm_counts[sig_atac$peak_id, ]
atac_heatmap_data <- log2(atac_heatmap_data + 1)

# Create column annotation
col_anno <- data.frame(
  Condition = c("CMP", "CMP", "CFUE", "CFUE"),
  row.names = colnames(atac_heatmap_data)
)

# Define colors
ann_colors <- list(
  Condition = c(CMP = "#377EB8", CFUE = "#E41A1C")
)

# Create heatmap with customized settings
pheatmap(
  atac_heatmap_data,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,  # Disable row clustering to maintain the LogFC-based order
  cluster_cols = TRUE,   # Enable column clustering for samples
  labels_row = sig_atac$gene_symbol,  # Use gene symbols as row labels
  fontsize_row = 8,
  fontsize_col = 10,
  main = "ATAC-seq Differential Peaks Heatmap",
  filename = "ATAC_diff_peaks_heatmap.png",
  width = 10,
  height = 12
)

# Save as PDF
pheatmap(
  atac_heatmap_data,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,  # Disable row clustering
  cluster_cols = TRUE,   # Enable column clustering for samples
  labels_row = sig_atac$gene_symbol,  # Use gene symbols as row labels
  fontsize_row = 8,
  fontsize_col = 10,
  main = "ATAC-seq Differential Peaks Heatmap",
  filename = "ATAC_diff_peaks_heatmap.pdf",
  width = 10,
  height = 12
)

# ---- ATAC-seq Volcano Plot ----

# Add status column to the results dataframe for visualization
atac_results_df$status <- "Not Sig"
atac_results_df$status[atac_results_df$padj < 0.05 & !is.na(atac_results_df$padj) & atac_results_df$log2FoldChange > 0] <- "Up"
atac_results_df$status[atac_results_df$padj < 0.05 & !is.na(atac_results_df$padj) & atac_results_df$log2FoldChange < 0] <- "Down"

# Define colors
colors <- c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")

# Display the number of significant differential peaks
cat("Differential accessible regions (FDR < 0.05):", sum(atac_results_df$status != "Not Sig"), "\n")
cat("Upregulated in CFUE:", sum(atac_results_df$status == "Up"), "\n")
cat("Downregulated in CFUE:", sum(atac_results_df$status == "Down"), "\n")

# Select the top 20 most significant peaks for labeling (based on p-value and fold change)
top_peaks <- atac_results_df %>%
  filter(padj < 0.05 & !is.na(padj) & !is.na(gene_symbol)) %>%
  mutate(abs_log2FC = abs(log2FoldChange)) %>%
  arrange(padj, desc(abs_log2FC)) %>%
  head(20)

# Create volcano plot
volcano_plot <- ggplot(atac_results_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_color_manual(values = colors) +
  geom_text_repel(data = top_peaks, 
                  aes(label = gene_symbol),
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) +
  labs(title = "Volcano Plot: CFUE vs CMP ATAC-seq Peaks",
       x = "Log2 Fold Change (CFUE/CMP)",
       y = "-log10(Adjusted p-value)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())

# Save volcano plot
ggsave("ATAC_volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)
ggsave("ATAC_volcano_plot.pdf", volcano_plot, width = 10, height = 8)

# Output completion information
cat("ATAC-seq heatmap and volcano plot have been generated.\n")

