# RNA-ATAC Integration Analysis using DESeq2
# This script integrates RNA-seq and ATAC-seq data using DESeq2

# Load necessary packages
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(grid)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# ---- PART 1: RNA-seq Differential Expression Analysis using DESeq2 ----

# Read RNA-seq counts data - updated file name
rna_counts <- read.csv("12_hsc_cmp_count_data.csv", row.names = 1)
cat("RNA-seq data dimensions:", dim(rna_counts), "\n")

# Create sample information for RNA-seq
rna_coldata <- data.frame(
  condition = factor(c("CMP", "CMP", "HSC", "HSC")),
  row.names = colnames(rna_counts)
)

# Check if RNA counts are integers
if(all(round(as.matrix(rna_counts)) == as.matrix(rna_counts))) {
  cat("RNA counts are already integers, using directly with DESeq2\n")
  rna_counts_int <- rna_counts
} else {
  cat("RNA counts contain non-integers, rounding for DESeq2 compatibility\n")
  # Round the counts to integers as DESeq2 requires integer counts
  rna_counts_int <- round(rna_counts)
}

# Create DESeq2 dataset for RNA-seq
dds_rna <- DESeqDataSetFromMatrix(countData = rna_counts_int,
                                  colData = rna_coldata,
                                  design = ~ condition)

# Filter low count genes
keep_rna <- rowSums(counts(dds_rna) >= 10) >= 2
dds_rna <- dds_rna[keep_rna,]
cat("Genes retained after filtering:", nrow(dds_rna), "\n")

# Run DESeq2 analysis for RNA-seq
dds_rna <- DESeq(dds_rna)

# Get results for HSC vs CMP comparison
rna_results <- results(dds_rna, contrast = c("condition", "HSC", "CMP"))

# Order results by adjusted p-value
rna_results <- rna_results[order(rna_results$padj),]

# Convert results to data frame for easier manipulation
rna_results_df <- as.data.frame(rna_results)
rna_results_df$gene_id <- rownames(rna_results_df)

# Print number of DEGs
deg_count <- sum(rna_results_df$padj < 0.05, na.rm = TRUE)
cat("Differentially expressed genes (FDR < 0.05):", deg_count, "\n")

# ---- PART 2: ATAC-seq Differential Accessibility Analysis using DESeq2 ----

# Read ATAC-seq counts data 
atac_counts <- read.csv("ATAC_counts_matrix_signal_based.csv", row.names = 1)
cat("ATAC-seq data dimensions:", dim(atac_counts), "\n")

# Create sample information for ATAC-seq
atac_coldata <- data.frame(
  condition = factor(c("CMP", "CMP", "HSC", "HSC")),
  row.names = colnames(atac_counts)
)

# Check if ATAC counts are integers
if(all(round(as.matrix(atac_counts)) == as.matrix(atac_counts))) {
  cat("ATAC counts are already integers, using directly with DESeq2\n")
  atac_counts_int <- atac_counts
} else {
  cat("ATAC counts contain non-integers, rounding for DESeq2 compatibility\n")
  # Round the counts to integers as DESeq2 requires integer counts
  atac_counts_int <- round(atac_counts)
}

# Create DESeq2 dataset for ATAC-seq
dds_atac <- DESeqDataSetFromMatrix(countData = atac_counts_int,
                                   colData = atac_coldata,
                                   design = ~ condition)

# Filter low count peaks
keep_atac <- rowSums(counts(dds_atac) >= 10) >= 2
dds_atac <- dds_atac[keep_atac,]
cat("ATAC peaks retained after filtering:", nrow(dds_atac), "\n")

# Run DESeq2 analysis for ATAC-seq
dds_atac <- DESeq(dds_atac)

# Get results for HSC vs CMP comparison
atac_results <- results(dds_atac, contrast = c("condition", "HSC", "CMP"))

# Order results by adjusted p-value
atac_results <- atac_results[order(atac_results$padj),]

# Convert results to data frame for easier manipulation
atac_results_df <- as.data.frame(atac_results)
atac_results_df$peak_id <- rownames(atac_results_df)

# Print number of DARs
dar_count <- sum(atac_results_df$padj < 0.05, na.rm = TRUE)
cat("Differentially accessible regions (FDR < 0.05):", dar_count, "\n")

# ---- PART 3: Integration of RNA and ATAC data by mapping peaks to genes ----
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)

# Function to convert peak IDs to GRanges
peakids_to_granges <- function(peak_ids) {
  # Extract chromosome, start, and end from peak IDs (chr_start_end format)
  # This handles format like "chr11_84823610_84823759"
  peak_parts <- strsplit(as.character(peak_ids), "_")
  
  # Create data frame with extracted information
  peaks_df <- data.frame(
    chr = sapply(peak_parts, function(x) x[1]),
    start = as.numeric(sapply(peak_parts, function(x) x[2])),
    end = as.numeric(sapply(peak_parts, function(x) x[3])),
    name = peak_ids
  )
  
  # Check extraction worked correctly
  cat("Example peak parsing:\n")
  print(head(peaks_df, 3))
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = peaks_df$chr,
    ranges = IRanges(start = peaks_df$start, end = peaks_df$end),
    name = peaks_df$name
  )
  
  return(gr)
}

# Create GRanges object from ATAC peak IDs
cat("Converting ATAC peak IDs to genomic ranges...\n")
atac_peaks_gr <- peakids_to_granges(rownames(atac_results_df))

# Use ChIPseeker to annotate peaks with gene information
cat("Annotating ATAC peaks with gene information...\n")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
atac_peak_anno <- annotatePeak(atac_peaks_gr, 
                               tssRegion = c(-3000, 3000),
                               TxDb = txdb, 
                               annoDb = "org.Mm.eg.db")

# Convert annotation to data frame
anno_df <- as.data.frame(atac_peak_anno)

# Create mapping from peak ID to ENSEMBL gene ID
peak_to_gene <- data.frame(
  peak_id = anno_df$name,
  entrez_id = anno_df$geneId,
  ensembl_id = NA,
  symbol = anno_df$SYMBOL,
  distance_to_tss = anno_df$distanceToTSS,
  annotation = anno_df$annotation,
  stringsAsFactors = FALSE
)

# Create mapping between Entrez ID and Ensembl ID
# Map Entrez IDs to Ensembl IDs
cat("Creating mapping between Entrez IDs and Ensembl IDs...\n")
entrez_ids <- unique(na.omit(peak_to_gene$entrez_id))

if(length(entrez_ids) > 0) {
  # Use org.Mm.eg.db to map Entrez IDs to Ensembl IDs
  ensembl_map <- AnnotationDbi::select(org.Mm.eg.db, 
                                       keys = entrez_ids,
                                       columns = c("ENSEMBL"), 
                                       keytype = "ENTREZID")
  
  # Create a mapping dictionary
  ensembl_dict <- setNames(ensembl_map$ENSEMBL, ensembl_map$ENTREZID)
  
  # Update the peak_to_gene data frame with Ensembl IDs
  peak_to_gene$ensembl_id <- ensembl_dict[peak_to_gene$entrez_id]
  
  # For peaks with no gene assignment, keep NA
  cat("Total peaks with Ensembl gene assignments:", sum(!is.na(peak_to_gene$ensembl_id)), 
      "out of", nrow(peak_to_gene), "\n")
} else {
  cat("No valid Entrez IDs found for mapping to Ensembl IDs\n")
}

# Keep only peaks with valid Ensembl gene assignments
peak_to_gene_valid <- peak_to_gene[!is.na(peak_to_gene$ensembl_id), ]

# Prioritize promoter peaks over distal peaks when multiple peaks map to the same gene
# This is important when multiple peaks are associated with the same gene
peak_to_gene_valid$is_promoter <- grepl("Promoter", peak_to_gene_valid$annotation)
peak_to_gene_valid$abs_distance <- abs(peak_to_gene_valid$distance_to_tss)

# Group by gene and sort by is_promoter (TRUE first) and then by distance to TSS
peak_to_gene_valid <- peak_to_gene_valid[order(peak_to_gene_valid$ensembl_id, 
                                               !peak_to_gene_valid$is_promoter,
                                               peak_to_gene_valid$abs_distance), ]

# Keep the top peak (closest to promoter) for each gene
peak_to_gene_unique <- peak_to_gene_valid[!duplicated(peak_to_gene_valid$ensembl_id), ]

# Add gene symbols to atac_results_df
atac_results_df$ensembl_id <- NA
atac_results_df$gene_symbol <- NA

# Create dictionaries for mapping
peak_to_ensembl <- setNames(peak_to_gene_unique$ensembl_id, peak_to_gene_unique$peak_id)
peak_to_symbol <- setNames(peak_to_gene_unique$symbol, peak_to_gene_unique$peak_id)

# Update atac_results_df with Ensembl IDs and gene symbols
matched_indices <- match(atac_results_df$peak_id, names(peak_to_ensembl))
atac_results_df$ensembl_id <- peak_to_ensembl[atac_results_df$peak_id]
atac_results_df$gene_symbol <- peak_to_symbol[atac_results_df$peak_id]

# Save the peak-to-gene mapping
write.csv(peak_to_gene, "ATAC_peak_to_gene_mapping.csv", row.names = FALSE)

# Get significant DEGs and DARs
sig_degs <- rownames(rna_results_df[rna_results_df$padj < 0.05 & !is.na(rna_results_df$padj), ])
sig_dars <- rownames(atac_results_df[atac_results_df$padj < 0.05 & !is.na(atac_results_df$padj), ])

# Get genes associated with significant DARs
sig_dar_genes <- na.omit(unique(atac_results_df$ensembl_id[atac_results_df$padj < 0.05 & !is.na(atac_results_df$padj)]))

# ---- Create Correlation Plot ----

# Display sample RNA-seq IDs
cat("Sample RNA-seq row names (first 5):\n")
print(head(rownames(rna_results_df), 5))

# Display sample ATAC-seq mapped Ensembl IDs
cat("\nSample ATAC-seq mapped Ensembl IDs (first 5 non-NA):\n")
print(head(na.omit(atac_results_df$ensembl_id), 5))

# Filter for only those with ENSEMBL gene IDs mapped from ATAC-seq
mapped_atac <- atac_results_df %>% 
  filter(!is.na(ensembl_id))

# Strip version numbers from RNA-seq Ensembl IDs in row names
rna_ids_no_version <- sub("\\.[0-9]+$", "", rownames(rna_results_df))
# Create a mapping from unversioned to versioned IDs
rna_id_map <- setNames(rownames(rna_results_df), rna_ids_no_version)

# Strip version numbers from ATAC-seq mapped Ensembl IDs if they exist
mapped_atac$ensembl_id_no_version <- sub("\\.[0-9]+$", "", mapped_atac$ensembl_id)

# Find matching genes between ATAC and RNA
common_genes <- intersect(mapped_atac$ensembl_id_no_version, rna_ids_no_version)
cat("Found", length(common_genes), "genes in common between RNA-seq and ATAC-seq datasets after removing version numbers\n")

# Filter data for matching genes
matched_atac <- mapped_atac[mapped_atac$ensembl_id_no_version %in% common_genes, ]

# Create data frame for correlation plot
correlation_data <- data.frame(
  ensembl_id = matched_atac$ensembl_id,
  ensembl_id_no_version = matched_atac$ensembl_id_no_version,
  gene_symbol = matched_atac$gene_symbol,
  peak_id = matched_atac$peak_id,
  ATAC_logFC = matched_atac$log2FoldChange,
  ATAC_padj = matched_atac$padj
)

# Add RNA data - need to match RNA row names using the unversioned IDs
correlation_data$RNA_logFC <- NA
correlation_data$RNA_padj <- NA

for (i in 1:nrow(correlation_data)) {
  unversioned_id <- correlation_data$ensembl_id_no_version[i]
  if (unversioned_id %in% rna_ids_no_version) {
    versioned_id <- rna_id_map[unversioned_id]
    correlation_data$RNA_logFC[i] <- rna_results_df[versioned_id, "log2FoldChange"]
    correlation_data$RNA_padj[i] <- rna_results_df[versioned_id, "padj"]
  }
}

# Filter for significant changes in both datasets
sig_corr_data <- correlation_data %>%
  filter(ATAC_padj < 0.05 & RNA_padj < 0.05 & !is.na(RNA_padj) & !is.na(ATAC_padj))

# Report the number of significant correlations
cat("Found", nrow(sig_corr_data), "genes with significant changes in both RNA-seq and ATAC-seq data\n")

# If no significant correlations were found, skip the correlation analysis
if(nrow(sig_corr_data) == 0) {
  cat("No genes with significant changes in both datasets. Skipping correlation analysis.\n")
} else {
  # Calculate correlation
  correlation <- cor(sig_corr_data$RNA_logFC, sig_corr_data$ATAC_logFC, 
                     method = "pearson", use = "complete.obs")
  r_squared <- correlation^2
  
  # Perform correlation test
  cor_test <- cor.test(sig_corr_data$RNA_logFC, sig_corr_data$ATAC_logFC)
  p_value <- cor_test$p.value
  
  # Create a quadrant-based color scheme
  sig_corr_data$quadrant <- "Not Significant"
  sig_corr_data$quadrant[sig_corr_data$RNA_logFC > 0 & sig_corr_data$ATAC_logFC > 0] <- "Up in Both"
  sig_corr_data$quadrant[sig_corr_data$RNA_logFC < 0 & sig_corr_data$ATAC_logFC < 0] <- "Down in Both"
  sig_corr_data$quadrant[sig_corr_data$RNA_logFC > 0 & sig_corr_data$ATAC_logFC < 0] <- "Up RNA, Down ATAC"
  sig_corr_data$quadrant[sig_corr_data$RNA_logFC < 0 & sig_corr_data$ATAC_logFC > 0] <- "Down RNA, Up ATAC"
  
  # Set color palette for quadrants
  quad_colors <- c("Up in Both" = "#E41A1C", 
                   "Down in Both" = "#377EB8", 
                   "Up RNA, Down ATAC" = "#4DAF4A", 
                   "Down RNA, Up ATAC" = "#984EA3")
  
  # Create scatter plot
  scatter_plot <- ggplot(sig_corr_data, aes(x = RNA_logFC, y = ATAC_logFC, color = quadrant)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
    scale_color_manual(values = quad_colors) +
    labs(
      x = expression(RNA-seq~log[2](fold~change)),
      y = expression(ATAC-seq~log[2](fold~change)),
      title = "Correlation between RNA-seq and ATAC-seq Changes",
      subtitle = paste("Pearson r =", round(correlation, 2), 
                       ", R² =", round(r_squared, 2), 
                       ", p =", format(p_value, digits = 3, scientific = TRUE))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  # Label top genes (highest absolute fold changes in both datasets)
  top_genes <- sig_corr_data %>%
    mutate(abs_rna_fc = abs(RNA_logFC),
           abs_atac_fc = abs(ATAC_logFC),
           combined_score = abs_rna_fc + abs_atac_fc) %>%
    arrange(desc(combined_score)) %>%
    head(15)
  
  scatter_plot <- scatter_plot +
    geom_text_repel(
      data = top_genes,
      aes(label = gene_symbol),
      box.padding = 0.5,
      point.padding = 0.2,
      size = 3,
      segment.color = "grey50"
    )
  
  # Save the scatter plot
  ggsave("RNA_ATAC_correlation_DESeq2.png", scatter_plot, width = 10, height = 8, dpi = 300)
  ggsave("RNA_ATAC_correlation_DESeq2.pdf", scatter_plot, width = 10, height = 8)
}

# Create Venn Diagram ----

# For a more realistic integration, we'll use our mapped genes
# Get list of genes with significant RNA expression changes (with version numbers)
rna_sig_genes_with_version <- rownames(rna_results_df)[rna_results_df$padj < 0.05 & !is.na(rna_results_df$padj)]

# Remove version numbers from RNA gene IDs
rna_sig_genes <- sub("\\.[0-9]+$", "", rna_sig_genes_with_version)

# Get list of genes with significant ATAC accessibility changes
atac_sig_genes_raw <- na.omit(unique(atac_results_df$ensembl_id[atac_results_df$padj < 0.05 & !is.na(atac_results_df$padj)]))

# Remove version numbers from ATAC gene IDs (if they exist)
atac_sig_genes <- sub("\\.[0-9]+$", "", atac_sig_genes_raw)

# Calculate overlaps using version-free IDs
atac_only <- setdiff(atac_sig_genes, rna_sig_genes)
rna_only <- setdiff(rna_sig_genes, atac_sig_genes)
overlap <- intersect(atac_sig_genes, rna_sig_genes)

# Log number of genes in each category
cat("Genes with significant RNA expression changes:", length(rna_sig_genes), "\n")
cat("Genes with significant ATAC accessibility changes:", length(atac_sig_genes), "\n")
cat("Genes with both RNA and ATAC changes:", length(overlap), "\n")

# Create Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(atac_sig_genes),
  area2 = length(rna_sig_genes),
  cross.area = length(overlap),
  category = c("ATAC", "RNA DEGs"),
  fill = c("#377EB8", "#E41A1C"),
  alpha = c(0.5, 0.5),
  cat.col = c("black", "black"),
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.05, 0.05),
  cat.pos = c(0, 0),
  euler.d = FALSE,
  scaled = FALSE
)

# Save Venn diagram
png("RNA_ATAC_venn_DESeq2.png", width = 800, height = 800, res = 120)
grid.draw(venn.plot)
dev.off()

# Save lists of genes in each category - use the version-free IDs for consistency
write.csv(data.frame(ensembl_id = overlap), "genes_with_both_RNA_ATAC_changes.csv", row.names = FALSE)
write.csv(data.frame(ensembl_id = rna_only), "genes_with_RNA_only_changes.csv", row.names = FALSE) 
write.csv(data.frame(ensembl_id = atac_only), "genes_with_ATAC_only_changes.csv", row.names = FALSE)

# ---- PART 4: Create a heatmap of concordant changes ----
# Get genes that change in the same direction in both datasets
concordant_data <- correlation_data %>%
  filter((RNA_logFC > 0 & ATAC_logFC > 0) | (RNA_logFC < 0 & ATAC_logFC < 0)) %>%
  filter(RNA_padj < 0.05 & ATAC_padj < 0.05) %>%
  arrange(desc(abs(RNA_logFC) * abs(ATAC_logFC))) %>%
  head(50)  # Top 50 concordant genes

# Check if we have concordant genes
if(nrow(concordant_data) == 0) {
  cat("No concordant genes found with significant changes in both datasets. Skipping heatmap.\n")
} else {
  # Get normalized counts
  rna_norm_counts <- counts(dds_rna, normalized = TRUE)
  atac_norm_counts <- counts(dds_atac, normalized = TRUE)
  
  # Try to extract RNA data for concordant genes
  rna_ids <- concordant_data$ensembl_id
  rna_data_exists <- all(rna_ids %in% rownames(rna_norm_counts))
  
  if(!rna_data_exists) {
    cat("Cannot find all concordant genes in RNA count data. Trying alternative mapping...\n")
    
    # Try to handle different ID formats in RNA data
    if(any(grepl("ENSMUSG", rownames(rna_norm_counts)))) {
      # Extract Ensembl IDs from RNA row names if they're in a different format
      rna_ensembl_ids <- sub("^(ENSMUSG[0-9]+).*$", "\\1", rownames(rna_norm_counts))
      
      # Create a mapping from extracted Ensembl IDs to row names
      rna_id_to_row <- setNames(rownames(rna_norm_counts), rna_ensembl_ids)
      
      # Get row names in RNA data that match the concordant genes
      matched_rna_rows <- rna_id_to_row[rna_ids]
      
      # Check if we found matches
      if(all(!is.na(matched_rna_rows))) {
        cat("Found alternative mapping for RNA data.\n")
        rna_data_exists <- TRUE
        
        # Extract RNA data using matched rows
        rna_heatmap_data <- log2(rna_norm_counts[matched_rna_rows, ] + 1)
      } else {
        cat("Still cannot find all concordant genes in RNA count data. Skipping heatmap.\n")
      }
    } else {
      cat("Cannot determine ID format in RNA count data. Skipping heatmap.\n")
    }
  } else {
    # Extract RNA data directly
    rna_heatmap_data <- log2(rna_norm_counts[rna_ids, ] + 1)
  }
  
  # If we have RNA data, continue with ATAC data and create heatmap
  if(exists("rna_heatmap_data")) {
    # Find the matching peak IDs for each gene in the concordant set
    atac_peak_ids <- concordant_data$peak_id
    
    # Check if all peak IDs exist in ATAC count data
    atac_data_exists <- all(atac_peak_ids %in% rownames(atac_norm_counts))
    
    if(!atac_data_exists) {
      cat("Cannot find all concordant peaks in ATAC count data. Skipping heatmap.\n")
    } else {
      # Extract ATAC data
      atac_heatmap_data <- log2(atac_norm_counts[atac_peak_ids, ] + 1)
      
      # Z-score normalize the data
      rna_z <- t(scale(t(rna_heatmap_data)))
      atac_z <- t(scale(t(atac_heatmap_data)))
      
      # Combine the datasets
      combined_data <- rbind(
        rna_z,
        atac_z
      )
      
      # Create annotation for the heatmap
      sample_anno <- data.frame(
        Condition = c("CMP", "CMP", "HSC", "HSC"),
        row.names = colnames(combined_data)
      )
      
      # For row annotation, create labels with gene symbols
      row_labels <- c(
        paste0("RNA:", concordant_data$gene_symbol),
        paste0("ATAC:", concordant_data$gene_symbol)
      )
      rownames(combined_data) <- row_labels
      
      row_anno <- data.frame(
        DataType = c(rep("RNA-seq", nrow(rna_z)), rep("ATAC-seq", nrow(atac_z))),
        Direction = c(
          ifelse(concordant_data$RNA_logFC > 0, "Up", "Down"),
          ifelse(concordant_data$ATAC_logFC > 0, "Up", "Down")
        ),
        row.names = row_labels
      )
      
      # Define colors
      ann_colors <- list(
        Condition = c(CMP = "#377EB8", HSC = "#E41A1C"),
        DataType = c("RNA-seq" = "#4DAF4A", "ATAC-seq" = "#984EA3"),
        Direction = c("Up" = "red", "Down" = "blue")
      )
      
      # Create the heatmap
      pheatmap(
        combined_data,
        annotation_col = sample_anno,
        annotation_row = row_anno,
        annotation_colors = ann_colors,
        show_rownames = TRUE,
        show_colnames = TRUE,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        fontsize_row = 8,
        filename = "RNA_ATAC_concordant_heatmap_DESeq2.png",
        width = 10,
        height = 14
      )
      
      # Save as PDF as well
      pheatmap(
        combined_data,
        annotation_col = sample_anno,
        annotation_row = row_anno,
        annotation_colors = ann_colors,
        show_rownames = TRUE,
        show_colnames = TRUE,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        fontsize_row = 8,
        filename = "RNA_ATAC_concordant_heatmap_DESeq2.pdf",
        width = 10,
        height = 14
      )
    }
  }
}