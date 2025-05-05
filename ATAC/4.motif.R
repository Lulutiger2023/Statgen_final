library(Biostrings)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MotifDb)
library(seqLogo)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Import peak files identified by DESeq2
up_peaks <- import.bed("upregulated_peaks_DESeq2.bed")
down_peaks <- import.bed("downregulated_peaks_DESeq2.bed")

# Extract genomic sequences from mouse reference genome
up_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, up_peaks)
down_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, down_peaks)

# Assign sequence names
names(up_sequences) <- paste0("up_peak_", 1:length(up_sequences))
names(down_sequences) <- paste0("down_peak_", 1:length(down_sequences))

# Write sequences to FASTA files
writeXStringSet(up_sequences, "upregulated_peak_sequences.fa")
writeXStringSet(down_sequences, "downregulated_peak_sequences.fa")

# Query MotifDb for hematopoietic transcription factors
mouse_motifs <- query(MotifDb, andStrings="Mmusculus")
gata1_motifs <- query(MotifDb, andStrings=c("Mmusculus", "GATA1"))
tal1_motifs <- query(MotifDb, andStrings=c("Mmusculus", "TAL1"))
klf1_motifs <- query(MotifDb, andStrings=c("Mmusculus", "KLF1"))
runx1_motifs <- query(MotifDb, andStrings=c("Mmusculus", "RUNX1"))

# Define simplified binding motifs for key transcription factors
simple_motifs <- list(
  GATA1 = "GATAA",      # GATA basic binding site
  GATA2 = "TGATA",      # Classical GATA binding site
  EBOX1 = "CACGTG",     # Classical E-box
  EBOX2 = "CAGCTG",     # Alternative E-box
  KLF1 = "CACCC",       # KLF1 binding site
  KLF2 = "CCACCC",      # Extended KLF binding site  
  RUNX1 = "TGTGGT",     # RUNX1 binding site
  RUNX2 = "ACCACA"      # RUNX1 reverse complement
)

# Function to search for motif patterns in sequences
search_motif_patterns <- function(sequences, patterns) {
  results <- data.frame(
    Motif = names(patterns),
    Pattern = unlist(patterns),
    Count = 0,
    Percentage = 0
  )
  
  for (i in 1:length(patterns)) {
    motif_name <- names(patterns)[i]
    pattern <- patterns[[i]]
    
    # Create DNA string pattern
    dna_pattern <- DNAString(pattern)
    
    # Search for exact matches
    total_matches <- 0
    sequences_with_motif <- 0
    
    for (j in 1:length(sequences)) {
      curr_matches <- countPattern(dna_pattern, sequences[[j]], max.mismatch=0)
      total_matches <- total_matches + curr_matches
      if (curr_matches > 0) {
        sequences_with_motif <- sequences_with_motif + 1
      }
    }
    
    percentage <- 100 * sequences_with_motif / length(sequences)
    
    results$Count[i] <- total_matches
    results$Percentage[i] <- percentage
  }
  
  return(results)
}

# Search for motifs in up- and down-regulated peaks
up_motif_results <- search_motif_patterns(up_sequences, simple_motifs)
down_motif_results <- search_motif_patterns(down_sequences, simple_motifs)

# Compare motif distribution
motif_comparison <- data.frame(
  Motif = up_motif_results$Motif,
  Pattern = up_motif_results$Pattern,
  Up_Count = up_motif_results$Count,
  Up_Percentage = up_motif_results$Percentage,
  Down_Count = down_motif_results$Count,
  Down_Percentage = down_motif_results$Percentage,
  Fold_Change = log2((up_motif_results$Percentage + 0.1) / (down_motif_results$Percentage + 0.1))
)

# Visualize results - Fold change bar plot
motif_comparison$Direction <- ifelse(motif_comparison$Fold_Change > 0, 
                                    "Enriched in Upregulated Peaks", 
                                    "Enriched in Downregulated Peaks")

p2 <- ggplot(motif_comparison, aes(x = reorder(Motif, Fold_Change), 
                                  y = Fold_Change, 
                                  fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(title = "Motif Enrichment Fold Change (log2) in Up vs Down Regulated Peaks",
       x = "Motif",
       y = "Fold Change (log2)") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("Enriched in Upregulated Peaks" = "red", 
                              "Enriched in Downregulated Peaks" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave("motif_fold_change.pdf", p2, width = 8, height = 6)