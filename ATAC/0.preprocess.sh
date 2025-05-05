#!/bin/bash
# Preprocessing script for ATAC-seq BED files to BigWig format

# Function to process a single BED file
process_bed_to_bw() {
    local file=$1
    local prefix=$2
    
    # Step 1: Extract columns 1-4 and sort
    cut -f 1,2,3,4 "$file" | sort -k 1,1 -k 2,2n > "${prefix}.bedGraph"
    
    # Step 2: Adjust coordinates (BED to BedGraph coordinate conversion)
    cat "${prefix}.bedGraph" | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > "${prefix}.bedGraph.tmp1"
    
    # Step 3: Sort and merge overlapping regions with mean calculation
    bedtools sort -i "${prefix}.bedGraph.tmp1" | bedtools merge -i - -c 4 -o mean > "${prefix}.bedGraph.tmp2"
    
    # Step 4: Final formatting
    cat "${prefix}.bedGraph.tmp2" | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > "${prefix}.merge.bedGraph"
    
    # Step 5: Convert to BigWig
    bedGraphToBigWig "${prefix}.merge.bedGraph" mm10.chrom.sizes "${prefix}.bw"
    
    # Clean up temporary files
    rm "${prefix}.bedGraph.tmp1" "${prefix}.bedGraph.tmp2"
}

# Process each sample
process_bed_to_bw "ENCFF832UUS.bed" "CMP_rep1"
process_bed_to_bw "ENCFF343PTQ.bed" "CMP_rep2"
process_bed_to_bw "ENCFF796ZSB.bed" "CFUE_rep1"
process_bed_to_bw "ENCFF599ZDJ.bed" "CFUE_rep2"

# Create a modified version of the merged peaks file with chr_start_end IDs
awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3,$5,$6}' all_replicated_merged.bed > all_replicated_merged.withpkname.bed

# Generate average signal over merged peaks with the modified peak IDs
bigWigAverageOverBed CMP_rep1.bw all_replicated_merged.withpkname.bed CMP_rep1_signal.tab
bigWigAverageOverBed CMP_rep2.bw all_replicated_merged.withpkname.bed CMP_rep2_signal.tab
bigWigAverageOverBed CFUE_rep1.bw all_replicated_merged.withpkname.bed CFUE_rep1_signal.tab
bigWigAverageOverBed CFUE_rep2.bw all_replicated_merged.withpkname.bed CFUE_rep2_signal.tab

echo "Preprocessing complete. Signal files generated for DESeq2 analysis with chr_start_end format peak IDs."