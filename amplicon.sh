#!/bin/bash

#------------------------------------------------------------------------------
# IMPORTANT: Run the script with 'source' to ensure medaka gets activated properly.
# For example: source amplicon.sh <input_file> <forward_primer> <reverse_primer> <min_amplicon_size> [max_amplicon_size]
# e.g. source amplicon.sh fastq_runid_cd2409afdc364972d4e64953023ac179690563d1_0.fastq CCCTTAGATRRYCTGGGCTGC CGTGTTACGACTTCTTC 250 350
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check input parameters  
#------------------------------------------------------------------------------  
# The script requires:
#   1) FASTQ file
#   2) Forward primer
#   3) Reverse primer
#   4) Min amplicon size
#   5) [Optional] Max amplicon size (if omitted, exact size = min size)
if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    echo "Usage: source $0 <input_file> <forward_primer> <reverse_primer> <min_amplicon_size> [max_amplicon_size]"
    return 1
fi

input_file="$1"
forward_primer="$2"
reverse_primer="$3"
min_amplicon_size="$4"
if [ "$#" -eq 5 ]; then
    max_amplicon_size="$5"
else
    max_amplicon_size="$4"
fi

#------------------------------------------------------------------------------
# Define cutadapt parameters
#------------------------------------------------------------------------------
# The error rate is set based on Reis et al., 2023.
error_rate=0.1

#------------------------------------------------------------------------------
# Step 1: Calculate minimum overlap values based on primer lengths
#------------------------------------------------------------------------------
# This is used by cutadapt to ensure that the adapter is recognized even if only a part
# of the primer sequence is present in the read.
min_overlap_forward=$(echo -n "$forward_primer" | wc -c)
# Subtract 1 because wc -c counts the newline character.
min_overlap_forward=$((min_overlap_forward - 1))

min_overlap_reverse=$(echo -n "$reverse_primer" | wc -c)
min_overlap_reverse=$((min_overlap_reverse - 1))

# Create the adapter input strings for cutadapt
forward_input="$forward_primer;min_overlap=$min_overlap_forward"
reverse_input="$reverse_primer;min_overlap=$min_overlap_reverse"

#------------------------------------------------------------------------------
# Step 2: Trim the reads using cutadapt
#------------------------------------------------------------------------------
# The reads are trimmed using two orientations:
#
# (a) Forward orientation:
#     Uses the forward primer as the adapter and the reverse primer as the anchored 
#     adapter.
echo "Trimming reads (forward orientation)..."
cutadapt -a "$forward_input" -g "$reverse_input" \
         --error-rate "$error_rate" --trimmed-only \
         -o trimmed_reads.fastq --cores 10 "$input_file" > adapter_summary.txt

# (b) Reverse complement orientation:
#     Reverse the FASTQ file using seqkit and then apply the trimmed adapter options.
echo "Trimming reads (reverse complement orientation)..."
cutadapt -a "$reverse_input" -g "$forward_input" \
         --error-rate "$error_rate" --trimmed-only \
         -o trimmed_reads_RC.fastq --cores 10 \
         <(seqkit seq --seq-type DNA -r -p "$input_file") > adapter_RC_summary.txt

# To avoid duplicate read IDs when merging, add "RC" to headers in the reverse complement file.
sed '1~4 s/^@/@RC/' trimmed_reads_RC.fastq > trimmed_reads_RC1.fastq
rm trimmed_reads_RC.fastq

#------------------------------------------------------------------------------
# Step 3: Concatenate trimmed reads and perform quality filtering
#------------------------------------------------------------------------------
# Merge the trimmed reads from both orientations into a single FASTQ file.
cat trimmed_reads.fastq trimmed_reads_RC1.fastq > concatenated_reads.fastq

# Filter the reads by trimming low-quality ends and discarding reads shorter than 50 nucleotides.
fastq_quality_trimmer -t 10 -l 50 -i concatenated_reads.fastq -o filtered_reads.fastq

#------------------------------------------------------------------------------
#  Step 4: Filter reads based on the amplicon size window  
#------------------------------------------------------------------------------  
echo "Filtering reads between ${min_amplicon_size}-${max_amplicon_size} bp"
seqkit seq -m "$min_amplicon_size" -M "$max_amplicon_size" \
           -i filtered_reads.fastq -o final_reads.fastq

#------------------------------------------------------------------------------
# Step 5: Calculate and save basic metrics
#------------------------------------------------------------------------------
# Count total reads by dividing the total number of lines in the FASTQ file by 4.
total_reads=$(( $(wc -l < final_reads.fastq) / 4 ))
echo "Total reads: $total_reads" >> metrics_results.txt

#------------------------------------------------------------------------------
# Step 6: Consensus building and polishing using NGSpeciesID and medaka
#------------------------------------------------------------------------------
echo "Running NGSpeciesID consensus and polishing..."
# Activate the NGSpeciesID conda environment
conda activate NGSpeciesID
# Activate medaka (requires sourcing)
. ~/medaka/bin/activate
NGSpeciesID --t 12 --ont --consensus --sample_size 300 --medaka \
            --fastq final_reads.fastq --outfolder consensus_results
conda deactivate

# Concatenate all consensus sequences from NGSpeciesID output into one final FASTA file.
cat consensus_results/medaka_cl_id*/consensus.fasta > final_consensus.fasta

echo "Pipeline completed. Final consensus file: final_consensus.fasta"
