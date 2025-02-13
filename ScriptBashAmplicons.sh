#!/bin/bash

# Run the script with 'source'

# Define the input file and the list of primer pairs to use
##"V1F CACCAGGTTGATTCTGCCTGAC 530R CCGCGGCKGCTGGCAC" "18SF CGGTAATTCCAGCTCCAAT 18SR TGACTTGCGCTTACTAGGG" "WL1 GCGCTACCTGGTTGATCCTGCC EukP3 GACGGGCGGTGTGTAC" "574f CGGTAAYTCCAGCTCYV UNonMetDB CTTTAARTTTCASYCTTGCG" "Am_RC GTGATGCCCTTAGAYRTCTTGG iLSUr ACCTGTCTCACGACGGTCTAAAC" "Kineto_RC AGCTGTAGGTGAACCTGCAGAAGGATC iLSUr ACCTGTCTCACGACGGTCTAAAC" "Greg_RC GAAGAAGTCGTAACACG iLSUr ACCTGTCTCACGACGGTCTAAAC" "V5_RC GAGTGGAKTGTGCTGNTTA iLSUr ACCTGTCTCACGACGGTCTAAAC" "V1F CACCAGGTTGATTCTGCCTGAC V5R TAANCAGCACAMTCCACTC" "AmF GAATTGACGGAAGGGCACAC AmR CCAAGAYRTCTAAGGGCATCAC" "KinetoF CTGCCAGTAGTCATATGCTTGTTTCAAGGA KinetoR GATCCTTCTGCAGGTTCACCTACAGCT" 
input_file="/home/edouard/SSU-Insects-17May2023/no_sample/20230517_1959_MN39555_FAW51311_74e0d713/basecall_pass_gpu_010623_q7/barcoder_SSU170523_mid_barcodes/barcode08/Barcode08_cat.fastq"
primer_pairs=("GregF CCCTTAGATRRYCTGGGCTGC GregR CGTGTTACGACTTCTTC") # Modify with you own primer pairs

# Define the cutadapt parameters
error_rate=0.1 #(value from Reis et al., 2023)

# Trim the reads for each primer pair
for pair in "${primer_pairs[@]}"; do
    forward_primer=$(echo "$pair" | awk '{print $2}')
    reverse_primer=$(echo "$pair" | awk '{print $4}')
    ## Calculate the min_overlap value based on the forward primer length
    min_overlap_forward=$(echo "$forward_primer" | wc -c) # count the number of characters in the forward primer
    min_overlap_forward=$((min_overlap_forward - 1)) # subtract 1 because wc -c counts the newline character as well
    ## Calculate the min_overlap value based on the reverse primer length
    min_overlap_reverse=$(echo "$reverse_primer" | wc -c) # count the number of characters in the reverse primer
    min_overlap_reverse=$((min_overlap_reverse - 1)) # subtract 1 because wc -c counts the newline character as well
    ## Create the input for the -a option with the adapter and the min_overlap value
    forward_input="$forward_primer;min_overlap=$min_overlap_forward"
    reverse_input="$reverse_primer;min_overlap=$min_overlap_reverse"
    ## Trim the reads based on the forward primer
    cutadapt -a "$forward_input" -g "$reverse_input" --error-rate "$error_rate" --trimmed-only -o trimmed_reads_${pair// /_}.fastq --cores 10 "$input_file" > adapter_summary.txt
    ## Trim the reads based on their reverse complement (reverse primer becomes the forward primer)
    cutadapt -a "$reverse_input" -g "$forward_input" --error-rate "$error_rate" --trimmed-only -o trimmed_reads_RC_${pair// /_}.fastq --cores 10 <(seqkit seq --seq-type DNA -r -p "$input_file") > adapter_RC_summary.txt
    ## Add 'RC' after the first @ in the trimmed_reads_RC_ files using sed (to avoid ID duplicates when the files are cat together)
    sed '1~4 s/^@/@RC/' trimmed_reads_RC_${pair// /_}.fastq > trimmed_reads_RC1_${pair// /_}.fastq
    rm trimmed_reads_RC_*.fastq
done

# Concatenate the trimmed read files for each primer pair and filter reads based on length
## All the working primers need to be present here
for pair in "${primer_pairs[@]}"; do
    ## Concatenate the two trimmed read files for the current primer pair
    cat trimmed_reads_${pair// /_}.fastq trimmed_reads_RC1_${pair// /_}.fastq > concatenated_reads_${pair// /_}.fastq
    fastq_quality_trimmer -t 10 -l 50 -i concatenated_reads_${pair// /_}.fastq -o filtered_reads_${pair// /_}.fastq # -t 10: Quality threshold - nucleotides with lower quality than 10 will be trimmed (from the end of the sequence). -l 50: Minimum length - sequences shorter than 50 (after trimming) will be discarded.
    ## Filter reads for length
    ### Check if the current primer pair is in the list of primer pairs that you want to use
    if [[ "$pair" == "V1F CACCAGGTTGATTCTGCCTGAC 530R CCGCGGCKGCTGGCAC" ]]; then
        ### Define length filtering thresholds for V1F/530R primer pair
        V1F_size=410
        V1F_min_length=380
        V1F_max_length=490
        ### Filter the reads based on their length
        seqkit seq -m "$V1F_min_length" -M "$V1F_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "18SF CGGTAATTCCAGCTCCAAT 18SR TGACTTGCGCTTACTAGGG" ]]; then
        ### Define length filtering thresholds for 18SF/18SR primer pair
        #### A variable cannot start with a number so use _18S instead of 18S
        _18S_size=950
        _18S_min_length=920
        _18S_max_length=1180
        ### Filter the reads based on their length
        seqkit seq -m "$_18S_min_length" -M "$_18S_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "WL1 GCGCTACCTGGTTGATCCTGCC EukP3 GACGGGCGGTGTGTAC" ]]; then
        ### Define length filtering thresholds for WL1/EukP3 primer pair
        WL1_size=1630
        WL1_min_length=1530
        WL1_max_length=1650
        ### Filter the reads based on their length
        seqkit seq -m "$WL1_min_length" -M "$WL1_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "GregF CCCTTAGATRRYCTGGGCTGC GregR CGTGTTACGACTTCTTC" ]]; then
        ### Define length filtering thresholds for GregF/GregR primer pair
        Greg_size=300
        Greg_min_length=250
        Greg_max_length=350
        ### Filter the reads based on their length
        seqkit seq -m "$Greg_min_length" -M "$Greg_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "574f CGGTAAYTCCAGCTCYV UNonMetDB CTTTAARTTTCASYCTTGCG" ]]; then
        ### Define length filtering thresholds for 574f/UNonMetDB primer pair
        _574_size=570
        _574_min_length=474
        _574_max_length=674
        ### Filter the reads based on their length
        seqkit seq -m "$_574_min_length" -M "$_574_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "Greg_RC GAAGAAGTCGTAACACG iLSUr ACCTGTCTCACGACGGTCTAAAC" ]]; then
        ### Define length filtering thresholds for Greg_RC primer pair
        Greg_RC_size=3350
        Greg_RC_min_length=3000
        Greg_RC_max_length=3700
        ### Filter reads based on length
        seqkit seq -m "$Greg_RC_min_length" -M "$Greg_RC_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "Am_RC GTGATGCCCTTAGAYRTCTTGG iLSUr ACCTGTCTCACGACGGTCTAAAC" ]]; then
        ### Define length filtering thresholds for Am_RC primer pair
        Am_RC_size=4200
        Am_RC_min_length=4000
        Am_RC_max_length=4400
        ### Filter reads based on length
        seqkit seq -m "$Am_RC_min_length" -M "$Am_RC_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "Kineto_RC AGCTGTAGGTGAACCTGCAGAAGGATC iLSUr ACCTGTCTCACGACGGTCTAAAC" ]]; then
        ### Define length filtering thresholds for Kineto_RC primer pair
        Kineto_RC_size=5000
        Kineto_RC_min_length=4800
        Kineto_RC_max_length=5200
        ### Filter reads based on length
        seqkit seq -m "$Kineto_RC_min_length" -M "$Kineto_RC_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "KinetoF CTGCCAGTAGTCATATGCTTGTTTCAAGGA KinetoR GATCCTTCTGCAGGTTCACCTACAGCT" ]]; then
        ### Define length filtering thresholds for KinetoF primer pair
        KinetoF_size=2000
        KinetoF_min_length=1800
        KinetoF_max_length=2200
        ### Filter reads based on length
        seqkit seq -m "$KinetoF_min_length" -M "$KinetoF_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "AmF GAATTGACGGAAGGGCACAC AmR CCAAGAYRTCTAAGGGCATCAC" ]]; then
        ### Define length filtering thresholds for AmF primer pair
        AmF_size=400
        AmF_min_length=200
        AmF_max_length=600
        ### Filter reads based on length
        seqkit seq -m "$AmF_min_length" -M "$AmF_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "V1F CACCAGGTTGATTCTGCCTGAC V5R TAANCAGCACAMTCCACTC" ]]; then
    ### Define length filtering thresholds for V1F/V5R primer pair
        V1FV5R_size=700
        V1FV5R_min_length=500
        V1FV5R_max_length=900
        ### Filter the reads based on their length
        seqkit seq -m "$V1FV5R_min_length" -M "$V1FV5R_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    elif [[ "$pair" == "V5_RC GAGTGGAKTGTGCTGNTTA iLSUr ACCTGTCTCACGACGGTCTAAAC" ]]; then
    ### Define length filtering thresholds for V5_RC/iLSUr primer pair
        V5R_RC_size=2800
        V5R_RC_min_length=2700
        V5R_RC_max_length=2900
        ### Filter the reads based on their length
        seqkit seq -m "$V5R_RC_min_length" -M "$V5R_RC_max_length" -i filtered_reads_${pair// /_}.fastq -o final_reads_${pair// /_}.fastq
    fi 
done

# Calculate metrics for concatenated reads and save results to a text file
for pair in "${primer_pairs[@]}"; do
    forward_primer=$(echo "$pair" | awk '{print $2}')
    reverse_primer=$(echo "$pair" | awk '{print $4}')
    total_reads=$(( $(wc -l < final_reads_${pair// /_}.fastq) / 4 ))
    echo "Metrics for primer pair $pair:" >> metrics_results.txt  # Append instead of overwrite
    echo "Total reads: $total_reads" >> metrics_results.txt
    # Calculate sample size for the following NGSpeciesID cmd, as 50% of total reads
    #sample_size=$(( total_reads / 2 ))
done

# Use of NGSpeciesID for reads consensus and polishing with medaka
conda activate NGSpeciesID
. ~/medaka/bin/activate # Activate medaka
for pair in "${primer_pairs[@]}"; do
    NGSpeciesID --t 12 --ont --consensus --sample_size 300  --medaka --fastq final_reads_${pair// /_}.fastq --outfolder ${pair// /_}_consensus 
    ## Subsampling value recommended in 'Rapid in situ...' Pomerantz at al., 2022, it creates an abundance_cutoff when medaka is forming draft consensus with >= 30 (10.0% of 300 reads) 
done
conda deactivate
# Concatenate all the consensus and polished sequences (medaka) into one fasta file
for pair in "${primer_pairs[@]}"; do
    cat ${pair// /_}_consensus/medaka_cl_id*/consensus.fasta > cat_${pair// /_}_consensus.fasta
done

# SCRIPT END

# QIIME OPTION

# Load Qiime2 for reads processing
conda activate qiime2-2023.5

# Import the data in qiime
for pair in "${primer_pairs[@]}"; do
    ## Create a manifest file with the path to the final read file
    echo -e "sample-id\tabsolute-filepath" > manifest_${pair// /_}.txt
    echo -e "Barcode01\t$(pwd)/final_reads_${pair// /_}.fastq" >> manifest_${pair// /_}.txt # Modify with Barcode number 
    ## Import the filtered FASTQ file into Qiime2 as a single-end sequence
    ### When importing sequences with a manifest file in QIIME, the QIIME artifact contains all the sequences and sample ID.
    qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path manifest_${pair// /_}.txt --output-path "Barcode01_${pair// /_}.qza" --input-format SingleEndFastqManifestPhred33V2
done

# Dereplicate the sequences
for pair in "${primer_pairs[@]}"; do
qiime vsearch dereplicate-sequences --i-sequences Barcode01_${pair// /_}.qza --o-dereplicated-table table_${pair// /_}.qza --o-dereplicated-sequences rep-seqs_${pair// /_}.qza # Use the imported QIIME artifacts as input
done

# Cluster the sequences de novo
for pair in "${primer_pairs[@]}"; do
qiime vsearch cluster-features-de-novo --i-sequences rep-seqs_${pair// /_}.qza --i-table table_${pair// /_}.qza --o-clustered-table table-dn_${pair// /_}.qza --o-clustered-sequences rep-seqs-dn_${pair// /_}.qza --p-perc-identity 0.89 #(value from ECOP conf, otherwise 95% is used in Reis et al., 2023) 
done

# Perform chimera detection on the consensus sequences
for pair in "${primer_pairs[@]}"; do
qiime vsearch uchime-denovo --i-sequences rep-seqs-dn_${pair// /_}.qza --i-table table-dn_${pair// /_}.qza --o-chimeras chimeras_${pair// /_}.qza --o-nonchimeras nonchimeras_${pair// /_}.qza --o-stats chimeras_stats_${pair// /_}.qza
done
#qiime metadata tabulate --m-input-file chimeras_stats.qza --o-visualization chimeras_stats.qzv

# Convert the Qiime2 artifact to a fasta file
for pair in "${primer_pairs[@]}"; do
qiime tools export --input-path nonchimeras_${pair// /_}.qza --output-path nonchimeras_exported_${pair// /_} # Add the pair of primers to the input and output paths
done

conda deactivate

# The script for extracting the amplicons ends here, the next steps are for phylogenetic analyses and should be done according to your different organism groups

# Align the dereplicated sequences (NGSpeciesID or qiime outputs) with MAFFT, curate with bmge, find the best model and make a tree with iqtree
mafft --thread -1 --auto input.fasta > mafft_output.fasta
bmge -i mafft_output.fasta -t DNA -of bmge_output.fasta
#or
trimal -in mafft_output.fasta -out trimal_output.fasta -automated1 # Use a heuristic selection of the automatic method based on similarity statistics. (Optimized for Maximum Likelihood phylogenetic tree reconstruction).
iqtree2 -s trimal_output.fasta -m MFP -bb 1000 -T AUTO --prefix
