from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import os

# Load sequences from fasta file
sequences = list(SeqIO.parse("sequences.fasta", "fasta"))

# Create local BLAST database
makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file="sequences.fasta")
stdout, stderr = makeblastdb_cline()

# Function to perform local BLAST and return percentage identity
def blast(seq1, seq2):
    SeqIO.write(seq1, "query.fasta", "fasta")  # Write seq1 to query.fasta
    SeqIO.write(seq2, "subject.fasta", "fasta")  # Write seq2 to subject.fasta

    # Perform local BLAST
    blastn_cline = NcbiblastnCommandline(query="query.fasta", subject="subject.fasta", outfmt=5, out="blast.xml")
    stdout, stderr = blastn_cline()

    # Parse BLAST result
    result_handle = open("blast.xml")
    blast_record = NCBIXML.read(result_handle)  # Use read() as there is only one record now
    max_identity = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            identity = hsp.identities / len(hsp.match) * 100
            if identity > max_identity:
                max_identity = identity
    result_handle.close()
    os.remove("query.fasta")
    os.remove("subject.fasta")
    os.remove("blast.xml")
    return max_identity

# Create DataFrame to store percentage identities
df = pd.DataFrame(index=[seq.id for seq in sequences], columns=[seq.id for seq in sequences])

# Perform BLAST for each pair of sequences and store percentage identity
for i in range(len(sequences)):
    for j in range(i+1, len(sequences)):  # Start from i+1 to avoid self-comparison
        identity = blast(sequences[i], sequences[j])
        df.loc[sequences[i].id, sequences[j].id] = identity
        df.loc[sequences[j].id, sequences[i].id] = identity

# Save DataFrame to CSV
df.to_csv('blast_results.csv')
