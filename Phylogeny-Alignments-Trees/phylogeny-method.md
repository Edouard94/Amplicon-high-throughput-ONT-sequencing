
The Basic Local Alignment Search Tool (BLAST) was employed for each retrieved Nanopore sequence against the non-redundant GenBank database to identify homologous nucleotide sequences (Altschul et al., 1997). The sequences from the top 50 BLASTn hits, based on the lowest e-value, were downloaded for each queried sequence. These 50 nucleotide hit sequences were aligned using MAFFT v7 (Katoh & Standley, 2013) before being trimmed using BMGE (Criscuolo & Gribaldo, 2010). For Gregarine Nanopore SSU sequences, this process included an additional set of 100 SSU sequences from Boisard et al. (2022), and LSU sequences were added to the LSU datasets. Maximum likelihood analyses were performed using FastTree version 2.1.11 with the -gtr and -gamma options (Price et al., 2009), and phylogenetic trees were visualized using PRESTO in the NGPhylogeny.fr online workflow (Lemoine et al., 2019).

Subsequently, hit sequences that were phylogenetically divergent were manually removed from the dataset for each tree. Datasets were then merged according to parasite phyla, and duplicate sequences were removed. The workflow was re-run using the selected and combined sequences to check and manually prune the phylogenetic trees by phyla. The final dataset for each parasite phylum was aligned in MAFFT version 7 (Katoh & Standley, 2013), trimmed with BMGE (Criscuolo & Gribaldo, 2010), and rooted maximum likelihood trees with 1,000 bootstraps were constructed in FastTree version 2.1.11 using the -gtr option (Price et al., 2009). The phylogenetic trees were visualized and produced using FigTree version 1.4.4 (Rambaut, 2018).

### References:
Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W. & Lipman D. J. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Research, 25: 3389–3402.

Boisard, J. et al. (2022) ‘Marine gregarine genomes reveal the breadth of apicomplexan diversity with a partially conserved glideosome machinery’, BMC Genomics, 23(1), p. 485. Available at: https://doi.org/10.1186/s12864-022-08700-8.

Criscuolo, A. & Gribaldo, S. (2010). BMGE (Block Mapping and Gathering with Entropy): a new software for selection of phylogenetic informative regions from multiple sequence alignments. BMC Evol Biol 10, 210. 

Katoh, K. & Standley, D.M. 2013. MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. Molecular Biology and Evolution 30: 772–780.

Lemoine, F., Correia, D., Lefort, V., Doppelt-Azeroual, O., Mareuil, F., Cohen-Boulakia, S., et al. 2019. NGPhylogeny.fr: new generation phylogenetic services for non-specialists. Nucleic Acids Research 47: W260–W265.

Price, M. N., Dehal, P. S. & Arkin, A. P. (2009). FastTree: computing large minimum evolution trees with profiles instead of a distance matrix. Molecular Biology and Evolution, 26, 1641–1650.

Rambaut, A. (2018). FigTree-v1.4.4. http://treebioedacuk/software/figtree/ (accessed April 5, 2023).
