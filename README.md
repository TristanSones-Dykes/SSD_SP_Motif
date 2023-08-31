# Ssd-1 Signal Recognition Particle targeting analysis using Signal Peptide, Hydrophobicity, and Binding Motif techniques

## This repo contains:
 - R scripts defining functions for Signal Peptide (SP) detection using [SignalP6](https://services.healthtech.dtu.dk/services/SignalP-6.0/) locally, for analysing hydrophobicity, and binding motif detection using Biostrings.
 - Python script for web scraping SP and Transmembrane Domain (TM) detetion using [Phobius](https://phobius.sbc.su.se).
 - Rmd documents running analyses and generating reports.
 - Saved SignalP and Phobius outputs for S. cerevisiae, S. pombe, N. crassa, A. fumigatus, A. aspergillus, Cryptococcus neoformans H99, and C. albicans.

## Main pipeline design

![Flow chart diagram of analysis pipelines](https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/exported%20image%20no%20background.png)

It takes in two FASTA files – the transcriptome and proteome of the target species – and detects SPs and TMs. It then uses these and hydrophobicity index scores to split this detected set into two groups, one that uses the SRP and one that doesn't.
Combining these groupings with those from primary Ssd-1 binding motif (CNYTCNYT) detection, and the results are these distribution and population differences:

<img src="[image1.png](https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/plots/plot%201.png)https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/plots/plot%201.png" width="50%"/> <img src="[image2.png](https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/plots/plot%202.png)https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/plots/plot%202.png" width="50%"/> 
