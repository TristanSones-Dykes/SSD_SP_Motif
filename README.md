# Ssd-1 Signal Recognition Particle targeting analysis using Signal Peptide, Hydrophobicity, and Binding Motif techniques

## This repo contains:
 - R scripts defining functions for Signal Peptide (SP) detection using [SignalP6](https://services.healthtech.dtu.dk/services/SignalP-6.0/) locally, for analysing hydrophobicity, and binding motif detection using Biostrings.
 - Python script for web scraping SP and Transmembrane Domain (TM) detetion using [Phobius](https://phobius.sbc.su.se).
 - Rmd documents running analyses and generating reports.
 - Saved SignalP and Phobius outputs for S. cerevisiae, S. pombe, N. crassa, A. fumigatus, A. aspergillus, Cryptococcus neoformans H99, and C. albicans.

## Main pipeline design

![Flow chart diagram of analysis pipelines](https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/exported%20image%20no%20background.png)

It takes in two FASTA files – the transcriptome and proteome of the target species – and detects SPs and TMs. It then uses these and hydrophobicity index scores to split this detected set into two groups, one that uses the SRP and one that doesn't.
Combining these groupings with those from primary Ssd-1 binding motif (CNYTCNYT) detection, and the results (for S. cerevisiae) are these distribution and population differences:

<div align="center">
<img src="https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/plots/plot%201.png" width=45%/> <img src="https://github.com/TristanSones-Dykes/SSD_SP_Motif/blob/master/plots/plot%202.png" width=45%/> 
</div>

This shows that the distribution of Ssd-1 targeted proteins is much more bimodal than non-targeted and a chi-squared test confirmed that the distribution of population distributions between the motif and non-motif groups are significantly dependent (p < 0.01).

## Installation instructions

There are a few steps to installation:
 - Install environment using conda/anaconda
 - Download and install SignalP
 - Add R libraries

#### Environment

To setup the environment, you need to have conda/anaconda installed and then use it to create an environment using the `requirements.txt` file at the top of the repo. I use [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Once you have it installed you can create the environment (once your current directory is the local copy of the repo) using:
```
conda create --name ssd_env --file requirements.txt
```
Feel free to replace `ssd_env` with whatever you want to call it.

You activate it using:
```
conda activate ssd_env
```

If it is too slow or seems to hang, you might want to consider using [Mamba](https://anaconda.org/conda-forge/mamba) (it is essentially conda written in C++ so is far faster).

#### SignalP installation

To install SignalP6, you need to fill in [this form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast) with your academic details and follow [their installation instructions](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md).

#### R libraries
When you open the project in a vscode terminal or RStudio, it should automatically call the renv using the .Rprofile. However, if it doesn't, run this and it should work.
```
if (!require("renv", quietly = TRUE))
    install.packages("renv")

renv::restore()
```

## How to use
You can generate the Rmd documents using `rmarkdown::render("<file_name>")` and can use them as a usage guide for the functions defined. If you want to use a set of the functions in another doc or project, you can call them using `source("<file_name>")`

For the species' data that is already in the repo, SignalP and Phobius results have been saved so they are called very quickly (as long as you use the correct file names), but new inputs will have to run the analyses (locally for SignalP, as long as you have installed it, or calling it remotely for Phobius) and this can take quite a lot of time.
