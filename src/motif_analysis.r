library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(ggseqlogo)


# ---- Description ----
# This script contains functions for finding motifs in an RNA/DNA sequence


#function that takes a 5'UTR FASTA and saves a FASTA of 
#those that have CNYTCNYT motifs n number of times in the first l nucleotides
find_motifs <- function(input_path, output_path, rna_length = 1000L, n = 2, l = 100L, motif = "CNYTCNYT") {
    # check L and rna_length are integers
    if (!is.integer(l)) {
        l <- as.integer(l)
    }
    if (!is.integer(rna_length)) {
        rna_length <- as.integer(rna_length)
    }

    # load fasta file as DNA string set
    input_RNA <- readDNAStringSet(input_path)

    # remove truncated sequences, which we don't need.
    input_RNA <- input_RNA[width(input_RNA) == 1000]

    # extract IDs
    input_ids <- names(input_RNA)

    # find motifs
    motif_matches <- tibble(ID = input_ids,
                            count = vcountPattern(pattern = DNAString(motif),
                                    subject = subseq(input_RNA, start = (rna_length - l) + 1, end = rna_length),
                                    fixed = "subject")) %>%
                    filter(count >= n)

    # save fasta
    writeXStringSet(input_RNA[motif_matches$ID], output_path)
}
#find_motifs("data/RNA/Asp_Ni_999.fasta", "results/out1.fasta", rna_length = 1000L, n = 2, l = 100L, motif = "CNYTCNYT")
