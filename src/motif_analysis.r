library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(ggseqlogo)


# ---- Description ----
# This script contains functions for finding motifs in an RNA/DNA sequence


# function that takes a 5'UTR FASTA returns a dataframe
# with the number of motifs found in each sequence
find_motifs <- function(input_path, motif = "CNYTCNYT") {
    # load fasta file as DNA string set
    input_RNA <- readDNAStringSet(input_path)

    # remove truncated sequences, which we don't need.
    input_RNA <- input_RNA[width(input_RNA) == 1000]

    # extract IDs
    input_ids <- names(input_RNA)

    # find motifs
    motif_matches <- tibble(seqid = input_ids,
                            count_up_1000 = vcountPattern(pattern = DNAString(motif),
                                    subject = input_RNA,
                                    fixed = "subject"),
                            count_up_200 = vcountPattern(pattern = DNAString(motif),
                                    subject = subseq(input_RNA, start = 801L, end = 1000L),
                                    fixed = "subject"),
                            count_up_100 = vcountPattern(pattern = DNAString(motif),
                                    subject = subseq(input_RNA, start = 901L, end = 1000L),
                                    fixed = "subject"))

    # return IDs
    return(motif_matches)
}