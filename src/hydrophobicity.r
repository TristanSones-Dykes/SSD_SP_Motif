library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(idpr)
library(rvest)
library(reticulate)


# ---- Description ----
# This script contains functions for calculating and plotting 
# the hydrophobicity of a protein sequence(s) on a kyte-doolittle scale.
# It also includes functions for finding the hydropathy window of a protein
# using the phobius online tool (does not require it to be installed)


source_python(here("src", "phobius.py"))
KD <- data.frame(V1 = KDNorm$V1, V2 = round((KDNorm$V2 * 9) - 4.5, 1))


#variables for testing
#paper_df <- read_csv(here("data", "Proteins", "tslil_paper", "tslil_paper_protein_annotations.csv"))
#accession_numbers <- paper_df %>%
#    pull(`Accession Number`) %>%
#    writeLines(here("data", "Proteins", "tslil_paper", "tslil_paper_yeast_accession.txt"))
#yeast_secretome <- readAAStringSet(here("data", "Proteins", "tslil_paper", "tslil_protein_strings.fasta"))
#example_protein <- yeast_secretome[1]
#scaledHydropathyLocal(example_protein, window = 9, scale = "Kyte-Doolittle", plotResults = TRUE)

# ---- Functions ----


# function that calculates the local hydropathy of a protein
# from idpr package, modified to use non normalised weights
hydropathy_local <- function(sequence, window = 9, scale = "Kyte-Doolittle") {
    seqVector <- sequenceCheck(sequence = sequence, nonstandardResidues = c('X'), method = "stop", 
        outputType = "vector", suppressOutputMessage = TRUE)
    if ((window%%2) == 0) {
        stop("Window must be an odd number")
    }

    names(seqVector) <- NULL
    seqLength <- length(seqVector)
    numberResiduesAnalyzed <- seqLength - (window - 1)
    positionVector <- ((window - 1)/2 + 1):(seqLength - (window - 
        1)/2)
    centerResidueVector <- seqVector[positionVector]
    windowVector <- rep(NA, numberResiduesAnalyzed)
    scoreVector <- rep(NA, numberResiduesAnalyzed)
    for (i in seq_len(numberResiduesAnalyzed)) {
        sequenceWindow <- seqVector[i:(i + (window - 1))]
        windowVector[i] <- paste0(sequenceWindow, collapse = "")
        windowValues <- KD$V2[match(sequenceWindow, KD$V1)]
        scoreVector[i] <- sum(windowValues)/window
    }
    windowDF <- data.frame(Position = positionVector, Window = windowVector, 
        CenterResidue = centerResidueVector, WindowHydropathy = scoreVector)
        
    return(windowDF)
}
#hydropathy_local(example_protein, 9)



# function that calculates the mean hydropathy of a protein
mean_hydropathy <- function(sequence) {
    seqCharacterVector <- sequenceCheck(sequence = sequence, method = "stop", 
        outputType = "vector", suppressOutputMessage = TRUE)
    
    seqLength <- length(seqCharacterVector)
    scoreVector <- KD$V2[match(seqCharacterVector, KDNorm$V1)]
    mean(scoreVector)
}
#mean_hydropathy(example_protein)



# function that uses phobius to find the hydropathy window of a protein
find_hydropathy_windows <- function(protein_AAStringSet, isString = FALSE) {

    # if it is a string, write it to a file
    if (isString) {
        writeXStringSet(protein_AAStringSet, file = "temp.fasta")
        protein_AAStringSet <- "temp.fasta"
    }
    
    # call phobius on file
    hydropathy_windows <- phobius(protein_AAStringSet)

    # if it is a string, delete the file
    if (isString) {
        file.remove("temp.fasta")
    }

    return(hydropathy_windows)
}
#find_hydropathy_windows(example_protein, isString = TRUE)
#find_hydropathy_windows(here("data", "Proteins", "tslil_paper", "ten_tslil_protein_strings.fasta"))



# function that takes in a protein AA string and returns a plot of the hydropathy window
plot_hydropathy_window <- function(protein_AAStringSet, window_size = 9) {
    # use phobius to find the hydropathy window
    hydropathy_window <- c(find_hydropathy_windows(protein_AAStringSet, isString = TRUE)[1,])

    hydropathy_df <- hydropathy_local(protein_AAStringSet)
    max_hydropathy_row <- hydropathy_df %>%
        filter(WindowHydropathy == max(WindowHydropathy))

    ggplot(hydropathy_df, aes(x = Position, y = WindowHydropathy)) + 
        geom_path() + 
        geom_hline(yintercept = mean_hydropathy(protein_AAStringSet), linetype = "dashed", colour = "grey") +
        geom_vline(xintercept = hydropathy_window$start, linetype = "dashed", colour = "red") +
        geom_vline(xintercept = hydropathy_window$end, linetype = "dashed", colour = "red") +
        geom_point(data = max_hydropathy_row, aes(x = Position, y = WindowHydropathy), colour = "black") +
        geom_text(data = max_hydropathy_row, aes(x = Position, y = WindowHydropathy, label = WindowHydropathy), vjust = -1) +
        labs(x = "Position", y = "Hydropathy", title = "Hydropathy of protein", colour = "Location")
}    
#plot_hydropathy_window(example_protein, 9)



# calculates the compund hydropathy score of a protein
# by multiplying the maximum hydropathy score by the length of the window
compound_hydropathy_score <- function(protein_AAStringSet, window_size = 9) {
    hydropathy_window <- c(find_hydropathy_windows(protein_AAStringSet, isString = TRUE)[1,])
    max_hydropathy_row <- hydropathy_local(protein_AAStringSet) %>%
        filter(WindowHydropathy == max(WindowHydropathy))
    
    return(max_hydropathy_row$WindowHydropathy * ((hydropathy_window$end - hydropathy_window$start) + 1))
}
#compound_hydropathy_score(example_protein, 9)

r_phobius <- function(protein_AA_path) {
    if (exists("phobius_output")) {
        print("phobius_output already exists, checking if the same...")
        if (length(phobius_output) == length(readAAStringSet(protein_AA_path))) {
            print("same length, not running again")
            return(phobius_output)
        }
    }

    # call phobius on file
    hydropathy_windows <- phobius(protein_AA_path) %>%
        mutate(type = case_when(
            type == "TRANSMEM" ~ "TM",
            type == "SIGNAL" ~ "SP",
            type == "NONE" ~ "OTHER"
        ))
    colnames(hydropathy_windows) <- c("seqid", "phobius_start", "phobius_end", "phobius_type")

    return(hydropathy_windows)
}
