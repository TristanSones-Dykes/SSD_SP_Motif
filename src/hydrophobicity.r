library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(idpr)
library(rvest)


KD <- data.frame(V1 = KDNorm$V1, V2 = round((KDNorm$V2 * 9) - 4.5, 1))
paper_df <- read_csv(here("data", "Proteins", "tslil_paper_protein_annotations.csv"))

accession_numbers <- paper_df %>%
    pull(`Accession Number`) %>%
    writeLines(here("data", "Proteins", "tslil_paper_yeast_accession.txt"))

yeast_secretome <- readAAStringSet(here("data", "Proteins", "tslil_protein_strings.fasta"))

example_protein <- yeast_secretome[1]
scaledHydropathyLocal(example_protein, window = 9, scale = "Kyte-Doolittle", plotResults = TRUE)

# function that calculates the local hydropathy of a protein
# from idpr package, modified to use non normalised weights
hydropathy_local <- function(sequence, window = 9, scale = "Kyte-Doolittle") {
    seqVector <- sequenceCheck(sequence = sequence, method = "stop", 
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
hydropathy_local(example_protein, 9)

mean_hydropathy <- function(sequence) {
    seqCharacterVector <- sequenceCheck(sequence = sequence, method = "stop", 
        outputType = "vector", suppressOutputMessage = TRUE)
    
    seqLength <- length(seqCharacterVector)
    scoreVector <- KD$V2[match(seqCharacterVector, KDNorm$V1)]
    mean(scoreVector)
}
mean_hydropathy(example_protein)

# function that takes in a protein AA string and returns a list
# with the hydropathy dataframe, the left and right limits of the window
find_hydropathy_window <- function(protein_AAStringSet, window_size = 9) {
    # calculate local and mean hydropathy
    hydropathy_df <- hydropathy_local(protein_AAStringSet, window = window_size, scale = "Kyte-Doolittle")
    limit <- meanScaledHydropathy(protein_AAStringSet)

    #pull out the maximum hydropathy row and values
    max_hydropathy_row <- hydropathy_df %>%
        filter(WindowHydropathy == max(WindowHydropathy))
    maximum <- max_hydropathy_row$WindowHydropathy
    centre <- max_hydropathy_row$Position

    left_len <- 1
    right_len <- 1

    # work backwards
    while ((hydropathy_df[which(hydropathy_df$Position == centre - left_len),]$WindowHydropathy > limit)
           && (centre - left_len >= 0)) {
        left_len <- left_len + 1
    }
    # work forwards
    while ((hydropathy_df[which(hydropathy_df$Position == centre + right_len),]$WindowHydropathy > limit)
           && (centre + right_len <= nrow(hydropathy_df))) {
        right_len <- right_len + 1
    }

    return(list(dat = hydropathy_df %>%
               mutate(isWindow = case_when(
                     (Position < centre - left_len) ~ "Outside",
                     (Position > centre + right_len) ~ "Outside",
                     TRUE ~ "Inside",
                )),
             left = centre - left_len,
             right = centre + right_len))
}
example_hydropathy_window <- find_hydropathy_window(example_protein, 9)
example_hydropathy_window

# function that takes in a protein AA string and returns a plot of the hydropathy window
plot_hydropathy_window <- function(protein_AAStringSet, window_size = 9) {
    hydropathy_window <- find_hydropathy_window(protein_AAStringSet, window_size)
    max_hydropathy_row <- hydropathy_window$dat %>%
        filter(WindowHydropathy == max(WindowHydropathy))

    ggplot(hydropathy_window$dat, aes(x = Position, y = WindowHydropathy, colour = isWindow)) + 
        geom_path(aes(group = 1)) + 
        geom_hline(yintercept = mean_hydropathy(protein_AAStringSet), linetype = "dashed", colour = "grey") +
        geom_vline(xintercept = hydropathy_window$left, linetype = "dashed", colour = "red") +
        geom_vline(xintercept = hydropathy_window$right, linetype = "dashed", colour = "red") +
        geom_point(data = max_hydropathy_row, aes(x = Position, y = WindowHydropathy), colour = "black") +
        geom_text(data = max_hydropathy_row, aes(x = Position, y = WindowHydropathy, label = WindowHydropathy), vjust = -1) +
        labs(x = "Position", y = "Hydropathy", title = "Hydropathy of protein", colour = "Location")
}    
plot_hydropathy_window(example_protein, 9)

# calculates the compund hydropathy score of a protein
# by multiplying the maximum hydropathy score by the length of the window
compound_hydropathy_score <- function(protein_AAStringSet, window_size = 9) {
    hydropathy_window <- find_hydropathy_window(protein_AAStringSet, window_size)
    max_hydropathy_row <- hydropathy_window$dat %>%
        filter(WindowHydropathy == max(WindowHydropathy))

    
    return(max_hydropathy_row$WindowHydropathy * (hydropathy_window$right - hydropathy_window$left))
}
compound_hydropathy_score(example_protein, 9)

string <- "ACDEFGHIKLMNPQRSTVWY"
hydro_df <- scaledHydropathyGlobal(string, scale = "Kyte-Doolittle")
hydro_df$mine <- mean_hydropathy(string)

# write first 10 proteins to file
writeXStringSet(yeast_secretome[1:10], file = here("data", "Proteins", "ten_tslil_protein_strings.fasta"))
