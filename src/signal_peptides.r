library(Biostrings)

signalp <- function(input_path, organism = "euk", mode = "fast", format = "none") {
    # extract file name from path, replace .fasta with _out
    file_name <- gsub(".fasta", "", basename(input_path))

    # create output path
    out_path <- paste0("results/", file_name)

    # check if output path exists, if it does, exit function
    if (!file.exists(out_path)) {
        # call signalp on file
        system(paste("signalp6 -ff", input_path, "-org", organism, "-od", out_path, "-fmt", format, "--mode", mode))
    }

    # read signalp output fasta
    signalp_output <- readAAStringSet(paste0(out_path, "/processed_entries.fasta"))

    #return gene IDs from stringset
    return(names(signalp_output))
}