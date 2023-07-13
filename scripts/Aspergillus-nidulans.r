library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(ggseqlogo)

# define fasta filename
As_up1000_file <- paste(here::here(), 
                        "motif_conservation",
                        "data", 
                        "Asp_Ni.fasta",
                        sep= "/")

# load fasta file as DNA string set
As_up1000_DSS <- readDNAStringSet(As_up1000_file)

# remove truncated sequences, which we don't need.
As_up1000_DSS <- As_up1000_DSS[width(As_up1000_DSS) == 1000]

# assign just the ORF id to the name, but keep all the info. 
As_namesinfo <- 
    tibble(everything = names(As_up1000_DSS)) %>%
    tidyr::separate(everything,into = c("id","strain","info","seqtype","seq","length"), sep = " \\| ") %>%
    dplyr::select(id, info,seq)

# print the DNA string set to check it came out ok
As_up1000_DSS

CNYTCNYT_count_df <- tibble(id = As_namesinfo$id,
                            info = As_namesinfo$info,
                            count_up1000 = vcountPattern(pattern = DNAString("CNYTCNYT"),
                                 subject = As_up1000_DSS,
                                 fixed = "subject"),
                            count_up200 = vcountPattern(pattern = DNAString("CNYTCNYT"),
                                 subject = subseq(As_up1000_DSS,start = 801L, end = 1000L),
                                 fixed = "subject"),
                            count_up100 = vcountPattern(pattern = DNAString("CNYTCNYT"),
                                 subject = subseq(As_up1000_DSS,start = 901L, end = 1000L),
                                 fixed = "subject")
)
summary(CNYTCNYT_count_df)
arrange(CNYTCNYT_count_df,desc(count_up100)) %>%
    head(n = 20)


# plots of counts to 1000bp upstream
CNYT_up1000_reducedmean <- CNYTCNYT_count_df %>%
    pull(count_up1000) %>%
    mean()
CNYT_up1000_reducedvar <- CNYTCNYT_count_df %>%
    pull(count_up1000) %>%
    var()
ggplot(data = CNYTCNYT_count_df, aes(x = count_up1000)) +
    geom_bar(width = 1, aes(y = ..prop..), stat = "count", group = 1) +
    geom_line(data = tibble(),
              aes(x = 0:20, y = dpois(0:20,lambda=CNYT_up1000_reducedmean)),
              colour = "blue") + 
    scale_x_continuous(limits = c(-0.5,16), oob = scales::squish)


# plots of counts to 100bp upstream
CNYT_up100_reducedmean <- CNYTCNYT_count_df %>%
    pull(count_up100) %>%
    mean()
CNYT_up100_reducedvar <- CNYTCNYT_count_df %>%
    pull(count_up100) %>%
    var()
ggplot(data = CNYTCNYT_count_df, aes(x = count_up100)) +
    geom_bar(width = 1, aes(y = ..prop..), stat = "count", group = 1) +
    geom_line(data = tibble(),
              aes(x = 0:20, y = dpois(0:20,lambda=CNYT_up100_reducedmean)),
              colour = "blue") + 
    scale_x_continuous(limits = c(-0.5,16), oob = scales::squish)


# saving ids of sequences with >= 1 motif to file
id_file_CNYTCNYT1_up100 <- paste(here::here(),
                                 "motif_conservation",
                                 "results",
                                 "An_id_list_CNYTCNYT1_up100.txt",
                                 sep = "/")
CNYTCNYT_count_df %>%
    filter(count_up100 >= 1) %>%
    pull(id) %>%
    write_lines(file = id_file_CNYTCNYT1_up100)
# saving ids of sequences with >= 2 motifs to file
id_file_CNYTCNYT2_up100 <- paste(here::here(),
                                 "motif_conservation",
                                 "results",
                                 "An_id_list_CNYTCNYT2_up100.txt",
                                 sep = "/")
CNYTCNYT_count_df %>%
    filter(count_up100 >= 2) %>%
    pull(id) %>%
    write_lines(file = id_file_CNYTCNYT2_up100)


# Create list of sequences with >= 1 motif in the first 100bp upstream
CNYTCNYT1_up100_DSS <- As_up1000_DSS[As_namesinfo$id %in% read_lines(id_file_CNYTCNYT1_up100)]
# Create list of sequences with >= 2 motifs in the first 100bp upstream
CNYTCNYT2_up100_DSS <- As_up1000_DSS[As_namesinfo$id %in% read_lines(id_file_CNYTCNYT2_up100)]


# Cut off the first 100bp of each sequence
CNYTCNYT1_up100_DSS <- subseq(CNYTCNYT1_up100_DSS, start = 901L)
# Generate a consensus sequence for the sequences with >= 1 motif
CNYTCNYT1_up100_consensus <- consensusMatrix(CNYTCNYT1_up100_DSS)


# Create a heatmap plot for the consensus matrix
heatmap(CNYTCNYT1_up100_consensus,
        Colv = NA,  # Turn off column clustering
        Rowv = NA,  # Turn off row clustering
        scale = "none",  # Use the original values without scaling
        col = colorRampPalette(c("white", "blue"))(100),  # Define the color palette
        xlab = "Position",
        ylab = "Base",
        main = "Consensus Matrix Heatmap")


# plot sequence logo for 100bp upstream
CNYTCNYT_vec <- as.character(CNYTCNYT1_up100_DSS, use.names = TRUE)
ggseqlogo(CNYTCNYT_vec)
# plot sequence logo for 15bp upstream
CNYTCNYT15_up100_DSS <- subseq(CNYTCNYT1_up100_DSS, start = 86L)
ggseqlogo(as.character(CNYTCNYT15_up100_DSS, use.names = TRUE))
# plot sequence logo for 15bp upstream for sequences with >= 2 motifs
CNYTCNYT15_up100_DSS <- subseq(CNYTCNYT2_up100_DSS, start = 986L)
ggseqlogo(as.character(CNYTCNYT15_up100_DSS, use.names = TRUE))

ASP_N_GO <- read_csv("motif_conservation/data/ASP_N.csv") %>%
    mutate(From = "ASP_Ni")
ASP_F_GO <- read_csv("motif_conservation/data/ASP_F.csv") %>%
    mutate(From = "ASP_Fu")

overlap <- intersect(ASP_N_GO$Name, ASP_F_GO$Name)

#plot bar chart of "Pct of bgd" for each GO name in overlap from both ASP_N_GO and ASP_F_GO
combined <- rbind(ASP_F_GO, ASP_N_GO)
combined <- combined[combined$Name %in% overlap,]

ggplot(combined, aes(x = reorder(Name, -`Pct of bgd`, FUN = min), y = `Pct of bgd`, group = From, fill = From)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "GO Name", y = "Pct of bgd", fill = "Type")


library(idpr)

# load protein FASTA file
AN8881_file <- readBStringSet("motif_conservation/data/proteins/A0A1U8QRJ5.fasta.txt")

# charge-hydropathy plot
AN8881_vec <- as.character(AN8881_file, use.names = TRUE)
chargeHydropathyPlot(AN8881_vec)

scaledHydropathyLocal(AN8881_vec, plotResults = TRUE, window = 9)

chargeCalculationGlobal(AN8881_vec, plotResults = TRUE)