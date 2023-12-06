library(here)
# grab functions from src
source(here("src", "signal_peptides.r"))
source(here("src", "hydrophobicity.r"))
source(here("src", "motif_analysis.r"))
library(mixtools)

###
# The purpose of this script is to extract all protein IDs from all species in the proteins file
# that are divided into two groups, cleaved and non cleaved signal peptides.
# This is using phobius labels.
#
# This is also where the figure data is generated for the writeup (stored in results).
###


# attach path to protein file names
protein_paths <- base::Map(paste, here("data", "proteins", "pub"), list.files(here("data", "proteins", "pub")), sep = "/")
species_names <- gsub(".fasta", "", list.files(here("data", "proteins", "pub")))

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, r_phobius)

for (i in 1:length(phobius_results)) {
    phobius_results[[i]] <- phobius_results[[i]] %>% 
        filter(phobius_end != 0) %>%
        mutate(window_length = phobius_end - phobius_start) %>% 
        mutate(species = species_names[i])
}

# join and reset row names
phobius_df <- do.call(rbind, phobius_results)
rownames(phobius_df) <- NULL

species_order <- list("S_Cerevisiae" = 1,
              "C_Albicans" = 2,
              "N_Crassa" = 3,
              "P_Oryzae" = 4,
              "Z_Tritici" = 5,
              "A_Fumigatus" = 6,
              "S_Pombe" = 7,
              "P_Graminis" = 8,
              "U_Maydis" = 9,
              "C_Neoformans" = 10,
              "R_Delemar" = 11,
              "B_Dendrobatitis" = 12)
phobius_df$species <- factor(phobius_df$species, levels = names(species_order))

plot <- ggplot(phobius_df, aes(x = window_length, colour = phobius_type)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    facet_wrap(~species, scales = "free_y", ncol = 1) + 
    labs(x = "Window Length (AA)", y = "Density") + 
    scale_x_continuous(breaks = seq(0, 35, 10))

# save 
ggsave(here("results", "figures", "phobius_window_length.jpg"), plot, width = 15, height = 50, units = "cm")

# bar plot of above categories
ggplot(phobius_df, aes(x = phobius_type, fill = phobius_type)) + 
    geom_bar(stat="count") + 
    facet_wrap(~species, scales = "free_y") + 
    labs(x = "Window Length (AA)", y = "Count")


# write each phobius_type group of each species to a text file
for (i in 1:length(phobius_results)) {
    SP <- phobius_results[[i]] %>% 
        filter(phobius_type == "SP") %>% 
        pull(seqid)
    TM <- phobius_results[[i]] %>% 
        filter(phobius_type == "TM") %>% 
        pull(seqid)

    write.table(SP, file = paste(here("results", "proteins"), paste(species_names[i], "SP.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(TM, file = paste(here("results", "proteins"), paste(species_names[i], "TM.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# find GMM with k = 2 for SC_first_60
SC_first_60 <- phobius_df %>% 
    filter(species == "SC_first_60")

model <- normalmixEM(SC_first_60$window_length, k = 2)
intersection <- uniroot(function(x) {
    dnorm(x, mean = model$mu[1], sd = model$sigma[1]) - dnorm(x, mean = model$mu[2], sd = model$sigma[2])
}, c(min(model$mu), max(model$mu)))$root

ggplot(SC_first_60, aes(x = window_length)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    stat_function(fun = dnorm, args = list(mean = model$mu[1], sd = model$sigma[1]), colour = "blue") + 
    stat_function(fun = dnorm, args = list(mean = model$mu[2], sd = model$sigma[2]), colour = "green") + 
    labs(x = "Window Length (AA)", y = "Count") + 
    geom_vline(xintercept = intersection, colour = "red")

# use intersection to classify SP and TM
SP <- SC_first_60 %>%
    filter(window_length <= intersection) %>% 
    pull(seqid)
TM <- SC_first_60 %>%
    filter(window_length > intersection) %>% 
    pull(seqid)

write.table(SP, file = paste(here("results", "proteins"), paste("SC_window_length", "SP.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TM, file = paste(here("results", "proteins"), paste("SC_window_length", "TM.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)


screened_non_srp <- read_tsv(here("data", "SC_screened.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)

verified_srp <- read_tsv(here("data", "SC_SRP.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)

verified_non_srp <- read_tsv(here("data", "SC_non_SRP.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)


labelled_df <- SC_first_60 %>% 
    mutate(verified = case_when(seqid %in% verified_non_srp ~ "non SRP",
                           seqid %in% screened_non_srp ~ "non SRP",
                           seqid %in% verified_srp ~ "SRP",
                           TRUE ~ "unlabelled")) %>% 
    filter(verified != "unlabelled")

# make contingency table
contingency_table <- labelled_df %>% 
    group_by(verified) %>% 
    summarise(SP = sum(window_length <= intersection),
              TM = sum(window_length > intersection))

# chisq test
chisq.test(contingency_table[,2:3])

phobius_proteins_SC <- SC_first_60 %>%
    pull(seqid)

SC_AA <- proteins[[6]]
writeXStringSet(SC_AA[phobius_proteins_SC], file = here("results", "proteins", "SC_phobius.fasta"))

# export data to results for figures
scale_name <- "Rose"
cols <- c("aa", scale_name)
rose <- scales[, cols]
colnames(rose) <- c("V1", "V2")

protein_path <- here("data", "Proteins", "SC_first_60.fasta")
AA_stringset <- readAAStringSet(protein_path)

phobius_output <- r_phobius(protein_path)
signalp_output <- signalp(protein_path)
combined_df <- signalp_output %>%
    full_join(phobius_output, by = "seqid")

rose_df <- add_compound_hydropathy_score(combined_df, AA_stringset, useSignalP = FALSE, scale = KD, include_max = TRUE) %>% 
    drop_na(compound_hydropathy) %>% 
    mutate(window_length = window_end - window_start)

write_csv(rose_df, here("results", "figures", "SC_first_60.csv"))
