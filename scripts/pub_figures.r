library(here)
library(tidyverse)
library(cowplot)
library(ggExtra)
library(ggplotify)
library(vcd)
theme_set(theme_cowplot())

###
# The purpose of this script is to generate all figures for the writeup
# It uses post-processed data in results so that it can be run without
# signalP and phobius related dependencies.
###

# Load data
rose_df <- read_csv(here("results", "figures", "SC_first_60.csv"))

screened_non_srp <- read_tsv(here("data", "SC_screened.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)
verified_srp <- read_tsv(here("data", "SC_SRP.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)
verified_non_srp <- read_tsv(here("data", "SC_non_SRP.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)

labelled_df <- rose_df %>% 
    mutate(`Experimental label` = case_when(seqid %in% verified_non_srp ~ "Cleaved SP",
                           seqid %in% screened_non_srp ~ "Cleaved SP",
                           seqid %in% verified_srp ~ "Non-cleaved SP",
                           TRUE ~ "unlabelled"))

verified_df <- labelled_df %>% 
    filter(`Experimental label` != "unlabelled")

# plot x axis (0, 35)
lower <- 5 
upper <- 35
x_delim <- seq(lower, upper, 10)
x_minor <- seq(lower, upper, 5)
lims <- c(lower, upper)

# Figure 1A - Image of cleaved vs non cleaved SP

# placeholder
# https://www.mdpi.com/1422-0067/22/21/11871
Fig_1A <- ggplot(data.frame(x = rnorm(100)), aes(x = x)) + 
    geom_histogram() + 
    ggtitle("Placeholder")

# Figure 1B - Plot of contingency table of SP/TM regions found by phobius
# and those verified experimentally

# make contingency table
contingency_table <- verified_df %>% 
    group_by(`Experimental label`) %>% 
    summarise(SP = sum(window_type == "TM"),
              TM = sum(window_type == "SP"))
contingency_table <- as.table(as.matrix(contingency_table[,2:3]))

names(dimnames(contingency_table)) <- c("Experimental label", "Phobius label")
rownames(contingency_table) <- c("Cleaved SP", "Non-cleaved SP")

# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
Fig_1B <- as.grob(~vcd::mosaic(contingency_table, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))


# FIGURE 1C - Scatter plot of only experimentally verified proteins
# with histograms on axis

base <- ggplot(verified_df, aes(x = window_length, y = max_hydropathy, colour = `Experimental label`, group = `Experimental label`)) +
    geom_point() + 
    xlab("Window length (aa)") + 
    ylab("Max hydropathy (Rose)") + 
    ggtitle("Lab verified SP/TM proteins") + 
    scale_x_continuous(breaks = x_delim, limits = lims, minor_breaks = x_minor)

Fig_1C <- ggMarginal(base, type = "histogram", groupColour = TRUE, groupFill = TRUE)


# Figure 1D - Scatter plot of all proteins with SP/TM regions
# found by phobius

base <- ggplot(labelled_df, aes(x = window_length, y = max_hydropathy, colour = `Experimental label`, group = `Experimental label`)) +
    geom_point() + 
    xlab("Window length (aa)") + 
    ylab("Max hydropathy (Rose)") + 
    ggtitle("Phobius detected SP/TM regions") + 
    scale_x_continuous(breaks = x_delim, limits = lims, minor_breaks = x_minor)

Fig_1D <- ggMarginal(base, type = "histogram", groupColour = TRUE, groupFill = TRUE)


# combine figures into 2x2 grid
Fig_1 <- plot_grid(Fig_1A, Fig_1B, Fig_1C, Fig_1D, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

# save figure
ggsave(here("results", "figures", "Fig_1.jpg"), Fig_1, width = 15, height = 10, dpi = 300)
