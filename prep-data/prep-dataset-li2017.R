### * Description

# Script to import into R the dataset from the article:

# Li, Lei, Clark J. Nelson, Josua Trösch, Ian Castleden, Shaobai Huang, and
# A. Harvey Millar. “Protein Degradation Rate in Arabidopsis Thaliana Leaf
# Growth and Development.” The Plant Cell 29, no. 2 (February 1, 2017):
# 207–28. https://doi.org/10.1105/tpc.16.00768.

# The data from this article is available on the Dryad depository:
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.q3h85

# The dataset citation is:
# Li, Lei, Clark J. Nelson, Josua Troesch, Ian Castleden, Shaobai Huang, and
# A. Harvey Millar. “Data from: Protein Degradation Rate in Arabidopsis
# Thaliana Leaf Growth and Development.” Dryad,
# 2018. https://doi.org/10.5061/DRYAD.Q3H85.


# Note that this script does not process the xlsx files from the Dryad zip
# archive directly. Instead, another script from this folder
# (`prep-dataset-li2017-from-xlsx.R`) is responsible for converting the xlsx
# files from Dryad into the appropriate tsv files that the present script
# processes.

# This is done so that the isotracer repository does not include the large xlsx
# files from the Dryad archive, but the more compact tsv files instead.


# In this study, sampling is done per leaf (l3, l5, l7), per day (t0, t1, t3,
# t5) and per replicate (3 replicates per sampling point).

# I understand that the actual samples used to generate total protein
# abundance, protein relatives abundances and protein labelling measurements
# are not exactly the same samples (a bit in the same way that we don't use the
# same invertebrate individuals to measure biomass and d15N in the Trinidad
# stream study).

# Note: proteins average abundances are provided in Sup. Table 7A of the paper,
# based on the PaxDB database.

### * Setup

suppressMessages(library(tidyverse))
suppressMessages(library(here))

### * Load the tsv files

# Function to adjust file paths
path <- function(f) {
    file.path(here::here(), "prep-data", f)
}

# Load the tsv files
design <- read_tsv(path("li-2017_design.tsv"),
                   col_types = cols(
                       sample = col_character(),
                       time_day = col_double(),
                       leaf_id = col_character()
                   ))
proteins <- read_tsv(path("li-2017_proteins.tsv"),
                     col_types = cols(
                         prot_id = col_character(),
                         description = col_character()
                     ))
prot_abundance <- read_tsv(path("li-2017_rel-abundances.tsv"),
                           col_types = cols(
                               prot_id = col_character(),
                               sample = col_character(),
                               ratio = col_double(),
                               n_quant = col_double(),
                               sd = col_double()
                           ))
prot_labelling <- read_tsv(path("li-2017_labelling.tsv"),
                           col_types = cols(
                               prot_id = col_character(),
                               sample = col_character(),
                               labeled_fraction = col_double(),
                               n_quant = col_double(),
                               sd = col_double()
                           ))

### * Merge design, relative abundance and labelling data

full <- bind_rows(
    prot_abundance %>% select(prot_id, sample, rel_abundance = ratio),
    prot_labelling %>% select(prot_id, sample, labeled_fraction)
    ) %>%
    left_join(design, by = "sample")

# Remove rows where there is no data
full <- full[!(is.na(full$rel_abundance) &
               is.na(full$labeled_fraction)), ]

### * Prepare a data count summary table (number of data points per protein)

n_abundance <- prot_abundance %>%
    select(prot_id, ratio) %>%
    group_by(prot_id) %>%
    summarize(n_abundance_data = sum(!is.na(ratio)), .groups = "drop")

n_labelling <- prot_labelling %>%
    select(prot_id, labeled_fraction) %>%
    group_by(prot_id) %>%
    summarize(n_labelling_data = sum(!is.na(labeled_fraction)), .groups = "drop")

prot_data_summary <- proteins %>%
    select(prot_id) %>%
    left_join(n_abundance, by = "prot_id") %>%
    left_join(n_labelling, by = "prot_id")
prot_data_summary$n_abundance_data[is.na(prot_data_summary$n_abundance_data)] <- 0
prot_data_summary$n_labelling_data[is.na(prot_data_summary$n_labelling_data)] <- 0

### * Save Li 2017 rda files

# Function to adjust file paths
pathRda <- function(f) {
    file.path(here::here(), "data", f)
}

# Save data

li2017 <- full
li2017_prots <- proteins
li2017_counts <- prot_data_summary

save(li2017, file = pathRda("li2017.rda"), version = 3, compress = "xz")
save(li2017_prots, file = pathRda("li2017_prots.rda"), version = 3, compress = "xz")
save(li2017_counts, file = pathRda("li2017_counts.rda"), version = 3, compress = "xz")
