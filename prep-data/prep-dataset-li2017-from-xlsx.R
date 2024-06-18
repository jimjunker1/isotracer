### * Description

# Script to process the protein abundance data and labeled fraction data from
# the article:

# Li, Lei, Clark J. Nelson, Josua Trösch, Ian Castleden, Shaobai Huang, and
# A. Harvey Millar. “Protein Degradation Rate in Arabidopsis Thaliana Leaf
# Growth and Development.” The Plant Cell 29, no. 2 (February 1, 2017):
# 207–28. https://doi.org/10.1105/tpc.16.00768.

# The data from this article is available on the Dryad depository:
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.q3h85

# The dataset was downloaded on 2020-08-14 (file
# `doi_10.5061_dryad.q3h85__v1.zip`, md5sum 75c0e09f389c771bd195206440b45710).

# The dataset citation is:
# Li, Lei, Clark J. Nelson, Josua Troesch, Ian Castleden, Shaobai Huang, and
# A. Harvey Millar. “Data from: Protein Degradation Rate in Arabidopsis
# Thaliana Leaf Growth and Development.” Dryad,
# 2018. https://doi.org/10.5061/DRYAD.Q3H85.

# The License of the dataset on Dryad is CC0 1.0 (Public Domain).
# Link: https://creativecommons.org/publicdomain/zero/1.0/

# This script assumes that the xlsx files from doi_10.5061_dryad.q3h85__v1.zip
# are present in the `./prep-data` folder. Those xlsx files are not included in
# the isotracer repository because of their large size. However, the product of
# this R script (a set of csv files) is included in the repository and
# version-controlled.

# Notes about the xlsx files:
# - Dataset 1 contains information about total protein abundance in the leaves
#   (this uses different samples from the 15N-labeled leaves).
# - Dataset 9 contains the abundances for each protein (relative to a common 15N
#   reference sample).
# - Dataset 14 contains information about labeled fractions for each protein,
#   for each sample (with SD per protein).


# In short: to reproduce the conversion from xlsx to csv files, download and
# unzip doi_10.5061_dryad.q3h85__v1.zip using the links above and run this
# script with the xlsx files in the `./prep-data` folder.

### * Setup

suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(readxl))

PROT_ABUNDANCE_XLSX_FILE <- file.path(here::here(), "prep-data",
  "TPC2016-00768-LSBR1_Supplemental_Data_Set_9.xlsx")
PROT_LABELLING_XLSX_FILE <- file.path(here::here(), "prep-data",
  "TPC2016-00768-LSBR1_Supplemental_Data_Set_14.xlsx")

### * Build the design table (mapping sample name to sampling design)

design <- tibble(
    sample = c("T0L3R1", "T0L3R2", "T0L3R3", "T0L5R1", "T0L5R2", "T0L5R3",
               "T0L7R1", "T0L7R2", "T0L7R3", "T1L3R1", "T1L3R2", "T1L3R3",
               "T1L5R1", "T1L5R2", "T1L5R3", "T1L7R1", "T1L7R2", "T1L7R3",
               "T3L3R1", "T3L3R2", "T3L3R3", "T3L5R1", "T3L5R2", "T3L5R3",
               "T3L7R1", "T3L7R2", "T3L7R3", "T5L3R1", "T5L3R2", "T5L3R3",
               "T5L5R1", "T5L5R2", "T5L5R3", "T5L7R1", "T5L7R2", "T5L7R3"),
    time_day = c(rep(c(0, 1, 3, 5), each = 9)),
    leaf_id = rep(rep(c("leaf_3", "leaf_5", "leaf_7"), each = 3), 4))

### * Prepare protein abundance data

### ** Documentation

# XLSX file: "TPC2016-00768-LSBR1_Supplemental_Data_Set_9.xlsx"

# The data is extracted from the spreadsheet "SD9B Protein Abundance Summary".

# The first line describes the content: "3632 proteins (p>0.95) from Arabidopsis
# leaf 3, leaf 5 or leaf 7 samples with a quantified protein abundance reported
# as a ratio to equivalent 15N reference peptide abundance.  AGI and
# description are from TAIR (Arabidopsis.org). See Supplemental Figure 2B for a
# flow diagram of the total experiment to determine protein abundance fold
# change relative to a 15N reference in 36 individual leaf samples. Each sample
# is identified by the timepoint (T), leaf number (L) and biological replicate
# (R), #quanted is the number of peptides quantified to make the ratio
# shown,and SD is the standard deviation across the quantified peptides for a
# protein."

### ** Load data

# Check that the expected sheet is present
sheets <- excel_sheets(PROT_ABUNDANCE_XLSX_FILE)
stopifnot("SD9B Protein Abundance Summary" %in% sheets)

# Load data and check that expected columns are present
d <- read_excel(PROT_ABUNDANCE_XLSX_FILE, "SD9B Protein Abundance Summary",
                skip = 1)
expected_columns <- c("AGI", "description", "...3", "T0L3R1", "#Quanted...5", "SD...6", 
  "T0L3R2", "#Quanted...8", "SD...9", "T0L3R3", "#Quanted...11", 
  "SD...12", "T0L5R1", "#Quanted...14", "SD...15", "T0L5R2", "#Quanted...17", 
  "SD...18", "T0L5R3", "#Quanted...20", "SD...21", "T0L7R1", "#Quanted...23", 
  "SD...24", "T0L7R2", "#Quanted...26", "SD...27", "T0L7R3", "#Quanted...29", 
  "SD...30", "T1L3R1", "#Quanted...32", "SD...33", "T1L3R2", "#Quanted...35", 
  "SD...36", "T1L3R3", "#Quanted...38", "SD...39", "T1L5R1", "#Quanted...41", 
  "SD...42", "T1L5R2", "#Quanted...44", "SD...45", "T1L5R3", "#Quanted...47", 
  "SD...48", "T1L7R1", "#Quanted...50", "SD...51", "T1L7R2", "#Quanted...53", 
  "SD...54", "T1L7R3", "#Quanted...56", "SD...57", "T3L3R1", "#Quanted...59", 
  "SD...60", "T3L3R2", "#Quanted...62", "SD...63", "T3L3R3", "#Quanted...65", 
  "SD...66", "T3L5R1", "#Quanted...68", "SD...69", "T3L5R2", "#Quanted...71", 
  "SD...72", "T3L5R3", "#Quanted...74", "SD...75", "T3L7R1", "#Quanted...77", 
  "SD...78", "T3L7R2", "#Quanted...80", "SD...81", "T3L7R3", "#Quanted...83", 
  "SD...84", "T5L3R1", "#Quanted...86", "SD...87", "T5L3R2", "#Quanted...89", 
  "SD...90", "T5L3R3", "#Quanted...92", "SD...93", "T5L5R1", "#Quanted...95", 
  "SD...96", "T5L5R2", "#Quanted...98", "SD...99", "T5L5R3", "#Quanted...101", 
  "SD...102", "T5L7R1", "#Quanted...104", "SD...105", "T5L7R2", 
  "#Quanted...107", "SD...108", "T5L7R3", "#Quanted...110", "SD...111", 
  "Brep NO", "Pep NO")
stopifnot(all(colnames(d) == expected_columns))

### ** Process time points

# Extract data one time point at a time
time_points <- colnames(d)[grepl("T[0-9]L[0-9]R[0-9]", colnames(d))]
data_tp <- list()
for (tp in time_points) {
    ratio_i <- which(colnames(d) == tp)
    n_quanted_i <- ratio_i + 1
    sd_i <- ratio_i + 2
    data <- cbind(d[, c("AGI", "description")],
                  d[, c(ratio_i, n_quanted_i, sd_i)])
    colnames(data) <- c("prot_id", "description", "ratio", "n_quant", "sd")
    data[["sample"]] <- tp
    data_tp[[tp]] <- as_tibble(data) %>%
        select(prot_id, description, sample, ratio, n_quant, sd)
}

prot_abundance <- bind_rows(data_tp) %>%
    arrange(prot_id)

### ** Checks

# Check that all sample names are known
stopifnot(all(prot_abundance$sample %in% design$sample))

### * Prepare protein labelling data

### ** Documentation

# XLSX file: "TPC2016-00768-LSBR1_Supplemental_Data_Set_14.xlsx"

# The data is extracted from a series of spreadsheets stored in the xlsx
# file. The documentation is contained in the spreadsheet "SD14D LPF in 27
# samples" in the same xlsx file.

# This is the content of this documentation spreadsheet:

# "Quantification of 27 individual sample peptide light to heavy ratios for
# determination of the labelled peptide fraction (LPF) for peptides by 15N
# progressively labelling in leaf samples. See Supplemental Figure 2 for flow
# diagram of the sampling and the mass spectrometry experiment. The legend
# provides explanation of each column.
#
# Protein         Arabidopsis Protein Identifier (AGI)
# Prot Prob.      TPP generated protein probability for the protein identification
# Unique Peps     No of unique peptides for this protein identificed
# Desc            TAIR Description of this protein
# Pep             Amino acid sequence of the peptide used for quantification
# LPF Pep         The labeled protein fraction measure for the peptide
# Pep enrichment  The 15N enrichment level of the heavy peptide fraction
# LPF Prot        The labeled protein fraction measure for the protein
# SD              standard deviation for the protein level LPF average
# #Quanted        Number of peptide quantitations used for the protein average
# 
# T1    Time 1, 22 day old plant leaf, 1 day following shift to 15N
# T3    Time 3, 24 day old plant leaf,  3 day following shift to 15N
# T5    Time five, 26 day old plant leaf,  5 days following shift to 15N
# L3    Leaf 3
# L5    Leaf 5
# L7    Leaf 7
# R1    replicate 1
# R2    replicate 2
# R3    replicate 3"
# 

### ** Load data

# Check that the expected worksheets are present
expected_sheets <- c("T1L3R1", "T1L3R2", "T1L3R3", "T1L5R1", "T1L5R2",
                     "T1L5R3", "1DL7R1", "1DL7R2", "1DL7R3", "T3L3R1",
                     "T3L3R2", "T3L3R3", "T3L5R1", "T3L5R2", "T3L5R3",
                     "T3L7R1", "T3L7R2", "T3L7R3", "T5L3R1", "T5L3R2",
                     "T5L3R3", "T5L5R1", "T5L5R2", "T5L5R3", "T5L7R1",
                     "T5L7R2", "T5L7R3")
actual_sheets <- excel_sheets(PROT_LABELLING_XLSX_FILE)
stopifnot(all(expected_sheets %in% actual_sheets))

# Load sheets one by one
data_sheets <- list()
expected_cols <- c("Protein", "Prot Prob.", "Unique Peps", "Desc", "Pep",
                   "LPF Pep", "Pep enrichment", "LPF Prot", "SD", "#Quanted")
for (sheet in expected_sheets) {
    x <- read_excel(PROT_LABELLING_XLSX_FILE, sheet)
    stopifnot(all(colnames(x) == expected_cols))
    x <- unique(x[, c("Protein", "Desc", "LPF Prot", "SD", "#Quanted")])
    colnames(x) <- c("prot_id", "description", "labeled_fraction", "sd", "n_quant")
    stopifnot(length(table(table(x$prot_id))) == 1)
    x[["sample"]] <- sheet
    data_sheets[[sheet]] <- x[, c("prot_id", "sample", "labeled_fraction", "n_quant", "sd", "description")]
}

prot_labelling <- bind_rows(data_sheets) %>%
    arrange(prot_id)

### ** Fix sample names

# We assume that some spreadsheet names were a mistake and correct them as
# follows:
corrections <- c("1DL7R1" = "T1L7R1", "1DL7R2" = "T1L7R2", "1DL7R3" = "T1L7R3")
i <- prot_labelling[["sample"]] %in% names(corrections)
prot_labelling[["sample"]][i] <- corrections[prot_labelling[["sample"]][i]]

### * Prepare tidy tables

### ** Separate tables into tidy tables

proteins <- bind_rows(
    prot_abundance %>% select(prot_id, description) %>% unique(),
    prot_labelling %>% select(prot_id, description) %>% unique()) %>%
    unique() %>%
    arrange(prot_id, description)

prot_abundance <- prot_abundance %>%
    select(- description) %>%
    arrange(prot_id, sample)

prot_labelling <- prot_labelling %>%
    select(- description) %>%
    arrange(prot_id, sample)

### * Save the tidy tables

# Function to adjust file paths
path <- function(f) {
    file.path(here::here(), "prep-data", f)
}

# Write tsv files
write_tsv(design, path("li-2017_design.tsv"))
write_tsv(proteins, path("li-2017_proteins.tsv"))
write_tsv(prot_abundance, path("li-2017_rel-abundances.tsv"))
write_tsv(prot_labelling, path("li-2017_labelling.tsv"))
