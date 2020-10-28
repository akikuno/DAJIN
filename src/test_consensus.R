################################################################################
#! Install required packages
################################################################################

options(repos = 'https://cloud.r-project.org/')
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# file_que_mids <- ".DAJIN_temp/consensus/temp/allele_id_barcode14_target_2"
# file_mutation_loc <- ".DAJIN_temp/clustering/temp/possible_true_mut_barcode14_target"

#===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que_mids <- args[1]
file_mutation_loc <- args[2]

#===========================================================
#? Inputs
#===========================================================

df_que_mids <- read_csv(file_que_mids,
    col_names = FALSE,
    col_types = cols())

df_mutation_loc <- read_tsv(file_mutation_loc,
    col_names = c("loc"),
    col_types = cols())

#===========================================================
#? Outputs
#===========================================================

output_suffix <- file_que_mids %>% str_remove(".*allele_id_")

################################################################################
#! MIDS max count
################################################################################

max_count <- function(x) {
    x %>% table() %>% which.max() %>% names()
}

df_que_max <-
    map_dfc(df_que_mids, max_count) %>%
    pivot_longer(cols = everything(), names_to = "loc", values_to = "mut") %>%
    mutate(loc = row_number())

#===========================================================
#? Summarize mutation position
#===========================================================

df_mut <-
    df_que_max %>%
    filter(mut != "M") %>%
    mutate(insnum = case_when(
        str_detect(mut, pattern = "S|D") ~ "0",
        TRUE ~ mut,
    )) %>%
    mutate(mut = case_when(
        !str_detect(mut, pattern = "S|D") ~ "I",
        TRUE ~ mut
    ))

################################################################################
#! Sequence error detection
################################################################################

df_mut <-
    df_mut %>%
    inner_join(df_mutation_loc, by = "loc")

################################################################################
#! Output results
################################################################################

write_delim(df_mut,
    sprintf(".DAJIN_temp/consensus/temp/mutation_%s", output_suffix),
    delim = " ",
    col_names = FALSE)
