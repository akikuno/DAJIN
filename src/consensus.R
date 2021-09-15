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
colnames(df_que_mids) <- seq(ncol(df_que_mids))
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

df_mut <-
    df_que_mids %>%
    select(df_mutation_loc$loc) %>%
    pivot_longer(cols = everything(), names_to = "loc", values_to = "mut") %>%
    mutate(loc = as.double(loc)) %>%
    # convert insertion numbers to "I"
    mutate(mut_I = if_else(mut %in% c("M", "D", "S"), mut, "I")) %>%
    group_by(loc) %>%
    mutate(mut_I_max = max_count(mut_I)) %>%
    filter(mut_I_max != "M") %>%
    filter(mut_I == mut_I_max) %>%
    # calculate max freq of insertion numbers
    mutate(mut_max = max_count(mut)) %>%
    select(loc = loc, mut = mut_I, insnum = mut_max) %>%
    distinct() %>%
    mutate(insnum = if_else(insnum %in% c("D", "S"), "0", insnum)) %>%
    arrange(loc)

################################################################################
#! Output results
################################################################################

write_delim(df_mut,
    sprintf(".DAJIN_temp/consensus/temp/mutation_%s", output_suffix),
    delim = " ",
    col_names = FALSE)
