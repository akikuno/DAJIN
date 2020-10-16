
################################################################################
#! Install required packages
################################################################################

options(repos='https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, parallel)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# file_mids <- ".DAJIN_temp/clustering/temp/tmp_MIDS_barcode43_wt"
# file_mask <- ".DAJIN_temp/clustering/temp/tmp_mask_barcode43_wt"
# threads <- 14L

#===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_mids <- args[1]
file_mask <- args[2]
file_output <- args[3]
threads <- as.integer(args[4])

#===========================================================
#? Inputs
#===========================================================

df_mids <- read_csv(file_mids,
    col_names = FALSE,
    col_types = cols(.default = "c"))

df_mask <- read_csv(file_mask,
    col_names = "mask",
    col_types = cols())

################################################################################
#! Sequence error detection
# Sequence error is 90% or less for M or 5% or more for each IDS item.
################################################################################

detect_seqerror <- function(x){
    x %>%
    mutate(MIDS = case_when(
        Var1 %in% 1:9 ~ "I",
        str_detect(Var1, "[a-z]") ~ "I",
        TRUE ~ as.character(Var1)
    )) %>%
    group_by(MIDS) %>%
    summarize(MIDS_SUM = sum(Freq), .groups = "drop_last") %>%
    mutate(Percentage = MIDS_SUM / sum(MIDS_SUM) * 100) %>%
    filter((MIDS == "M" & Percentage < 90) | MIDS != "M" & Percentage > 5) %>%
    count %>%
    mutate(n = if_else(n == 0, 1, 2))
}

num_seqerror <- df_mids %>%
    mclapply(table, mc.cores = threads) %>%
    mclapply(as.data.frame, mc.cores = threads) %>%
    mclapply(detect_seqerror, mc.cores = threads) %>%
    unlist

logic_mask <- df_mask %>% pull(mask) %>% str_detect("[acgt]")
num_seqerror[logic_mask] <- 2

write_csv(as_tibble(num_seqerror), file_output, col_names = FALSE)