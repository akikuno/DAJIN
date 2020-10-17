################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, parallel, furrr)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# control <- "barcode43"
# file_control_mids <- sprintf(".DAJIN_temp/clustering/temp/tmp_MIDS_%s_wt", control)
# threads <- 14L
# plan(multisession, workers = threads)

# ===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_control_mids <- args[1]
threads <- as.integer(args[2])
plan(multisession, workers = threads)

#===========================================================
#? Inputs
#===========================================================

df_control_mids <- read_csv(file_control_mids,
    col_names = FALSE,
    col_types = cols())
colnames(df_control_mids) <- 1:ncol(df_control_mids)

################################################################################
#! MIDS scoring
################################################################################

df_control_freq_wt <-
    df_control_mids %>%
    pivot_longer(col = everything(), names_to = "loc", values_to = "MIDS") %>%
    nest(nest = c(MIDS)) %>%
    mutate(control_freq = future_map(nest,
        ~ .x %>% count(MIDS) %>% mutate(Freq = n / sum(n) * 100) %>% select(-n)
    )) %>%
    mutate(loc = as.double(loc)) %>%
    select(loc, control_freq)

################################################################################
#! Save results
################################################################################

saveRDS(df_control_freq_wt, ".DAJIN_temp/clustering/temp/df_control_freq_wt.RDS")

# percentage_MIDS <- function(x) {
#     x %>%
#     table(dnn = list("value")) %>%
#     as.data.frame(responseName = "freq") %>%
#     mutate(MIDS = case_when(
#         value %in% 1:9 ~ "I",
#         str_detect(value, "[a-z]") ~ "I",
#         TRUE ~ as.character(value)
#     )) %>%
#     group_by(MIDS) %>%
#     summarize(MIDS_SUM = sum(freq), .groups = "drop_last") %>%
#     mutate(Percentage = MIDS_SUM / sum(MIDS_SUM) * 100) %>%
#     select(-MIDS_SUM)
# }

# df_control_per <-
#     df_control_mids %>%
#     future_map(percentage_MIDS)

################################################################################
#! Save results
################################################################################

# saveRDS(df_control_per, ".DAJIN_temp/clustering/temp/df_control_per_wt.rds")