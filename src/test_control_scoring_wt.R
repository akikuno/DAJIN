################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, parallel)

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
    mutate(control_freq = mclapply(nest,
        function(x)
            x %>% count(MIDS) %>% mutate(Freq = n / sum(n) * 100) %>% select(-n),
        mc.cores = threads)) %>%
    mutate(loc = as.double(loc)) %>%
    select(loc, control_freq)

################################################################################
#! Save results
################################################################################

saveRDS(df_control_freq_wt, ".DAJIN_temp/clustering/temp/df_control_freq_wt.RDS")
