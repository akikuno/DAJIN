################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, parallel, vroom)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# control <- "barcode42"
# file_control_mids <- sprintf(".DAJIN_temp/clustering/temp/tmp_MIDS_%s_wt", control)
# threads <- 14L

# ===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_control_mids <- args[1]
threads <- as.integer(args[2])

#===========================================================
#? Inputs
#===========================================================

df_control_mids <- vroom(file_control_mids,
    col_names = FALSE,
    col_types = cols(),
    num_threads = threads)
colnames(df_control_mids) <- seq_len(ncol(df_control_mids))

################################################################################
#! MIDS scoring
################################################################################

df_control_freq_wt <-
    df_control_mids %>%
    pivot_longer(col = everything(), names_to = "loc", values_to = "MIDS") %>%
    nest(nest = c(MIDS)) %>%
    mutate(control_freq = mclapply(nest,
        function(x)
            x %>% count(MIDS) %>% mutate(freq = n / sum(n) * 100),
        mc.cores = threads)) %>%
    mutate(loc = as.double(loc), mut = 0) %>%
    select(loc, control_freq, mut)

################################################################################
#! Save results
################################################################################

saveRDS(df_control_freq_wt, ".DAJIN_temp/clustering/temp/df_control_freq_wt.RDS")
