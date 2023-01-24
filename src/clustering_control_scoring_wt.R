################################################################################
# ! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(readr, stringr, tibble, dplyr, tidyr, purrr, parallel, tidyfast)

################################################################################
# ! I/O naming
################################################################################

# ===========================================================
# ? Auguments
# ===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_control_mids <- args[1]
threads <- as.integer(args[2])

# ===========================================================
# ? Inputs
# ===========================================================

df_control_mids <- read_csv(file_control_mids,
    col_names = FALSE,
    col_types = cols(),
    num_threads = threads
)
colnames(df_control_mids) <- seq_len(ncol(df_control_mids))

################################################################################
# ! MIDS scoring
################################################################################

df_control_freq_wt <-
    df_control_mids %>%
    dt_pivot_longer(names_to = "loc", values_to = "MIDS") %>%
    nest(nest = c(MIDS)) %>%
    mutate(control_freq = mclapply(nest,
        function(x) {
              x %>%
                  count(MIDS) %>%
                  mutate(freq = n / sum(n) * 100)
          },
        mc.cores = threads
    )) %>%
    mutate(loc = as.double(loc), mut = 0) %>%
    select(loc, control_freq, mut)

################################################################################
# ! Save results
################################################################################

saveRDS(df_control_freq_wt, ".DAJIN_temp/clustering/temp/df_control_freq_wt.RDS")
