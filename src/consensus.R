################################################################################
# Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(readr, stringr, tibble, dplyr, tidyr, purrr)

################################################################################
# I/O naming
################################################################################

# ===========================================================
# Auguments
# ===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que_mids <- args[1]
file_mutation_loc <- args[2]

# ===========================================================
# Inputs
# ===========================================================

df_que_mids <- read_csv(file_que_mids,
    col_names = FALSE,
    col_types = cols()
)

df_mutation_loc <- read_tsv(file_mutation_loc,
    col_names = c("loc"),
    col_types = cols()
)

# ===========================================================
# Outputs
# ===========================================================

output_suffix <- file_que_mids %>% str_remove(".*allele_id_")

################################################################################
# MIDS max count
################################################################################

max_count <- function(x) {
    y <- table(x)
    if (names(which.max(y)) != "S" && "S" %in% names(y)) {
        y["M"] <- y["M"] + y["S"]
    }
    y %>%
        which.max() %>%
        names()
}

colnames(df_que_mids) <- seq_len(ncol(df_que_mids))
df_que_max <-
    map_dfc(df_que_mids, max_count) %>%
    pivot_longer(cols = everything(), names_to = "loc", values_to = "mut") %>%
    mutate(loc = as.integer(loc))

# ===========================================================
# Summarize mutation position
# ===========================================================

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
# Sequence error detection
################################################################################

df_mut <-
    df_mut %>%
    inner_join(df_mutation_loc, by = "loc")

################################################################################
# ! Output results
################################################################################

write_delim(df_mut,
    sprintf(".DAJIN_temp/consensus/temp/mutation_%s", output_suffix),
    delim = " ",
    col_names = FALSE
)
