################################################################################
# ! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(readr, stringr, tibble, dplyr, tidyr, purrr)

################################################################################
# ! I/O naming
################################################################################

# ===========================================================
# ? Auguments
# ===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
threads <- as.integer(args[2])

# ===========================================================
# ? Inputs
# ===========================================================

df_wt <- readRDS(".DAJIN_temp/clustering/temp/df_control_freq_wt.RDS")

df_label <- read_csv(
    file_name,
    col_names = c("loc", "mut"),
    col_types = cols()
)

# ===========================================================
# ? Outputs
# ===========================================================

output_label <-
    file_name %>%
    str_remove_all(".*control_score_") %>%
    str_remove(".csv$")

################################################################################
# ! Detect allele type
################################################################################

logic_point_mutation <- FALSE
logic_inversion <- FALSE
logic_insertion <- FALSE
logic_deletion <- FALSE

if (sum(df_label$mut) == 1) {
    logic_point_mutation <- TRUE
} else if (any(df_label$mut == 2)) {
    logic_inversion <- TRUE
} else if (length(df_label$loc) > length(df_wt$loc)) {
    logic_insertion <- TRUE
} else {
    logic_deletion <- TRUE
}

# ===========================================================
# ? Point mutation
# ===========================================================

if (logic_point_mutation) {
    output_label <- "wt"

    mut_loc <- df_label %>%
        filter(mut == 1) %>%
        pull(loc)

    df_control_freq_label <- df_label

    for (num_label_loc in df_label$loc) {
        if (num_label_loc != mut_loc) {
            df_control_freq_label$mut[num_label_loc] <- df_wt$control_freq[num_label_loc]
        }
    }

    df_control_freq_label <-
        df_control_freq_label %>%
        rename(control_freq = mut) %>%
        mutate(mut = df_label$mut)
}

# ===========================================================
# ? Inversion
# ===========================================================

if (logic_inversion) {
    invert_loc <-
        df_label %>%
        filter(mut == 2) %>%
        arrange(desc(loc)) %>%
        pull(loc)

    df_label$loc[df_label$mut == 2] <- invert_loc

    df_control_freq_label <-
        df_wt %>%
        mutate(loc = df_label$loc, mut = df_label$mut) %>%
        arrange(loc)
}

# ===========================================================
# ? Insertion
# ===========================================================

if (logic_insertion) {
    df_control_freq_label <- df_label

    num_ref_loc <- 1
    for (num_label_loc in df_label$loc) {
        if (df_label$mut[num_label_loc] == 0) {
            df_control_freq_label$mut[num_label_loc] <- df_wt$control_freq[num_ref_loc]
            num_ref_loc <- num_ref_loc + 1
        }
    }
    df_control_freq_label <-
        df_control_freq_label %>%
        rename(control_freq = mut) %>%
        mutate(mut = df_label$mut)
}

# ===========================================================
# ? Deletion
# ===========================================================

if (logic_deletion) {
    df_control_freq_label <-
        df_wt[df_label$mut == 0, ] %>%
        mutate(loc = row_number(), mut = 0)
}

################################################################################
# ! Save results
################################################################################

saveRDS(
    df_control_freq_label,
    sprintf(".DAJIN_temp/clustering/temp/df_control_freq_%s.RDS", output_label)
)
