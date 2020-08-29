################################################################################
#! Install required packages
################################################################################

options(repos='https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# file_que <- ".DAJIN_temp/consensus/temp/allele_id_barcode14_target_4"
# file_control <- ".DAJIN_temp/clustering/temp/control_score_target"
# barcode <- "barcode14"

#===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que <- args[1]
file_control <- args[2]
barcode <- args[3]

#===========================================================
#? Inputs
#===========================================================

df_que <- read_csv(file_que,
    col_names = FALSE,
    col_types = cols(.default = "c"))

df_control <- read_csv(sprintf("%s_%s", file_control, barcode),
    col_names = c("score"),
    col_types = cols())

#===========================================================
#? Outputs
#===========================================================

output_suffix <- file_que %>% str_remove(".*allele_id_")

################################################################################
#! Sequence error compensation
################################################################################

df_que[, pull(df_control) == 2] <- "M"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! MIDS count
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

max_count <- function(x){
    x %>% table() %>% which.max() %>% names()
}

df_que_max <- lapply(df_que, max_count) %>% unlist

#===========================================================
#? Summarize mutation position
#===========================================================

df_mut <- tibble(loc = which(df_que_max != "M"),
    mut = df_que_max[df_que_max != "M"] %>% as.character) %>%
    mutate(insnum = as.character(mut))

df_mut$insnum[str_detect(df_mut$insnum, pattern = "S|D")] <- 0
df_mut$mut[!str_detect(df_mut$mut, pattern = "S|D")] <- "I"

################################################################################
#! Output results
################################################################################

write_delim(df_mut,
    sprintf(".DAJIN_temp/consensus/temp/mutation_%s", output_suffix),
    delim = " ",
    col_names = FALSE)