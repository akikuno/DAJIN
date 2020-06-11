#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Install required packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos='https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Argument
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# df_que <- read_csv(".DAJIN_temp/clustering/temp/allele_id_barcode23_target_2", col_names = FALSE, col_types = cols(.default = "c"))
# control <- read_csv(".DAJIN_temp/clustering/temp/control_score_target", col_names = FALSE, col_types = cols())
# cluster <- 1
# suffix <- ".DAJIN_temp/clustering/temp/allele_id_barcode23_target_2" %>% str_remove(".*allele_id_")

args <- commandArgs(trailingOnly = TRUE)
df_que <- read_csv(args[1], col_names = FALSE, col_types = cols(.default = "c"))
control <- read_csv(args[2], col_names = FALSE, col_types = cols())
cluster <- args[3]
suffix <- args[1] %>% str_remove(".*allele_id_")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! MIDS count
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

max_count <- function(x){
    x %>% table() %>% which.max() %>% names()
}
df_que <- apply(df_que, 2, max_count)

df_que[pull(control) == 2] <- "M"

df_mut <- tibble(loc = which(df_que != "M"),
    mut = df_que[df_que != "M"] %>% as.character()) %>%
    mutate(insnum = as.character(mut))

df_mut$insnum[str_detect(df_mut$insnum, pattern = "S|D")] <- 0
df_mut$mut[!str_detect(df_mut$mut, pattern = "S|D")] <- "I"

write_delim(df_mut, sprintf(".DAJIN_temp/clustering/temp/mutation_%s", suffix),
    delim=" ", col_names = FALSE)