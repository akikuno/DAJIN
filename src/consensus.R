#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Install required packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos='https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Argument
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# que_name <- ".DAJIN_temp/consensus/temp/allele_id_barcode23_target_3"
# df_que <- read_csv(que_name, col_names = FALSE, col_types = cols(.default = "c"))
# df_control <- read_csv(".DAJIN_temp/clustering/temp/control_score_target", col_names = c("score"), col_types = cols())
# cluster <- 1
# suffix <- que_name %>% str_remove(".*allele_id_")

args <- commandArgs(trailingOnly = TRUE)
df_que <- read_csv(args[1], col_names = FALSE, col_types = cols(.default = "c"))
df_control <- read_csv(args[2], col_names = c("score"), col_types = cols())
cluster <- args[3]
suffix <- args[1] %>% str_remove(".*allele_id_")

#===========================================================
#? Match or 0 at sequence error loci
#===========================================================

#--------------------------------------
#* Insertion
#--------------------------------------

tmp_inserr <- df_que %>%
    select(which(df_control$score == 100)) %>%
    lapply(function(x) x %>% table %>% which.max %>% names) %>%
    unlist %>%
    str_detect("M")

tmp_inserr[tmp_inserr == TRUE] <- 2
tmp_inserr[tmp_inserr == FALSE] <- 1

df_control$score[df_control$score == 100] <- tmp_inserr

#--------------------------------------
#* Input sequence error (M/0)
#--------------------------------------

df_que[pull(df_control) == 2] <- "M"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! MIDS count
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

max_count <- function(x){
    x %>% table() %>% which.max() %>% names()
}
df_que <- apply(df_que, 2, max_count)

#===========================================================
#? Summarize mutation position
#===========================================================

df_mut <- tibble(loc = which(df_que != "M"),
    mut = df_que[df_que != "M"] %>% as.character()) %>%
    mutate(insnum = as.character(mut))

df_mut$insnum[str_detect(df_mut$insnum, pattern = "S|D")] <- 0
df_mut$mut[!str_detect(df_mut$mut, pattern = "S|D")] <- "I"

write_delim(df_mut, sprintf(".DAJIN_temp/consensus/temp/mutation_%s", suffix),
    delim=" ", col_names = FALSE)