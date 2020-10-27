################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
options(dplyr.summarise.inform = FALSE)
options(future.globals.maxSize = Inf)
options(warn = -1)

if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
if (!requireNamespace("reticulate", quietly = T)) install.packages("reticulate")
pacman::p_load(tidyverse, parallel, furrr)

DAJIN_Python <- reticulate:::conda_list()$python %>%
    str_subset("DAJIN/bin/python")
Sys.setenv(RETICULATE_PYTHON = DAJIN_Python)
reticulate::use_condaenv("DAJIN")

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# barcode <- "barcode23"
# allele <- "wt"

# if(allele == "abnormal") control_allele <- "wt"
# if(allele != "abnormal") control_allele <- allele
# file_que_mids <- sprintf(".DAJIN_temp/clustering/temp/query_score_%s_%s", barcode, allele)
# file_que_label <- sprintf(".DAJIN_temp/clustering/temp/query_labels_%s_%s", barcode, allele)
# file_control_score <- sprintf(".DAJIN_temp/clustering/temp/df_control_freq_%s.RDS", control_allele)
# threads <- 12L
# plan(multiprocess, workers = threads)

# ===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que_mids <- args[1]
file_que_label <- args[2]
file_control_score <- args[3]
threads <- as.integer(args[4])
plan(multiprocess, workers = threads)

#===========================================================
#? Inputs
#===========================================================

df_que_mids <- read_csv(file_que_mids,
    col_names = FALSE,
    col_types = cols())
colnames(df_que_mids) <- seq_len(ncol(df_que_mids))

df_que_label <- read_csv(file_que_label,
    col_names = c("id", "strand", "barcode"),
    col_types = cols())

df_control_score <- readRDS(file_control_score)

#===========================================================
#? Outputs
#===========================================================

output_suffix <-
    str_remove(file_que_label, ".*labels_")

################################################################################
#! MIDS scoring
################################################################################

df_que_score <-
    df_que_mids %>%
    pivot_longer(col = everything(), names_to = "loc", values_to = "MIDS") %>%
    group_by(loc) %>%
    nest(nest = c(MIDS)) %>%
    mutate(que_freq = mclapply(nest,
        function(x)
            x %>% count(MIDS) %>% mutate(freq = n / sum(n) * 100),
        mc.cores = threads)) %>%
    mutate(loc = as.double(loc)) %>%
    select(loc, que_freq)

################################################################################
#! MIDS subtraction
################################################################################

tmp <-
    inner_join(df_que_score, df_control_score, by = "loc")

list_mids_score <-
    future_map2(tmp$que_freq, tmp$control_freq, function(x, y) {
    if (y == 1) {
        x %>%
        rename(score = freq) %>%
        mutate(score = replace_na(score, 0))
    } else {
        full_join(x, y, by = "MIDS", suffix = c("_x", "_y")) %>%
        mutate(score = freq_x - freq_y) %>%
        select(-contains("freq")) %>%
        mutate(score = replace_na(score, 0))
    }
})
rm(tmp)

################################################################################
#! Score each reads
################################################################################

df_score <-
    future_map2_dfc(df_que_mids, list_mids_score,
    function(x, y) {
        tmp1 <- x %>% as_tibble() %>% set_names("MIDS")
        tmp2 <- y
        left_join(tmp1, tmp2, by = "MIDS") %>%
            mutate(score = replace_na(score, 0)) %>%
            pull(score)
    })

df_score[, colSums(df_score) == 0] <- 10^-100

################################################################################
#! PCA
################################################################################

prcomp_result <- prcomp(df_score, scale = FALSE)

num_components <- 1:10

prcomp_loading <-
    sweep(prcomp_result$rotation, 2, prcomp_result$sdev, FUN = "*")[, 1:10] %>%
    as.data.frame()

df_coord <- prcomp_result$x[, num_components] %>% as_tibble
num_prop_variance <- summary(prcomp_result)$importance[2, num_components]

output_pca <- map2_dfc(df_coord, num_prop_variance, ~ .x * .y)

rm(prcomp_result)

################################################################################
#! Clustering
################################################################################

input_hdbscan <- output_pca

joblib <- reticulate::import("joblib")
h <- reticulate::import("hdbscan")

min_cluster_sizes <-
    seq(nrow(input_hdbscan)/20, nrow(input_hdbscan)/2, length = 10) %>%
    as.integer %>%
    `+`(2) %>%
    unique

hd <- function(x) {
    cl <- h$HDBSCAN(min_samples = 1L, min_cluster_size = as.integer(x),
        memory = joblib$Memory(cachedir = ".DAJIN_temp/clustering/temp", verbose = 0))
    cl$fit_predict(input_hdbscan) %>% table %>% length
}

#===========================================================
#? Clustering with multile cluster sizes
#? to find the most frequent cluster numbers
#===========================================================

int_cluster_nums <-
    mclapply(min_cluster_sizes, hd,
    mc.cores = as.integer(threads)) %>%
    unlist %>%
    .[. != 1]

#===========================================================
#? Extract cluster size with the smallest cluster size
#? and the most frequent cluster numbers
#===========================================================

int_cluster_nums_opt <-
    int_cluster_nums %>%
    as_tibble %>%
    mutate(id = row_number()) %>%
    add_count(value, name = "count") %>%
    slice_max(count) %>%
    slice_min(id) %>%
    pull(id)

if(length(int_cluster_nums_opt) == 0)
    int_cluster_nums_opt <- which.max(min_cluster_sizes)

#===========================================================
#? Clustering with optimized cluster sizes
#===========================================================

cl <- h$HDBSCAN(min_samples = 1L,
    min_cluster_size = as.integer(min_cluster_sizes[int_cluster_nums_opt]),
    memory = joblib$Memory(cachedir = ".DAJIN_temp/clustering/temp", verbose = 0))

int_hdbscan_clusters <- cl$fit_predict(input_hdbscan) + 1

################################################################################
#! Output results
################################################################################

write_csv(tibble(cl = int_hdbscan_clusters),
    sprintf(".DAJIN_temp/clustering/temp/int_hdbscan_clusters_%s", output_suffix),
    col_names = F
)

saveRDS(df_score,
    sprintf(".DAJIN_temp/clustering/temp/df_score_%s.RDS", output_suffix)
)
