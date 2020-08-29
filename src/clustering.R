################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
if (!requireNamespace("reticulate", quietly = T)) install.packages("reticulate")
pacman::p_load(tidyverse, parallel)

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

# file_que <- ".DAJIN_temp/clustering/temp/query_score_barcode02_inversion"
# file_label <- ".DAJIN_temp/clustering/temp/query_labels_barcode02_inversion"
# file_control <- ".DAJIN_temp/clustering/temp/control_score_inversion"
# threads <- 14L

#===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que <- args[1]
file_label <- args[2]
file_control <- args[3]
threads <- as.integer(args[4])

#===========================================================
#? Inputs
#===========================================================

df_que <- read_csv(file_que,
    col_names = FALSE,
    col_types = cols())

df_label <- read_csv(file_label,
    col_names = c("id", "strand", "barcode"),
    col_types = cols())

df_control <- read_csv(file_control,
    col_names = c("score"),
    col_types = cols())

#===========================================================
#? Outputs
#===========================================================

output_suffix <- file_label %>% str_remove(".*labels_")

################################################################################
#! Pre-processing
################################################################################

#===========================================================
#? MIDS scoring
#===========================================================

mids <- c("M", "I", "D", "S", "=")
mids_score <- matrix(NA, length(mids), ncol(df_que))
rownames(mids_score) <- mids

for (row in seq_along(mids)) {
    mids_score[row, ] <-
    mclapply(seq_along(colnames(df_que)) %>% as.list,
        function(x) {
            pull(df_que[, x]) %>% str_count(pattern = mids[row]) %>% sum()
            },
        mc.cores = as.integer(threads)) %>%
        unlist
}

mids_score["M", ] <- 0
mids_score["D", ] <- mids_score["D", ] * -1

df_score <- mclapply(seq_along(colnames(df_que)),
    function(x) pull(df_que[, x]) %>% mids_score[, x][.],
    mc.cores = as.integer(threads)) %>%
    as.data.frame

rm(mids_score)

#===========================================================
#? Match or 0 at sequence error loci
#===========================================================

#--------------------------------------
#* Insertion
#--------------------------------------

tmp_score100 <- which(df_control$score == 100)

if (length(tmp_score100) != 0) {
    df_control$score[seq(tmp_score100[1], tail(tmp_score100, 1))] <- 100

    tmp_ins_error <- df_que %>%
    select(which(df_control$score == 100)) %>%
    mclapply(
        function(x) x %>% table %>% which.max %>% names,
        mc.cores = as.integer(threads)) %>%
    unlist %>%
    str_replace("M", "2") %>%
    str_replace("[^M2]", "1") %>%
    as.numeric

    df_control$score[df_control$score == 100] <- tmp_ins_error
    rm(tmp_ins_error)
}
rm(tmp_score100)

#--------------------------------------
#* Input sequence error (M/0)
#--------------------------------------

df_que <- df_que[, !pull(df_control) == 2]

df_score[, pull(df_control) == 2] <- 0
df_score <- df_score[, colSums(df_score) != 0]

################################################################################
#! PCA
################################################################################

pca_res <- prcomp(df_score, scale = FALSE)

components <- 1:10
output_pca <- pca_res$x[, components]
output_pca_importance <- summary(pca_res)$importance[2, components]
for (i in components) {
    output_pca[, i] <- output_pca[, i] * output_pca_importance[i]
}
rm(pca_res)

################################################################################
#! Clustering
################################################################################

joblib <- reticulate::import("joblib")
h <- reticulate::import("hdbscan")

min_cluster_sizes <- seq(nrow(output_pca)/20, nrow(output_pca)/2, length = 10) %>%
    unique
min_cluster_sizes <- as.integer(min_cluster_sizes + 2)

hd <- function(x) {
    cl <- h$HDBSCAN(min_samples = 1L, min_cluster_size = as.integer(x),
        memory = joblib$Memory(cachedir = ".DAJIN_temp/clustering/temp", verbose = 0))
    cl$fit_predict(output_pca) %>% table %>% length
}

cl_nums <- mclapply(min_cluster_sizes, hd,
    mc.cores = as.integer(threads)) %>%
    unlist

cl_num_opt <- cl_nums %>% table()

if ((cl_num_opt %>% names != 1) %>% sum > 0) {
    cl_num_opt <- cl_num_opt[names(cl_num_opt) != 1] %>%
        which.max() %>%
        names()
} else {
    cl_num_opt <- cl_num_opt %>%
        which.max() %>%
        names()
}
cl_num_opt <- which(cl_nums == cl_num_opt) %>% max()

cl <- h$HDBSCAN(min_samples = 1L,
    min_cluster_size = as.integer(min_cluster_sizes[cl_num_opt]),
    memory = joblib$Memory(cachedir = ".DAJIN_temp/clustering/temp", verbose = 0))

hdbscan_cl <- cl$fit_predict(output_pca) + 1

################################################################################
#! Extract mutation frequency scores in each cluster
################################################################################

df_cluster <- tibble(loc = integer(), cluster = integer(), score = double())

tmp_df_score <- df_score %>% colSums / nrow(df_score)
for (i in unique(hdbscan_cl)) {
    tmp_score <- df_score[hdbscan_cl == i, ] %>%
        colSums / sum(hdbscan_cl == i)

    tmp_score <- tmp_df_score - tmp_score

    tmp_df <- tibble(
        loc = seq_along(colnames(df_score)),
        cluster = i,
        score = tmp_score
    )
    df_cluster <- df_cluster %>% bind_rows(tmp_df)
}
rm(df_score)

#### TEST <<<<<<<<<<<<
# ggplot(df_cluster, aes(x=loc, y=score)) +
# geom_point() +
# facet_wrap(~cluster, nrow=1)
#### TEST >>>>>>>>>>>>

################################################################################
#! Cosine similarity to merge similar clusters
# if two sequences are simillar, merge them
################################################################################

calc_cosine_sim <- function(a, b) {
    crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))
    }

cossim_merged_cl <- hdbscan_cl

cluster <- df_cluster$cluster %>%
    unique()

if (length(cluster) > 1) {
    df_cossim <- NULL
    cl_combn <- combn(cluster, 2)

    i <- 1
    for (i in seq(ncol(cl_combn))) {
        score_1 <- df_cluster %>%
            filter(cluster == cl_combn[1, i]) %>%
            pull(score)
        #
        score_2 <- df_cluster %>%
            filter(cluster == cl_combn[2, i]) %>%
            pull(score)

        cos_sim <- calc_cosine_sim(score_1, score_2)

        df_ <- tibble(
            one = cl_combn[1, i],
            two = cl_combn[2, i],
            score = cos_sim
        )
        df_cossim <- bind_rows(df_cossim, df_)
    }
    df_cossim_extracted <- df_cossim %>% filter(score > 0.6)

    if (nrow(df_cossim_extracted) != 0) {
        for (i in seq_along(rownames(df_cossim_extracted))) {
            pattern_ <- df_cossim_extracted[i, ]$one
            query_ <- df_cossim_extracted[i, ]$two
            cossim_merged_cl[cossim_merged_cl == pattern_] <- query_
        }
        for (i in seq_along(rownames(df_cossim_extracted))) {
            pattern_ <- df_cossim_extracted[i, ]$two
            query_ <- df_cossim_extracted[i, ]$one
            cossim_merged_cl[cossim_merged_cl == pattern_] <- query_
        }
    }
}

pattern_ <- cossim_merged_cl %>%
    unique() %>%
    sort()
query_ <- cossim_merged_cl %>%
    unique() %>%
    order() %>%
    sort()

cossim_merged_cl <- lapply(pattern_,
    function(x) which(x == cossim_merged_cl)) %>%
    set_names(query_) %>%
    map2_df(., query_, function(x, y) {
        tibble(pattern = x,
        query = rep(y, length(x))
        )}) %>%
    arrange(pattern) %>%
    pull(query)

################################################################################
#! Sequence identity to merge similar clusters
# if two sequences are the same, merge them
################################################################################

if (length(query_) > 1) {
    df_cossim <- NULL
    cl_combn <- combn(query_, 2)

    tmp_seq <- mclapply(cossim_merged_cl %>% unique %>% sort,
            function(x) {
                df_que[cossim_merged_cl == x, ] %>%
                lapply(function(x) x %>% table %>% which.max %>% names) %>%
                unlist %>%
                str_c(collapse = "")
            },
            mc.cores = as.integer(threads))

    df_cossim <- NULL
    for (i in seq(ncol(cl_combn))) {
            df_ <- tibble(
                one = cl_combn[1, i],
                two = cl_combn[2, i],
                score = identical(
                    tmp_seq[[cl_combn[1, i]]],
                    tmp_seq[[cl_combn[2, i]]]
                    )
            )
            df_cossim <- bind_rows(df_cossim, df_)
    }

    df_cossim_extracted <- df_cossim %>% filter(score == TRUE)

    if (nrow(df_cossim_extracted) != 0) {
        for (i in seq_along(rownames(df_cossim_extracted))) {
            pattern_ <- df_cossim_extracted[i, ]$one
            query_ <- df_cossim_extracted[i, ]$two
            cossim_merged_cl[cossim_merged_cl == pattern_] <- query_
        }
    }
}

pattern_ <- cossim_merged_cl %>%
    unique() %>%
    sort()
query_ <- cossim_merged_cl %>%
    unique() %>%
    order() %>%
    sort()

cossim_merged_cl <- lapply(pattern_,
    function(x) which(x == cossim_merged_cl)) %>%
    set_names(query_) %>%
    map2_df(., query_, function(x, y) {
        tibble(pattern = x,
        query = rep(y, length(x))
        )}) %>%
    arrange(pattern) %>%
    pull(query)

################################################################################
#! Remove reads with strand specific mutation
################################################################################

tmp_cl <-
    tibble(cl = cossim_merged_cl, strand = df_label$strand) %>%
    group_by(cl) %>%
    count(strand) %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq < 0.95 & freq > 0.05) %>%
    pull(cl) %>%
    unique

retain_reads <- cossim_merged_cl %in% tmp_cl

################################################################################
#! Output results
################################################################################

result <- tibble(read_id = df_label$id[retain_reads],
    cluster = cossim_merged_cl[retain_reads])

write_tsv(result,
    sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
    col_names = F
)

barcode <- output_suffix %>% str_remove("_.*$")

write_tsv(df_control,
    sprintf("%s_%s", file_control, barcode),
    col_names = F
)