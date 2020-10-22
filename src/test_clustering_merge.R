################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
options(dplyr.summarise.inform = FALSE)
options(future.globals.maxSize = Inf)
options(warn = -1)

if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, parallel, furrr)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# barcode <- "barcode02"
# allele <- "wt"

# if (allele == "abnormal") control_allele <- "wt"
# if (allele != "abnormal") control_allele <- allele
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

df_cluster <-
    read_csv(
        sprintf(".DAJIN_temp/clustering/temp/df_cluster_%s", output_suffix),
        col_names = c("loc", "cluster", "score"),
        col_types = cols())

int_hdbscan_clusters <-
    read_csv(
        sprintf(".DAJIN_temp/clustering/temp/int_hdbscan_clusters_%s", output_suffix),
        col_names = c("cl"),
        col_types = cols()) %>%
        pull(cl)

prcomp_loading <-
    read_csv(
        sprintf(".DAJIN_temp/clustering/temp/prcomp_loading_%s", output_suffix),
        col_names = FALSE,
        col_types = cols())

################################################################################
#! Merge clusters
################################################################################

#===========================================================
#? Cosine similarity to merge similar clusters
# if two sequences are simillar, merge them
#===========================================================

merged_clusters <- int_hdbscan_clusters

calc_cosine_sim <- function(x, y) {
    crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
    }

cluster <- unique(merged_clusters)

if (length(cluster) > 1) {
    cl_combn <- combn(cluster, 2)
    df_cossim <-
        seq_along(cluster) %>%
        sort %>%
        map_dfr(function(x) {
            score_1 <- df_cluster %>%
                filter(cluster == cl_combn[1, x]) %>%
                pull(score)
            score_2 <- df_cluster %>%
                filter(cluster == cl_combn[2, x]) %>%
                pull(score)

            cos_sim <- calc_cosine_sim(score_1, score_2)

            tibble(
                one = cl_combn[1, x],
                two = cl_combn[2, x],
                score = cos_sim
            )
        })
    df_cossim_extracted <- df_cossim %>% filter(score > 0.5)

    if (nrow(df_cossim_extracted) > 0) {
        for (i in seq_along(rownames(df_cossim_extracted))) {
            pattern_ <- df_cossim_extracted[i, ]$one
            query_ <- df_cossim_extracted[i, ]$two
            merged_clusters[merged_clusters == pattern_] <- query_
        }
        for (i in seq_along(rownames(df_cossim_extracted))) {
            pattern_ <- df_cossim_extracted[i, ]$two
            query_ <- df_cossim_extracted[i, ]$one
            merged_clusters[merged_clusters == pattern_] <- query_
        }
    }
}

#===========================================================
#? Merge clusters with strand specific mutation into a major cluster
#===========================================================

tmp_tb_strand <-
    tibble(cl = merged_clusters, strand = df_que_label$strand)

tmp_dual_adaptor <-
    tmp_tb_strand %>%
    count(strand) %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq < 0.90 & freq > 0.10)

if (nrow(tmp_dual_adaptor) > 0) {
    tmp_biased_cl <-
        tmp_tb_strand %>%
        count(cl, strand) %>%
        mutate(freq = n / sum(n)) %>%
        mutate(bias = if_else(freq > 0.90, TRUE, FALSE)) %>%
        filter(bias == TRUE) %>%
        pull(cl) %>%
        unique

    tmp_cl_max <-
        tmp_tb_strand %>%
        count(cl) %>%
        filter(!(cl %in% tmp_biased_cl)) %>%
        slice_max(n, n = 1) %>%
        pull(cl) %>%
        as.integer()

    merged_clusters <-
        merged_clusters %>%
        if_else(. %in% tmp_biased_cl, tmp_cl_max, .)
}

#===========================================================
#? Merge clusters with the same mutations
#===========================================================

loading_score <-
    prcomp_loading %>%
    mutate(loc = row_number()) %>%
    pivot_longer(-loc, names_to = "PC", values_to = "score") %>%
    group_by(loc) %>%
    summarize(score = sum(score)) %>%
    select(score)

add_rnorm <-
    rnorm(nrow(loading_score) * 100, mean = 0, sd = 0.1) %>%
    as_tibble() %>%
    rename(score = value)

hotelling_mut <-
    bind_rows(loading_score, add_rnorm) %>%
    mutate(loc = row_number()) %>%
    summarize(score = score,
        mean = mean(score),
        var = mean((score - mean(score))^2)) %>%
    summarize(anomaly_score = (score - mean)^2 / var) %>%
    mutate(loc = row_number(), threshold = qchisq(0.99, 1)) %>%
    filter(anomaly_score > threshold) %>%
    select(loc)

if (nrow(hotelling_mut) > 0) {
    possible_true_mut <-
        inner_join(hotelling_mut, df_control_score, by = "loc") %>%
        unnest(control_freq) %>%
        filter(MIDS != "M") %>%
        group_by(loc) %>%
        slice_max(freq, n = 1) %>%
        filter(freq < 10) %>%
        pull(loc) %>%
        unique()
} else {
    possible_true_mut <- as.integer()
}

if (sum(df_control_score$mut) == 1) {
    possible_true_mut <-
        c(possible_true_mut, which(df_control_score$mut == 1)) %>%
        unique() %>%
        sort()
}

if (length(unique(merged_clusters)) > 1 && length(possible_true_mut) > 0) {

    max_mids <-
        future_map_chr(possible_true_mut, function(x) {
        df_que_mids[, x] %>%
        rename(MIDS = colnames(.)) %>%
        count(MIDS) %>%
        mutate(freq = n / sum(n) * 100) %>%
        slice_max(freq, n = 1) %>%
        slice_sample(MIDS, n = 1) %>%
        pull(MIDS)
    })

    retain_seq_consensus <-
        mclapply(merged_clusters %>% unique, function(cl) {
            map2_chr(possible_true_mut, max_mids, function(x, y) {
                df_que_mids[merged_clusters == cl, x] %>%
                rename(MIDS = colnames(.)) %>%
                count(MIDS) %>%
                mutate(freq = n / sum(n) * 100) %>%
                slice_max(freq, n = 1) %>%
                slice_sample(MIDS, n = 1) %>%
                mutate(MIDS = if_else(freq < 90, y, MIDS)) %>%
                pull(MIDS)
                }) %>%
            str_c(collapse = "")
        }, mc.cores = as.integer(threads)) %>%
        set_names(merged_clusters %>% unique)

    cl_combn <- combn(unique(merged_clusters), 2)

    df_consensus <-
        seq_along(unique(merged_clusters)) %>%
        sort %>%
        map_dfr(function(x) {
            tmp_seq1 <- retain_seq_consensus[[as.character(cl_combn[1, x])]]
            tmp_seq2 <- retain_seq_consensus[[as.character(cl_combn[2, x])]]

            tibble(
                one = cl_combn[1, x],
                two = cl_combn[2, x],
                score = identical(tmp_seq1, tmp_seq2)
            )
        })
    df_consensus_extracted <- df_consensus %>% filter(score == TRUE)

    if (nrow(df_consensus_extracted) != 0) {
        for (i in seq_along(rownames(df_consensus_extracted))) {
            pattern_ <- df_consensus_extracted[i, ]$one
            query_ <- df_consensus_extracted[i, ]$two
            merged_clusters[merged_clusters == pattern_] <- query_
        }
    }
} else {
    merged_clusters <-
        rep(1, length(merged_clusters))
}

################################################################################
#! Output results
################################################################################

pattern_ <- merged_clusters %>%
    unique() %>%
    sort()
query_ <- merged_clusters %>%
    unique() %>%
    seq_along()

merged_clusters <-
    merged_clusters %>%
    map(function(cl) {
        map2(pattern_, query_, function(x, y) {
                if_else(cl %in% x, as.integer(y), NULL)
        })
    }) %>%
    unlist %>%
    .[!is.na(.)]


df_readid_cluster <-
    tibble(read_id = df_que_label$id,
    cluster = merged_clusters)

write_tsv(df_readid_cluster,
    sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
    col_names = F
)

write_tsv(tibble(loc = possible_true_mut),
    sprintf(".DAJIN_temp/clustering/temp/possible_true_mut_%s", output_suffix),
    col_names = F
)
