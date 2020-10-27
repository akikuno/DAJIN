################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
options(dplyr.summarise.inform = FALSE)
options(future.globals.maxSize = Inf)
options(warn = -1)

if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, furrr, vroom)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# barcode <- "barcode47"
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

df_que_mids <- vroom(file_que_mids,
    col_names = FALSE,
    col_types = cols(),
    num_threads = threads)
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

int_hdbscan_clusters <-
    read_csv(
        sprintf(".DAJIN_temp/clustering/temp/int_hdbscan_clusters_%s", output_suffix),
        col_names = c("cl"),
        col_types = cols()) %>%
        pull(cl)

df_score <- readRDS(sprintf(".DAJIN_temp/clustering/temp/df_score_%s.RDS", output_suffix))

################################################################################
#! Merge clusters
################################################################################

merged_clusters <- int_hdbscan_clusters

#===========================================================
#? Merge clusters with the same mutations
#===========================================================

sequence_error <-
    df_control_score %>%
    unnest(control_freq) %>%
    filter(MIDS != "M" & freq > 10) %>%
    pull(loc) %>%
    unique()

df_score[, sequence_error] <- 10^-100

add_rnorm <-
    rnorm(ncol(df_que_mids) * 100, mean = 0, sd = 5) %>%
    as_tibble() %>%
    rename(score = value)

hotelling <-
    map_dfr(unique(merged_clusters), ~
        df_score[merged_clusters == .x, ] %>%
        colSums() %>%
        scale() %>%
        as_tibble %>%
        rename(score = colnames(.)) %>%
        bind_rows(add_rnorm) %>%
        mutate(loc = row_number()) %>%
        summarize(score = score,
            mean = mean(score),
            var = mean((score - mean(score))^2)) %>%
        summarize(anomaly_score = (score - mean)^2 / var) %>%
        mutate(loc = row_number(), threshold = qchisq(0.95, 1)) %>%
        filter(anomaly_score > threshold & loc <= ncol(df_que_mids)) %>%
        select(loc) %>%
        mutate(cl = .x)
    )

hotelling_common <-
    hotelling %>%
    add_count(loc, name = "count_loc") %>%
    mutate(common = if_else(
        count_loc == length(unique(merged_clusters)),
        TRUE,
        FALSE)) %>%
    select(loc, common) %>%
    distinct()

possible_true_mut <- hotelling_common$loc %>% sort

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
        }) %>%
    set_names(possible_true_mut)

    # common
    tmp_ <-
        hotelling_common %>%
        filter(common == TRUE) %>%
        pull(loc)

    consensus_common <-
        map(unique(merged_clusters), function(cl) {
            map(tmp_, function(loc) {
                df_que_mids[merged_clusters == cl, loc] %>%
                rename(MIDS = colnames(.)) %>%
                count(MIDS) %>%
                mutate(freq = n / sum(n) * 100) %>%
                slice_max(freq, n = 1) %>%
                slice_sample(MIDS, n = 1) %>%
                mutate(MIDS = if_else(freq > 75, MIDS, max_mids[names(max_mids) == loc])) %>%
                pull(MIDS)
            }) %>%
            unlist %>%
            str_c(collapse = "")
        }) %>%
        set_names(unique(merged_clusters))


    # uncommon
    tmp_ <-
        hotelling_common %>%
        filter(common == FALSE) %>%
        pull(loc)

    consensus_uncommon <-
        map(unique(merged_clusters), function(cl) {
            map(tmp_, function(y) {
                control <-
                    df_control_score %>%
                    filter(loc == y) %>%
                    unnest(control_freq) %>%
                    mutate(MIDS = if_else(mut == 1, "M", NULL)) %>%
                    mutate(freq = if_else(mut == 1, 100, NULL))

                df_que_mids[merged_clusters == cl, y] %>%
                rename(MIDS = colnames(.)) %>%
                count(MIDS) %>%
                mutate(freq = n / sum(n) * 100) %>%
                full_join(., control, by = "MIDS", suffix = c("_x", "_y")) %>%
                slice_max(freq_x, n = 1) %>%
                filter(MIDS != "M" & freq_y < 5) %>%
                slice_sample(MIDS, n = 1) %>%
                mutate(MIDS = if_else(freq_x > 90, MIDS, NULL)) %>%
                pull(MIDS)
            }) %>%
            unlist %>%
            str_c(collapse = "")
        }) %>%
        set_names(unique(merged_clusters))

    retain_seq_consensus <-
        map2(consensus_common, consensus_uncommon, function(x, y) {
            str_c(x, y)
        })

    cl_combn <- combn(unique(merged_clusters), 2)

    df_consensus <-
        seq(ncol(cl_combn)) %>%
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
        slice_sample(n = 1) %>%
        pull(cl)

    merged_clusters <-
        merged_clusters %>%
        if_else(. %in% tmp_biased_cl, tmp_cl_max, .)
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
