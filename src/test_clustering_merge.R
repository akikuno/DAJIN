################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
options(dplyr.summarise.inform = FALSE)
options(warn = -1)

if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, parallel)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# barcode <- "barcode04"
# allele <- "target"

# if (allele == "abnormal") control_allele <- "wt"
# if (allele != "abnormal") control_allele <- allele
# file_que_mids <- sprintf(".DAJIN_temp/clustering/temp/query_score_%s_%s", barcode, allele)
# file_que_label <- sprintf(".DAJIN_temp/clustering/temp/query_labels_%s_%s", barcode, allele)
# file_control_score <- sprintf(".DAJIN_temp/clustering/temp/df_control_freq_%s.RDS", control_allele)
# threads <- 12L

# ===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que_mids <- args[1]
file_que_label <- args[2]
file_control_score <- args[3]
threads <- as.integer(args[4])

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

calc_cosine_sim <- function(a, b) {
    crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))
    }

merged_clusters <- int_hdbscan_clusters

cluster <-
    df_cluster$cluster %>%
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
    df_cossim_extracted <- df_cossim %>% filter(score > 0.5)

    if (nrow(df_cossim_extracted) != 0) {
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

logic_dual <-
    tmp_tb_strand %>%
    count(strand) %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq < 0.90 & freq > 0.10) %>%
    nrow %>%
    `>`(0)

if (logic_dual) {
    tmp_biased_cl <-
        tmp_tb_strand %>%
        group_by(cl) %>%
        count(strand) %>%
        mutate(freq = n / sum(n)) %>%
        mutate(bias = if_else(freq > 0.90, TRUE, FALSE)) %>%
        filter(bias == TRUE) %>%
        pull(cl) %>%
        unique

    tmp_cl_max <-
        tmp_tb_strand %>%
        group_by(cl) %>%
        count() %>%
        ungroup(cl) %>%
        filter(!(cl %in% tmp_biased_cl)) %>%
        slice_max(n, n = 1) %>%
        pull(cl)

    merged_clusters <-
        ifelse(merged_clusters %in% tmp_biased_cl, tmp_cl_max, merged_clusters)
}

#===========================================================
#? Merge clusters with the same mutations
#===========================================================
#* TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 大型欠損では異常塩基>>正常塩基になり, 異常塩基の同定数が少なくなってしまうため
# サンプル正常サンプルを大量に投下して異常塩基の数を相対的に減らして
# ホテリングの異常スコアを際立たせる.
abs_score <-
    prcomp_loading %>%
    mutate(loc = row_number()) %>%
    pivot_longer(-loc, names_to = "PC", values_to = "score") %>%
    group_by(loc) %>%
    summarize(score = sum(score)) %>%
    select(score)

add_rnorm <-
    rnorm(nrow(prcomp_loading) * 10, mean = 0, sd = 0.1) %>%
    as_tibble() %>%
    rename(score = value)

hotelling_mut <-
    bind_rows(abs_score, add_rnorm) %>%
    mutate(loc = row_number()) %>%
    summarize(score = score,
        mean = mean(score),
        var = mean((score - mean(score))^2)) %>%
    summarize(anomaly_score = (score - mean)^2 / var) %>%
    mutate(loc = row_number(), threshold = qchisq(0.99, 1)) %>%
    filter(anomaly_score > threshold) %>%
    select(loc)

    # mutate(loc = row_number()) %>%
    # ggplot(aes(x = loc, y = anomaly_score)) + geom_point()
#* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# hotelling_mut <-
#     prcomp_loading %>%
#     mutate(loc = row_number()) %>%
#     pivot_longer(-loc, names_to = "PC", values_to = "score") %>%
#     group_by(loc) %>%
#     summarize(score = sum(abs(score))) %>%
#     summarize(score = score,
#         mean = mean(score),
#         var = mean((score - mean(score))^2)) %>%
#     summarize(anomaly_score = (score - mean)^2 / var) %>%
#     mutate(loc = row_number(), threshold = qchisq(0.99, 1)) %>%
#     filter(anomaly_score > threshold) %>%
#     select(loc)

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
        unique()
}

cl_nums <- length(unique(merged_clusters))
shared_true_mut <- as.integer()

if (cl_nums > 1 && length(possible_true_mut) > 0) {
    shared_true_mut <-
        map_dfr(merged_clusters %>% unique, function(x) {
            df_que_mids[merged_clusters == x, possible_true_mut] %>%
            pivot_longer(col = everything(),
                names_to = "loc",
                values_to = "MIDS") %>%
            group_by(loc) %>%
            count(MIDS) %>%
            mutate(freq = n / sum(n) * 100) %>%
            group_by(loc) %>%
            slice_max(n, n = 1) %>%
            mutate(cl = x)
        }) %>%
        filter(freq > 75) %>%
        select(loc, cl, MIDS) %>%
        distinct(loc, MIDS, .keep_all = TRUE) %>%
        arrange(loc) %>%
        summarise(MIDS = MIDS, cl = cl, n = n()) %>%
        filter(n > 1) %>%
        pull(loc) %>%
        as.integer() %>%
        unique
}

if (length(shared_true_mut) > 0) {

    retain_seq_consensus <-
        mclapply(merged_clusters %>% unique,
            function(x) {
                df_que_mids[merged_clusters == x, ] %>%
                select(all_of(shared_true_mut)) %>%
                lapply(function(x) x %>% table %>% which.max %>% names) %>%
                unlist %>%
                str_c(collapse = "")
            },
            mc.cores = as.integer(threads)) %>%
        set_names(merged_clusters %>% unique)

    query_ <- merged_clusters %>% unique
    if (length(query_) > 1) {
        df_consensus <- NULL
        cl_combn <- combn(query_, 2)

        for (i in seq(ncol(cl_combn))) {
                df_ <- tibble(
                    one = cl_combn[1, i],
                    two = cl_combn[2, i],
                    score = identical(
                        retain_seq_consensus[[as.character(cl_combn[1, i])]],
                        retain_seq_consensus[[as.character(cl_combn[2, i])]]
                        )
                )
                df_consensus <- bind_rows(df_consensus, df_)
        }

        df_consensus_extracted <- df_consensus %>% filter(score == TRUE)

        if (nrow(df_consensus_extracted) != 0) {
            for (i in seq_along(rownames(df_consensus_extracted))) {
                pattern_ <- df_consensus_extracted[i, ]$one
                query_ <- df_consensus_extracted[i, ]$two
                merged_clusters[merged_clusters == pattern_] <- query_
            }
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
    order() %>%
    sort()

merged_clusters <-
    lapply(pattern_,
        function(x) which(x == merged_clusters)) %>%
        set_names(query_) %>%
        map2_df(., query_, function(x, y) {
            tibble(pattern = x,
            query = rep(y, length(x))
            )}) %>%
    arrange(pattern) %>%
    pull(query)

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
