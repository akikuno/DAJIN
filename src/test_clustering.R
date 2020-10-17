################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
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

# barcode <- "barcode19"
# control <- "barcode43"
# allele <- "flox_deletion"

# if(allele == "abnormal") control_allele <- "wt"
# if(allele != "abnormal") control_allele <- allele
# file_que_mids <- sprintf(".DAJIN_temp/clustering/temp/query_score_%s_%s", barcode, allele)
# file_que_label <- sprintf(".DAJIN_temp/clustering/temp/query_labels_%s_%s", barcode, allele)
# file_control_score <- sprintf(".DAJIN_temp/clustering/temp/df_control_freq_%s.RDS", control_allele)
# threads <- 14L
# plan(multisession, workers = threads)

# ===========================================================
#? Auguments
#===========================================================

args <- commandArgs(trailingOnly = TRUE)
file_que_mids <- args[1]
file_que_label <- args[2]
file_control_score <- args[3]
threads <- as.integer(args[4])
plan(multisession, workers = threads)

#===========================================================
#? Inputs
#===========================================================

df_que_mids <- read_csv(file_que_mids,
    col_names = FALSE,
    col_types = cols())
colnames(df_que_mids) <- 1:ncol(df_que_mids)

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
            x %>% count(MIDS) %>% mutate(Freq = n / sum(n) * 100) %>% select(-n),
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
        rename(score = Freq) %>%
        mutate(score = if_else(MIDS == "M", 0, score)) %>%
        mutate(score = replace_na(score, 0))
    } else {
        full_join(x, y, by = "MIDS", suffix = c("_x", "_y")) %>%
        mutate(score = Freq_x - Freq_y) %>%
        select(-contains("Freq")) %>%
        mutate(score = if_else(MIDS == "M", 0, score)) %>%
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

df_score <- df_score[, colSums(df_score) != 0]

################################################################################
#! PCA
################################################################################

prcomp_result <- prcomp(df_score, scale = FALSE)

num_components <- 1:10
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

#===========================================================
#? Clustering with optimized cluster sizes
#===========================================================

cl <- h$HDBSCAN(min_samples = 1L,
    min_cluster_size = as.integer(min_cluster_sizes[int_cluster_nums_opt]),
    memory = joblib$Memory(cachedir = ".DAJIN_temp/clustering/temp", verbose = 0))

int_hdbscan_clusters <- cl$fit_predict(input_hdbscan) + 1

################################################################################
#! Extract mutation frequency scores in each cluster
################################################################################

df_cluster <- tibble(loc = integer(), cluster = integer(), score = double())

tmp_df_score <- df_score %>% colSums / nrow(df_score)
for (i in unique(int_hdbscan_clusters)) {
    tmp_score <- df_score[int_hdbscan_clusters == i, ] %>%
        colSums / sum(int_hdbscan_clusters == i)

    tmp_score <- tmp_df_score - tmp_score

    tmp_df <- tibble(
        loc = seq_along(colnames(df_score)),
        cluster = i,
        score = tmp_score
    )
    df_cluster <- df_cluster %>% bind_rows(tmp_df)
}
rm(tmp_df, tmp_df_score, tmp_score)

################################################################################
#! Cosine similarity to merge similar clusters
# if two sequences are simillar, merge them
################################################################################

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

pattern_ <- merged_clusters %>%
    unique() %>%
    sort()
query_ <- merged_clusters %>%
    unique() %>%
    order() %>%
    sort()

merged_clusters <- lapply(pattern_,
    function(x) which(x == merged_clusters)) %>%
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

seq_consensus <- mclapply(merged_clusters %>% unique %>% sort,
        function(x) {
            df_que_mids[merged_clusters == x, ] %>%
            lapply(function(x) x %>% table %>% which.max %>% names) %>%
            unlist %>%
            str_c(collapse = "")
        },
        mc.cores = as.integer(threads))

if (length(query_) > 1) {
    df_consensus <- NULL
    cl_combn <- combn(query_, 2)

    for (i in seq(ncol(cl_combn))) {
            df_ <- tibble(
                one = cl_combn[1, i],
                two = cl_combn[2, i],
                score = identical(
                    seq_consensus[[cl_combn[1, i]]],
                    seq_consensus[[cl_combn[2, i]]]
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

pattern_ <- merged_clusters %>%
    unique() %>%
    sort()
query_ <- merged_clusters %>%
    unique() %>%
    order() %>%
    sort()

merged_clusters <- lapply(pattern_,
    function(x) which(x == merged_clusters)) %>%
    set_names(query_) %>%
    map2_df(., query_, function(x, y) {
        tibble(pattern = x,
        query = rep(y, length(x))
        )}) %>%
    arrange(pattern) %>%
    pull(query)

################################################################################
#! Merge clusters
################################################################################

tmp_tibble <-
    tibble(cl = merged_clusters, strand = df_que_label$strand)

#===========================================================
#? Define "dual-end" or "single-end" adaptor
#===========================================================

logic_dual <-
    tmp_tibble %>%
    count(strand) %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq < 0.90 & freq > 0.10) %>%
    nrow %>%
    `>`(0)

#===========================================================
#? Merge clusters with strand specific mutation into a major cluster
#===========================================================

if (logic_dual) {
    tmp_biased_cl <-
        tmp_tibble %>%
        group_by(cl) %>%
        count(strand) %>%
        mutate(freq = n / sum(n)) %>%
        mutate(bias = if_else(freq > 0.90, TRUE, FALSE)) %>%
        filter(bias == TRUE) %>%
        pull(cl) %>%
        unique

    tmp_cl_max <-
        tmp_tibble %>%
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
#? Merge clusters with fuzzy mutations
#===========================================================

possible_true_mut <-
    inner_join(df_que_score, df_control_score, by = "loc") %>%
    unnest(que_freq, control_freq) %>%
    filter(MIDS != "M") %>%
    mutate(score = Freq - Freq1) %>%
    filter(score > 5) %>%
    select(loc, MIDS)

retain_seq_consensus <-
    seq_consensus[merged_clusters %>% unique %>% sort]

retain_seq_consensus <-
    map(retain_seq_consensus, function(x) {
        map_chr(possible_true_mut$loc, function(y) str_sub(x, start = y, end = y))
    })

query_ <- merged_clusters %>% unique %>% sort
if (length(query_) > 1) {
    df_consensus <- NULL
    cl_combn <- combn(query_, 2)

    for (i in seq(ncol(cl_combn))) {
            df_ <- tibble(
                one = cl_combn[1, i],
                two = cl_combn[2, i],
                score = identical(
                    retain_seq_consensus[[cl_combn[1, i]]],
                    retain_seq_consensus[[cl_combn[2, i]]]
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

pattern_ <- merged_clusters %>%
    unique() %>%
    sort()
query_ <- merged_clusters %>%
    unique() %>%
    order() %>%
    sort()

merged_clusters <- lapply(pattern_,
    function(x) which(x == merged_clusters)) %>%
    set_names(query_) %>%
    map2_df(., query_, function(x, y) {
        tibble(pattern = x,
        query = rep(y, length(x))
        )}) %>%
    arrange(pattern) %>%
    pull(query)

################################################################################
#! Format df_readid_cluster
################################################################################

df_readid_cluster <-
    tibble(read_id = df_que_label$id,
    cluster = merged_clusters)

################################################################################
#! Generate df_mutation_score
################################################################################

df_mutation_score <-
    future_map_dfr(list_mids_score,
    function(x) {
        x %>%
        summarize(sum = sum(score)) %>%
        mutate(sum = replace_na(sum, 0))
        }
    ) %>%
    mutate(num = row_number())

################################################################################
#! Output results
################################################################################

write_tsv(df_readid_cluster,
    sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
    col_names = F
)

write_tsv(df_mutation_score,
    sprintf(".DAJIN_temp/clustering/temp/control_score_%s", output_suffix),
    col_names = F
)