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

# barcode <- "barcode04"
# allele <- "target"

# if(allele == "abnormal") control_allele <- "wt"
# if(allele != "abnormal") control_allele <- allele
# file_que_mids <- sprintf(".DAJIN_temp/clustering/temp/query_score_%s_%s", barcode, allele)
# file_que_label <- sprintf(".DAJIN_temp/clustering/temp/query_labels_%s_%s", barcode, allele)
# file_control_score <- sprintf(".DAJIN_temp/clustering/temp/df_control_freq_%s.RDS", control_allele)
# threads <- 2L
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

# tmp_df_score <- df_score %>% colSums / nrow(df_score)
for (i in unique(int_hdbscan_clusters)) {
    tmp_score <-
        df_score[int_hdbscan_clusters == i, ] %>%
        colSums / sum(int_hdbscan_clusters == i)

    # tmp_score <- tmp_df_score - tmp_score

    tmp_df <- tibble(
        loc = seq_along(colnames(df_score)),
        cluster = i,
        score = tmp_score
    )
    df_cluster <- df_cluster %>% bind_rows(tmp_df)
}
rm(tmp_df, tmp_df_score, tmp_score)

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
#? Merge clusters with solid/fuzzy mutations
#? if there is a possible true mutation locus
#===========================================================

tmp_que_score <- df_que_score %>% unnest(que_freq)
tmp_control_score <- df_control_score %>% unnest(control_freq)

tmp_possible_true_mut <-
    full_join(tmp_que_score, tmp_control_score, by = c("loc","MIDS"), suffix = c("_x", "_y")) %>%
    mutate(Freq_x = replace_na(Freq_x, 0)) %>%
    mutate(Freq_y = replace_na(Freq_y, 0)) %>%
    filter(MIDS != "M") %>%
    mutate(score = Freq_x - Freq_y) %>%
    filter((Freq_y < 5 & score > 5)) %>%
    filter(!(Freq_y == 0 & Freq_x < 20)) %>%
    pull(loc)

if(length(tmp_possible_true_mut) > 0) {
    possible_true_mut <-
        df_que_mids[, tmp_possible_true_mut] %>%
        cbind(cl = merged_clusters) %>%
        pivot_longer(-cl, names_to = "loc", values_to = "MIDS") %>%
        group_by(cl, loc) %>%
        count(MIDS) %>%
        mutate(Freq = n / sum(n) * 100) %>%
        filter(MIDS != "M") %>%
        filter(Freq > 75) %>%
        ungroup() %>%
        select(loc, MIDS) %>%
        unique
} else {
    possible_true_mut <- tibble()
}

if(nrow(possible_true_mut) > 0) {
    solid_clusters <-
        map_lgl(merged_clusters %>% unique, function(x){
            df_que_mids[merged_clusters == x, possible_true_mut$loc] %>%
            pivot_longer(col = everything(), names_to = "loc", values_to = "MIDS") %>%
            group_by(loc) %>%
            count(MIDS) %>%
            mutate(Freq = n / sum(n) * 100) %>%
            mutate(tf = if_else(Freq > 75, TRUE, FALSE)) %>%
            pull(tf) %>%
            all()
        })

    solid_cluster_numbers <- unique(merged_clusters)[solid_clusters]
    fuzzy_cluster_numbers <- unique(merged_clusters)[!solid_clusters]

    #--------------------------------------
    #* Merge clusters with solid mutations
    #* if their mutations are the same
    #--------------------------------------

    retain_seq_consensus <-
        mclapply(solid_cluster_numbers,
            function(x) {
                df_que_mids[merged_clusters == x, ] %>%
                lapply(function(x) x %>% table %>% which.max %>% names) %>%
                unlist %>%
                str_c(collapse = "")
            },
            mc.cores = as.integer(threads)) %>%
        map(function(x) {
            map_chr(possible_true_mut$loc, function(y) str_sub(x, start = y, end = y))
        }) %>%
        set_names(solid_cluster_numbers)

    query_ <- solid_cluster_numbers %>% unique
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

    #--------------------------------------
    #* Merge clusters with fuzzy mutations into a major cluster
    #--------------------------------------

    tmp_cl_max <-
        merged_clusters %>%
        as_tibble() %>%
        count(value) %>%
        filter(!(value %in% fuzzy_cluster_numbers)) %>%
        slice_max(n, n = 1) %>%
        pull(value)

    merged_clusters <-
        ifelse(merged_clusters %in% fuzzy_cluster_numbers, tmp_cl_max, merged_clusters)
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

write_tsv(possible_true_mut,
    sprintf(".DAJIN_temp/clustering/temp/possible_true_mut_%s", output_suffix),
    col_names = F
)