################################################################################
#! Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan, parallel)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# file_que <- ".DAJIN_temp/clustering/temp/query_score_barcode01_target"
# file_label <- ".DAJIN_temp/clustering/temp/query_labels_barcode01_target"
# file_control <- ".DAJIN_temp/clustering/temp/control_score_target"
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

df_que <- read_csv(file_que, col_names = FALSE, col_types = cols())
df_label <- read_csv(file_label, col_names = c("id", "label"), col_types = cols())
df_control <- read_csv(file_control, col_names = c("score"), col_types = cols())

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
        function(x) pull(df_que[, x]) %>% str_count(pattern = mids[row]) %>% sum(),
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

tmp <- which(df_control$score == 100)

if (length(tmp) != 0){
    df_control$score[seq(tmp[1], tail(tmp, 1))] <- 100

    tmp_inserr <- df_que %>%
        select(which(df_control$score == 100)) %>%
        mclapply(
            function(x) x %>% table %>% which.max %>% names,
            mc.cores = as.integer(threads)) %>%
        unlist %>%
        str_detect("M")

    tmp_inserr[tmp_inserr == TRUE] <- 2
    tmp_inserr[tmp_inserr == FALSE] <- 1

    df_control$score[df_control$score == 100] <- tmp_inserr
}
rm(tmp, tmp_inserr)
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

if (nrow(output_pca) < 500) {
    cl_sizes <- seq(10, nrow(output_pca), length = 10) %>% as.integer
} else {
    cl_sizes <- seq(25, 500, length = 10) %>% as.integer
}

cl_nums <- mclapply(cl_sizes,
    function(x) hdbscan(output_pca, minPts = x)$cluster %>% table %>% length,
    mc.cores = as.integer(1)) %>%
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

hdbscan_cl <- hdbscan(output_pca, minPts = cl_sizes[cl_num_opt])$cluster + 1

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

# ggplot(df_cluster, aes(x=loc, y=score)) +
# geom_point() +
# facet_wrap(~cluster, nrow=1)

################################################################################
#! Cosine similarity to merge similar clusters
# if two sequences are simillar, merge them
################################################################################

calc_cosine_sim <- function(a, b) crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))

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

for (i in seq_along(pattern_)) {
    cossim_merged_cl[cossim_merged_cl == pattern_[i]] <- query_[i]
}

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
                score = identical(tmp_seq[[cl_combn[1, i]]], tmp_seq[[cl_combn[2, i]]])
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

for (i in seq_along(pattern_)) {
    cossim_merged_cl[cossim_merged_cl == pattern_[i]] <- query_[i]
}

################################################################################
#! Output results
################################################################################

result <- tibble(read_id = df_label$id, cossim_merged_cl)

write_tsv(result,
    sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
    col_names = F
)

barcode <- output_suffix %>% str_remove("_.*$")

write_tsv(df_control,
    sprintf(".DAJIN_temp/clustering/temp/control_score_target_%s", barcode),
    col_names = F
)