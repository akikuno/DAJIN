# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cloud.r-project.org/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

args <- commandArgs(trailingOnly = TRUE)
df_que <- read_csv(args[1], col_names = FALSE, col_types = cols())
df_label <- read_csv(args[2], col_names = c("id", "label"), col_types = cols())
control <- read_csv(args[3], col_names = FALSE, col_types = cols())

output_suffix <- args[2] %>% str_remove(".*labels_")

mids_score <- matrix(NA, 4, ncol(df_que))
rownames(mids_score) <- c("M", "I", "D", "S")

for (row in 1:4) {
    mutation <- rownames(mids_score)[row]
    for (col in 1 : ncol(df_que)) {
        mids_score[row, col] <-
            str_count(pull(df_que[, col]), pattern = mutation) %>%
            sum()
    }
}
mids_score[1, ] <- 0
mids_score["D", ] <- mids_score["D", ] * -1

for (col in 1:ncol(df_que)) {
    res <- mids_score[, col]
    df_que[, col] <- pull(df_que[, col]) %>% res[.]
}

df_que[, pull(control) == 2] <- 0
df_que <- df_que[, colSums(df_que) != 0]

rm(mids_score, control)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# input: df_que
# --------------------------------------------------
pca_res <- prcomp(df_que, scale. = F)

components <- 1:10
output_pca <- pca_res$x[, components]
output_pca_importance <- summary(pca_res)$importance[2, components]
for (i in components) {
    output_pca[, i] <- output_pca[, i] * output_pca_importance[i]
}
rm(pca_res)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if (nrow(output_pca) < 250) {
    cl_sizes <- seq(10, nrow(output_pca), length = 10)
} else {
    cl_sizes <- seq(25, 250, length = 10)
}

cl_nums <- c()
i <- 1
for (i in seq_along(cl_sizes)) {
    cl <- hdbscan(output_pca, minPts = cl_sizes[i])
    cl_nums[i] <- cl$cluster %>%
        table() %>%
        length()
}

cl_num_opt <- cl_nums %>% table()

if (length(cl_num_opt[names(cl_num_opt) != 1]) > 0) {
    cl_num_opt <- cl_num_opt[names(cl_num_opt) != 1] %>%
        which.max() %>%
        names()
} else {
    cl_num_opt <- cl_num_opt %>%
        which.max() %>%
        names()
}
cl_num_opt <- which(cl_nums == cl_num_opt) %>% max()

cl <- hdbscan(output_pca, minPts = cl_sizes[cl_num_opt])
hdbscan_cl <- cl$cluster + 1

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract feature nucleotide position
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# input: df_que, hdbscan_cl
# --------------------------------------------------

zero_to_one <- function(x) (x - min(x)) / (max(x) - min(x))

df_cluster <- tibble(loc = integer(), cluster = integer(), score = double())
for (i in unique(hdbscan_cl)) {
    tmp_score <- df_que[hdbscan_cl == i, ] %>%
        colSums() %>%
        zero_to_one()

    tmp_df <- tibble(
        loc = seq_len(ncol(df_que)),
        cluster = i,
        score = tmp_score
    )
    df_cluster <- df_cluster %>% bind_rows(tmp_df)
}
rm(df_que)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cosine similarity to merge similar clusters
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# input: df_cluster, hdbscan_cl
output_cl <- hdbscan_cl
# --------------------------------------------------
cosine_sim <- function(a, b) crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))

cluster <- df_cluster %>%
    group_by(cluster) %>%
    count() %>%
    pull(cluster)

if (length(cluster) > 1) {
    cl_combn <- combn(cluster, 2)
    df_cossim <- NULL
    i <- 1
    for (i in seq(ncol(cl_combn))) {
        df_1 <- df_cluster %>%
            filter(cluster == cl_combn[1, i]) %>%
            pull(score)
        #
        df_2 <- df_cluster %>%
            filter(cluster == cl_combn[2, i]) %>%
            pull(score)

        df_ <- tibble(
            one = cl_combn[1, i],
            two = cl_combn[2, i],
            score = cosine_sim(df_1, df_2)
        )
        df_cossim <- bind_rows(df_cossim, df_)
    }
    df_cossim_filtered <- df_cossim %>% filter(score > 0.95)
    #
    if (nrow(df_cossim_filtered) != 0) {
        for (i in seq_len(nrow(df_cossim_filtered))) {
            pattern_ <- df_cossim_filtered[i, ]$one
            query_ <- df_cossim_filtered[i, ]$two
            output_cl[output_cl == pattern_] <- query_
        }
    }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pattern_ <- output_cl %>%
    unique() %>%
    sort()
query_ <- output_cl %>%
    unique() %>%
    order() %>%
    sort()

for (i in seq_along(pattern_)) {
    output_cl[output_cl == pattern_[i]] <- query_[i]
}

result <- tibble(read_id = df_label$id, output_cl)
write_tsv(result,
    sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
    col_names = F
)