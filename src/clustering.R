# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan, vroom, umap)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# suffix <- "barcode14_target_HOGE"
# que <- paste0(".DAJIN_temp/clustering//query_score_", sprintf("%s",suffix))
# label <- paste0(".DAJIN_temp/clustering//query_labels_", sprintf("%s", suffix))
# df_que <- vroom(que, col_names = F, col_types = cols())
# df_label <- vroom(label, col_names = c("id", "label"), col_types = cols())
# output_suffix <- label %>% str_remove(".*labels_")

args <- commandArgs(trailingOnly = TRUE)
df_que <- vroom(args[1], col_names = F, col_types = cols())
df_label <- vroom(args[2], col_names = c("id", "label"), col_types = cols())
output_suffix <- args[2] %>% str_remove(".*labels_")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Format input
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df_que <- t(df_que)
df_que[is.na(df_que)] <- 0

df_que_rm0 <- df_que[, colSums(df_que) != 0]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# print("Dimension reduction by PCA...")

input_pca <- df_que_rm0
# //////////////////////////////////////////////////////////
pca_res <- prcomp(input_pca, scale. = F)

components <- seq(1, 10)
output_pca <- pca_res$x[, components]
output_pca_importance <- summary(pca_res)$importance[2, components]
for (i in components) output_pca[, i] <- output_pca[, i] * output_pca_importance[i]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# print("Clustering by HDBSCAN...")
input_hdbscan <- output_pca
# //////////////////////////////////////////////////////////
if (nrow(input_hdbscan) < 250) {
    cl_sizes <- seq(10, nrow(input_hdbscan), length = 10)
} else {
    cl_sizes <- seq(25, 250, length = 10)
}

cl_nums <- c()
i <- 1
for (i in seq_along(cl_sizes)) {
    cl <- hdbscan(input_hdbscan, minPts = cl_sizes[i])
    cl_nums[i] <- cl$cluster %>%
        table() %>%
        length()
}
# cl_nums %>% table()
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

cl <- hdbscan(input_hdbscan, minPts = cl_sizes[cl_num_opt])
output_hdbscan <- cl$cluster + 1

# output_hdbscan %>% table()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pca_plot <- as_tibble(output_pca) %>%
    bind_cols(cluster = factor(output_hdbscan)) %>%
    select(PC1, PC2, cluster)

g <- ggplot(
    data = pca_plot,
    aes(
        x = PC1, y = PC2,
        color = cluster
    )
) +
    geom_point(size = 3) +
    theme_bw(base_size = 20) # +

ggsave(
    plot = g,
    filename = sprintf(".DAJIN_temp/clustering/pca_%s.png", output_suffix),
    width = 10, height = 8
)

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Extract feature nucleotide position
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_que <- df_que_rm0
input_cl <- output_hdbscan

zero_to_one <- function(x) (x - min(x)) / (max(x) - min(x))

df_cluster <- tibble(loc = integer(), cluster = integer(), score = double())
for (i in unique(input_cl)) {
    tmp_score <- input_que[input_cl == i, ] %>%
        colSums() # %>%
    # zero_to_one()

    tmp_df <- tibble(
        loc = seq_len(ncol(input_que)),
        cluster = i,
        score = tmp_score
    )
    df_cluster <- df_cluster %>% bind_rows(tmp_df)
}

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Cosine similarity to merge similar clusters
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# print("Merge similar clusters...")
input_cossim <- df_cluster
output_cl <- output_hdbscan
# //////////////////////////////////////////////////////////
cosine_sim <- function(a, b) crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))

cluster <- input_cossim %>%
    group_by(cluster) %>%
    count() %>%
    select(cluster)

if (nrow(cluster) > 1) {
    cl_combn <- combn(cluster$cluster, 2)
    df_cossim <- NULL
    i <- 1
    for (i in seq_along(1:ncol(cl_combn))) {
        df_1 <- input_cossim %>%
            filter(cluster == cl_combn[1, i]) %>%
            select(score)
        #
        df_2 <- input_cossim %>%
            filter(cluster == cl_combn[2, i]) %>%
            select(score)

        df_ <- tibble(
            one = cl_combn[1, i],
            two = cl_combn[2, i],
            score = cosine_sim(df_1$score, df_2$score)
        )
        df_cossim <- bind_rows(df_cossim, df_)
    }
    df_cossim_filtered <- df_cossim %>% filter(score > 0.90)
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

for (i in seq_along(pattern_)) output_cl[output_cl == pattern_[i]] <- query_[i]

result <- tibble(read_id = df_label$id, output_cl)
write_tsv(result,
    sprintf(".DAJIN_temp/clustering/hdbscan_%s", output_suffix),
    col_names = F
)
