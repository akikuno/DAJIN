# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos='https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# suffix <- "barcode08_normal"
# que <- paste0(".DAJIN_temp/clustering/temp/query_score_", sprintf("%s",suffix))
# label <- paste0(".DAJIN_temp/clustering/temp/query_labels_", sprintf("%s", suffix))
# df_que <- read_csv(que, col_names = F, col_types = cols())
# df_que <- t(df_que)
# df_que[is.na(df_que)] <- 0
# df_label <- read_csv(label, col_names = c("id", "label"), col_types = cols())
# output_suffix <- label %>% str_remove(".*labels_")

# df_score <- read_csv("tmp_score", col_names = F, col_types = cols())
# df_weight <- read_csv("tmp_weight", col_names = F, col_types = cols())
# df_label <- read_csv("tmp_label", col_names = c("id", "label"), col_types = cols())

# dim(df_score)
# dim(df_weight)
# dim(df_label)

# df_que <- read_csv("test_MIDS.csv", col_names = FALSE, col_types = cols())
# df_label <- read_csv("tmp_label", col_names = c("id", "label"), col_types = cols())
# control <- read_csv("test_control", col_names = FALSE, col_types = cols())

args <- commandArgs(trailingOnly = TRUE)
df_que <- read_csv(args[1], col_names = FALSE, col_types = cols())
df_label <- read_csv(args[2], col_names = c("id", "label"), col_types = cols())
control <- read_csv(args[3], col_names = FALSE, col_types = cols())

output_suffix <- args[2] %>% str_remove(".*labels_")

score_MIDS <- matrix(NA,4,ncol(df_que))
rownames(score_MIDS) <- c("M","I", "D", "S")

for (row in 1:4){
    mutation = rownames(score_MIDS)[row]
    for(col in 1:ncol(df_que)){
        score_MIDS[row,col] <-
            str_count(pull(df_que[,col]),pattern = mutation) %>%
            sum()
    }
}
score_MIDS[1,] <- 0
score_MIDS["D",] <- score_MIDS["D",]*-1

for(col in 1:ncol(df_que)){
    res <- score_MIDS[,col]
    df_que[,col] <- pull(df_que[,col]) %>% res[.]
}

df_que[,pull(control)==2] <- 0
df_que <- df_que[, colSums(df_que) != 0]

rm(score_MIDS, control)

# args <- commandArgs(trailingOnly = TRUE)
# df_que <- read_csv(args[1], col_names = F, col_types = cols(), delim = ",")
# df_label <- read_csv(args[2], col_names = c("id", "label"), col_types = cols())
# output_suffix <- args[2] %>% str_remove(".*labels_")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Format input
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# df_que_new <- df_score * t(df_weight[1,])

# test1 <- df_res[1:2,50:75] %>% as.matrix
# test2 <- df_que[1:2,50:75] %>% as.matrix
# rownames(test1) <- rownames(test2)
# colnames(test1) <- colnames(test2)
# identical(test1,test2)
# test1
# test2

# df[,52] %>% table()

# df_que_rm0 <- df_que[, colSums(df_que) != 0]

# test1 <- df_que_rm0 %>% as.matrix
# test2 <- df_que_rm02 %>% as.matrix
# rownames(test1) <- rownames(test2)
# colnames(test1) <- colnames(test2)
# identical(test1,test2)
# test1
# test2

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# input: df_que
# --------------------------------------------------
pca_res <- prcomp(df_que, scale. = F)

components <- 1:10
output_pca <- pca_res$x[, components]
output_pca_importance <- summary(pca_res)$importance[2, components]
for (i in components){
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

cl <- hdbscan(output_pca, minPts = cl_sizes[cl_num_opt])
hdbscan_cl <- cl$cluster + 1

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Extract feature nucleotide position
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Cosine similarity to merge similar clusters
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

for (i in seq_along(pattern_)){
    output_cl[output_cl == pattern_[i]] <- query_[i]
}

result <- tibble(read_id = df_label$id, output_cl)
write_tsv(result,
    sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
    col_names = F
)



# # TEST >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# hdbscan_cl %>% table()
# pca_plot <- as_tibble(output_pca) %>%
#     bind_cols(cluster = factor(hdbscan_cl)) %>%
#     select(PC1, PC2, cluster)

# g <- ggplot(
#     data = pca_plot,
#     aes(
#         x = PC1, y = PC2,
#         color = cluster
#     )
# ) +
#     geom_point(size = 3) +
#     theme_bw(base_size = 20) # +

# ggsave(
#     plot = g,
#     filename = "tmp_pca1.png",
#     width = 10, height = 8
# )


