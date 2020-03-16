# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan, vroom)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

args <- commandArgs(trailingOnly = TRUE)
df_que <- vroom(args[1], col_names = F)
df_label <- vroom(args[2], col_names = c("id", "label"))
output_suffix <- args[2] %>% str_remove(".*labels_")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Format input
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

df_que <- t(df_que)
df_que[is.na(df_que)] <- 0

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_pca <- df_que
# output: output_pca # PC1 and PC2
# //////////////////////////////////////////////////////////
pca_res <- prcomp(input_pca, scale. = F)
output_pca <- pca_res$x[, 1:2]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_hdbscan <- output_pca
# output: output_hdbscan # Cluster + 1
# //////////////////////////////////////////////////////////
cl_nums <- c()
cl_sizes <- seq(25, 250, by = 25)
i <- 1
for (i in seq_along(cl_sizes)) {
  cl <- hdbscan(input_hdbscan, minPts = cl_sizes[i])
#   print(cl_size)
#   print(cl$cluster %>% table())
  cl_nums[i] <- cl$cluster %>%
    table() %>%
    length()
}
# cl_nums %>% table()
cl_num_opt <- cl_nums %>%
  table() %>%
  which.max() %>%
  names() #
cl_num_opt <- which(cl_nums == cl_num_opt) %>% min()
# cl_num_opt
cl_num_opt <- cl_sizes[cl_num_opt]
cl <- hdbscan(input_hdbscan, minPts = cl_num_opt)
output_hdbscan <- cl$cluster + 1

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
g <- ggplot(
    data = as.data.frame(output_pca),
    aes(
        x = PC1, y = PC2,
        color = factor(output_hdbscan)
    )
) +
    geom_point(size = 3) +
    theme_bw(base_size = 20) +
    theme(legend.position = "none")
ggsave(plot = g, filename = sprintf("pca_%s.png", output_suffix), width = 10, height = 8)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract feature nucleotide position
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_que <- df_que
input_cl <- output_hdbscan
# //////////////////////////////////////////////////////////

zero_to_one <- function(x) (x - min(x)) / (max(x) - min(x))

df_cluster <- tibble(loc = integer(), cluster = integer(), score = double())
for (i in unique(input_cl)) {
  tmp_score <- input_que[input_cl == i, ] %>%
    colSums() %>%
    zero_to_one()

  tmp_df <- tibble(
    loc = seq_len(ncol(input_que)),
    cluster = i,
    score = tmp_score
  )
  df_cluster <- df_cluster %>% bind_rows(tmp_df)
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cosine similarity to merge similar clusters
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
    i <- 8
    for (i in seq_along(1:ncol(cl_combn))) {
        df_1 <- input_cossim %>%
            filter(cluster == cl_combn[1, i]) %>%
            select(score)
        #
        df_2 <- input_cossim %>%
            filter(cluster == cl_combn[2, i]) %>%
            select(score)
        #
        df_ <- tibble(
            one = cl_combn[1, i],
            two = cl_combn[2, i],
            score = cosine_sim(df_1$score, df_2$score)
        )
        df_cossim <- bind_rows(df_cossim, df_)
    }
    df_cossim <- df_cossim %>% filter(score > 0.80)
    #
    if (nrow(df_cossim) != 0) {
        for (i in 1:nrow(df_cossim)) {
            pattern_ <- df_cossim[i, ]$two
            query_ <- df_cossim[i, ]$one
            output_cl[output_cl == pattern_] <- query_
        }
    }
}

pattern_ <- output_cl %>%
    unique() %>%
    sort()
query_ <- output_cl %>%
    unique() %>%
    order() %>%
    sort()

for(i in seq_along(pattern_)) output_cl[output_cl == pattern_[i]] <- query_[i]
output_cl %>% unique()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
result <- tibble(read_id = df_label$id, output_cl)
write_tsv(result, sprintf(".tmp_/clustering_id_%s", output_suffix), col_names = F)

