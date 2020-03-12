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
df_ref <- vroom(args[1], col_names = F)
df_que <- vroom(args[2], col_names = F)
df_label <- vroom(args[3], col_names = c("id", "label"))
output_suffix <- args[3] %>% str_remove(".*labels_")

# df_ref <- vroom(".tmp_/clustering_score_control_target_merged_target", col_names = F)
# df_que <- vroom(".tmp_/clustering_score_target_merged_target", col_names = F)
# df_label <- vroom(".tmp_/clustering_labels_target_merged_target", col_names = c("id", "label"))
# output_suffix <- ".tmp_/clustering_labels_target_merged_target" %>% str_remove(".*labels_")

df_ref <- t(df_ref)
df_que <- t(df_que)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_ <- df_que
output_ <- "pca_res_reduction"
# //////////////////////////////////////////////////////////

pca_res <- prcomp(input_, scale. = F)
assign(output_, pca_res$x[, 1:5])

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_ <- pca_res_reduction
output_ <- "cl$cluster"
# //////////////////////////////////////////////////////////
cl_nums <- c()
i <- 1
for (cl_num in seq(2, 500, by = 25)) {
  cl <- hdbscan(pca_res_reduction, minPts = cl_num)
  cl_nums[i] <- cl$cluster %>%
    table() %>%
    length()
  i <- i + 1
}
# cl_nums %>% table()
cl_num_opt <- cl_nums %>%
  table() %>%
  which.max() %>%
  names() #
cl_num_opt <- which(cl_nums == cl_num_opt) %>% min()
# cl_num_opt
cl_num_opt <- seq(2, 500, by = 25)[cl_num_opt]
cl <- hdbscan(pca_res_reduction, minPts = cl_num_opt)
# cl$cluster %>% table()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract feature nucleotide position
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_ref <- df_ref
input_que <- df_que
output_df <- NULL
# //////////////////////////////////////////////////////////
tmp <- input_ref %>%
  as_tibble() %>%
  summarise_all(sum, na.rm = TRUE) %>%
  t()
png("test.png")
plot(tmp[,1])
dev.off()


zero_to_one <- function(x) (x - min(x)) / (max(x) - min(x))

input_ref_norm <- input_ref %>%
  as_tibble() %>%
  summarise_all(sum, na.rm = TRUE) %>%
  zero_to_one()
#
# subtract from control ------------------------------------
df_sum_after <- tibble(
  loc = seq_len(ncol(input_ref)),
  cluster = "cont", score = t(input_ref_norm[1, ])
)

for (i in unique(cl$cluster)) {
  tmp <- input_que[cl$cluster == i, ] %>%
    colSums() %>%
    zero_to_one()
  #
  tmp <- tmp - input_ref_norm
  tmp[tmp < 0] <- 0
  tmp_df <- tibble(
    loc = seq_len(ncol(input_ref)),
    cluster = sprintf("allele%d", i + 1),
    score = t(tmp[1, ])
  )
  df_sum_after <- df_sum_after %>% bind_rows(tmp_df)
}
# df_sum_after %>% group_by(cluster) %>% count()

df_ref_blacklist_1 <- df_sum_after %>%
  filter(cluster == "cont") %>%
  filter(score == 1)

df_ref_blacklist_2 <- df_sum_after %>%
  filter(cluster == "cont") %>%
  filter(score != 1) %>%
  filter(score > quantile(score, 0.95))
# plot(df_sum_after[df_sum_after$cluster == "cont", ]$score)

df_ref_blacklist <- bind_rows(df_ref_blacklist_1, df_ref_blacklist_2)

output_df <- df_sum_after %>%
  filter(!(loc %in% df_ref_blacklist$loc))

# output_df %>% group_by(cluster) %>% count()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cosine similarity to merge similar clusters
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_ <- output_df
output_df <- NULL
output_cl <- sprintf("allele%d", cl$cluster + 1)
# //////////////////////////////////////////////////////////
cosine_sim <- function(a, b) crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))

cluster <- input_ %>%
  filter(cluster != "cont") %>%
  group_by(cluster) %>%
  count() %>%
  select(cluster)

if (nrow(cluster) > 1) {
  cl_combn <- combn(cluster$cluster, 2)
  df_cossim <- NULL
  i <- 1
  for (i in seq_len(ncol(cl_combn))) {
    df_1 <- input_ %>%
      filter(cluster == cl_combn[1, i]) %>%
      select(score)
    #
    df_2 <- input_ %>%
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
  df_cossim <- df_cossim %>% filter(score > 0.95)
  #
  for (i in seq_len(nrow(df_cossim))) {
    pattern_ <- df_cossim[i, ]$two
    query_ <- df_cossim[i, ]$one
    input_$cluster <- input_$cluster %>% str_replace(pattern_, query_)
    output_cl <- output_cl %>% str_replace(pattern_, query_)
  }
}

output_df <- input_ %>%
  filter(cluster != "cont") %>%
  group_by(cluster) %>%
  filter(score > quantile(score, 0.99)) %>%
  filter(!(loc %in% df_ref_blacklist$loc)) %>%
  select(loc, cluster) %>%
  distinct(loc) %>%
  arrange(cluster, loc)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
g <- ggplot(
  data = as.data.frame(pca_res_reduction[, c(1, 2)]),
  aes(
    x = PC1, y = PC2,
    color = factor(cl$cluster + 1)
  )
) +
  geom_point(size = 3) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")
ggsave(plot = g, filename = "test_pca.png", width = 10, height = 8)
#
write_tsv(output_df, sprintf(".tmp_/clustering_features_%s",output_suffix), col_names = F)
#
result <- tibble(read_id = df_label$id, output_cl)
write_tsv(result, sprintf(".tmp_/clustering_id_%s", output_suffix), col_names = F)
