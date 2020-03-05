# setwd("/mnt/ssd_2tb/ssd_folder/mizuno_sensei_nanopore/190921-4th-run/190921-4th/stx2_clustering")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

df_cont <- read_csv("test_cont.csv", col_names = F)
df <- read_csv("test.csv", col_names = F)
label <- read_tsv("test_labels.txt", col_names = c("id", "label"))

df_cont <- t(df_cont)
df <- t(df)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
data <- data.frame(df) %>% bind_cols(label = label$id)
# data <- data[, colSums(data %>% select(-label)) > 0]

tmp <- data %>%
  select(-label) %>%
  as.matrix()
# which(apply(tmp, 2, var) == 0)
# which(apply(tmp, 2, sum) == 0)
pca_res <- prcomp(tmp, scale. = F)
pca_res_reduction <- pca_res$x[, 1:5]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cl_nums <- c()
i <- 1
for (cl_num in seq(2, 500, by = 25)) {
  print(cl_num)
  cl <- hdbscan(pca_res_reduction, minPts = cl_num)
  cl_nums[i] <- cl$cluster %>%
    table() %>%
    length()
  i <- i + 1
}
cl_nums %>% table()
cl_num_opt <- cl_nums %>%
  table() %>%
  which.max() %>%
  names()
cl_num_opt <- which(cl_nums == cl_num_opt) %>% min()
cl_num_opt
cl_num_opt <- seq(2, 500, by = 25)[cl_num_opt]
cl <- hdbscan(pca_res_reduction, minPts = cl_num_opt)
cl$cluster %>% table()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract feature nucleotide position
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
zero_to_one <- function(x) (x - min(x)) / (max(x) - min(x))

df_cont_norm <- df_cont %>%
  as_tibble() %>%
  summarise_all(sum) %>%
  zero_to_one()
#
df_sum_before <- tibble(loc = c(1:ncol(df_cont)), cluster = "cont", score = t(df_cont_norm[1, ]))

# data <- data %>% bind_cols(cluster = cl$cluster + 1)
for (i in unique(cl$cluster + 1)) {
  tmp <- data[cl$cluster + 1 == i, ] %>%
    select(-label) %>%
    colSums() %>%
    zero_to_one()
  #
  tmp_df <- tibble(loc = c(1:ncol(df_cont)), cluster = sprintf("cluster%d", i), score = tmp)
  df_sum_before <- df_sum_before %>% bind_rows(tmp_df)
}
df_sum_before <- df_sum_before %>% mutate(size = ifelse(cluster == "cont", 5, 1))
g1 <- ggplot(df_sum_before, aes(x = loc, y = score)) +
  geom_point(aes(color = cluster, size = size))
g1

# subtract from control ------------------------------------
df_sum_after <- tibble(loc = c(1:ncol(df_cont)), cluster = "cont", score = t(df_cont_norm[1, ]))
for (i in unique(cl$cluster + 1)) {
  tmp <- data[cl$cluster + 1 == i, ] %>%
    select(-label) %>%
    colSums() %>%
    zero_to_one()
  #
  tmp <- tmp - df_cont_norm
  tmp[tmp < 0] <- 0
  tmp_df <- tibble(loc = c(1:ncol(df_cont)), cluster = sprintf("cl%d", i), score = t(tmp[1, ]))
  df_sum_after <- df_sum_after %>% bind_rows(tmp_df)
}

df_sum_after <- df_sum_after %>% mutate(size = ifelse(cluster == "cont", 3, 2))
# g2 <- ggplot(df_sum_after, aes(x = loc, y = score)) +
#   geom_point(aes(color = cluster, size = size))
# g2

# g3 <- ggplot(df_sum_after %>% filter(cluster != "cont"), aes(x = loc, y = score)) +
#   geom_point(aes(color = cluster, size = size))
# g3


# df_sum_before %>% filter(loc == 549)
# df_sum_after %>% filter(loc == 549)

df_cont_blacklist <- df_sum_after %>%
  filter(cluster == "cont") %>%
  filter(score > quantile(score, 0.99))

df_filtered <- df_sum_after %>%
  group_by(cluster) %>%
  filter(score > quantile(score, 0.99)) %>%
  filter(!(loc %in% df_cont_blacklist$loc))

# df_c2 <- df_sum_after %>%
#   filter(cluster == "cluster2") %>%
#   filter(score > quantile(score, 0.99)) %>%
#   filter(!(loc %in% df_cont_blacklist$loc))
# plot(df_c2$score)

# df_c3 <- df_sum_after %>%
#   filter(cluster == "cluster3") %>%
#   filter(score > quantile(score, 0.99)) %>%
#   filter(!(loc %in% df_cont_blacklist$loc))
# plot(df_c3$score)

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

ggsave(plot = g, filename = "test.png", width = 10, height = 8)
#
result <- df_filtered %>%
  select(loc, cluster) %>%
  arrange(cluster)
write_tsv(result, "test_position.txt", col_names = F)
#
result <- tibble(read_id = label$id, sprintf("cl%d", cl$cluster + 1)) %>% arrange(read_id)
write_tsv(result, "test_result.txt", col_names = F)
