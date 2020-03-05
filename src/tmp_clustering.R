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

df <- read_csv("test.csv", col_names = F)
label <- read_tsv("test_labels.txt", col_names = c("id", "labels"))

df <- t(df)
read_id <- label$id

# data <- df[,colSums(df)>10]
data <- data.frame(df) %>% bind_cols(labels = label$id)
# data <- data[,colSums(data %>% select(-labels))>10]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tmp <- data %>%
  select(-labels) %>%
  as.matrix()
pca_res <- prcomp(tmp)
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


result <- tibble(read_id, cl$cluster + 1) %>% arrange(read_id)
write_tsv(result, "test_result.txt", col_names = F)
