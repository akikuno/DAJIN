# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan, vroom, umap)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# suffix <- "barcode02_wt"
# que <- paste0(".DAJIN_temp/clustering//4_score_", sprintf("%s",suffix), "_9797")
# label <- paste0(".DAJIN_temp/clustering//2_labels_", sprintf("%s", suffix), "_9797")
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
df_que <- as.data.frame(df_que)
# dim(df_que)
# png("hoge1.png")
# plot(as.integer(df_que[1, ]))
# dev.off()
# png("hoge2.png")
# plot(as.integer(df_que[2, ]))
# dev.off()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Filter non-informative information
# 0付近の値
# 全てのリードに共通している値
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# tmp_sum <- as.data.frame(df_que)
# apply(tmp_sum, 2, function(x) quantile(x, 0.6)) %>% summary()
# apply(tmp_sum, 2, function(x) quantile(x, 0.5)) %>% summary()
# apply(tmp_sum, 2, function(x) quantile(x, 0.4)) %>% summary()
# summary(tmp_sum[1,])

# summary(tmp_sum)
# png("hoge_sum.png")
# plot(tmp_sum)
# dev.off()

tmp_sum <- colSums(df_que)
tmp_sum %>% table() %>% which.max() %>% names()
plus <- nrow(df_que) * sum(tmp_sum > 0) / length(tmp_sum)
minus <- -1 * nrow(df_que) * sum(tmp_sum < 0) / length(tmp_sum)

tmp <- (tmp_sum > plus | tmp_sum < minus)

# test <- tmp_sum %>% table()
# test[1811]
# sum(test==1811)
# quan_min <- quantile(colSums(df_que), 0.05)
# quan_mid <- c(quantile(colSums(df_que), 0.45), quantile(colSums(df_que), 0.55))
# quan_max <- quantile(colSums(df_que), 0.95)
# tmp <- (colSums(df_que) > quan_min) * (colSums(df_que) < quan_max)
# tmp <- tmp * (colSums(df_que) < quan_mid[1]) | (colSums(df_que) > quan_mid[2])
sum(tmp)
#
df_filtered <- df_que[, tmp]
# df_filtered <- df_que[, tmp_sum!=0]
dim(df_filtered)
# png("fuga1.png")
# plot(as.integer(df_filtered[1, ]))
# dev.off()
# png("fuga2.png")
# plot(as.integer(df_filtered[2, ]))
# dev.off()

# df_que <- df_que[, colSums(df_que) > summary(colSums(df_que))[3]]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# print("Dimension reduction by PCA...")
input_pca <- df_que
input_pca <- df_filtered
# output: output_pca # PC1 and PC2
# //////////////////////////////////////////////////////////
# pca_res <- prcomp(input_pca, scale. = F)
# output_pca <- pca_res$x[, 1:2]

pca_res <- umap(input_pca)
output_pca <- pca_res$layout
colnames(output_pca) <- c("PC1", "PC2")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# print("Clustering by HDBSCAN...")
input_hdbscan <- output_pca
# output: output_hdbscan # Cluster + 1
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
#   print(cl_size)
#   print(cl$cluster %>% table())
  cl_nums[i] <- cl$cluster %>%
    table() %>%
    length()
}
# cl_nums %>% table()
cl_num_opt <- cl_nums %>% table()
# もし最頻クラスタ数が1つだった場合は次に頻度の高いクラスタ数にする。
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
# クラスタ数が1つしか無い場合には1を代入する
# cl_num_opt
cl <- hdbscan(input_hdbscan, minPts = cl_sizes[cl_num_opt])
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
    theme_bw(base_size = 20) # +
    # theme(legend.position = "none")
ggsave(plot = g,
    filename = sprintf(".DAJIN_temp/clustering/pca_%s.png", output_suffix),
    width = 10, height = 8)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract feature nucleotide position
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_que <- df_filtered
input_cl <- output_hdbscan
# Stx2 #26 abnormalのケースでは、
# 0-1正規化をかますとCosine類似度が0.95ほどになり、理想的なクラスタができないため
# 正規化しないことにする。
# //////////////////////////////////////////////////////////

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
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cosine similarity to merge similar clusters
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Merge similar clusters...")
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
        #
        # png("hoge_score1.png")
        # plot(df_1$score)
        # dev.off()
        # png("hoge_score2.png")
        # plot(df_2$score)
        # dev.off()

        df_ <- tibble(
            one = cl_combn[1, i],
            two = cl_combn[2, i],
            score = cosine_sim(df_1$score, df_2$score)
        )
        df_cossim <- bind_rows(df_cossim, df_)
    }
    df_cossim <- df_cossim %>% filter(score > 0.90)
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
# output_cl %>% unique()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
result <- tibble(read_id = df_label$id, output_cl)
write_tsv(result,
    sprintf(".DAJIN_temp/clustering/hdbscan_%s", output_suffix),
    col_names = F)

