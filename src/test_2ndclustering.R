# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, dbscan, vroom)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

data <- vroom("test_tmp_cluster", col_names = F)
colnames(data) <- c("coodinate", "value", "mutation", "size", "cluster")
data_error <- data %>% filter(cluster == "control") %>% filter(mutation!="Match") %>% select(coodinate)

data_error$coodinate
# data[]
data %>% filter(cluster != "control") %>% filter(coodinate %in% data_error$coodinate) %>% mutate(mutation = "Match", size = 1)
data %>%
    group_by(cluster, mutation) %>%
    count()
g <- ggplot(data, aes(x = coodinate, y = value, color = mutation)) +
    geom_point(aes(size = size)) +
    scale_color_manual(values = c("blue", "red", "gray", "green")) +
    theme_bw(base_size = 20) +
    facet_wrap(~cluster, ncol = 1)
g
ggsave("test_barcode29_alleleall.png", g, width = 8, height = 4)

data <- vroom("test_tmp_allele3", col_names = F)
colnames(data) <- c("loc", "value", "mutation", "size")
data %>%
    group_by(mutation) %>%
    count()
g <- ggplot(data, aes(x = loc, y = value, color = mutation)) +
    geom_point(aes(size = size)) +
    scale_color_manual(values = c("blue", "red", "gray", "green")) +
    theme_bw()

ggsave("test_barcode29_allele3.png", g, width = 8, height = 1)

