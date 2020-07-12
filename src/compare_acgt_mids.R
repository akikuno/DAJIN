library(tidyverse)

df <- read_csv("compare_acgt_mids.csv",
    col_names=c("umap1","umap2","id", "conv","bfaf"),
    col_types=cols())

df$id <- df$id %>%
    str_replace("abtest_wt_","") %>%
    str_replace("_simulated","") %>%
    str_replace("nega","")

df <- df %>% filter(!str_detect(id, "ins"))

p <- ggplot(df, aes(x=umap1, y=umap2, fill=id)) +
    geom_point(aes(color = id, shape = id, alpha = 0.5)) +
    theme_bw() +
    facet_wrap(~ conv + bfaf)

ggsave("umap_compare_acgt_mids.png", p)
ggsave("umap_compare_acgt_mids.pdf", p)
