library(tidyverse)

df_lof <- read_csv("results_lof_add.csv")
colnames(df_lof) <- c("barcodeID", "outliers", "score", "metric", "iter")

df_lof$barcodeID = df_lof$barcodeID %>% str_remove("abtest_target_") %>% str_remove("_simulated")

df_lof %>% filter(outliers=="abnormal") %>% filter(barcodeID=="negacont") %>% filter(metric=="jaccard")
df_lof %>% filter(metric=="jaccard")

df_lof %>% filter(outliers=="abnormal") %>%
    ggplot(aes(x=barcodeID, y=score)) +
    geom_point() +
    ylim(0,1000) +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(width=0.2, height=0.05, alpha=0.35, size=0.75) +
    facet_wrap(~ metric)