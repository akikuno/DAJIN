library(tidyverse)

df_lof <- read_csv("results_lof_neighbor.csv",
    c("barcodeID", "outliers", "score", "neighbors", "iter"))

df_lof$barcodeID = df_lof$barcodeID %>%
    str_remove("abtest_target_") %>%
    str_remove("_simulated") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("negacont", "cont")

df_lof %>% filter(outliers=="abnormal") %>%
    ggplot(aes(x=barcodeID, y=score)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.2, height=0.05, alpha=0.35, size=0.75) +
    ylim(0,1000) +
    labs(x = NULL, y = "Numbers of the detected abnormal reads") +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~ neighbors, ncol = 4)