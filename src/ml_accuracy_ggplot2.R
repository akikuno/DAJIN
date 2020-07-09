library(tidyverse)

df <- read_csv("accuracy_anomaly_detection.csv",
    col_names=c("id","value","L2", "MIDS"))

df$id <- df$id %>%
    str_replace("abtest_wt_","") %>%
    str_replace("_simulated","")

p <- ggplot(df, aes(x=id, y=value)) +
    geom_boxplot() +
    # geom_point() +
    facet_wrap(~ L2 + MIDS)

ggsave("accuracy_anomaly_detection.png", p, width = 12)