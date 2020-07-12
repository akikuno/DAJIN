library(tidyverse)

df <- read_csv("accuracy_anomaly_detection.csv",
    col_names=c("id","Accuracy","L2", "MIDS"),
    col_types=cols())

df$id <- df$id %>%
    str_replace("abtest_wt_","") %>%
    str_replace("_simulated","")

df <- df %>% filter(!str_detect(id, "ins"))

p <- ggplot(df, aes(x=id, y=Accuracy)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.2, height=0.05, alpha=0.35, size=0.75) +
    theme_bw() +
    facet_wrap(~ L2 + MIDS)

ggsave("accuracy_anomaly_detection.png", p, width = 12)
ggsave("accuracy_anomaly_detection.pdf", p, width = 12)