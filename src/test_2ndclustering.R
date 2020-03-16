# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, vroom)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
args <- commandArgs(trailingOnly = TRUE)
data <- vroom(args[1], col_names = F)
output_suffix <- args[1] %>% str_remove(".*plot_")
# data <- vroom::vroom(".tmp_/clustering_plot_barcode04_target", col_names = F); output_suffix <- "fuga"
colnames(data) <- c("coodinate", "value", "mutation", "cluster")

myColors <- c("gray50", "red", "blue", "green")
names(myColors) <- c("Match","Insertion","Deletion","Substitution")
colScale <- scale_colour_manual(name = "mutation", values = myColors)

# data %>%
#     group_by(cluster, mutation) %>%
#     count()
g <- ggplot(data, aes(x = coodinate, y = value, color = mutation)) +
    geom_point(data=(data %>% filter(mutation=="Match")), aes(size = 1)) +
    geom_point(data=(data %>% filter(mutation!="Match")), aes(size = 2)) +
    colScale +
    guides(size = FALSE) +
    ggtitle(output_suffix) +
    xlab("Genomic Coodinate") +
    ylab("") +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    axis.text.y = element_blank()) +
    facet_wrap( ~ cluster, ncol = 1)
# g
ggsave(sprintf("%s.png", output_suffix), g, width = 15, height = 4)

