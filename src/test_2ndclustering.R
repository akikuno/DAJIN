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
data <- vroom(args[1], col_names = F, col_types = cols())
cutsites <- vroom(args[2], col_names = F, col_types = cols(), delim = ",")
output_suffix <- args[1] # %>% str_remove(".*plot_")

# data <- vroom(".DAJIN_temp/clustering//5_plot_barcode02_wt_9797", col_names = F); output_suffix <- "fuga"
# cutsites <- vroom(".DAJIN_temp/clustering/tmp_mutation_barcode02_wt_9797", col_names = F, delim=",")
colnames(data) <- c("coodinate", "value", "mutation", "cluster")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Color pallet
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
my_colors <- c("gray50", "red", "blue", "green", "pink")
names(my_colors) <- c("Match",
    "Insertion",
    "Deletion",
    "Substitution",
    "Target")
col_scale <- scale_colour_manual(name = "mutation", values = my_colors)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
g <- ggplot(data, aes(x = coodinate, y = value, color = mutation)) +
    geom_point(aes(size = 1))

for (i in cutsites) {
    g <- g + geom_vline(
        xintercept = i, linetype = "dotted",
        color = "black", size = 1
    )
}

if(nrow(data %>% filter(mutation != "Match")) > 0) {
    g <- g + geom_point(
        data = (data %>% filter(mutation != "Match")),
        aes(size = 2))
}

g <- g +
    col_scale +
    guides(size = FALSE) +
    ggtitle(output_suffix) +
    xlab("Genomic Coodinate") +
    ylab("") +
    theme_bw(base_size = 20) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        axis.text.y = element_blank()
    ) +
    facet_wrap(~cluster, ncol = 1)
# g
ggsave(sprintf("%s.png", output_suffix), g, width = 10, height = 4)
# ggsave(sprintf("%s.svg", output_suffix), g, width = 15, height = 4)
