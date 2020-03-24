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
colnames(data) <- c("coordinate", "M", "I", "D", "S", "cluster")

cutsites <- vroom(args[2], col_names = F, col_types = cols(), delim = ",")

output_suffix <- args[1]
print(output_suffix)

g_data <- data %>%
    select(-cluster) %>%
    pivot_longer(-coordinate, names_to = "mut", values_to = "score")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Color pallet
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

my_colors <- c("white", "#FC4E07", "#56B4E9", "#009E73")
names(my_colors) <- c(
    "M",
    "I",
    "D",
    "S"
)
col_scale <- scale_colour_manual(name = "mut", values = my_colors)

g_data$mut <- factor(g_data$mut, levels = c("M", "I", "D", "S"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Stacked plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

g <- ggplot(g_data, aes(x = coordinate, y = score, fill = mut))

for (i in cutsites) {
    g <- g + geom_point(
        x = i, y = 1.05, shape = 25, size = 1, color = "red", fill = "red"
    )
}

g <- g +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(name = "mut", values = my_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(y = "percent") +
    theme_light() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()
    )


ggsave(sprintf("%s.png", output_suffix), g, width = 5, height = 1, dpi = 2400)



# test import -------------------------------------------
# data <- vroom(".DAJIN_temp/clustering//plot_barcode26_abnormal_5018", col_names = F); output_suffix <- "fuga"
# cutsites <- vroom(".DAJIN_temp/clustering/tmp_mutation_barcode26_abnormal_5018", col_names = F, delim=",")
# data <- vroom("tmp", col_names = F); output_suffix <- "fuga"
# colnames(data) <- c("coordinate", "M", "I", "D", "S", "cluster")

# colnames(data) <- c("mutation", "coordinate", "M", "I", "D", "S")
# test -------------------------------------------


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# for (i in cutsites) {
#     g <- g + geom_segment(
#         x = i, xend = i,
#         y = -0.1, yend = 0,
#     arrow = arrow(length = unit(1, "npc"))
#     )
# }
# for (i in cutsites) {
#     g <- g + geom_vline(
#         xintercept = i, linetype = "dotted",
#         color = "black", size = 0.5, alpha = 0.5
#     )
# }

# g <- ggplot(data, aes(x = coodinate, y = value, color = mutation)) +
#     geom_point(aes(size = 1))

# if(nrow(data %>% filter(mutation != "Match")) > 0) {
#     g <- g + geom_point(
#         data = (data %>% filter(mutation != "Match")),
#         aes(size = 2))
# }

# for (i in cutsites) {
#     g <- g + geom_vline(
#         xintercept = i, linetype = "dotted",
#         color = "black", size = 1
#     )
# }

# g <- g +
#     col_scale +
#     guides(size = FALSE) +
#     # ggtitle(output_suffix) +
#     xlab("Genomic Coodinate") +
#     ylab("") +
#     theme_bw(base_size = 20) +
#     theme(
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#     ) +
#     facet_wrap(~cluster)
# # g
# ggsave(sprintf("%s.png", output_suffix), g, width = 10, height = 3)
# ggsave(sprintf("%s.svg", output_suffix), g, width = 15, height = 4)


# g <- ggplot(data, aes(coordinate, genome))
# g <- g + geom_point(alpha = data$M, color = "gray50", size = data$M * 5, shape = 15)
# g <- g + geom_point(alpha = data$S, color = "green", size = data$S * 10, shape = 15)
# g <- g + geom_point(alpha = data$I, color = "red", size = data$I * 10, shape = 15)
# g <- g + geom_point(alpha = data$D, color = "blue", size = data$D * 10, shape = 15)

# for (i in cutsites) {
#     g <- g + geom_vline(
#         xintercept = i, linetype = "dotted",
#         color = "black", size = 1
#     )
# }
# ggsave(sprintf("%s.png", output_suffix), g, width = 10, height = 2)


# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Import data
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# args <- commandArgs(trailingOnly = TRUE)
# data <- vroom(args[1], col_names = F, col_types = cols())
# cutsites <- vroom(args[2], col_names = F, col_types = cols(), delim = ",")
# output_suffix <- args[1]
# # output_suffix <- str_c(c(args[1], args[3]), collapse = "_")
# print(output_suffix)
# # data <- vroom(".DAJIN_temp/clustering//plot_barcode26_abnormal_5018", col_names = F); output_suffix <- "fuga"
# # cutsites <- vroom(".DAJIN_temp/clustering/tmp_mutation_barcode26_abnormal_5018", col_names = F, delim=",")
# colnames(data) <- c("coodinate", "value", "mutation", "cluster", "insertion_num")

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Color pallet
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# my_colors <- c("gray50", "red", "blue", "green", "pink")
# names(my_colors) <- c("Match",
#     "Insertion",
#     "Deletion",
#     "Substitution",
#     "Target")
# col_scale <- scale_colour_manual(name = "mutation", values = my_colors)

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Plot
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# g <- ggplot(data, aes(x = coodinate, y = value, color = mutation)) +
#     geom_point(aes(size = 1))

# if(nrow(data %>% filter(mutation != "Match")) > 0) {
#     g <- g + geom_point(
#         data = (data %>% filter(mutation != "Match")),
#         aes(size = 2))
# }

# for (i in cutsites) {
#     g <- g + geom_vline(
#         xintercept = i, linetype = "dotted",
#         color = "black", size = 1
#     )
# }

# g <- g +
#     col_scale +
#     guides(size = FALSE) +
#     # ggtitle(output_suffix) +
#     xlab("Genomic Coodinate") +
#     ylab("") +
#     theme_bw(base_size = 20) +
#     theme(
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#     ) +
#     facet_wrap(~cluster)
# # g
# ggsave(sprintf("%s.png", output_suffix), g, width = 10, height = 3)
# # ggsave(sprintf("%s.svg", output_suffix), g, width = 15, height = 4)
