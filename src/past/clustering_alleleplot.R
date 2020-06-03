# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install required packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, vroom)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
args <- commandArgs(trailingOnly = TRUE)

data <- vroom(args[1], col_names = F, col_types = cols())
cutsites <- vroom(args[2], col_names = F, delim = ",", col_types = cols())
colnames(data) <- c("coord", "mut", "cluster", "ins_num")

# tmp <- ".DAJIN_temp/tmp_barcode04_target_HOGE_3"
# data <- vroom(".DAJIN_temp/tmp_barcode04_target_HOGE_3", col_names = F)
# cutsites <- vroom(".DAJIN_temp/clustering/tmp_mutation_barcode04_target_HOGE", col_names = F, delim=",")
# colnames(data) <- c("coord", "mut", "cluster", "ins_num")

output_suffix <- args[1] %>%
    str_remove(".*tmp_") %>%
    str_remove("HOGE_|FUGA_|FOO_")
output_dir <- "DAJIN_Report/alleletypes/"
# print(output_suffix)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Color pallet
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

my_colors <- c("gray80", "#FC4E07", "#56B4E9", "#009E73", "#FFA100")
names(my_colors) <- c("M", "I", "D", "S", "T")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

g_data <- data %>%
    select(-cluster, -ins_num) %>%
    mutate(value = "@")
g_data$mut <- factor(g_data$mut, levels = c("M", "I", "D", "S", "T"))

g <- ggplot(g_data, aes(x = coord, y = value, color = mut)) +
    geom_point(color = "gray80")

for (i in cutsites) {
    g <- g + geom_vline(
        xintercept = i,
        linetype = "dotted",
        color = "black",
    )
}

g <- g + geom_point(data = (g_data %>% filter(mut == "T")))

if(nrow(g_data %>% filter(!str_detect(mut, "M|T"))) > 0) {
    g <- g + geom_point(
        data = (g_data %>% filter(!str_detect(mut, "M|T"))),
        aes(size = 2))
}
g <- g +
    scale_color_manual(name = "mut", values = my_colors) +
    theme_light() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
    )
ggsave(sprintf("%s/%s.png",output_dir, output_suffix), g, width = 3, height = 1, dpi = 350)
