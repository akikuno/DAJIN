#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Install required packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

options(repos = 'https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, RColorBrewer)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Format
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

df <- read_csv(".DAJIN_temp/details/Details.csv", col_types = cols())
df <-
    df %>%
    rename_with(~ str_replace_all(.x, " ", "_")) %>%
    mutate(
        Allele_type = case_when(
            Design == "+" ~ "intact target",
            Allele_type == "wt" & Indel == "-" ~ "intact wt",
            TRUE ~ as.character(Allele_type)
        )
    )


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Plot
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#==========================================================
#? Color setting
#==========================================================

color <- c(
        "abnormal" = "#C0C0C0",
        "intact wt" = "#3CB371",
        "wt" = "#77D9A8",
        "intact target" = "#ff4b00",
        "target" = "#F6AA00"
        )

allele_others <- df$Allele_type %>%
        str_remove("abnormal|wt|target") %>%
        str_remove("intact ") %>%
        unique() %>%
        .[-1]

if (length(allele_others) > 0) {
    color_brewer <- c(
        brewer.pal(n = 8, "Pastel2"),
        brewer.pal(n = 8, "Set1"),
        brewer.pal(n = 8, "Set2"),
        brewer.pal(n = 12, "Set3")
        )

    color_names <- allele_others %>%
        seq_along() %>%
        color_brewer[.]

    color_others <- allele_others %>%
        set_names(color_names, .)

    color <- c(color, color_others)
}

#==========================================================
#? Plot
#==========================================================

p <-
    ggplot(df, aes(x = Sample, y = `%_of_reads`, fill = Allele_type)) +
    geom_col(position = position_stack(), color = "black", size = 0.5) +
    scale_fill_manual(
        name = NULL,
        values = color) +
    labs(x = NULL, y = "Percentage of reads") +
    theme_bw(base_size = 20) +
    theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Save figure
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

width <-
    df$Sample %>% unique %>% length %>% `*`(0.5) %>% `+`(4)

ggsave(p, filename = ".DAJIN_temp/details/Details.pdf",
    width = width, height = 7)
