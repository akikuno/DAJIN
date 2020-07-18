#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Install required packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

options(repos = 'https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Format
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

df <- read_csv(".DAJIN_temp/details/Details.csv", col_types = cols())
df <- df %>%
    rename_with(~ str_replace_all(.x, " ", "_")) %>%
    mutate(
        Allele_type = case_when(
            Design == "+" ~ "intact target",
            Allele_type == "wt" && Indel == "-" ~ "intact wt",
            TRUE ~ as.character(Allele_type)
        )
    )

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Plot
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

p <- ggplot(df, aes(x = Sample, y = `%_of_reads`, fill = Allele_type)) +
    geom_col(position = position_stack()) +
    scale_fill_manual(
        name = "Allele type",
        values = c(
            "abnormal" = "#C0C0C0",
            "intact wt" = "#3CB371",
            "wt" = "#96DBB4",
            "intact target" = "#FB6E52",
            "target" = "#FDB8AB"
            )) +
    labs(x = NULL, y = "% of reads") +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Save figure
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ggsave(p, filename = ".DAJIN_temp/details/Details.pdf",width = 12, height = 7)
