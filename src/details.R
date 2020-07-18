#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Install required packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = 'https://cloud.r-project.org/')
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, gghighlight)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! TEST
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

df <- read_csv(".DAJIN_temp/details/Details.csv")
df <- df %>% rename_with(~ str_replace(.x, " ", "_"))

p <- ggplot(df, aes(x = Sample, y = `%_of_reads`, fill = Allele_type)) +
    geom_col(position = position_stack()) +
    gghighlight(Design == '+')
p