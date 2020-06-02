# ============================================================================
# Environment setting ----
# ============================================================================
if(!requireNamespace("pacman", quietly = T)) install.packages("pacman")
if(!requireNamespace("reticulate", quietly = T)) install.packages("reticulate")
reticulate::use_condaenv("DAJIN")
pacman::p_load(tidyverse, keras, reticulate, caret, e1071, gridExtra, svglite, progress)
# reticulate::py_config()

np <- import("numpy")
pd <- import("pandas")
np_utils <- import("keras.utils.np_utils")
model_selection <- import("sklearn.model_selection")
tf <- import("tensorflow")
# tf_logger <- tf$get_logger()
# tf_logger$setLevel("INFO")
# tf$autograph$set_verbosity(1)
# logging <- import("logging")
# tf_logger$setLevel(logging$ERROR)

# os <- import("os")
# os$environ['TF_CPP_MIN_LOG_LEVEL'] <- '3' 

regularizers <- import("keras.regularizers")
ftools <- import("functools")
# <<<<<<<<<<

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
# file_name <- ".DAJIN_temp/data/DAJIN.txt"

output_dir <- ".DAJIN_temp/results"
output_npz <- file_name %>% str_replace(".txt*$", ".npz")
output_model <- file_name %>% str_replace(".txt*$", ".h5")
output_file_name <- file_name %>% str_remove(".*/") %>% str_remove(".txt.*$")

df_anomaly <- read_tsv(".DAJIN_temp/data/DAJIN_anomaly_classification.txt",
  col_names = c("seqID", "barcodeID", "abnormal_prediction"),
col_types = cols()
)
labels_index <- read_tsv(
  ".DAJIN_temp/data/DAJIN_anomaly_classification_labels.txt",
  col_names = "label", col_types = cols()
)


# ============================================================================
# # Load One hot data and model
# ============================================================================
np_load <- ftools$partial(np$load, allow_pickle = TRUE)
npz <- np_load(output_npz)
X_real <- npz["X_real"]

model <- tf$keras$models$load_model(output_model)

# ============================================================================
# # Prediction
# ============================================================================
print("Predict allele types...")
predict <- model %>% predict(X_real)
predict <- apply(predict, 1, function(x) which.max(x))

df_result <- df_anomaly %>%
  mutate(predict = as.character(predict)) %>%
  mutate(predict = if_else(abnormal_prediction == "abnormal", "abnormal", predict)) %>%
  select(barcodeID, seqID, predict) %>%
  arrange(seqID)

for(i in seq_along(labels_index$label)) {
  df_result$predict <- df_result$predict %>%
    str_replace(as.character(i), as.character(labels_index$label[i])) %>%
    str_remove("_simulated.*")
}

# ============================================================================
# ## Output result
# ============================================================================
write_tsv(df_result, ".DAJIN_temp/data/DAJIN_prediction_result.txt", col_names = F)

df_plot <- df_result %>%
  group_by(barcodeID, predict) %>%
  count(name = "value")

g_bar <- ggplot(df_plot, aes(x = barcodeID, y = value, fill = predict)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sample ID", y = "DAJIN Predicted Allele Proportion") +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(as.factor(df_plot$barcodeID)))) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )
ggsave(sprintf("%s/DAJIN_allele_prediction.jpg", output_dir), g_bar, width = 10, height = 10)
