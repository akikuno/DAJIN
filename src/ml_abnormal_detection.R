# One-hot encording in R
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
model_selection <- import("sklearn.model_selection")
tf <- import("tensorflow")
np_utils <- import("keras.utils.np_utils")
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

df <- read_tsv(file_name, col_names = c("ID", "SEQ", "LABEL"), col_types = cols())

df$SEQ <- str_c("MIDS=", df$SEQ, sep = "")
df_sim <- df %>% filter(str_detect(LABEL,"simulate"))
df_real <- df %>% filter(!str_detect(LABEL,"simulate"))

# ============================================================================
# # one-hot encoding ----
# ============================================================================
one_hot_encode <- function(seq){
  seq_len <- seq[1] %>% str_count()
  array_oh <- array(0, dim = c(length(seq), seq_len, 4))
  #
  list_seq <- seq %>%
    str_split("")
  #
  pb <- progress_bar$new(
      format = "  One-hot encording... [:bar] :percent eta: :eta",
      total = length(list_seq), clear = FALSE, width = 100)
  #
  for(i in seq_along(list_seq)){
    pb$tick()
    x <- list_seq[[i]]
    for(j in seq_along(x)){
      if(x[j] == "M") array_oh[i,j,] <- c(1,0,0,0)
      else if(x[j] == "I") array_oh[i,j,] <- c(0,1,0,0)
      else if(x[j] == "D") array_oh[i,j,] <- c(0,0,1,0)
      else if(x[j] == "S") array_oh[i,j,] <- c(0,0,0,1)
    }
  }
  return(array_oh)
}
# # ////
X_all <- one_hot_encode(df$SEQ)
#
X_sim <- X_all[str_detect(df$LABEL, "simulate"), , ]
X_real <- X_all[!str_detect(df$LABEL, "simulate"), , ]
#
# # ////
# np_load = ftools$partial(np$load, allow_pickle=TRUE)
# npz = np_load(output_npz)

# X_all = npz["X_all"]
# X_sim = npz["X_sim"]
# X_real = npz["X_real"]

labels_categorical <- pd$factorize(df_sim$LABEL)[[1]] %>%
  np_utils$to_categorical()
labels_id <- pd$factorize(df_sim$LABEL)[[2]] %>% str_remove("_aligned_reads")

train_test <- model_selection$train_test_split(
  X_sim, labels_categorical,
  test_size = 0.2, shuffle = TRUE
)

X_train <- train_test[[1]]
X_test <- train_test[[2]]
Y_train <- train_test[[3]]
Y_test <- train_test[[4]]

# ============================================================================
# model construction ---- 
# ============================================================================
# os <- import("os")
# os$environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

model <- keras_model_sequential()
model %>%
  layer_conv_1d(
    filters = 32, kernel_size = 32, activation = "relu",
    input_shape = c(dim(X_train)[2], 4),
    name = "1st_Conv1D"
  ) %>%
  layer_max_pooling_1d(pool_size = 4, name = "1st_MaxPooling1D") %>%
  layer_conv_1d(
    filters = 32, kernel_size = 16,
    activation = "relu", name = "2nd_Conv1D"
  ) %>%
  layer_max_pooling_1d(pool_size = 4, name = "2nd_MaxPooling1D") %>%
  layer_conv_1d(
    filters = 32, kernel_size = 8,
    activation = "relu", name = "3rd_Conv1D"
  ) %>%
  layer_max_pooling_1d(pool_size = 4, name = "3rd_MaxPooling1D") %>%
  layer_conv_1d(
    filters = 32, kernel_size = 4,
    activation = "relu", name = "4th_Conv1D"
  ) %>%
  layer_max_pooling_1d(pool_size = 4, name = "4th_MaxPooling1D") %>%
  layer_flatten(name = "flatten") %>%
  layer_dense(64, activation = "relu", name = "1st_FC")

# L2-norm
alpha <- c(0.1) # hyperparameter
model %>% layer_dense(64,
  activation = "linear",
  activity_regularizer = regularizers$l2(alpha), name = "L2-softmax"
)
#
model %>% layer_dense(length(labels_id), activation = "softmax", name = "final_layer")
# summary(model)

model %>% compile(
  optimizer = "adam",
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)
# <<<<<<<<<<<<<

print("Model fitting...")
suppressMessages(
    history <- model %>% fit(X_train, Y_train,
    epochs = 20, verbose = 0,
    validation_split = 0.2, shuffle = TRUE
  )
)
# plot(history)

# ============================================================================
# Abnormal allele detection ---------------
# ============================================================================

cosine_sim <- function(a, b) crossprod(a, b) / sqrt(crossprod(a) * crossprod(b))

get_cosine_score <- function(model, train, test){
  model_ <- keras_model(model$get_layer(index = 0L)$input,
                        model$get_layer(index = -2L)$output)
  # print(summary(model_))
  train_vector <- model_ %>% predict(train, verbose = 0)
  test_vector <- model_ %>% predict(test, verbose = 0)
  cosine_sim(train_vector[1,], test_vector[1,])
  # train_list <- split(train_vector, rep(1:nrow(train_vector), each = ncol(train_vector)))
  # score <- parallel::mclapply(train_list, function(x) {
  #   apply(x, 1, function(y) cosine_sim(test_vector[i, ], y)) %>% max()
  # })
  score <- rep(0, nrow(test_vector))
  pb <- progress_bar$new(
    format = "  Cosine similarity... [:bar] :percent eta: :eta",
    total = length(score), clear = FALSE, width = 100
  )
  for(i in seq_along(score)){
    pb$tick()
    score[i] <- apply(train_vector, 1, function(x) cosine_sim(test_vector[i, ], x)) %>% max()
    }
  return(score)
}
# cos_max <- function(train_vector, test_vector){
#   score <- apply(train_vector, 1, function(x) cosine_sim(test_vector[i, ], x)) %>% max()
#   return(score)
#   }

cos_all <- get_cosine_score(model, X_train[0:500, , ], X_all)
#
df_all <- tibble(barcode = df$LABEL, cos_similarity = cos_all)
df_all <- df_all %>% mutate(label = if_else(str_detect(barcode, "simulate"), "simulated", "real")) # -----------------------------
cos_threshold <- df_all %>%
  filter(label == "simulated") %>%
  summarize(quantile(cos_similarity, 0.001)) %>%
  unlist()

# ============================================================================
# Output the results ---------------
# ============================================================================

# Save one-hot encording ...
np$savez_compressed(output_npz,
  X_all = X_all,
  X_sim = X_sim,
  X_real = X_real
)
# Save the model ...
model$save(output_model)
# Save Anomaly annotation ...
df_anomaly <- df_all %>%
  filter(label == "real") %>%
  mutate(abnormal_prediction = if_else(cos_similarity > cos_threshold, "normal", "abnormal")) %>%
  bind_cols(., df_real) %>%
  select(ID, barcode, abnormal_prediction)
#
write_tsv(df_anomaly, ".DAJIN_temp/data/DAJIN_anomaly_classification.txt", col_names = F)
# Save labels
write_tsv(data.frame(labels_id), ".DAJIN_temp/data/DAJIN_anomaly_classification_labels.txt", col_names = F)

# confusion matrix for classification -----------------------
predicted_labels <- model %>% predict_classes(X_test)
df_id <- data.frame(label = seq_along(labels_id) - 1, id = labels_id)
df_pred <- data.frame(label = predicted_labels)
predicted_labels <- inner_join(df_id, df_pred, by = "label") %>%
  .$id %>%
  as.factor()

df_true <- apply(Y_test, 1, which.max) %>% data.frame(label = . - 1)
true_labels <- inner_join(df_id, df_true, by = "label") %>%
  .$id %>%
  as.factor()

cm <- confusionMatrix(predicted_labels, true_labels)
cm <- data.frame(cm$table)

plotTable <- cm %>%
  mutate(goodbad = ifelse(cm$Prediction == cm$Reference, "correct", "incorrect")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq / sum(Freq))
# plotTable

g_cm <- ggplot(plotTable, aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", prop)), vjust = .5, fontface = "bold", alpha = 1) +
  scale_fill_manual(values = c(correct = "#FF4B00", incorrect = "#5290FF")) +
  xlim(rev(levels(cm$Reference))) +
  theme_bw(base_size = 15, base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
ggsave(sprintf("%s/confusion_matrix.jpg", output_dir), g_cm, width = 8, height = 6)

# Cosine similarity...
g_cos <- ggplot(df_all, aes(x = barcode, y = cos_similarity)) +
  geom_boxplot(aes(color = label)) +
  geom_hline(yintercept = cos_threshold, linetype = "dashed") +
  expand_limits(y = 0) +
  labs(x = "Sample ID", y = "Cosine similarity") +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(as.factor(df_all$barcode)))) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )
ggsave(sprintf("%s/cosine_similarity.jpg", output_dir), g_cos, width = 10, height = 10)
