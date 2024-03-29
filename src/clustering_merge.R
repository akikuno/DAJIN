################################################################################
# Install required packages
################################################################################

options(repos = "https://cloud.r-project.org/")
options(readr.show_progress = FALSE)
options(dplyr.summarise.inform = FALSE)
options(future.globals.maxSize = Inf)
options(warn = -1)

if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(readr, stringr, tibble, dplyr, tidyr, purrr, furrr)

################################################################################
# I/O naming
################################################################################

# ===========================================================
# Auguments
# ===========================================================

args <- commandArgs(trailingOnly = TRUE)
query_score <- args[1]
query_label <- args[2]
control_RDS <- args[3]
threads <- as.integer(args[4])
plan(multicore, workers = threads)

# ===========================================================
# Inputs
# ===========================================================

df_que_mids <- read_csv(query_score,
  col_names = FALSE,
  col_types = cols(),
  num_threads = threads
)
colnames(df_que_mids) <- seq_len(ncol(df_que_mids))

df_que_label <- read_csv(query_label,
  col_names = c("id", "strand", "barcode"),
  col_types = cols()
)

df_control_score <- readRDS(control_RDS)

# ===========================================================
# Outputs
# ===========================================================

output_suffix <-
  str_remove(query_label, ".*labels_")

int_hdbscan_clusters <-
  read_csv(
    sprintf(".DAJIN_temp/clustering/temp/int_hdbscan_clusters_%s", output_suffix),
    col_names = c("cl"),
    col_types = cols()
  ) %>%
  pull(cl)

df_score <- readRDS(sprintf(".DAJIN_temp/clustering/temp/df_score_%s.RDS", output_suffix))

################################################################################
# Extract mutations by Hotelling's T2 statistics
################################################################################

merged_clusters <- int_hdbscan_clusters

# ===========================================================
# Replace sequence error
# ===========================================================

sequence_error <-
  df_control_score %>%
  filter(mut == 0) %>%
  unnest(control_freq) %>%
  filter(MIDS != "M" & freq > 10) %>%
  pull(loc) %>%
  unique()

df_score[, sequence_error] <- 10^-10000

# ===========================================================
# Extract mutations at each clusters
# ===========================================================

hotelling <-
  future_map_dfr(unique(merged_clusters), function(x) {
    tmp_scale <-
      df_score[merged_clusters == x, ] %>%
      colSums() %>%
      scale() %>%
      .[c(1:100, tail(seq(ncol(df_score)), 100))]
    tmp_mean <- mean(tmp_scale)
    tmp_sd <- sd(tmp_scale)
    add_rnorm <-
      rnorm(ncol(df_que_mids) * 100, mean = tmp_mean, sd = tmp_sd) %>%
      as_tibble() %>%
      rename(score = value)

    df_score[merged_clusters == x, ] %>%
      colSums() %>%
      scale() %>%
      as_tibble() %>%
      rename(score = colnames(.)) %>%
      bind_rows(add_rnorm) %>%
      summarize(
        score = score,
        mean = mean(score),
        var = mean((score - mean(score))^2)
      ) %>%
      summarize(anomaly_score = (score - mean)^2 / var) %>%
      mutate(loc = row_number(), threshold = qchisq(0.95, 1)) %>%
      filter(anomaly_score > threshold & loc <= ncol(df_que_mids)) %>%
      select(loc) %>%
      mutate(cl = x)
  })

# ===========================================================
# Exclude fuzzy mutations
# ===========================================================

tmp_possible_true_mut <- hotelling$loc %>%
  sort() %>%
  unique()

lgl_possible_true_mut <-
  future_map_lgl(tmp_possible_true_mut, function(loc) {
    tmp_ <-
      map_lgl(unique(merged_clusters), function(cl) {
        df_que_mids[merged_clusters == cl, loc] %>%
          table(dnn = "MIDS") %>%
          as_tibble() %>%
          mutate(freq = n / sum(n) * 100) %>%
          slice_max(freq, n = 1) %>%
          head(n = 1) %>%
          mutate(lgl = if_else(MIDS != "M" & freq > 50, TRUE, FALSE)) %>%
          pull(lgl)
      })
    if_else(any(tmp_), TRUE, FALSE)
  })

possible_true_mut <-
  tmp_possible_true_mut[lgl_possible_true_mut]

# Force to add a point mutation location
if (sum(df_control_score$mut) == 1) {
  possible_true_mut <-
    append(possible_true_mut, which(df_control_score$mut == 1)) %>%
    unique()
}

# ===========================================================
# Divide mutation into common and uncommon
# ===========================================================

hotelling_common <-
  hotelling %>%
  filter(loc %in% possible_true_mut) %>%
  add_count(loc, name = "count_loc") %>%
  mutate(common = if_else(
    count_loc == length(unique(merged_clusters)),
    TRUE,
    FALSE
  )) %>%
  select(loc, common) %>%
  distinct()

################################################################################
# Merge clusters
################################################################################

# ===========================================================
# Merge clusters with the same mutation profiles
# ===========================================================

if (length(unique(merged_clusters)) > 1 & length(possible_true_mut) > 0) {
  max_mids <-
    future_map_chr(possible_true_mut, function(x) {
      df_que_mids[, x] %>%
        table(dnn = "MIDS") %>%
        as_tibble() %>%
        mutate(freq = n / sum(n) * 100) %>%
        slice_max(freq, n = 1) %>%
        head(n = 1) %>%
        pull(MIDS)
    }) %>%
    set_names(possible_true_mut)

  # common
  tmp_ <-
    hotelling_common %>%
    filter(common == TRUE) %>%
    pull(loc)

  consensus_common <-
    map(unique(merged_clusters), function(cl) {
      future_map(tmp_, function(loc) {
        df_que_mids[merged_clusters == cl, loc] %>%
          table(dnn = "MIDS") %>%
          as_tibble() %>%
          mutate(freq = n / sum(n) * 100) %>%
          slice_max(freq, n = 1) %>%
          head(n = 1) %>%
          mutate(MIDS = if_else(freq > 75, MIDS, max_mids[names(max_mids) == loc])) %>%
          pull(MIDS)
      }) %>%
        unlist() %>%
        str_c(collapse = "")
    }) %>%
    set_names(unique(merged_clusters))


  # uncommon
  tmp_ <-
    hotelling_common %>%
    filter(common == FALSE) %>%
    pull(loc)

  consensus_uncommon <-
    map(unique(merged_clusters), function(cl) {
      future_map(tmp_, function(y) {
        control <-
          df_control_score %>%
          filter(loc == y) %>%
          unnest(control_freq)
        control <-
          if (any(str_detect(colnames(control), "MIDS"))) {
            control %>%
              mutate(MIDS = case_when(mut == 1 ~ "M", TRUE ~ MIDS)) %>%
              mutate(freq = case_when(mut == 1 ~ 100, TRUE ~ freq))
          } else {
            control %>%
              mutate(MIDS = "M", freq = 0) %>%
              bind_rows(tibble(MIDS = c("I", "D", "S"), freq = c(0, 0, 0)))
          }

        df_que_mids[merged_clusters == cl, y] %>%
          table(dnn = "MIDS") %>%
          as_tibble() %>%
          mutate(freq = n / sum(n) * 100) %>%
          full_join(., control, by = "MIDS", suffix = c("_x", "_y")) %>%
          slice_max(freq_x, n = 1) %>%
          filter(MIDS != "M" & freq_y < 5) %>%
          head(n = 1) %>%
          mutate(MIDS = if_else(freq_x > 90, MIDS, NULL)) %>%
          pull(MIDS)
      }) %>%
        unlist() %>%
        .[!is.na(.)] %>%
        str_c(collapse = "")
    }) %>%
    set_names(unique(merged_clusters))

  retain_seq_consensus <-
    map2(consensus_common, consensus_uncommon, function(x, y) {
      str_c(x, y)
    })

  cl_combn <- combn(unique(merged_clusters), 2)

  df_consensus <-
    seq(ncol(cl_combn)) %>%
    map_dfr(function(x) {
      tmp_seq1 <- retain_seq_consensus[[as.character(cl_combn[1, x])]]
      tmp_seq2 <- retain_seq_consensus[[as.character(cl_combn[2, x])]]

      tibble(
        one = cl_combn[1, x],
        two = cl_combn[2, x],
        score = identical(tmp_seq1, tmp_seq2)
      )
    })
  df_consensus_extracted <- df_consensus %>% filter(score == TRUE)

  if (nrow(df_consensus_extracted) != 0) {
    for (i in seq_along(rownames(df_consensus_extracted))) {
      pattern_ <- df_consensus_extracted[i, ]$one
      query_ <- df_consensus_extracted[i, ]$two
      merged_clusters[merged_clusters == pattern_] <- query_
    }
  }
} else {
  merged_clusters <-
    rep(1, length(merged_clusters))
}

# ===========================================================
# Merge clusters with Cosine similarity
# ===========================================================

if (length(unique(merged_clusters)) > 1 & length(possible_true_mut) > 1) {
  calc_cosine_sim <- function(x, y) {
    crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
  }

  cl_combn <- combn(unique(merged_clusters), 2)

  df_cossim <-
    future_map_dfr(seq(ncol(cl_combn)), function(x) {
      score_1 <-
        df_score[merged_clusters == cl_combn[1, x], ] %>%
        colSums()
      score_2 <-
        df_score[merged_clusters == cl_combn[2, x], ] %>%
        colSums()

      tibble(
        one = cl_combn[1, x],
        two = cl_combn[2, x],
        score = calc_cosine_sim(score_1, score_2)
      )
    })

  cossim_thres <- if_else(str_detect(query_score, "abnormal"), 0.95, 0.80)
  df_cossim_extracted <- df_cossim %>% filter(score > cossim_thres)

  if (nrow(df_cossim_extracted) != 0) {
    for (i in seq_along(rownames(df_cossim_extracted))) {
      pattern_ <- df_cossim_extracted[i, ]$one
      query_ <- df_cossim_extracted[i, ]$two
      merged_clusters[merged_clusters == pattern_] <- query_
    }
    for (i in seq_along(rownames(df_cossim_extracted))) {
      pattern_ <- df_cossim_extracted[i, ]$two
      query_ <- df_cossim_extracted[i, ]$one
      merged_clusters[merged_clusters == pattern_] <- query_
    }
  }
}

# ===========================================================
# Merge clusters with strand specific mutation into a major cluster
# ===========================================================

tmp_tb_strand <-
  tibble(cl = merged_clusters, strand = df_que_label$strand)

tmp_dual_adaptor <-
  tmp_tb_strand %>%
  count(strand) %>%
  mutate(freq = n / sum(n)) %>%
  filter(freq < 0.90 & freq > 0.10)

if (nrow(tmp_dual_adaptor) > 0) {
  tmp_biased_cl <-
    tmp_tb_strand %>%
    count(cl, strand) %>%
    mutate(freq = n / sum(n)) %>%
    mutate(bias = if_else(freq > 0.90, TRUE, FALSE)) %>%
    filter(bias == TRUE) %>%
    pull(cl) %>%
    unique()

  tmp_cl_max <-
    tmp_tb_strand %>%
    count(cl) %>%
    filter(!(cl %in% tmp_biased_cl)) %>%
    slice_max(n, n = 1) %>%
    slice_sample(n = 1) %>%
    pull(cl)

  merged_clusters <-
    merged_clusters %>%
    if_else(. %in% tmp_biased_cl, tmp_cl_max, .)
}

################################################################################
# Output results
################################################################################

merged_clusters <- as.factor(merged_clusters) %>% as.integer()

df_readid_cluster <-
  tibble(
    read_id = df_que_label$id,
    cluster = merged_clusters
  )

write_tsv(df_readid_cluster,
  sprintf(".DAJIN_temp/clustering/temp/hdbscan_%s", output_suffix),
  col_names = F
)

write_tsv(tibble(loc = possible_true_mut),
  sprintf(".DAJIN_temp/clustering/temp/possible_true_mut_%s", output_suffix),
  col_names = F
)
