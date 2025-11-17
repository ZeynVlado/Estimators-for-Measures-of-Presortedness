

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# CONFIG

ROOT_IN  <- "C:/Users/markg/CLionProjects/DisorderMetrics/src/experiment_data_output"
ROOT_OUT <- "C:/Users/markg/CLionProjects/DisorderMetrics/src/plots/scatter_plot"

DIR_PERFAMILY        <- file.path(ROOT_OUT, "perFamily")
DIR_RANDOM_COMBINED  <- file.path(ROOT_OUT, "random_combined")
DIR_RUN_COMBINED     <- file.path(ROOT_OUT, "run_array_combined")
DIR_ALL_FAMILIES     <- file.path(ROOT_OUT, "all_families")

METRIC_COLS <- c(
  "inv_norm","runs_norm","rem_norm","osc_norm","dis_norm","ham_norm",
  "run_entr_norm","val_freq_entr_norm","distinct_norm","distinct_entr_norm"
)

METRIC_MAP  <- c(
  inv_norm            = "inv",
  runs_norm           = "runs",
  rem_norm            = "rem",
  osc_norm            = "osc",
  dis_norm            = "dis",
  ham_norm            = "ham",
  run_entr_norm       = "run_entr",
  val_freq_entr_norm  = "valfreq_entr",
  distinct_norm       = "distinct",
  distinct_entr_norm  = "distinct_entr"
)

INV_MAP <- setNames(names(METRIC_MAP), METRIC_MAP)

POINT_ALPHA <- 0.7
POINT_SIZE  <- 1.6
IMG_WIDTH   <- 1800  # px
IMG_HEIGHT  <- 1200  # px
IMG_DPI     <- 200

GREEN <- "#39FF14"
RED   <- "red"


# HELPERS

norm_path <- function(p) gsub("\\\\", "/", normalizePath(p, winslash = "/", mustWork = FALSE))
dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

find_arrays_csv_up <- function(start_dir) {
  cur <- norm_path(start_dir)
  repeat {
    cand <- file.path(cur, "arrays_metrics.csv")
    if (file.exists(cand)) return(cand)
    parent <- norm_path(dirname(cur))
    if (identical(parent, cur)) break
    cur <- parent
  }
  NA_character_
}

read_arrays_csv <- function(csv_path) {
  message("  [read_arrays_csv] ", csv_path)
  df <- read.csv(csv_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  expected <- c("n", METRIC_COLS)
  miss <- setdiff(expected, names(df))
  if (length(miss)) stop("arrays_metrics.csv missing: ", paste(miss, collapse=", "), "\nFile: ", csv_path)
  df[expected]
}

read_sample_csv <- function(csv_path) {
  message("  [read_sample_csv] ", csv_path)
  df <- read.csv(csv_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  expected <- c("n", METRIC_COLS)
  miss <- setdiff(expected, names(df))
  if (length(miss)) stop("sample_metrics.csv missing: ", paste(miss, collapse=", "), "\nFile: ", csv_path)
  df[expected]
}

list_leaf_sample_dirs <- function() {
  all <- list.files(ROOT_IN, pattern = "^sample_metrics\\.csv$", recursive = TRUE, full.names = TRUE)
  unique(norm_path(dirname(all)))
}

rewrite_rel_for_perFamily <- function(rel) {
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  if (length(parts) >= 3 && parts[1] == "run_array" && grepl("^[0-9]+$", parts[2])) {
    parts <- parts[-2]  # удалить N
  }
  paste(parts, collapse = "/")
}

collect_pairs_for_leaf <- function(sample_dir) {
  sample_csv <- file.path(sample_dir, "sample_metrics.csv")
  if (!file.exists(sample_csv)) {
    message("  [collect_pairs_for_leaf] sample_metrics.csv NOT FOUND: ", sample_dir)
    return(NULL)
  }
  arrays_csv <- find_arrays_csv_up(sample_dir)
  if (is.na(arrays_csv) || !file.exists(arrays_csv)) {
    message("  [collect_pairs_for_leaf] arrays_metrics.csv NOT FOUND for: ", sample_dir)
    return(NULL)
  }
  
  df_s <- read_sample_csv(sample_csv)
  df_a <- read_arrays_csv(arrays_csv)
  
  n_pairs <- min(nrow(df_s), nrow(df_a))
  if (n_pairs == 0) {
    message("  [collect_pairs_for_leaf] empty arrays/samples: ", sample_dir)
    return(NULL)
  }
  
  df_s <- df_s[seq_len(n_pairs), , drop = FALSE]
  df_a <- df_a[seq_len(n_pairs), , drop = FALSE]
  
  s_long <- df_s %>%
    mutate(.row = seq_len(n())) %>%
    pivot_longer(cols = all_of(METRIC_COLS), names_to = "metric_col", values_to = "sample_value")
  a_long <- df_a %>%
    mutate(.row = seq_len(n())) %>%
    pivot_longer(cols = all_of(METRIC_COLS), names_to = "metric_col", values_to = "arr_value")
  
  pairs <- s_long %>%
    inner_join(a_long, by = c(".row","metric_col")) %>%
    select(.row, metric_col, arr_value, sample_value)
  
  rel <- substring(norm_path(sample_dir), nchar(norm_path(ROOT_IN)) + 2L)
  rel_segs <- strsplit(rel, "/", fixed = TRUE)[[1]]
  family <- if (length(rel_segs)) rel_segs[1] else NA_character_
  
  pairs$n <- df_a$n[pairs$.row]
  # S/design/variant — из пути: .../samples/<S>/<design>/<variant>/
  segs <- strsplit(norm_path(sample_dir), "/")[[1]]
  idx  <- which(tolower(segs) == "samples")
  ssize   <- if (length(idx) && idx + 1 <= length(segs)) segs[idx + 1] else NA_character_
  design  <- if (length(idx) && idx + 2 <= length(segs)) segs[idx + 2] else NA_character_
  variant <- if (length(idx) && idx + 3 <= length(segs)) segs[idx + 3] else NA_character_
  
  pairs$ssize   <- ssize
  pairs$design  <- design
  pairs$variant <- variant
  pairs$family  <- family
  pairs$metric  <- unname(METRIC_MAP[pairs$metric_col])  # короткое имя метрики
  
  if (!is.na(family) && family == "run_array") {
    method <- if (length(rel_segs) >= 3 && grepl("^[0-9]+$", rel_segs[2])) rel_segs[3] else NA_character_
    pairs$run_method <- method
  } else {
    pairs$run_method <- NA_character_
  }
  
  dbg <- pairs %>% count(metric, name = "rows")
  if (nrow(dbg)) {
    msg <- paste(apply(dbg, 1, function(r) paste0(r[1], "=", r[2])), collapse = ", ")
    message("  [collect_pairs_for_leaf] metrics rows: ", msg)
  }
  
  pairs
}

compute_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  mx <- mean(x); my <- mean(y)
  vx <- var(x);  vy <- var(y)
  sx <- sqrt(vx); sy <- sqrt(vy)
  if (sx == 0 || sy == 0) return(NA_real_)
  rho <- suppressWarnings(cor(x, y, method = "pearson"))
  if (!is.finite(rho)) return(NA_real_)
  (2 * rho * sx * sy) / (vx + vy + (mx - my)^2)
}

compute_stats <- function(metric_short, x, y) {
  mean_arr <- mean(x); var_arr <- var(x); med_arr <- median(x); sd_arr <- sd(x)
  mean_sam <- mean(y); var_sam <- var(y); med_sam <- median(y); sd_sam <- sd(y)
  
  pearson  <- if (sd_arr > 0 && sd_sam > 0) suppressWarnings(cor(x, y, method="pearson"))  else NA_real_
  spearman <- if (sd_arr > 0 && sd_sam > 0) suppressWarnings(cor(x, y, method="spearman")) else NA_real_
  ccc_val  <- compute_ccc(x, y)
  sd_ratio <- if (sd_arr > 0) sd_sam / sd_arr else NA_real_
  bias     <- mean_sam - mean_arr
  
  data.frame(
    metric = metric_short,
    mean_arr = mean_arr,
    variance_arr = var_arr,
    median_arr = med_arr,
    standard_deviation_arr = sd_arr,
    mean_sample = mean_sam,
    variance_sample = var_sam,
    median_sample = med_sam,
    standard_deviation_sample = sd_sam,
    pearson_corr = pearson,
    spearman_corr = spearman,
    ccc = ccc_val,
    sd_ratio = sd_ratio,
    bias = bias,
    check.names = FALSE
  )
}

create_scatter_plot_xy <- function(x, y, metric_short, output_file) {
  plot_data <- data.frame(array_metric = x, sample_metric = y)
  plot_data <- na.omit(plot_data)
  message("    [plot] metric=", metric_short, "  n=", nrow(plot_data))
  if (nrow(plot_data) == 0) { message("    -> skip (no data): ", output_file); return(invisible(NULL)) }
  
  p <- ggplot(plot_data, aes(x = array_metric, y = sample_metric)) +
    geom_point(alpha = POINT_ALPHA, size = POINT_SIZE, color = GREEN) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = RED) +
    labs(
      title = paste("Scatter Plot:", toupper(metric_short), "Metric"),
      x = paste("Array", toupper(metric_short), "(normalized)"),
      y = paste("Sample", toupper(metric_short), "(normalized)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = RED),
      axis.title = element_text(color = "red"),
      axis.text  = element_text(color = "red")
    ) +
    xlim(0, 1) + ylim(0, 1) +
    coord_fixed(ratio = 1)
  
  dir_create(dirname(output_file))
  message("    [ggsave] -> ", norm_path(output_file))
  ggsave(output_file, p,
         width = IMG_WIDTH/IMG_DPI, height = IMG_HEIGHT/IMG_DPI,
         dpi = IMG_DPI, units = "in")
  message("    [saved PNG] ", norm_path(output_file))
}

write_stats_csv <- function(stats_df, out_csv) {
  dir_create(dirname(out_csv))
  header <- paste0('"', names(stats_df), '"', collapse = ",")
  vals <- mapply(function(v, isn) {
    if (isn) formatC(v, digits = 16, format = "g") else paste0('"', v, '"')
  }, stats_df[1, ], sapply(stats_df[1, ], is.numeric), SIMPLIFY = TRUE, USE.NAMES = FALSE)
  line <- paste(vals, collapse = ",")
  con <- file(out_csv, open = "wt", encoding = "UTF-8")
  writeLines(header, con); writeLines(line, con); close(con)
  message("    [saved CSV] ", norm_path(out_csv))
}

safe_bind_rows <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) return(NULL)
  do.call(rbind, lst)
}

# PER-FAMILY 

message(">>> Step 1/4: perFamily ...")
leaf_dirs <- list_leaf_sample_dirs()
message("  [perFamily] leaves found: ", length(leaf_dirs))

png_cnt <- 0L; csv_cnt <- 0L
t1 <- proc.time()

for (d in leaf_dirs) {
  message("  [leaf] ", d)
  pairs <- collect_pairs_for_leaf(d)
  if (is.null(pairs) || nrow(pairs) == 0) { message("   -> skip (no pairs)"); next }
  
  rel <- substring(norm_path(d), nchar(norm_path(ROOT_IN)) + 2L)
  rel_fixed <- rewrite_rel_for_perFamily(rel)
  out_dir  <- file.path(DIR_PERFAMILY, rel_fixed)
  
  for (metric_s in unique(pairs$metric)) {
    dfm <- pairs %>% filter(metric == !!metric_s)
    if (nrow(dfm) == 0) next
    
    out_png <- file.path(out_dir, paste0("scatter_plot)", metric_s, ".png"))
    out_csv <- file.path(out_dir, paste0(metric_s, ".csv"))
    
    message("   -> metric=", metric_s,
            "  rows=", nrow(dfm),
            "  out_dir=", norm_path(out_dir))
    x <- dfm$arr_value
    y <- dfm$sample_value
    
    st <- compute_stats(metric_s, x, y)
    create_scatter_plot_xy(x, y, metric_s, out_png)
    write_stats_csv(st, out_csv)
    png_cnt <- png_cnt + 1L; csv_cnt <- csv_cnt + 1L
  }
}
elapsed1 <- (proc.time() - t1)[["elapsed"]]
message("<<< perFamily DONE. PNG=", png_cnt, " CSV=", csv_cnt, " time=", round(elapsed1, 2), "s")


# RANDOM-COMBINED (perm + random + multiset)

message(">>> Step 2/4: random_combined ...")
leaf_perm  <- leaf_dirs[grepl("(^|/)permutation(/|$)",          leaf_dirs)]
leaf_rand  <- leaf_dirs[grepl("(^|/)random_array(/|$)",         leaf_dirs)]
leaf_multi <- leaf_dirs[grepl("(^|/)multiset_permutation(/|$)", leaf_dirs)]
message("  [random_combined] leaves: perm=", length(leaf_perm),
        " random=", length(leaf_rand), " multiset=", length(leaf_multi))

pairs_perm  <- safe_bind_rows(lapply(leaf_perm,  collect_pairs_for_leaf))
pairs_rand  <- safe_bind_rows(lapply(leaf_rand,  collect_pairs_for_leaf))
pairs_multi <- safe_bind_rows(lapply(leaf_multi, collect_pairs_for_leaf))

pairs_random_combined <- safe_bind_rows(list(pairs_perm, pairs_rand, pairs_multi))

if (!is.null(pairs_random_combined) && nrow(pairs_random_combined) > 0) {
  g <- pairs_random_combined %>%
    filter(metric %in% unname(METRIC_MAP)) %>%
    group_by(n, ssize, design, variant, metric) %>%   # <-- важная правка
    group_split()
  
  message("  [random_combined] groups: ", length(g))
  
  for (dfm in g) {
    nval    <- dfm$n[1]
    ssize   <- dfm$ssize[1]
    design  <- dfm$design[1]
    variant <- dfm$variant[1]
    metric_s <- dfm$metric[1]
    
    out_dir <- file.path(DIR_RANDOM_COMBINED,
                         as.character(nval), "samples", ssize, design, variant)
    out_png <- file.path(out_dir, paste0("scatter_plot)", metric_s, ".png"))
    out_csv <- file.path(out_dir, paste0(metric_s, ".csv"))
    
    message("   -> N=", nval, "  S=", ssize,
            "  design=", design, "  variant=", variant,
            "  metric=", metric_s,
            "  rows=", nrow(dfm),
            "  out_dir=", norm_path(out_dir))
    
    x <- dfm$arr_value
    y <- dfm$sample_value
    st <- compute_stats(metric_s, x, y)
    create_scatter_plot_xy(x, y, metric_s, out_png)
    write_stats_csv(st, out_csv)
  }
} else {
  message("  [random_combined] no data.")
}
message("<<< random_combined DONE.")

# RUN-ARRAY-COMBINED 

message(">>> Step 3/4: run_array_combined ...")
leaf_run <- leaf_dirs[grepl("(^|/)run_array(/|$)", leaf_dirs)]
message("  [run_array_combined] leaves(run_array): ", length(leaf_run))
pairs_run_all <- safe_bind_rows(lapply(leaf_run, collect_pairs_for_leaf))

if (!is.null(pairs_run_all) && nrow(pairs_run_all) > 0) {
  g <- pairs_run_all %>%
    filter(metric %in% unname(METRIC_MAP)) %>%
    group_by(n, ssize, design, variant, metric) %>%   # <-- важная правка
    group_split()
  
  message("  [run_array_combined] groups: ", length(g))
  
  for (dfm in g) {
    nval    <- dfm$n[1]
    ssize   <- dfm$ssize[1]
    design  <- dfm$design[1]
    variant <- dfm$variant[1]
    metric_s <- dfm$metric[1]
    
    out_dir <- file.path(DIR_RUN_COMBINED,
                         as.character(nval), "r=all", "samples", ssize, design, variant)
    out_png <- file.path(out_dir, paste0("scatter_plot)", metric_s, ".png"))
    out_csv <- file.path(out_dir, paste0(metric_s, ".csv"))
    
    message("   -> N=", nval, "  S=", ssize,
            "  design=", design, "  variant=", variant,
            "  metric=", metric_s,
            "  rows=", nrow(dfm),
            "  out_dir=", norm_path(out_dir))
    
    x <- dfm$arr_value
    y <- dfm$sample_value
    st <- compute_stats(metric_s, x, y)
    create_scatter_plot_xy(x, y, metric_s, out_png)
    write_stats_csv(st, out_csv)
  }
} else {
  message("  [run_array_combined] no data.")
}
message("<<< run_array_combined DONE.")

# ALL-FAMILIES 

message(">>> Step 4/4: all_families ...")
pairs_all <- safe_bind_rows(list(pairs_perm, pairs_rand, pairs_multi, pairs_run_all))

if (!is.null(pairs_all) && nrow(pairs_all) > 0) {
  g <- pairs_all %>%
    filter(metric %in% unname(METRIC_MAP)) %>%
    group_by(n, ssize, design, variant, metric) %>%    
    group_split()
  
  message("  [all_families] groups: ", length(g))
  
  for (dfm in g) {
    nval    <- dfm$n[1]
    ssize   <- dfm$ssize[1]
    design  <- dfm$design[1]
    variant <- dfm$variant[1]
    metric_s <- dfm$metric[1]
    
    out_dir <- file.path(DIR_ALL_FAMILIES,
                         as.character(nval), "samples", ssize, design, variant)
    out_png <- file.path(out_dir, paste0("scatter_plot)", metric_s, ".png"))
    out_csv <- file.path(out_dir, paste0(metric_s, ".csv"))
    
    message("   -> N=", nval, "  S=", ssize,
            "  design=", design, "  variant=", variant,
            "  metric=", metric_s,
            "  rows=", nrow(dfm),
            "  out_dir=", norm_path(out_dir))
    
    x <- dfm$arr_value
    y <- dfm$sample_value
    st <- compute_stats(metric_s, x, y)
    create_scatter_plot_xy(x, y, metric_s, out_png)
    write_stats_csv(st, out_csv)
  }
} else {
  message("  [all_families] no data.")
}
message("<<< all_families DONE.")
message(">>> ALL DONE.")






