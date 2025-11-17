
suppressPackageStartupMessages({
  library(ggplot2)
})

# CONFIG

ROOT_IN  <- "C:/Users/markg/CLionProjects/DisorderMetrics/src/experiment_data_input"
ROOT_OUT <- "C:/Users/markg/CLionProjects/DisorderMetrics/src/plots/coverage_plot"

IMG_WIDTH_PX  <- 1800
IMG_HEIGHT_PX <- 600
IMG_DPI       <- 200

GREEN <- "#39FF14"
RED   <- "red"


# HELPERS

norm_path <- function(p) gsub("\\\\", "/", normalizePath(p, winslash = "/", mustWork = FALSE))
dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

list_leaf_dirs <- function() {
  all <- list.files(ROOT_IN, pattern = "^sample_indices\\.csv$", recursive = TRUE, full.names = TRUE)
  unique(norm_path(dirname(all)))
}

read_sample_indices <- function(csv_path) {
  message("[read] ", csv_path)
  df <- read.csv(csv_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  num_cols <- sapply(df, is.numeric)
  if (!any(num_cols)) return(NULL)
  as.matrix(df[, num_cols, drop = FALSE])
}

# pick 3 random rows
pick_rows <- function(n_rows) {
  if (n_rows <= 0) return(integer(0))
  sample.int(n_rows, size = min(3L, n_rows), replace = FALSE)
}


plot_coverage_single_sample <- function(idx_vec, n, title_suffix) {
  idx_vec <- as.integer(idx_vec[is.finite(idx_vec)])
  idx_vec <- unique(idx_vec[idx_vec >= 1L & idx_vec <= n])
  
  cov <- integer(n)
  if (length(idx_vec)) cov[idx_vec] <- 1L
  
  df <- data.frame(index = seq_len(n), coverage = cov, check.names = FALSE)
  ggplot(df, aes(x = index, y = coverage)) +
    geom_col(fill = GREEN, width = 1) +
    labs(
      title = paste0("Sample Index Coverage — ", title_suffix),
      x = "Index",
      y = "Hit (0/1)"
    ) +
    theme_minimal() +
    theme(
      text        = element_text(color = RED),   
      plot.title  = element_text(hjust = 0.5, face = "bold"),
      axis.title  = element_text(),
      axis.text   = element_text()
    ) +
    ylim(0, 1)
}

# MAIN

leaf_dirs <- list_leaf_dirs()
message("Found sample leaves (INPUT): ", length(leaf_dirs))

for (d in leaf_dirs) {
  csv <- file.path(d, "sample_indices.csv")
  mat <- read_sample_indices(csv)
  if (is.null(mat)) { message("  skip (no numeric indices): ", csv); next }
  
  max_idx <- suppressWarnings(max(mat, na.rm = TRUE))
  if (!is.finite(max_idx) || max_idx < 1) { message("  skip (invalid indices): ", csv); next }
  n <- as.integer(max_idx)
  
  rel <- substring(norm_path(d), nchar(norm_path(ROOT_IN)) + 2L)
  out_dir <- norm_path(file.path(ROOT_OUT, rel))
  dir_create(out_dir)
  
  rows <- pick_rows(nrow(mat))
  if (length(rows) == 0) next
  message("  leaf: ", rel, "  -> rows picked: ", paste(rows, collapse = ", "))
  
  for (j in seq_along(rows)) {
    r <- rows[j]
    p <- plot_coverage_single_sample(mat[r, ], n, title_suffix = paste("row", r))
    out_png <- file.path(out_dir, paste0("coverage_plot_", j, ".png"))
    ggsave(out_png, p,
           width = IMG_WIDTH_PX/IMG_DPI, height = IMG_HEIGHT_PX/IMG_DPI,
           dpi = IMG_DPI, units = "in")
    message("    saved: ", out_png)
  }
}

message(">>> DONE (coverage plots: three random charts per leaf, red text style).")








