
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# CONFIG

ROOT_IN  <- "C:/Users/markg/CLionProjects/DisorderMetrics/src/experiment_data_input"
ROOT_OUT <- "C:/Users/markg/CLionProjects/DisorderMetrics/src/plots/disrtibution_plot"  

IMG_WIDTH_PANEL  <- 2200  
IMG_HEIGHT_PANEL <- 600 
IMG_DPI          <- 200

IMG_WIDTH_SINGLE  <- 900 
IMG_HEIGHT_SINGLE <- 600  

GREEN <- "#39FF14"
RED   <- "red"


# HELPERS

norm_path <- function(p) gsub("\\\\", "/", normalizePath(p, winslash = "/", mustWork = FALSE))
dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

list_arrays_csv <- function() {
  all <- list.files(ROOT_IN, pattern = "^arrays\\.csv$", recursive = TRUE, full.names = TRUE)
  unique(norm_path(all))
}

row_to_long <- function(df_row) {
  num_cols <- sapply(df_row, is.numeric)
  if (!any(num_cols)) return(data.frame(index = integer(0), value = numeric(0)))
  sub <- df_row[, num_cols, drop = FALSE]
  
  if ("n" %in% names(sub) && ncol(sub) >= 1) {
    sub <- sub[, setdiff(names(sub), "n"), drop = FALSE]
  }
  
  if (ncol(sub) == 0) return(data.frame(index = integer(0), value = numeric(0)))
  vals <- as.numeric(sub[1, , drop = TRUE])
  data.frame(index = seq_along(vals), value = vals, check.names = FALSE)
}

draw_single_barchart <- function(df_long, title_suffix = NULL) {
  ggplot(df_long, aes(x = index, y = value)) +
    geom_col(fill = GREEN, width = 1) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = RED),
      axis.title = element_text(color = "black"),
      axis.text  = element_text(color = "black")
    ) +
    labs(
      title = if (is.null(title_suffix)) "Distribution (array)" else paste("Distribution (array) —", title_suffix),
      x = "Index", y = "Value"
    )
}

save_panel_1x5 <- function(list_plots, out_png) {
  if (length(list_plots) == 0) return(invisible(NULL))
  ok <- requireNamespace("patchwork", quietly = TRUE)
  dir_create(dirname(out_png))
  if (ok) {
    panel <- Reduce(`|`, list_plots)
    ggsave(
      out_png, panel,
      width  = IMG_WIDTH_PANEL / IMG_DPI,
      height = IMG_HEIGHT_PANEL / IMG_DPI,
      dpi    = IMG_DPI, units = "in"
    )
    message("  [saved panel] ", norm_path(out_png))
  } else {
    warning("Package 'patchwork' not installed; panel PNG skipped: ", out_png)
  }
}

pick_rows_random <- function(n_rows, max_take = 5L) {
  if (n_rows <= 0) return(integer(0))
  sample.int(n_rows, size = min(max_take, n_rows), replace = FALSE)
}


# MAIN

arrays_files <- list_arrays_csv()
message("Found arrays.csv files: ", length(arrays_files))

for (csv in arrays_files) {
  rel <- substring(norm_path(csv), nchar(norm_path(ROOT_IN)) + 2L) # path/to/arrays.csv
  out_dir <- norm_path(file.path(ROOT_OUT, dirname(rel)))
  dir_create(out_dir)
  
  message("\n[processing] ", norm_path(csv))
  message("  -> out_dir: ", out_dir)
  
  df <- tryCatch(
    read.csv(csv, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) { message("  read.csv ERROR: ", e$message); NULL }
  )
  if (is.null(df) || nrow(df) == 0) { message("  -> skip: empty or unreadable"); next }
  
  rows <- pick_rows_random(nrow(df), max_take = 5L)
  if (length(rows) == 0) { message("  -> skip: no rows selected"); next }
  
  plots <- vector("list", length(rows))
  
  for (j in seq_along(rows)) {
    r <- rows[j]
    df_long <- row_to_long(df[r, , drop = FALSE])
    if (nrow(df_long) == 0) { message("    row ", r, ": skip (no numeric data)"); next }
    
    p <- draw_single_barchart(df_long, title_suffix = paste("row", r))
    plots[[j]] <- p
    
    single_png <- file.path(out_dir, paste0("distribution_plot_", j, ".png"))
    ggsave(
      single_png, p,
      width  = IMG_WIDTH_SINGLE / IMG_DPI,
      height = IMG_HEIGHT_SINGLE / IMG_DPI,
      dpi    = IMG_DPI, units = "in"
    )
    message("    [saved] ", norm_path(single_png))
  }
  
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) > 0) {
    panel_png <- file.path(out_dir, "distribution_plot_pannel.png")
    save_panel_1x5(plots, panel_png)
  } else {
    message("  -> panel skipped: no valid plots")
  }
}

message("\n>>> DONE (distribution plots, arrays only).")