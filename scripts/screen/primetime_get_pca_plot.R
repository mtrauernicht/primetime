#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

option_list <- list(
  make_option(c("--results"), type = "character"),
  make_option(c("--results-file"), type = "character", dest = "results_file"),
  make_option(c("--output"), type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$output) || (is.null(opt$results) && is.null(opt$results_file))) {
  stop("Missing required arguments: results and output")
}

result_files <- character(0)
if (!is.null(opt$results) && nzchar(opt$results)) {
  result_files <- c(result_files, unlist(strsplit(opt$results, ",")))
}
if (!is.null(opt$results_file) && nzchar(opt$results_file)) {
  if (!file.exists(opt$results_file)) stop("Results list file not found: ", opt$results_file)
  result_files <- c(result_files, unlist(strsplit(paste(readLines(opt$results_file, warn = FALSE), collapse = "\n"), "[[:space:]]+")))
}
result_files <- unique(trimws(result_files))
result_files <- result_files[result_files != ""]

empty_plot <- function(title, message_line) {
  ggplot() +
    annotate("text", x = 0, y = 0.2, label = title, fontface = "bold", size = 5) +
    annotate("text", x = 0, y = -0.2, label = message_line, size = 4) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
}

comparison_parts <- function(path) {
  pieces <- strsplit(sub("\\.txt$", "", basename(path)), "_vs_", fixed = TRUE)[[1]]
  if (length(pieces) != 2 || any(pieces == "")) return(NULL)
  list(contrast = pieces[1], reference = pieces[2])
}

select_spread_labels <- function(coords, x_col, y_col, label_col) {
  coords$display_label <- ""
  if (nrow(coords) == 0) return(coords)
  distances <- sqrt(
    (coords[[x_col]] - mean(coords[[x_col]], na.rm = TRUE))^2 +
    (coords[[y_col]] - mean(coords[[y_col]], na.rm = TRUE))^2
  )
  n_labels <- min(nrow(coords), max(10, ceiling(nrow(coords) * 0.6)), 60)
  top_indices <- order(distances, decreasing = TRUE)[seq_len(n_labels)]
  coords$display_label[top_indices] <- as.character(coords[[label_col]][top_indices])
  coords
}

plot_scores <- function(scores, x_col, y_col, title, x_label, y_label) {
  scores <- select_spread_labels(scores, x_col, y_col, "label")
  scores$status <- ifelse(grepl("control|ctrl|dmso", tolower(scores$label)), "Control", "Sample")
  ggplot(scores, aes_string(x = x_col, y = y_col, color = "status")) +
    geom_point(size = 3.5) +
    geom_text_repel(
      data = subset(scores, display_label != ""),
      aes(label = display_label),
      size = 2.8,
      max.overlaps = Inf,
      force = 1,
      box.padding = 0.3,
      point.padding = 0.1,
      segment.alpha = 0.5,
      min.segment.length = 0
    ) +
    scale_color_manual(values = c("Control" = "#f37f80", "Sample" = "#6495ed"), name = NULL) +
    labs(title = title, x = x_label, y = y_label) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

build_matrix <- function(rows) {
  if (length(rows) == 0) return(NULL)
  combined <- rbindlist(rows, fill = TRUE)
  combined <- combined[is.finite(value)]
  if (nrow(combined) == 0) return(NULL)
  if ("reference_condition" %in% names(combined)) {
    reference_counts <- unique(combined[, .(condition, reference_condition)])[, .(n_refs = .N), by = condition]
    combined <- reference_counts[combined, on = "condition"]
    combined[, condition := ifelse(n_refs > 1, paste0(condition, " vs ", reference_condition), condition)]
  }
  combined <- combined[, .(value = mean(value, na.rm = TRUE)), by = .(condition, tf)]
  wide <- dcast(combined, condition ~ tf, value.var = "value", fill = 0)
  if (nrow(wide) < 2 || ncol(wide) < 3) return(NULL)
  matrix_values <- as.matrix(wide[, -1, with = FALSE])
  rownames(matrix_values) <- wide$condition
  variable_columns <- apply(matrix_values, 2, function(values) sd(values, na.rm = TRUE) > 0)
  if (sum(variable_columns) < 2) return(NULL)
  matrix_values[, variable_columns, drop = FALSE]
}

pca_plot <- function(matrix_values, title) {
  if (is.null(matrix_values) || nrow(matrix_values) < 2 || ncol(matrix_values) < 2) {
    return(empty_plot(title, "Insufficient variable data for PCA"))
  }
  pca <- tryCatch(prcomp(matrix_values, scale. = TRUE, center = TRUE), error = function(e) NULL)
  if (is.null(pca) || ncol(pca$x) < 2) return(empty_plot(title, "PCA could not be computed"))
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$label <- rownames(scores)
  variance <- summary(pca)$importance[2, 1:2] * 100
  plot_scores(scores, "PC1", "PC2", title,
              sprintf("PC1 (%.1f%% variance)", variance[1]),
              sprintf("PC2 (%.1f%% variance)", variance[2]))
}

compute_umap <- function(matrix_values) {
  if (is.null(matrix_values) || nrow(matrix_values) < 3 || ncol(matrix_values) < 2) return(NULL)
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else getwd()
  project_root <- dirname(dirname(dirname(normalizePath(script_path))))
  python_bin <- file.path(project_root, ".venv", "bin", "python")
  if (!file.exists(python_bin)) python_bin <- Sys.which("python")
  if (!nzchar(python_bin)) return(NULL)

  input_file <- tempfile(fileext = ".tsv")
  output_file <- tempfile(fileext = ".tsv")
  write.table(data.frame(label = rownames(matrix_values), scale(matrix_values)), input_file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  python_code <- paste(
    "import csv, sys, numpy as np, umap",
    "rows = list(csv.reader(open(sys.argv[1]), delimiter='\\t'))",
    "labels = [row[0] for row in rows[1:]]",
    "matrix = np.array([[float(x) for x in row[1:]] for row in rows[1:]])",
    "coords = umap.UMAP(n_components=2, random_state=42, n_neighbors=max(2, min(5, len(labels)-1)), min_dist=0.3).fit_transform(matrix)",
    "writer = csv.writer(open(sys.argv[2], 'w', newline=''), delimiter='\\t')",
    "writer.writerow(['label', 'UMAP1', 'UMAP2'])",
    "writer.writerows([[label, x, y] for label, (x, y) in zip(labels, coords)])",
    sep = "\n"
  )
  python_file <- tempfile(fileext = ".py")
  writeLines(python_code, python_file)
  suppressWarnings(system2(python_bin, c(python_file, input_file, output_file), stdout = TRUE, stderr = TRUE))
  if (!file.exists(output_file) || file.info(output_file)$size == 0) return(NULL)
  tryCatch(read.delim(output_file, check.names = FALSE), error = function(e) NULL)
}

activity_rows <- list()
for (path in result_files) {
  if (!file.exists(path) || file.info(path)$size == 0) next
  parts <- comparison_parts(path)
  if (is.null(parts)) next
  data <- tryCatch(fread(path, header = TRUE, na.strings = c("", "NA")), error = function(e) NULL)
  if (is.null(data) || !"tf" %in% names(data)) next
  data <- data[!grepl("^RANDOM", tf)]
  for (condition in c(parts$reference, parts$contrast)) {
    if (condition %in% names(data)) {
      activity_rows[[length(activity_rows) + 1]] <- data.table(tf = data$tf, condition = condition, value = as.numeric(data[[condition]]))
    }
  }
}

activity_matrix <- build_matrix(activity_rows)

barcode_activity_path <- file.path(dirname(dirname(normalizePath(opt$output, mustWork = FALSE))), "tmp_primetime", "activity", "barcode_activity.txt")
barcode_matrix <- NULL
if (file.exists(barcode_activity_path) && file.info(barcode_activity_path)$size > 0) {
  barcode_data <- tryCatch(fread(barcode_activity_path), error = function(e) NULL)
  if (!is.null(barcode_data) && all(c("cDNA_sample", "tf") %in% names(barcode_data))) {
    value_column <- if ("log2_mean_RPM" %in% names(barcode_data)) "log2_mean_RPM" else if ("mean_RPM" %in% names(barcode_data)) "mean_RPM" else NULL
    if (!is.null(value_column)) {
      barcode_data <- barcode_data[!grepl("^RANDOM", tf)]
      barcode_data[, value := as.numeric(get(value_column))]
      if (value_column == "mean_RPM") barcode_data[, value := log2(value + 1)]
      barcode_rows <- list(data.table(tf = barcode_data$tf, condition = barcode_data$cDNA_sample, value = barcode_data$value))
      barcode_matrix <- build_matrix(barcode_rows)
    }
  }
}

umap_coords <- compute_umap(barcode_matrix)
umap_panel <- if (is.null(umap_coords)) {
  empty_plot("UMAP of Individual Sample Activities", "UMAP requires at least three samples and the Python umap package")
} else {
  names(umap_coords)[names(umap_coords) == "label"] <- "label"
  plot_scores(umap_coords, "UMAP1", "UMAP2", "UMAP of Individual Sample Activities", "UMAP 1", "UMAP 2")
}

plots <- list(
  umap_panel,
  pca_plot(barcode_matrix, "PCA of Individual Sample Activities"),
  pca_plot(activity_matrix, "PCA of Primetime-Computed Activities")
)

pdf(opt$output, width = 8, height = 6)
for (plot in plots) print(plot)
dev.off()
