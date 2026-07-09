#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggthemes)
})

option_list <- list(
  make_option(c("--results"), type = "character", help = "Comma-separated list of comparison result files"),
  make_option(c("--results-file"), type = "character", dest = "results_file", help = "Path to file with newline-separated comparison result files"),
  make_option(c("--bleedthrough"), type = "character", help = "Path to bleedthrough data file"),
  make_option(c("--output"), type = "character", help = "Output PDF path for UMAP plot")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if ((is.null(opt$results) && is.null(opt$results_file)) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else getwd()
project_root <- dirname(dirname(normalizePath(script_path)))
python_bin <- file.path(project_root, ".venv", "bin", "python")
if (!file.exists(python_bin)) {
  python_bin <- Sys.which("python")
}

write_placeholder_pdf <- function(path, message_line) {
  pdf(path, width = 8, height = 4)
  plot.new()
  text(0.5, 0.6, "UMAP of TF Activities Across Conditions", cex = 1.2, font = 2)
  text(0.5, 0.4, message_line, cex = 1)
  invisible(dev.off())
}

result_files <- character(0)
if (!is.null(opt$results)) {
  result_files <- c(result_files, unlist(strsplit(opt$results, ",")))
}
if (!is.null(opt$results_file)) {
  result_files <- c(result_files, trimws(readLines(opt$results_file, warn = FALSE)))
}
result_files <- unique(result_files[result_files != ""])

condition_activity_df <- data.frame()
reference_conditions <- character(0)
primetime_logfc_df <- data.frame()

compute_umap_coords <- function(input_matrix, labels, python_exe) {
  if (nrow(input_matrix) < 3 || ncol(input_matrix) < 2) {
    return(NULL)
  }

  scaled_matrix <- scale(input_matrix)
  scaled_matrix[is.na(scaled_matrix)] <- 0

  input_file <- tempfile(fileext = ".tsv")
  output_file <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(label = labels, as.data.frame(scaled_matrix, check.names = FALSE)),
    file = input_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  python_code <- c(
    "import csv, sys, numpy as np",
    "from pathlib import Path",
    "import umap",
    "input_path = Path(sys.argv[1])",
    "output_path = Path(sys.argv[2])",
    "with input_path.open() as fh:",
    "    reader = csv.reader(fh, delimiter='\\t')",
    "    header = next(reader)",
    "    rows = list(reader)",
    "labels = [row[0] for row in rows]",
    "matrix = np.array([[float(value) for value in row[1:]] for row in rows], dtype=float)",
    "n_obs = matrix.shape[0]",
    "n_neighbors = max(2, min(5, n_obs - 1))",
    "reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=n_neighbors, min_dist=0.3)",
    "coords = reducer.fit_transform(matrix)",
    "with output_path.open('w', newline='') as fh:",
    "    writer = csv.writer(fh, delimiter='\\t')",
    "    writer.writerow(['label', 'UMAP1', 'UMAP2'])",
    "    for label, (x, y) in zip(labels, coords):",
    "        writer.writerow([label, x, y])"
  )
  py_script <- tempfile(fileext = ".py")
  writeLines(python_code, py_script)

  suppressWarnings(system2(python_exe, c(py_script, input_file, output_file), stdout = TRUE, stderr = TRUE))
  if (!file.exists(output_file) || file.info(output_file)$size == 0) {
    return(NULL)
  }

  coords <- suppressWarnings(read.table(output_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  if (nrow(coords) == 0 || !all(c("label", "UMAP1", "UMAP2") %in% colnames(coords))) {
    return(NULL)
  }
  coords
}

for (this_file in result_files) {
  if (!file.exists(this_file) || file.info(this_file)$size == 0) {
    next
  }

  file_base <- sub("\\.txt$", "", basename(this_file))
  comparison_parts <- strsplit(file_base, "_vs_")[[1]]
  if (length(comparison_parts) != 2) {
    next
  }

  ref_condition <- comparison_parts[2]
  contrast_condition <- comparison_parts[1]
  comparison_df <- suppressWarnings(read.table(this_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
  if (!all(c("tf", ref_condition, contrast_condition) %in% colnames(comparison_df))) {
    next
  }

  if (all(c("tf", "logFC") %in% colnames(comparison_df))) {
    comparison_label <- paste0(contrast_condition, " vs ", ref_condition)
    comparison_df %>%
      select(tf, logFC) %>%
      mutate(
        logFC = as.numeric(logFC),
        comparison = comparison_label,
        reference_condition = ref_condition
      ) %>%
      filter(!grepl("^RANDOM", tf)) %>%
      as.data.frame() -> new_logfc_df
    primetime_logfc_df <- bind_rows(primetime_logfc_df, new_logfc_df)
  }

  comparison_df %>%
    select(tf, all_of(c(ref_condition, contrast_condition))) %>%
    pivot_longer(cols = all_of(c(ref_condition, contrast_condition)), names_to = "condition", values_to = "activity") %>%
    mutate(activity = as.numeric(activity)) %>%
    filter(!grepl("^RANDOM", tf)) %>%
    group_by(condition, tf) %>%
    summarise(activity = mean(activity, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame() -> new_df

  condition_activity_df <- bind_rows(condition_activity_df, new_df)
  reference_conditions <- c(reference_conditions, ref_condition)
}

plots_to_print <- list()

if (nrow(condition_activity_df) > 0) {
  condition_activity_df <- condition_activity_df %>%
    group_by(condition, tf) %>%
    summarise(activity = mean(activity, na.rm = TRUE), .groups = "drop")

  condition_wide <- condition_activity_df %>%
    pivot_wider(names_from = tf, values_from = activity) %>%
    arrange(condition)

  condition_names <- condition_wide$condition
  condition_matrix <- condition_wide %>% select(-condition) %>% as.matrix()
  rownames(condition_matrix) <- condition_names

  valid_columns <- apply(condition_matrix, 2, function(x) !all(is.na(x)) && sd(x, na.rm = TRUE) > 0)
  condition_matrix <- condition_matrix[, valid_columns, drop = FALSE]
  condition_matrix <- condition_matrix[complete.cases(condition_matrix), , drop = FALSE]

  coords <- compute_umap_coords(condition_matrix, rownames(condition_matrix), python_bin)
  if (!is.null(coords)) {
    coords$status <- ifelse(coords$label %in% unique(reference_conditions), "Reference", "Contrast")
    coords$label <- factor(coords$label, levels = coords$label[order(coords$UMAP1, coords$UMAP2)])

    p_conditions <- ggplot(coords, aes(x = UMAP1, y = UMAP2, color = status, label = label)) +
      geom_point(size = 3.5) +
      geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2) +
      scale_color_manual(values = c("Reference" = "#f37f80", "Contrast" = "#6495ed"), name = NULL) +
      theme_bw() +
      labs(
        title = "UMAP of Primetime Condition Activities",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    plots_to_print[[length(plots_to_print) + 1]] <- p_conditions
  }
}

if (nrow(primetime_logfc_df) > 0) {
  logfc_wide <- primetime_logfc_df %>%
    group_by(comparison, tf, reference_condition) %>%
    summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = tf, values_from = logFC, values_fill = 0) %>%
    arrange(comparison)

  comparison_labels <- logfc_wide$comparison
  logfc_matrix <- logfc_wide %>% select(-comparison, -reference_condition) %>% as.matrix()
  rownames(logfc_matrix) <- comparison_labels

  valid_logfc_cols <- apply(logfc_matrix, 2, function(x) !all(is.na(x)) && sd(x, na.rm = TRUE) > 0)
  logfc_matrix <- logfc_matrix[, valid_logfc_cols, drop = FALSE]
  logfc_matrix[is.na(logfc_matrix)] <- 0

  coords_logfc <- compute_umap_coords(logfc_matrix, rownames(logfc_matrix), python_bin)
  if (!is.null(coords_logfc)) {
    ref_lookup <- setNames(logfc_wide$reference_condition, logfc_wide$comparison)
    coords_logfc$reference_condition <- ref_lookup[coords_logfc$label]
    coords_logfc$status <- ifelse(grepl("control|ctrl|dmso", tolower(coords_logfc$reference_condition)), "vs Control", "Other reference")
    coords_logfc$label <- factor(coords_logfc$label, levels = coords_logfc$label[order(coords_logfc$UMAP1, coords_logfc$UMAP2)])

    p_logfc <- ggplot(coords_logfc, aes(x = UMAP1, y = UMAP2, color = status, label = label)) +
      geom_point(size = 3.5) +
      geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2) +
      scale_color_manual(values = c("vs Control" = "#f37f80", "Other reference" = "#6495ed"), name = NULL) +
      theme_bw() +
      labs(
        title = "UMAP of Primetime logFC Profiles (Comparison x TF)",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    plots_to_print[[length(plots_to_print) + 1]] <- p_logfc
  }
}

output_root <- dirname(dirname(normalizePath(opt$output, mustWork = FALSE)))
barcode_activity_path <- file.path(output_root, "tmp_primetime", "activity", "barcode_activity.txt")

if (file.exists(barcode_activity_path) && file.info(barcode_activity_path)$size > 0) {
  qc_df <- suppressWarnings(read.table(barcode_activity_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
  if (all(c("cDNA_sample", "tf") %in% colnames(qc_df))) {
    if (!"log2_mean_RPM" %in% colnames(qc_df) && "mean_RPM" %in% colnames(qc_df)) {
      qc_df$log2_mean_RPM <- log2(as.numeric(qc_df$mean_RPM) + 1)
    }
    if ("log2_mean_RPM" %in% colnames(qc_df)) {
      sample_activity_df <- qc_df %>%
        filter(!grepl("^RANDOM", tf)) %>%
        mutate(log2_mean_RPM = as.numeric(log2_mean_RPM)) %>%
        group_by(cDNA_sample, tf) %>%
        summarise(activity = mean(log2_mean_RPM, na.rm = TRUE), .groups = "drop")

      sample_wide <- sample_activity_df %>%
        pivot_wider(names_from = tf, values_from = activity) %>%
        arrange(cDNA_sample)

      sample_names <- sample_wide$cDNA_sample
      sample_matrix <- sample_wide %>% select(-cDNA_sample) %>% as.matrix()
      rownames(sample_matrix) <- sample_names

      valid_sample_cols <- apply(sample_matrix, 2, function(x) !all(is.na(x)) && sd(x, na.rm = TRUE) > 0)
      sample_matrix <- sample_matrix[, valid_sample_cols, drop = FALSE]
      sample_matrix <- sample_matrix[complete.cases(sample_matrix), , drop = FALSE]

      coords_samples <- compute_umap_coords(sample_matrix, rownames(sample_matrix), python_bin)
      if (!is.null(coords_samples)) {
        coords_samples$status <- ifelse(grepl("control|ctrl|dmso", tolower(coords_samples$label)), "Control", "Sample")
        coords_samples$label <- factor(coords_samples$label, levels = coords_samples$label[order(coords_samples$UMAP1, coords_samples$UMAP2)])

        p_samples <- ggplot(coords_samples, aes(x = UMAP1, y = UMAP2, color = status, label = label)) +
          geom_point(size = 3.5) +
          geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2) +
          scale_color_manual(values = c("Control" = "#f37f80", "Sample" = "#6495ed"), name = NULL) +
          theme_bw() +
          labs(
            title = "UMAP of Individual Sample Activities (QC Barcode Activity)",
            x = "UMAP 1",
            y = "UMAP 2"
          ) +
          theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
        plots_to_print[[length(plots_to_print) + 1]] <- p_samples
      }
    }
  }
}

if (length(plots_to_print) == 0) {
  write_placeholder_pdf(opt$output, "No valid activity matrix could be constructed")
  quit(status = 0)
}

pdf(opt$output, width = 8, height = 6)
for (this_plot in plots_to_print) {
  print(this_plot)
}
invisible(dev.off())
