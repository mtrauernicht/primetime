suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
    library(ggplot2)
    library(ggthemes)
    library(reshape2)
    library(pheatmap)
})

option_list <- list(
    make_option(c("--results"), type = "character", help = "Comma-separated list of results files"),
    make_option(c("--results-file"), type = "character", help = "File containing newline-separated results files"),
    make_option(c("--output"), type = "character", help = "Output file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Be tolerant to optparse naming differences for --results-file.
results_file_opt <- opt$results_file
if (is.null(results_file_opt) && !is.null(opt$results.file)) {
    results_file_opt <- opt$results.file
}
if (is.null(results_file_opt) && !is.null(opt$`results-file`)) {
    results_file_opt <- opt$`results-file`
}

if (is.null(opt$output) || (is.null(opt$results) && is.null(results_file_opt))) {
    stop("Missing required arguments")
}

result_files <- character(0)
if (!is.null(opt$results)) {
    result_files <- c(result_files, unlist(strsplit(opt$results, ",")))
}
if (!is.null(results_file_opt)) {
    filelist_lines <- readLines(results_file_opt, warn = FALSE)
    # Support both newline-separated and whitespace-separated path lists.
    filelist_tokens <- unlist(strsplit(paste(filelist_lines, collapse = "\n"), "[[:space:]]+"))
    result_files <- c(result_files, filelist_tokens)
}
result_files <- unique(trimws(result_files))
result_files <- result_files[result_files != ""]

all_comparisons_df <- data.frame()
processed_files <- 0L

extract_comparison_parts <- function(file_base) {
    vs_pos <- regexpr("_vs_", file_base, fixed = TRUE)[1]
    if (vs_pos <= 0) {
        return(NULL)
    }
    contrast <- substr(file_base, 1, vs_pos - 1)
    reference <- substr(file_base, vs_pos + 4, nchar(file_base))
    if (contrast == "" || reference == "") {
        return(NULL)
    }
    list(contrast = contrast, reference = reference)
}

for (this_comparison in result_files) {
    if (!file.exists(this_comparison) || file.info(this_comparison)$size == 0) {
        next
    }

    file_base <- sub("\\.txt$", "", basename(this_comparison))
    comparison_parts <- extract_comparison_parts(file_base)
    if (is.null(comparison_parts)) {
        next
    }

    comparison_df <- suppressWarnings(read.table(this_comparison, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))

    if (!all(c("tf", "logFC") %in% colnames(comparison_df))) {
        next
    }

    new_df <- comparison_df %>%
        select(tf, logFC, any_of("sig")) %>%
        mutate(
            logFC = as.numeric(logFC),
            sig = if ("sig" %in% colnames(.)) as.character(sig) else "NS",
            condition_raw = comparison_parts$contrast,
            reference_condition = comparison_parts$reference,
            comparison = paste0(comparison_parts$contrast, " vs ", comparison_parts$reference)
        ) %>%
        filter(!grepl("^RANDOM", tf)) %>%
        as.data.frame()

    all_comparisons_df <- bind_rows(all_comparisons_df, new_df)
    processed_files <- processed_files + 1L
}

if (nrow(all_comparisons_df) == 0) {
    message("Heatmap: processed 0 comparison files with usable fold-change columns out of ", length(result_files), " listed files")
    pdf(opt$output, width = 8, height = 4)
    plot.new()
    text(0.5, 0.5, "No comparable conditions found")
    invisible(dev.off())
    quit(status = 0)
}

message("Heatmap: processed ", processed_files, " comparison files with usable fold-change columns")

# Use plain condition names if each condition maps to one reference, otherwise use full comparison labels.
condition_ref_count <- all_comparisons_df %>%
    distinct(condition_raw, reference_condition) %>%
    count(condition_raw, name = "n_refs")

all_comparisons_df <- all_comparisons_df %>%
    left_join(condition_ref_count, by = "condition_raw") %>%
    mutate(condition = ifelse(n_refs > 1, comparison, condition_raw))

heatmap_wide <- all_comparisons_df %>%
    group_by(condition, tf) %>%
    summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = tf, values_from = logFC, values_fill = 0)

sig_wide <- all_comparisons_df %>%
    group_by(condition, tf) %>%
    summarise(is_sig = any(sig != "NS", na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = tf, values_from = is_sig, values_fill = FALSE)

if (nrow(heatmap_wide) == 0 || ncol(heatmap_wide) <= 1) {
    pdf(opt$output, width = 8, height = 4)
    plot.new()
    text(0.5, 0.5, "No fold-change values available for heatmap")
    invisible(dev.off())
    quit(status = 0)
}

matrix_values <- as.matrix(heatmap_wide %>% select(-condition))
rownames(matrix_values) <- heatmap_wide$condition
sig_matrix <- as.matrix(sig_wide %>% select(-condition))
rownames(sig_matrix) <- sig_wide$condition

# Ensure there are no NA/NaN values for pheatmap.
matrix_values[is.na(matrix_values)] <- 0

# Drop rows/columns with zero variance to avoid clustering artifacts.
if (nrow(matrix_values) > 1) {
    row_keep <- apply(matrix_values, 1, function(x) sd(x) > 0)
    if (any(row_keep)) {
        matrix_values <- matrix_values[row_keep, , drop = FALSE]
        sig_matrix <- sig_matrix[row_keep, , drop = FALSE]
    }
}
if (ncol(matrix_values) > 1) {
    col_keep <- apply(matrix_values, 2, function(x) sd(x) > 0)
    if (any(col_keep)) {
        matrix_values <- matrix_values[, col_keep, drop = FALSE]
        sig_matrix <- sig_matrix[, col_keep, drop = FALSE]
    }
}

if (nrow(matrix_values) == 0 || ncol(matrix_values) == 0) {
    pdf(opt$output, width = 8, height = 4)
    plot.new()
    text(0.5, 0.5, "No variable fold-change values available for heatmap")
    invisible(dev.off())
    quit(status = 0)
}

matrix_plot <- matrix_values
matrix_plot[!sig_matrix] <- NA_real_

# Compute dendrograms from the zero-filled matrix to avoid NA issues in hclust.
row_cluster <- FALSE
col_cluster <- FALSE
if (nrow(matrix_values) >= 2) {
    row_cluster <- hclust(dist(matrix_values), method = "ward.D2")
}
if (ncol(matrix_values) >= 2) {
    col_cluster <- hclust(dist(t(matrix_values)), method = "ward.D2")
}

pdf(
    opt$output,
    width = max(10, min(40, 4 + 0.16 * ncol(matrix_values))),
    height = max(6, min(30, 3 + 0.45 * nrow(matrix_values)))
)
pheatmap(
    matrix_plot,
    cluster_rows = row_cluster,
    cluster_cols = col_cluster,
    border_color = NA,
    color = colorRampPalette(c("#6495ed", "white", "#f37f80"))(100),
    breaks = seq(-1, 1, length.out = 101),
    na_col = "#eeedf1",
    main = "TF Fold-Change Heatmap Across Conditions",
    cellwidth = 6,
    cellheight = 6,
    fontsize_row = 6,
    fontsize_col = 6,
    angle_col = 90
)
invisible(dev.off())