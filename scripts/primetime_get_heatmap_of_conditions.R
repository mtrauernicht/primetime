suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
    library(ggplot2)
    library(ggthemes)
    library(reshape2)
    library(pheatmap)
    library(data.table)
})

option_list <- list(
    make_option(c("--results"), type = "character", help = "Comma-separated list of results files"),
    make_option(c("--results-file"), type = "character", help = "File containing newline-separated results files"),
    make_option(c("--read-count-summary"), type = "character", default = "", help = "TSV with per-sample read counts and read_count_lt_25000 column"),
    make_option(c("--sample-viability-summary"), type = "character", default = "", help = "TSV with per-sample mean viability values"),
    make_option(c("--sample_info"), type = "character", default = "", help = "Optional CSV with sample metadata columns such as sample and class"),
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


# Be tolerant to optparse naming differences for --read-count-summary.
read_count_summary_opt <- opt$read_count_summary
if (is.null(read_count_summary_opt) && !is.null(opt$read.count.summary)) {
    read_count_summary_opt <- opt$read.count.summary
}
if (is.null(read_count_summary_opt) && !is.null(opt$`read-count-summary`)) {
    read_count_summary_opt <- opt$`read-count-summary`
}

sample_viability_summary_opt <- opt$sample_viability_summary
if (is.null(sample_viability_summary_opt) && !is.null(opt$sample.viability.summary)) {
    sample_viability_summary_opt <- opt$sample.viability.summary
}
if (is.null(sample_viability_summary_opt) && !is.null(opt$`sample-viability-summary`)) {
    sample_viability_summary_opt <- opt$`sample-viability-summary`
}

if (is.null(opt$output) || (is.null(opt$results) && is.null(results_file_opt))) {
    stop("Missing required arguments")
}

sample_info_path <- opt$sample_info
if (is.null(sample_info_path) || is.na(sample_info_path)) sample_info_path <- ""

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

normalize_condition_key <- function(x) {
    key <- trimws(as.character(x))
    key <- gsub("[-[:space:]]+", "_", key)
    key <- gsub("[^A-Za-z0-9_]", "", key)
    tolower(key)
}

sample_info_lookup <- data.frame(
    sample_key = character(0),
    sample_target = character(0),
    sample_class = character(0)
)

if (nzchar(trimws(sample_info_path)) && file.exists(sample_info_path) && file.info(sample_info_path)$size > 0) {
    sample_info_df <- suppressWarnings(read.csv(sample_info_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "N/A")))
    if (nrow(sample_info_df) > 0) {
        sample_col <- intersect(c("sample", "condition", "sample_name", "name", "sample_id"), colnames(sample_info_df))[1]
        target_col <- intersect(c("sample_target", "target", "Sample_target", "Target", "gene_target"), colnames(sample_info_df))[1]
        class_col <- intersect(c("class", "sample_class", "compound_class", "Class"), colnames(sample_info_df))[1]

        if (!is.na(sample_col)) {
            sample_info_lookup <- sample_info_df %>%
                mutate(
                    sample_key = normalize_condition_key(.data[[sample_col]]),
                    sample_name = trimws(as.character(.data[[sample_col]])),
                    target_value = if (!is.na(target_col) && target_col %in% colnames(.)) trimws(as.character(.data[[target_col]])) else NA_character_,
                    sample_target = ifelse(
                        !is.na(target_value) & nzchar(target_value),
                        paste0(sample_name, "_", target_value),
                        sample_name
                    ),
                    sample_class = if (!is.na(class_col) && class_col %in% colnames(.)) trimws(as.character(.data[[class_col]])) else NA_character_
                ) %>%
                mutate(
                    sample_target = ifelse(is.na(sample_target) | trimws(sample_target) == "", sample_name, sample_target),
                    sample_class = ifelse(is.na(sample_class) | trimws(sample_class) == "", "unknown", sample_class)
                ) %>%
                group_by(sample_key) %>%
                summarise(
                    sample_target = first(sample_target[!is.na(sample_target) & nzchar(sample_target)], default = first(sample_key)),
                    sample_class = first(sample_class[!is.na(sample_class) & nzchar(sample_class)], default = "unknown"),
                    .groups = "drop"
                ) %>%
                as.data.frame()
        }
    }
}

sample_viability_lookup <- data.frame(sample_key = character(0), mean_viability = numeric(0))
if (!is.null(sample_viability_summary_opt) && nzchar(trimws(sample_viability_summary_opt)) && file.exists(sample_viability_summary_opt) && file.info(sample_viability_summary_opt)$size > 0) {
    sample_viability_df <- suppressWarnings(read.delim(sample_viability_summary_opt, stringsAsFactors = FALSE, check.names = FALSE))
    if (all(c("sample", "mean_viability") %in% colnames(sample_viability_df))) {
        sample_viability_lookup <- sample_viability_df %>%
            transmute(
                sample_key = normalize_condition_key(sample),
                mean_viability = suppressWarnings(as.numeric(mean_viability))
            ) %>%
            filter(is.finite(mean_viability)) %>%
            group_by(sample_key) %>%
            summarise(mean_viability = mean(mean_viability), .groups = "drop") %>%
            as.data.frame()
    }
}

pvalue_candidates <- c("p_adjusted", "adj.P.Val", "P.Value", "pvalue")

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
    select(tf, logFC, any_of(c("sig", pvalue_candidates))) %>%
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
    mutate(
        condition_key = normalize_condition_key(condition_raw),
        sample_key = condition_key
    ) %>%
    left_join(sample_info_lookup, by = c("sample_key" = "sample_key")) %>%
    left_join(sample_viability_lookup, by = c("sample_key" = "sample_key")) %>%
    mutate(
        sample_label = ifelse(!is.na(sample_target) & trimws(sample_target) != "", sample_target, condition_raw),
        condition = ifelse(n_refs > 1, comparison, sample_label),
        sample_class = ifelse(is.na(sample_class) | trimws(sample_class) == "", "unknown", sample_class),
        mean_viability = suppressWarnings(as.numeric(mean_viability))
    )

# Build a flexible p-value column for downstream filtering.
pvalue_column <- intersect(pvalue_candidates, colnames(all_comparisons_df))[1]
if (!is.na(pvalue_column) && nzchar(pvalue_column)) {
    all_comparisons_df$p_value_filter <- suppressWarnings(as.numeric(all_comparisons_df[[pvalue_column]]))
} else {
    all_comparisons_df$p_value_filter <- NA_real_
}

heatmap_wide <- all_comparisons_df %>%
    group_by(sample_label, tf) %>%
    summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = tf, values_from = logFC, values_fill = 0)

sig_wide <- all_comparisons_df %>%
    group_by(sample_label, tf) %>%
    summarise(is_sig = any(sig != "NS", na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = tf, values_from = is_sig, values_fill = FALSE)

if (nrow(heatmap_wide) == 0 || ncol(heatmap_wide) <= 1) {
    pdf(opt$output, width = 8, height = 4)
    plot.new()
    text(0.5, 0.5, "No fold-change values available for heatmap")
    invisible(dev.off())
    quit(status = 0)
}

matrix_values <- as.matrix(heatmap_wide %>% select(-sample_label))
rownames(matrix_values) <- heatmap_wide$sample_label
sig_matrix <- as.matrix(sig_wide %>% select(-sample_label))
rownames(sig_matrix) <- sig_wide$sample_label

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

row_condition_lookup <- all_comparisons_df %>%
    distinct(sample_label, condition, condition_raw, sample_class, mean_viability) %>%
    mutate(condition_key = normalize_condition_key(condition_raw))

low_read_by_condition <- data.frame(
    condition_raw = character(0),
    condition_key = character(0),
    read_count = numeric(0)
)

if (!is.null(read_count_summary_opt) && nzchar(trimws(read_count_summary_opt)) && file.exists(read_count_summary_opt) && file.info(read_count_summary_opt)$size > 0) {
    read_count_df <- suppressWarnings(read.delim(
        read_count_summary_opt,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        fill = TRUE,
        quote = "",
        comment.char = "",
        na.strings = c("", "NA")
    ))

    if (nrow(read_count_df) > 0) {
        if (!"condition" %in% colnames(read_count_df)) {
            if ("replicate" %in% colnames(read_count_df)) {
                read_count_df$condition <- sub("_[0-9]+$", "", as.character(read_count_df$replicate))
            } else if ("sample" %in% colnames(read_count_df)) {
                read_count_df$condition <- as.character(read_count_df$sample)
            }
        }

        if (all(c("condition", "total_read_count") %in% colnames(read_count_df))) {
            low_read_by_condition <- read_count_df %>%
                mutate(
                    condition_raw = trimws(as.character(condition)),
                    condition_key = normalize_condition_key(condition_raw),
                    total_read_count = suppressWarnings(as.numeric(total_read_count))
                ) %>%
                filter(!is.na(total_read_count), is.finite(total_read_count)) %>%
                group_by(condition_raw, condition_key) %>%
                summarise(read_count = min(total_read_count, na.rm = TRUE), .groups = "drop") %>%
                as.data.frame()
        }
    }
}

if (nrow(low_read_by_condition) > 0) {
    message("Heatmap annotations: conditions with read-count values from read_count_summary = ", nrow(low_read_by_condition))
}

row_annotation_df <- row_condition_lookup %>%
    left_join(low_read_by_condition %>% select(condition_key, read_count) %>% distinct(), by = "condition_key") %>%
    mutate(
        read_count = suppressWarnings(as.numeric(read_count)),
        sample_class = ifelse(is.na(sample_class) | trimws(sample_class) == "", "unknown", sample_class),
        mean_viability = suppressWarnings(as.numeric(mean_viability))
    ) %>%
    select(sample_label, sample_class, mean_viability, read_count) %>%
    distinct(sample_label, .keep_all = TRUE) %>%
    as.data.frame()

rownames(row_annotation_df) <- row_annotation_df$sample_label
row_annotation_df <- row_annotation_df[rownames(matrix_values), c("sample_class", "mean_viability", "read_count"), drop = FALSE]

# Some runs have no usable metadata for row labels. Guarantee finite annotation values so
# pheatmap can render a valid heatmap instead of failing on all-NA/Inf annotation ranges.
if (nrow(row_annotation_df) == 0) {
    row_annotation_df <- data.frame(
        sample_class = character(0),
        mean_viability = numeric(0),
        read_count = numeric(0),
        row.names = character(0)
    )
}
row_annotation_df$sample_class <- ifelse(is.na(row_annotation_df$sample_class) | trimws(row_annotation_df$sample_class) == "", "unknown", row_annotation_df$sample_class)
row_annotation_df$mean_viability <- suppressWarnings(as.numeric(row_annotation_df$mean_viability))
row_annotation_df$read_count <- suppressWarnings(as.numeric(row_annotation_df$read_count))
row_annotation_df$mean_viability[!is.finite(row_annotation_df$mean_viability)] <- 0
row_annotation_df$read_count[!is.finite(row_annotation_df$read_count)] <- 0

# Emphasize low-count samples: values above 50k are saturated to white in the annotation scale.
read_count_cap <- 50000
row_annotation_df$read_count <- pmin(read_count_cap, row_annotation_df$read_count)

finite_read_counts <- row_annotation_df$read_count[is.finite(row_annotation_df$read_count)]
if (length(finite_read_counts) > 0) {
    message("Heatmap annotations: read_count min on heatmap rows = ", round(min(finite_read_counts), 1))
    message("Heatmap annotations: read_count max on heatmap rows = ", round(max(finite_read_counts), 1))
} else {
    message("Heatmap annotations: read_count min/max unavailable; all read counts were missing or non-finite")
}
message("Heatmap annotations: read_count cap for coloring = ", read_count_cap)

sample_class_levels <- unique(as.character(row_annotation_df$sample_class))
annotation_colors <- list(
    sample_class = setNames(grDevices::hcl.colors(length(sample_class_levels), palette = "Set2"), sample_class_levels),
    mean_viability = colorRampPalette(c("#ffffff", "#2e7d32"))(100),
    read_count = colorRampPalette(c("#1680a9", "#ffffff"))(100)
)

plot_no_rows_page <- function(label_text) {
    plot.new()
    text(0.5, 0.5, label_text)
}

plot_filtered_heatmap <- function(matrix_source, annotation_source, keep_rows, title_text, break_limits = c(-2, 2), mask_nonsig = FALSE, sig_source = NULL) {
    keep_rows <- rownames(matrix_source)[rownames(matrix_source) %in% keep_rows]
    if (length(keep_rows) == 0) {
        plot_no_rows_page(paste0(title_text, "\nNo rows passed filter"))
        return(invisible(NULL))
    }

    matrix_sub <- matrix_source[keep_rows, , drop = FALSE]
    annotation_sub <- annotation_source[keep_rows, , drop = FALSE]

    if (mask_nonsig) {
        if (is.null(sig_source)) {
            stop("sig_source is required when mask_nonsig = TRUE")
        }
        sig_sub <- sig_source[keep_rows, , drop = FALSE]
        matrix_sub[!sig_sub] <- NA_real_
    }

    if (nrow(matrix_sub) == 0 || ncol(matrix_sub) == 0) {
        plot_no_rows_page(paste0(title_text, "\nNo matrix values available"))
        return(invisible(NULL))
    }

    row_cluster_sub <- FALSE
    col_cluster_sub <- FALSE
    matrix_for_clustering <- matrix_sub
    matrix_for_clustering[is.na(matrix_for_clustering)] <- 0
    if (nrow(matrix_sub) >= 2) {
        row_cluster_sub <- hclust(dist(matrix_for_clustering), method = "ward.D2")
    }
    if (ncol(matrix_sub) >= 2) {
        col_cluster_sub <- hclust(dist(t(matrix_for_clustering)), method = "ward.D2")
    }

    pheatmap(
        matrix_sub,
        cluster_rows = row_cluster_sub,
        cluster_cols = col_cluster_sub,
        annotation_row = annotation_sub,
        annotation_colors = annotation_colors,
        border_color = NA,
        color = colorRampPalette(c("#6495ed", "white", "#f37f80"))(100),
        breaks = seq(break_limits[1], break_limits[2], length.out = 101),
        na_col = "#eeedf1",
        main = title_text,
        cellwidth = 6,
        cellheight = 6,
        fontsize_row = 6,
        fontsize_col = 6,
        angle_col = 90
    )
}

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

# Heatmap 1: full fold-change heatmap across all conditions (all values, none masked).
pheatmap(
    matrix_values,
    cluster_rows = row_cluster,
    cluster_cols = col_cluster,
    annotation_row = row_annotation_df,
    annotation_colors = annotation_colors,
    border_color = NA,
    color = colorRampPalette(c("#6495ed", "white", "#f37f80"))(100),
    breaks = seq(-2, 2, length.out = 101),
    main = "TF Fold-Change Heatmap Across Conditions (All Values)",
    cellwidth = 6,
    cellheight = 6,
    fontsize_row = 6,
    fontsize_col = 6,
    angle_col = 90
)

# Heatmap 2: use all rows and median-center per TF.
rows_keep_center <- rownames(matrix_values)
if (length(rows_keep_center) == 0) {
    plot_no_rows_page("TF Fold-Change Heatmap (Median-centered per TF)\nNo rows available")
} else {
    matrix_values_centered_input <- matrix_values[rows_keep_center, , drop = FALSE]
    annotation_centered <- row_annotation_df[rows_keep_center, , drop = FALSE]

    tf_median_logfc <- apply(matrix_values_centered_input, 2, median, na.rm = TRUE)
    matrix_values_centered <- sweep(matrix_values_centered_input, 2, tf_median_logfc, FUN = "-")

    centered_limit <- 2

    row_cluster_centered <- FALSE
    col_cluster_centered <- FALSE
    if (nrow(matrix_values_centered) >= 2) {
        row_cluster_centered <- hclust(dist(matrix_values_centered), method = "ward.D2")
    }
    if (ncol(matrix_values_centered) >= 2) {
        col_cluster_centered <- hclust(dist(t(matrix_values_centered)), method = "ward.D2")
    }

    pheatmap(
        matrix_values_centered,
        cluster_rows = row_cluster_centered,
        cluster_cols = col_cluster_centered,
        annotation_row = annotation_centered,
        annotation_colors = annotation_colors,
        border_color = NA,
        color = colorRampPalette(c("#6495ed", "white", "#f37f80"))(100),
        breaks = seq(-centered_limit, centered_limit, length.out = 101),
        main = "TF Fold-Change Heatmap (Median-centered per TF)",
        cellwidth = 6,
        cellheight = 6,
        fontsize_row = 6,
        fontsize_col = 6,
        angle_col = 90
    )
}

# Heatmap 3: keep only samples with at least one strong/sig TF change.
condition_effect_filter <- all_comparisons_df %>%
    group_by(condition) %>%
    summarise(
        has_large_abs_logfc = any(is.finite(logFC) & abs(logFC) > 1, na.rm = TRUE),
        has_sig_pval = any(is.finite(p_value_filter) & p_value_filter < 0.05, na.rm = TRUE),
        keep_condition = has_large_abs_logfc | has_sig_pval,
        .groups = "drop"
    )

rows_keep_effect <- condition_effect_filter %>%
    filter(keep_condition) %>%
    pull(condition) %>%
    as.character()

plot_filtered_heatmap(
    matrix_source = matrix_values,
    annotation_source = row_annotation_df,
    keep_rows = rows_keep_effect,
    title_text = "TF Fold-Change Heatmap (Filtered: |logFC| > 1 or p-value < 0.05)",
    break_limits = c(-2, 2),
    mask_nonsig = FALSE
)

# Heatmap 4: same filtered set as Heatmap 3, but with non-significant TF values greyed out.
plot_filtered_heatmap(
    matrix_source = matrix_values,
    annotation_source = row_annotation_df,
    keep_rows = rows_keep_effect,
    title_text = "TF Fold-Change Heatmap (Filtered: |logFC| > 1 or p-value < 0.05, non-significant values greyed out)",
    break_limits = c(-2, 2),
    mask_nonsig = TRUE,
    sig_source = sig_matrix
)

# Heatmap 5: signed significance score, logFC * -log10(p-value).

## Compute a logFC x -log10(p-value) score for each TF in each condition.
score_df <- all_comparisons_df %>%
    mutate(
        signed_score = logFC * -log10(p_value_filter)
    ) %>%
    group_by(sample_label, tf) %>%
    summarise(signed_score = mean(signed_score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = tf, values_from = signed_score, values_fill = 0)


## Convert to matrix and plot the signed score heatmap.
score_matrix <- as.matrix(score_df %>% select(-sample_label))
rownames(score_matrix) <- score_df$sample_label
score_matrix[!is.finite(score_matrix)] <- 0

## Subset the annotation table to match the filtered rows, same as the other heatmaps.
score_annotation <- row_annotation_df[rownames(row_annotation_df) %in% rownames(score_matrix), , drop = FALSE]
score_annotation <- score_annotation[rownames(score_matrix), , drop = FALSE]

## Plot the signed score heatmap if we have a valid matrix.
if (!is.null(score_matrix) && nrow(score_matrix) > 0 && ncol(score_matrix) > 0) {
    score_row_cluster <- if (nrow(score_matrix) >= 2) hclust(dist(score_matrix), method = "ward.D2") else FALSE
    score_col_cluster <- if (ncol(score_matrix) >= 2) hclust(dist(t(score_matrix)), method = "ward.D2") else FALSE
    pheatmap(
        score_matrix,
        cluster_rows = score_row_cluster,
        cluster_cols = score_col_cluster,
        annotation_row = score_annotation,
        annotation_colors = annotation_colors,
        border_color = NA,
        color = colorRampPalette(c("#6495ed", "white", "#f37f80"))(100),
        breaks = seq(-4, 4, length.out = 101),
        na_col = "#eeedf1",
        main = "TF Signed Score Heatmap Across Conditions (logFC * -log10(p-value))",
        cellwidth = 6,
        cellheight = 6,
        fontsize_row = 6,
        fontsize_col = 6,
        angle_col = 90
    )
} else {
    plot.new()
    text(0.5, 0.5, "TF Signed Score Heatmap\nNo rows passed |score| >= 2 filter")
}

# Sixth heatmap: The same as before, but filtering for rows with at least one |score| >= 1.
## Keep only rows with at least one |score| >= 1.
score_threshold <- 1
rows_keep_score <- rownames(score_matrix)[apply(score_matrix, 1, function(x) any(abs(x) >= score_threshold))]
score_matrix <- score_matrix[rows_keep_score, , drop = FALSE]

## Plot the signed score heatmap if we have a valid matrix.
if (!is.null(score_matrix) && nrow(score_matrix) > 0 && ncol(score_matrix) > 0) {
    score_row_cluster <- if (nrow(score_matrix) >= 2) hclust(dist(score_matrix), method = "ward.D2") else FALSE
    score_col_cluster <- if (ncol(score_matrix) >= 2) hclust(dist(t(score_matrix)), method = "ward.D2") else FALSE
    pheatmap(
        score_matrix,
        cluster_rows = score_row_cluster,
        cluster_cols = score_col_cluster,
        annotation_row = score_annotation,
        annotation_colors = annotation_colors,
        border_color = NA,
        color = colorRampPalette(c("#6495ed", "white", "#f37f80"))(100),
        breaks = seq(-4, 4, length.out = 101),
        na_col = "#eeedf1",
        main = "TF Signed Score Heatmap Across Conditions (logFC * -log10(p-value)), |score| >= 1",
        cellwidth = 6,
        cellheight = 6,
        fontsize_row = 6,
        fontsize_col = 6,
        angle_col = 90
    )
} else {
    plot.new()
    text(0.5, 0.5, "TF Signed Score Heatmap\nNo rows passed |score| >= 1 filter")
}

invisible(dev.off())