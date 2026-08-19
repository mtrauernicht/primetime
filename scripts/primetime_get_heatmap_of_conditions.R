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
    make_option(c("--negative-correlation-stats"), type = "character", default = "", help = "TSV with comparison and correlation flag columns (supports legacy has_negative_correlation)"),
    make_option(c("--read-count-summary"), type = "character", default = "", help = "TSV with per-sample read counts and read_count_lt_25000 column"),
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

# Be tolerant to optparse naming differences for --negative-correlation-stats.
negative_correlation_stats_opt <- opt$negative_correlation_stats
if (is.null(negative_correlation_stats_opt) && !is.null(opt$negative.correlation.stats)) {
    negative_correlation_stats_opt <- opt$negative.correlation.stats
}
if (is.null(negative_correlation_stats_opt) && !is.null(opt$`negative-correlation-stats`)) {
    negative_correlation_stats_opt <- opt$`negative-correlation-stats`
}

# Be tolerant to optparse naming differences for --read-count-summary.
read_count_summary_opt <- opt$read_count_summary
if (is.null(read_count_summary_opt) && !is.null(opt$read.count.summary)) {
    read_count_summary_opt <- opt$read.count.summary
}
if (is.null(read_count_summary_opt) && !is.null(opt$`read-count-summary`)) {
    read_count_summary_opt <- opt$`read-count-summary`
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

# Build a flexible p-value column for downstream filtering.
pvalue_candidates <- c("p_adjusted", "adj.P.Val", "P.Value", "pvalue")
pvalue_column <- intersect(pvalue_candidates, colnames(all_comparisons_df))[1]
if (!is.na(pvalue_column) && nzchar(pvalue_column)) {
    all_comparisons_df$p_value_filter <- suppressWarnings(as.numeric(all_comparisons_df[[pvalue_column]]))
} else {
    all_comparisons_df$p_value_filter <- NA_real_
}

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

normalize_condition_key <- function(x) {
    key <- trimws(as.character(x))
    key <- gsub("[-[:space:]]+", "_", key)
    key <- gsub("[^A-Za-z0-9_]", "", key)
    tolower(key)
}

row_condition_lookup <- all_comparisons_df %>%
    distinct(condition, condition_raw) %>%
    mutate(condition_key = normalize_condition_key(condition_raw))

correlation_by_condition <- data.frame(
    condition_raw = character(0),
    condition_key = character(0),
    pDNA_correlation_strength = numeric(0),
    correlation_direction = character(0)
)

if (!is.null(negative_correlation_stats_opt) && nzchar(trimws(negative_correlation_stats_opt)) && file.exists(negative_correlation_stats_opt) && file.info(negative_correlation_stats_opt)$size > 0) {
    negative_corr_df <- suppressWarnings(read.delim(
        negative_correlation_stats_opt,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        fill = TRUE,
        quote = "",
        comment.char = "",
        na.strings = c("", "NA")
    ))

    if (nrow(negative_corr_df) > 0 && "comparison" %in% colnames(negative_corr_df)) {
        correlation_by_condition <- negative_corr_df %>%
            filter(!is.na(comparison), trimws(as.character(comparison)) != "") %>%
            mutate(
                condition_raw = trimws(sub("_[0-9]+$", "", comparison)),
                condition_key = normalize_condition_key(condition_raw),
                corr_numeric = suppressWarnings(as.numeric(corr)),
                slope_numeric = suppressWarnings(as.numeric(slope)),
                corr_p_adj_numeric = suppressWarnings(as.numeric(corr_p_two_sided_adj)),
                has_correlation = if ("has_correlation" %in% colnames(.)) {
                    tolower(trimws(as.character(has_correlation))) %in% c("true", "t", "1")
                } else {
                    !is.na(corr_p_adj_numeric) & corr_p_adj_numeric < 0.0001
                },
                has_negative_correlation = if ("has_negative_correlation" %in% colnames(.)) {
                    tolower(trimws(as.character(has_negative_correlation))) %in% c("true", "t", "1")
                } else {
                    has_correlation & !is.na(corr_numeric) & corr_numeric < 0
                },
                has_positive_correlation = if ("has_positive_correlation" %in% colnames(.)) {
                    tolower(trimws(as.character(has_positive_correlation))) %in% c("true", "t", "1")
                } else {
                    has_correlation & !is.na(slope_numeric) & !is.na(corr_numeric) & slope_numeric > 0.125 & corr_numeric > 0
                },
                has_correlation = has_correlation | has_negative_correlation | has_positive_correlation,
                pDNA_correlation_strength = ifelse(has_correlation & !is.na(corr_numeric), pmin(1, abs(corr_numeric)), ifelse(has_correlation, 1, 0)),
                correlation_direction = case_when(
                    has_positive_correlation ~ "positive_slope_gt_0.125",
                    has_negative_correlation ~ "negative",
                    has_correlation ~ "correlated_other",
                    TRUE ~ "none"
                )
            ) %>%
            group_by(condition_raw, condition_key) %>%
            summarise(
                pDNA_correlation_strength = max(pDNA_correlation_strength, na.rm = TRUE),
                has_positive_correlation = any(correlation_direction == "positive_slope_gt_0.125"),
                has_negative_correlation = any(correlation_direction == "negative"),
                has_other_correlation = any(correlation_direction == "correlated_other"),
                correlation_direction = case_when(
                    has_positive_correlation ~ "positive_slope_gt_0.125",
                    has_negative_correlation ~ "negative",
                    has_other_correlation ~ "correlated_other",
                    TRUE ~ "none"
                ),
                .groups = "drop"
            ) %>%
            select(condition_raw, condition_key, pDNA_correlation_strength, correlation_direction) %>%
            as.data.frame()
    }
}

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
    left_join(correlation_by_condition %>% select(condition_key, pDNA_correlation_strength, correlation_direction) %>% distinct(), by = "condition_key") %>%
    left_join(low_read_by_condition %>% select(condition_key, read_count) %>% distinct(), by = "condition_key") %>%
    mutate(
        pDNA_correlation_strength = ifelse(is.na(pDNA_correlation_strength), 0, pDNA_correlation_strength),
        pDNA_correlation_strength = pmax(0, pmin(1, pDNA_correlation_strength)),
        correlation_direction = ifelse(is.na(correlation_direction), "none", correlation_direction),
        read_count = suppressWarnings(as.numeric(read_count))
    ) %>%
    select(condition, pDNA_correlation_strength, correlation_direction, read_count) %>%
    distinct(condition, .keep_all = TRUE) %>%
    as.data.frame()

rownames(row_annotation_df) <- row_annotation_df$condition
row_annotation_df <- row_annotation_df[rownames(matrix_values), c("pDNA_correlation_strength", "correlation_direction", "read_count"), drop = FALSE]
row_annotation_df$correlation_direction <- factor(
    row_annotation_df$correlation_direction,
    levels = c("none", "correlated_other", "negative", "positive_slope_gt_0.125")
)

# Emphasize low-count samples: values above 50k are saturated to white in the annotation scale.
read_count_cap <- 50000
row_annotation_df$read_count <- ifelse(
    is.na(row_annotation_df$read_count),
    NA_real_,
    pmin(read_count_cap, row_annotation_df$read_count)
)

message("Heatmap annotations: read_count min on heatmap rows = ", round(min(row_annotation_df$read_count, na.rm = TRUE), 1))
message("Heatmap annotations: read_count max on heatmap rows = ", round(max(row_annotation_df$read_count, na.rm = TRUE), 1))
message("Heatmap annotations: read_count cap for coloring = ", read_count_cap)
message("Heatmap annotations: pDNA_correlation_strength max on heatmap rows = ", round(max(row_annotation_df$pDNA_correlation_strength, na.rm = TRUE), 3))

annotation_colors <- list(
    pDNA_correlation_strength = colorRampPalette(c("#ffffff", "#c0392b"))(100),
    correlation_direction = c(
        none = "#f0f0f0",
        correlated_other = "#f1c40f",
        negative = "#4c78a8",
        positive_slope_gt_0.125 = "#d62728"
    ),
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
    if (nrow(matrix_sub) >= 2) {
        row_cluster_sub <- hclust(dist(matrix_sub), method = "ward.D2")
    }
    if (ncol(matrix_sub) >= 2) {
        col_cluster_sub <- hclust(dist(t(matrix_sub)), method = "ward.D2")
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
pheatmap(
    matrix_plot,
    cluster_rows = row_cluster,
    cluster_cols = col_cluster,
    annotation_row = row_annotation_df,
    annotation_colors = annotation_colors,
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

# Also include a full fold-change heatmap where non-significant values are not greyed out.
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

# Third heatmap: remove low-read-count samples and samples with significant
# negative/positive control-residual correlation.
rows_keep_qc <- rownames(row_annotation_df)[
    !(row_annotation_df$correlation_direction %in% c("negative", "positive_slope_gt_0.125")) &
    (is.na(row_annotation_df$read_count) | row_annotation_df$read_count >= 25000)
]

plot_filtered_heatmap(
    matrix_source = matrix_values,
    annotation_source = row_annotation_df,
    keep_rows = rows_keep_qc,
    title_text = "TF Fold-Change Heatmap (Filtered: read_count >= 25k and no significant +/- correlation)",
    break_limits = c(-2, 2),
    mask_nonsig = FALSE
)

# Fourth heatmap: keep only samples with at least one strong/sig TF change.
condition_effect_filter <- all_comparisons_df %>%
    group_by(condition) %>%
    summarise(
        has_large_abs_logfc = any(is.finite(logFC) & abs(logFC) > 0.5, na.rm = TRUE),
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
    title_text = "TF Fold-Change Heatmap (Filtered: |logFC| > 0.5 or p-value < 0.05)",
    break_limits = c(-2, 2),
    mask_nonsig = FALSE
)

# Fifth heatmap: use all rows and median-center per TF.
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

invisible(dev.off())