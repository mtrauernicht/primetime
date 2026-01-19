suppressPackageStartupMessages({
    # No external plotting deps to avoid env failures
})

# Parse arguments without extra dependencies
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0 || idx == length(args)) return(NA_character_)
    args[idx + 1]
}
opt <- list(
    results = get_arg("--results"),
    output = get_arg("--output")
)
if (is.na(opt$results) || is.na(opt$output)) {
    stop("Missing required arguments --results and --output")
}

# Construct the df (skip missing/empty files gracefully)
all_rows <- list()
comparison_files <- strsplit(opt$results, split = ",")[[1]]

safe_read_results <- function(path) {
    if (is.null(path) || is.na(path) || !file.exists(path) || file.info(path)$size == 0) return(NULL)
    tryCatch(read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) NULL)
}

required_cols <- c("tf", "sig", "logFC")
row_idx <- 1

for (this_comparison in comparison_files) {
    df <- safe_read_results(this_comparison)
    if (is.null(df)) next

    # Skip files missing required columns or with non-numeric logFC
    if (!all(required_cols %in% names(df))) next
    df$logFC <- suppressWarnings(as.numeric(df$logFC))
    df <- df[is.finite(df$logFC), c("tf", "sig", "logFC")]
    if (nrow(df) == 0) next

    tmp <- strsplit(sub("\\.txt$", "", basename(this_comparison)), "_vs_")[[1]]
    if (length(tmp) < 2) next
    comparison_name <- sprintf("%s vs. %s", tmp[1], tmp[2])

    df$comparison <- comparison_name
    all_rows[[row_idx]] <- df
    row_idx <- row_idx + 1
}

if (length(all_rows) > 0) {
    all_comparisons_df <- do.call(rbind, all_rows)
} else {
    all_comparisons_df <- data.frame()
}

if (length(comparison_files) == 0) {
    pdf(opt$output, width = 8, height = 6)
    plot.new(); text(0.5, 0.5, "No input comparisons provided")
    invisible(dev.off())
    quit(status = 0)
}

if (nrow(all_comparisons_df) == 0) {
    pdf(opt$output, width = 8, height = 6)
    plot.new(); text(0.5, 0.5, "No significant TFs to plot")
    invisible(dev.off())
    quit(status = 0)
}

# Create matrix for pheatmap (only significant TFs)
heatmap_dt <- all_comparisons_df[all_comparisons_df$sig != "NS", c("tf", "comparison", "logFC")]
if (nrow(heatmap_dt) == 0) {
    pdf(opt$output, width = 8, height = 6)
    plot.new(); text(0.5, 0.5, "No significant TFs to plot")
    invisible(dev.off())
    quit(status = 0)
}

# Collapse duplicates by mean logFC per tf/comparison
agg <- aggregate(logFC ~ tf + comparison, data = heatmap_dt, FUN = mean)
tfs <- sort(unique(agg$tf))
comparisons <- sort(unique(agg$comparison))
heatmap_matrix <- matrix(0, nrow = length(tfs), ncol = length(comparisons),
                         dimnames = list(tfs, comparisons))
for (i in seq_len(nrow(agg))) {
    heatmap_matrix[agg$tf[i], agg$comparison[i]] <- agg$logFC[i]
}

if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
    pdf(opt$output, width = 8, height = 6)
    plot.new(); text(0.5, 0.5, "No significant TFs to plot")
    invisible(dev.off())
    quit(status = 0)
}

# Create annotation for significant changes
sig_rows <- all_comparisons_df[all_comparisons_df$sig != "NS", ]
if (nrow(sig_rows) > 0) {
    comps <- unique(sig_rows$comparison)
    total_changes <- setNames(numeric(length(comps)), comps)
    for (comp in comps) {
        total_changes[comp] <- sum(sig_rows$comparison == comp)
    }
    annotation_col <- data.frame(total_changes = total_changes[comparisons], row.names = comparisons)
} else {
    annotation_col <- data.frame(row.names = comparisons)
}

# Define color palette (diverging: blue for down, red for up)
color_palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)

# Label columns with total change counts if available
col_labels <- colnames(heatmap_matrix)
if (!is.null(annotation_col) && nrow(annotation_col) > 0 && "total_changes" %in% colnames(annotation_col)) {
    lab_with_counts <- paste0(col_labels, " (n=", annotation_col[col_labels, "total_changes"], ")")
    names(lab_with_counts) <- col_labels
    col_labels <- lab_with_counts[colnames(heatmap_matrix)]
}

# Plot heatmap using base graphics to avoid external dependencies
pdf(opt$output, width = 12, height = 10)
heatmap(
    heatmap_matrix,
    Colv = NA,
    scale = "none",
    col = color_palette,
    margins = c(10, 10),
    labRow = rownames(heatmap_matrix),
    labCol = col_labels,
    main = "TF Activity Changes Across Comparisons"
)
invisible(dev.off())