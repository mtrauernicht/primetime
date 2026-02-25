suppressPackageStartupMessages({
    library(pheatmap)
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
    results_file = get_arg("--results-file"),
    bleedthrough = get_arg("--bleedthrough"),
    output = get_arg("--output")
)
if ((is.na(opt$results) && is.na(opt$results_file)) || is.na(opt$output)) {
    stop("Missing required arguments: provide --results or --results-file, and --output")
}

# Construct the df (skip missing/empty files gracefully)
all_rows <- list()
comparison_files <- character(0)
if (!is.na(opt$results)) {
    comparison_files <- c(comparison_files, strsplit(opt$results, split = ",")[[1]])
}
if (!is.na(opt$results_file)) {
    if (!file.exists(opt$results_file)) {
        stop("Results list file does not exist: ", opt$results_file)
    }
    listed <- readLines(opt$results_file, warn = FALSE)
    listed <- trimws(listed)
    listed <- listed[listed != ""]
    comparison_files <- c(comparison_files, listed)
}
comparison_files <- unique(comparison_files)

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

# Exclude negative controls
heatmap_dt <- heatmap_dt[!grepl("Random", heatmap_dt$tf, ignore.case = TRUE), ]
if (nrow(heatmap_dt) == 0) {
    pdf(opt$output, width = 8, height = 6)
    plot.new(); text(0.5, 0.5, "No significant TFs to plot (after filtering controls)")
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
# Using the same colors as in lollipop plots, with darker extremes
color_palette <- colorRampPalette(c("#2468E5", "#6495ed", "#FFFFFF", "#f37f80", "#E01518"))(100)

# Set fixed limits at -2 and +2 for symmetric scale
color_limits <- c(-2, 2)
breaks <- seq(color_limits[1], color_limits[2], length.out = 101)

# Clamp values to the color limits
heatmap_matrix_clamped <- heatmap_matrix
heatmap_matrix_clamped[heatmap_matrix_clamped < color_limits[1]] <- color_limits[1]
heatmap_matrix_clamped[heatmap_matrix_clamped > color_limits[2]] <- color_limits[2]

# Get original column labels before modification
col_labels <- colnames(heatmap_matrix_clamped)

# Load bleedthrough data and add to annotation BEFORE modifying column names
if (!is.na(opt$bleedthrough) && file.exists(opt$bleedthrough)) {
    bleedthrough_df <- read.delim(opt$bleedthrough, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    message("Loaded bleedthrough data with ", nrow(bleedthrough_df), " rows")
    message("First few conditions: ", paste(head(bleedthrough_df$condition, 3), collapse = ", "))
    message("First few comparison names: ", paste(head(col_labels, 3), collapse = " | "))
    
    # For each comparison, get bleedthrough values for both reference and contrast
    contrast_bleedthrough <- numeric(length(col_labels))
    reference_bleedthrough <- numeric(length(col_labels))
    names(contrast_bleedthrough) <- col_labels
    names(reference_bleedthrough) <- col_labels
    
    for (i in seq_along(col_labels)) {
        # Try both comparison formats: "X vs. Y" and "X_vs_Y"
        comp_parts <- if (grepl(" vs\\. ", col_labels[i])) {
            strsplit(col_labels[i], " vs\\. ")[[1]]
        } else {
            strsplit(col_labels[i], "_vs_")[[1]]
        }
        
        if (length(comp_parts) == 2) {
            ref_cond <- comp_parts[1]  # Reference condition
            contrast_cond <- comp_parts[2]  # Contrast condition
            
            ref_bleed <- bleedthrough_df$bleedthrough[bleedthrough_df$condition == ref_cond]
            contrast_bleed <- bleedthrough_df$bleedthrough[bleedthrough_df$condition == contrast_cond]
            
            # Take median if multiple values (from multiple replicates)
            if (length(ref_bleed) > 0) {
                reference_bleedthrough[col_labels[i]] <- median(ref_bleed, na.rm = TRUE)
            } else {
                reference_bleedthrough[col_labels[i]] <- NA
            }
            
            if (length(contrast_bleed) > 0) {
                contrast_bleedthrough[col_labels[i]] <- median(contrast_bleed, na.rm = TRUE)
            } else {
                contrast_bleedthrough[col_labels[i]] <- NA
            }
        } else {
            reference_bleedthrough[col_labels[i]] <- NA
            contrast_bleedthrough[col_labels[i]] <- NA
        }
    }
    
    # Cap bleedthrough at 0.4 for visualization
    contrast_bleedthrough_capped <- pmin(contrast_bleedthrough, 0.4, na.rm = FALSE)
    
    message("Bleedthrough range: ", 
            round(min(contrast_bleedthrough, na.rm=TRUE), 3), " to ", 
            round(max(contrast_bleedthrough, na.rm=TRUE), 3))
    message("Median contrast bleedthrough: ", round(median(contrast_bleedthrough, na.rm=TRUE), 3))
    message("Median reference bleedthrough: ", round(median(reference_bleedthrough, na.rm=TRUE), 3))
    
    # Add continuous bleedthrough values to annotation (contrast condition)
    annotation_col$Bleedthrough_contrast <- contrast_bleedthrough_capped[rownames(annotation_col)]
    
    # Add reference bleedthrough as text annotation (will be displayed as labels)
    annotation_col$Bleedthrough_ref <- reference_bleedthrough[rownames(annotation_col)]
}

# Label columns with total change counts if available (do this AFTER bleedthrough annotation)
if (!is.null(annotation_col) && nrow(annotation_col) > 0 && "total_changes" %in% colnames(annotation_col)) {
    lab_with_counts <- paste0(col_labels, " (n=", annotation_col[col_labels, "total_changes"], ")")
    names(lab_with_counts) <- col_labels
    # Update both the matrix column names and annotation rownames
    colnames(heatmap_matrix_clamped) <- lab_with_counts
    rownames(annotation_col) <- lab_with_counts
}

# Define annotation colors
annotation_colors <- list()
if ("Bleedthrough_contrast" %in% colnames(annotation_col) || "Bleedthrough_ref" %in% colnames(annotation_col)) {
    # Use continuous gradient from white to orange for bleedthrough annotations (0 to 0.4)
    bleedthrough_gradient <- colorRampPalette(c("white", "#FF8C00"))(100)
    
    # For continuous annotations in pheatmap, we need to map each value to a color
    # Map values from 0-0.4 range to colors (normalizing to 0-100 index)
    if ("Bleedthrough_contrast" %in% colnames(annotation_col)) {
        # Normalize contrast values to 0-1, then to 0-99 (100 colors indexed 0-99)
        contrast_vals <- annotation_col$Bleedthrough_contrast
        contrast_normalized <- pmin(pmax(contrast_vals / 0.4, 0), 1)  # Normalize to 0-1, capped at 0.4
        contrast_indices <- round(contrast_normalized * 99) + 1  # Map to 1-100
        annotation_colors$Bleedthrough_contrast <- bleedthrough_gradient[contrast_indices]
    }
    
    if ("Bleedthrough_ref" %in% colnames(annotation_col)) {
        # Same for reference values
        ref_vals <- annotation_col$Bleedthrough_ref
        ref_normalized <- pmin(pmax(ref_vals / 0.4, 0), 1)
        ref_indices <- round(ref_normalized * 99) + 1
        annotation_colors$Bleedthrough_ref <- bleedthrough_gradient[ref_indices]
    }
}

# Plot heatmap with pheatmap
pheatmap(
    heatmap_matrix_clamped,
    color = color_palette,
    breaks = breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    fontsize = 10,
    border_color = NA,
    fontsize_row = 9,
    fontsize_col = 9,
    cellwidth = 9,
    cellheight = 9,
    angle_col = 90,
    main = "TF Activity Changes Across Comparisons",
    legend = TRUE,
    legend_breaks = seq(color_limits[1], color_limits[2], by = 1),
    legend_labels = as.character(seq(color_limits[1], color_limits[2], by = 1)),
    filename = opt$output
)