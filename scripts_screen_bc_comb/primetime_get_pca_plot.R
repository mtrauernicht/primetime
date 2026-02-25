#!/usr/bin/env Rscript
# ==============================================================================
# PCA plot for TF activity across all conditions
# ==============================================================================
# Generate PCA plot from comparative analysis outputs
# - Reference conditions colored red, contrasts in other colors
# - Label conditions with highest spread
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# Parse arguments
option_list <- list(
  make_option(c("--results"), type="character", 
              help="Comma-separated list of comparison result files"),
  make_option(c("--results-file"), type="character", dest="results_file",
              help="Path to file with newline-separated comparison result files"),
  make_option(c("--bleedthrough"), type="character",
              help="Path to bleedthrough data file"),
  make_option(c("--output"), type="character", 
              help="Output PDF path for PCA plot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ((is.null(opt$results) && is.null(opt$results_file)) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments")
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
safe_read <- function(path) {
  if (!file.exists(path)) {
    message("Skipping missing file: ", path)
    return(NULL)
  }
  if (file.info(path)$size == 0) {
    message("Skipping empty file: ", path)
    return(NULL)
  }
  tryCatch({
    fread(path, header=TRUE, na.strings=c("", "NA"))
  }, error = function(e) {
    message("Failed to read ", path, ": ", conditionMessage(e))
    NULL
  })
}

write_placeholder_pdf <- function(path, message_line) {
  pdf(path, width=8, height=4)
  plot.new()
  text(0.5, 0.6, "PCA of TF Activities Across Conditions", cex=1.2, font=2)
  text(0.5, 0.4, message_line, cex=1)
  dev.off()
  cat(message_line, "\n")
}

# Parse input files (allow either comma-separated arg or newline-separated file)
result_files <- character(0)
if (!is.null(opt$results)) {
  result_files <- c(result_files, unlist(strsplit(opt$results, ",")))
}
if (!is.null(opt$results_file)) {
  if (!file.exists(opt$results_file)) {
    stop("Results list file not found: ", opt$results_file)
  }
  listed <- trimws(readLines(opt$results_file, warn = FALSE))
  listed <- listed[listed != ""]
  result_files <- c(result_files, listed)
}
result_files <- unique(result_files)

# Read and merge all comparison data
all_data <- list()
conditions <- c()
reference_conditions <- c()

for (file in result_files) {
  dt <- safe_read(file)
  if (is.null(dt) || nrow(dt) == 0) {
    next
  }
  dt <- as.data.table(dt)
  
  # Extract conditions from filename: ref_vs_contrast.txt
  basename <- basename(file)
  basename <- gsub(".txt$", "", basename)
  parts <- strsplit(basename, "_vs_")[[1]]
  
  if (length(parts) == 2) {
    ref <- parts[1]
    contrast <- parts[2]
    
    # Check if the condition columns exist in the data
    if (ref %in% names(dt) && contrast %in% names(dt) && "tf" %in% names(dt)) {
      ref_vals <- suppressWarnings(as.numeric(dt[[ref]]))
      contrast_vals <- suppressWarnings(as.numeric(dt[[contrast]]))
      
      if (all(is.na(ref_vals)) || all(is.na(contrast_vals))) {
        next
      }

      ref_data <- data.table(
        tf = dt$tf,
        condition = ref,
        activity = ref_vals
      )
      contrast_data <- data.table(
        tf = dt$tf,
        condition = contrast,
        activity = contrast_vals
      )
      
      all_data[[length(all_data) + 1]] <- ref_data
      all_data[[length(all_data) + 1]] <- contrast_data
      
      conditions <- c(conditions, ref, contrast)
      reference_conditions <- c(reference_conditions, ref)
    }
  }
}

if (length(all_data) == 0) {
  write_placeholder_pdf(opt$output, "No valid comparison files to plot")
  quit(status=0)
}

# Combine all data
combined_dt <- rbindlist(all_data)

if (nrow(combined_dt) == 0) {
  write_placeholder_pdf(opt$output, "No activity values available for PCA")
  quit(status=0)
}

# Get unique conditions and references
unique_conditions <- unique(conditions)
unique_references <- unique(reference_conditions)

# Pivot to wide format for PCA (rows=TFs, cols=conditions)
wide_dt <- dcast(combined_dt, tf ~ condition, value.var="activity", fun.aggregate=mean)

# Remove TFs with missing data
wide_dt <- wide_dt[complete.cases(wide_dt)]

if (nrow(wide_dt) == 0 || ncol(wide_dt) <= 2) {
  write_placeholder_pdf(opt$output, "Insufficient data for PCA")
  quit(status=0)
}

# Extract matrix for PCA
tf_names <- wide_dt$tf
pca_matrix <- as.matrix(wide_dt[, -"tf"])
rownames(pca_matrix) <- tf_names

# Transpose for PCA (want conditions as observations)
pca_matrix_t <- t(pca_matrix)

if (nrow(pca_matrix_t) < 2 || ncol(pca_matrix_t) < 2) {
  write_placeholder_pdf(opt$output, "Need at least 2 conditions and TFs for PCA")
  quit(status=0)
}

# Perform PCA
pca_result <- tryCatch({
  prcomp(pca_matrix_t, scale.=TRUE, center=TRUE)
}, error = function(e) {
  write_placeholder_pdf(opt$output, paste("PCA failed:", conditionMessage(e)))
  quit(status=0)
})

# Extract PC scores
pca_scores <- as.data.table(pca_result$x[, 1:2])
pca_scores$condition <- rownames(pca_result$x)

# Add bleedthrough status to PCA scores
if (!is.null(opt$bleedthrough) && !is.na(opt$bleedthrough) && file.exists(opt$bleedthrough)) {
  bleedthrough_df <- fread(opt$bleedthrough)
  message("Loaded bleedthrough data with ", nrow(bleedthrough_df), " rows")
  
  # For each condition, determine bleedthrough status
  pca_scores$bleedthrough_status <- sapply(pca_scores$condition, function(cond) {
    bleed_val <- bleedthrough_df$bleedthrough[bleedthrough_df$condition == cond]
    if (length(bleed_val) > 0) {
      max_bleed <- max(bleed_val, na.rm = TRUE)
      if (max_bleed > 0.2) {
        return("High")
      } else if (max_bleed > 0.1) {
        return("Medium")
      } else {
        return("Low")
      }
    } else {
      return("Unknown")
    }
  })
  pca_scores$bleedthrough_status <- factor(pca_scores$bleedthrough_status,
                                            levels = c("Low", "Medium", "High", "Unknown"))
  message("Bleedthrough status: Low=", sum(pca_scores$bleedthrough_status == "Low"),
          ", Medium=", sum(pca_scores$bleedthrough_status == "Medium"),
          ", High=", sum(pca_scores$bleedthrough_status == "High"),
          ", Unknown=", sum(pca_scores$bleedthrough_status == "Unknown"))
} else {
  message("No bleedthrough data provided - using reference vs contrast coloring")
  # Default to reference vs contrast if no bleedthrough data
  pca_scores$bleedthrough_status <- factor(
    ifelse(pca_scores$condition %in% unique_references, "Reference", "Contrast"),
    levels = c("Reference", "Contrast")
  )
}

# Calculate spread from centroid
centroid <- c(mean(pca_scores$PC1), mean(pca_scores$PC2))
pca_scores$distance <- sqrt((pca_scores$PC1 - centroid[1])^2 + 
                             (pca_scores$PC2 - centroid[2])^2)

# Identify top conditions by spread (top 25% or at least 3)
n_label <- max(3, ceiling(nrow(pca_scores) * 0.25))
pca_scores$label <- ""
top_indices <- order(pca_scores$distance, decreasing=TRUE)[1:min(n_label, nrow(pca_scores))]
pca_scores$label[top_indices] <- pca_scores$condition[top_indices]

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, ]
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

# Extract and analyze loadings (contributions of each TF to each PC)
loadings_matrix <- pca_result$rotation[, 1:2]
loadings_dt <- as.data.table(loadings_matrix, keep.rownames = "tf")
setnames(loadings_dt, "PC1", "loading_PC1")
setnames(loadings_dt, "PC2", "loading_PC2")

# Add absolute values for ranking
loadings_dt[, abs_loading_PC1 := abs(loading_PC1)]
loadings_dt[, abs_loading_PC2 := abs(loading_PC2)]

# Get top 5 TFs for each PC
top_n <- 5
top_pc1 <- loadings_dt[order(-abs_loading_PC1)][1:top_n]
top_pc2 <- loadings_dt[order(-abs_loading_PC2)][1:top_n]

cat("\n==== PC1 Analysis (", pc1_var, "% variance) ====\n", sep="")
cat("Top ", top_n, " TFs influencing PC1:\n", sep="")
for (i in seq_len(nrow(top_pc1))) {
  direction <- if (top_pc1$loading_PC1[i] > 0) "positive" else "negative"
  cat(sprintf("  %d. %s: %.4f (%s)\n", i, top_pc1$tf[i], top_pc1$loading_PC1[i], direction))
}

cat("\n==== PC2 Analysis (", pc2_var, "% variance) ====\n", sep="")
cat("Top ", top_n, " TFs influencing PC2:\n", sep="")
for (i in seq_len(nrow(top_pc2))) {
  direction <- if (top_pc2$loading_PC2[i] > 0) "positive" else "negative"
  cat(sprintf("  %d. %s: %.4f (%s)\n", i, top_pc2$tf[i], top_pc2$loading_PC2[i], direction))
}

# Define colors (same as bleedthrough plots and heatmap)
corColors <- brewer.pal(n = 7, name = "RdYlBu")[2:6]
bleedthrough_colors <- c(
  "Low" = corColors[5],      # Blue (low bleedthrough, good)
  "Medium" = corColors[2],   # Yellowish (medium bleedthrough, caution)
  "High" = corColors[1],     # Orange/Red (high bleedthrough, bad)
  "Unknown" = "grey80",
  "Reference" = "red",       # Fallback if no bleedthrough data
  "Contrast" = "steelblue"   # Fallback if no bleedthrough data
)

# Create PCA plot
p <- ggplot(pca_scores, aes(x=PC1, y=PC2, color=bleedthrough_status, label=label)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=bleedthrough_colors, name="Bleedthrough") +
  geom_text_repel(
    data=subset(pca_scores, label != ""),
    size=3,
    box.padding=0.5,
    point.padding=0.3,
    max.overlaps=Inf
  ) +
  labs(
    title="PCA of TF Activities Across Conditions",
    x=paste0("PC1 (", pc1_var, "% variance)"),
    y=paste0("PC2 (", pc2_var, "% variance)")
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position="bottom",
    plot.title=element_text(hjust=0.5, face="bold")
  )

# Create loadings plot for PC1
p_pc1_loadings <- ggplot(top_pc1, aes(x=reorder(tf, loading_PC1), y=loading_PC1, fill=loading_PC1)) +
  geom_col() +
  scale_fill_gradient2(low="#1e3a8a", mid="white", high="#c0392b", midpoint=0) +
  coord_flip() +
  labs(title=paste0("Top TFs driving PC1 (", pc1_var, "% variance)"),
       x="TF", y="Loading") +
  theme_bw() +
  theme(legend.position="none", axis.title.y=element_blank())

# Create loadings plot for PC2
p_pc2_loadings <- ggplot(top_pc2, aes(x=reorder(tf, loading_PC2), y=loading_PC2, fill=loading_PC2)) +
  geom_col() +
  scale_fill_gradient2(low="#1e3a8a", mid="white", high="#c0392b", midpoint=0) +
  coord_flip() +
  labs(title=paste0("Top TFs driving PC2 (", pc2_var, "% variance)"),
       x="TF", y="Loading") +
  theme_bw() +
  theme(legend.position="none", axis.title.y=element_blank())

# Combine plots (PCA in center, loadings on sides)
library(patchwork)
p_combined <- p + (p_pc1_loadings / p_pc2_loadings) + plot_layout(widths=c(2, 1))

# Save combined plot
ggsave(opt$output, p_combined, width=14, height=8)

cat("PCA plot saved to:", opt$output, "\n")
cat("Total conditions:", nrow(pca_scores), "\n")
if ("bleedthrough_status" %in% names(pca_scores) && is.factor(pca_scores$bleedthrough_status)) {
  cat("Low bleedthrough:", sum(pca_scores$bleedthrough_status == "Low", na.rm=TRUE), "\n")
  cat("Medium bleedthrough:", sum(pca_scores$bleedthrough_status == "Medium", na.rm=TRUE), "\n")
  cat("High bleedthrough:", sum(pca_scores$bleedthrough_status == "High", na.rm=TRUE), "\n")
  cat("Unknown bleedthrough:", sum(pca_scores$bleedthrough_status == "Unknown", na.rm=TRUE), "\n")
}
