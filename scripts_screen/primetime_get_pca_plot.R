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
})

# Parse arguments
option_list <- list(
  make_option(c("--results"), type="character", 
              help="Comma-separated list of comparison result files"),
  make_option(c("--output"), type="character", 
              help="Output PDF path for PCA plot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$results) || is.null(opt$output)) {
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

# Parse input files
result_files <- unlist(strsplit(opt$results, ","))

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
    
    # Extract TF activity columns (assume corrected_activity columns exist)
    activity_cols <- grep("corrected_activity", names(dt), value=TRUE)
    
    if (length(activity_cols) >= 2 && "tf" %in% names(dt)) {
      # Assume first is reference, second is contrast
      ref_col <- activity_cols[1]
      contrast_col <- activity_cols[2]
      
      ref_vals <- suppressWarnings(as.numeric(dt[[ref_col]]))
      contrast_vals <- suppressWarnings(as.numeric(dt[[contrast_col]]))
      
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

# Add color indicator for reference vs contrast
pca_scores$type <- ifelse(pca_scores$condition %in% unique_references, 
                           "Reference", "Contrast")

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

# Create PCA plot
p <- ggplot(pca_scores, aes(x=PC1, y=PC2, color=type, label=label)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("Reference"="red", "Contrast"="steelblue")) +
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
    y=paste0("PC2 (", pc2_var, "% variance)"),
    color="Condition Type"
  ) +
  theme_bw() +
  theme(
    legend.position="bottom",
    plot.title=element_text(hjust=0.5, face="bold")
  )

# Save plot
ggsave(opt$output, p, width=10, height=8)

cat("PCA plot saved to:", opt$output, "\n")
cat("Total conditions:", nrow(pca_scores), "\n")
cat("Reference conditions:", sum(pca_scores$type == "Reference"), "\n")
cat("Contrast conditions:", sum(pca_scores$type == "Contrast"), "\n")
