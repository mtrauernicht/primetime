suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0 || idx == length(args)) return(NA_character_)
    args[idx + 1]
}

opt <- list(
    results_file = get_arg("--results-file"),
    barcode_dir = get_arg("--barcode-dir"),
    design = get_arg("--design"),
    well_map = get_arg("--well-map"),
    output_dir = get_arg("--output-dir")
)

# Validate inputs
required_args <- c("results_file", "barcode_dir", "well_map", "output_dir")
missing_args <- required_args[sapply(opt[required_args], is.na)]
if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

# Load results files - consolidate all into single dataframe with comparison info
message("Loading results files...")
comparison_files <- trimws(readLines(opt$results_file))
comparison_files <- comparison_files[comparison_files != ""]

all_results <- NULL
for (file in comparison_files) {
    if (!file.exists(file)) {
        warning("File not found: ", file)
        next
    }
    df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    comparison_name <- basename(file)
    comparison_name <- sub("\\.txt$", "", comparison_name)
    
    df$comparison_id <- comparison_name
    
    if (is.null(all_results)) {
        all_results <- df
    } else {
        # Bind only common columns
        common_cols <- intersect(colnames(all_results), colnames(df))
        all_results <- rbind(all_results[, common_cols], df[, common_cols])
    }
}

if (is.null(all_results) || nrow(all_results) == 0) {
    stop("No valid results files loaded")
}

message("Loaded results for ", length(unique(all_results$comparison_id)), " comparisons")
message("Total TF-comparison pairs: ", nrow(all_results))

# Load design file if provided (for sample -> condition mapping)
design <- NULL
if (!is.na(opt$design) && file.exists(opt$design)) {
    message("Loading design file...")
    design <- read.delim(opt$design, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
    message("No design file provided or file not found - will extract sample info from filenames")
}

# Load well map (maps well position -> sample name)
message("Loading well map...")
well_map_data <- read.csv(opt$well_map, header = FALSE, stringsAsFactors = FALSE, col.names = c("well", "sample"))
well_map <- setNames(well_map_data$sample, well_map_data$well)

# Load barcode annotation files
message("Loading barcode annotations...")
barcode_files <- list.files(opt$barcode_dir, pattern = "\\.cluster\\.annotated\\.txt$", full.names = TRUE)

all_barcodes <- list()
for (file in barcode_files) {
    sample_id <- basename(file)
    sample_id <- sub("\\.cluster\\.annotated\\.txt$", "", sample_id)
    
    tryCatch({
        df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        all_barcodes[[sample_id]] <- df
    }, error = function(e) {
        warning("Failed to load barcode file: ", file)
    })
}

message("Loaded barcode annotations for ", length(all_barcodes), " samples")

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Get unique TFs from results (exclude random controls)
all_tfs <- unique(all_results$tf)
if ("tf" %in% colnames(all_results)) {
    all_tfs <- all_tfs[!is.na(all_tfs)]
    all_tfs <- all_tfs[!grepl("^RANDOM", all_tfs)]
    all_tfs <- sort(all_tfs)
}

message("Found ", length(all_tfs), " unique TFs to process")

# Function to convert well position to plate coordinates
well_to_coords <- function(well) {
    # 384-well plate: rows A-P (16), columns 1-24
    row_char <- substr(well, 1, 1)
    col_num <- as.numeric(substr(well, 2, nchar(well)))
    
    rows <- LETTERS[1:16]
    row_idx <- match(row_char, rows)
    
    if (is.na(row_idx) || is.na(col_num) || col_num < 1 || col_num > 24) {
        return(list(row = NA, col = NA))
    }
    
    return(list(row = row_idx, col = col_num))
}

# Function to create 384-well plate visualization for a TF (multi-page PDF, one per plate)
create_plate_overview <- function(tf_name, all_results, barcode_list, well_mapping, design_data) {
    # 384-well plate layout: 16 rows (A-P), 24 columns (1-24)
    rows <- LETTERS[1:16]
    cols <- 1:24
    
    # Get all comparisons for this TF
    tf_results <- all_results[all_results$tf == tf_name, ]
    
    if (nrow(tf_results) == 0) {
        return(NULL)
    }
    
    # Create reverse lookup: sample_name -> wells
    sample_to_wells <- split(names(well_mapping), well_mapping)
    
    # Detect if multiple plates exist by checking if well_mapping has plate prefixes
    # Format could be: well -> "plate1_sample" or just well -> "sample"
    sample_names <- unique(well_mapping)
    
    # Check if samples have format "plateX_..." (where X is plate ID)
    has_plate_prefix <- any(grepl("^[^_]+_", sample_names))
    
    # Extract plate IDs if they exist
    if (has_plate_prefix) {
        plate_ids <- unique(sapply(strsplit(sample_names, "_"), function(x) x[1]))
        # Filter to keep only actual plate identifiers (not just first part of sample name)
        # Assume plate IDs are short alphanumeric strings
        plate_ids <- plate_ids[grepl("^(plate|[0-9]+|[A-Z][0-9]+)$", plate_ids, ignore.case = TRUE)]
        if (length(plate_ids) == 0) {
            plate_ids <- c("default")  # Fall back to single plate
        }
    } else {
        plate_ids <- c("default")  # Single plate
    }
    
    # Create a plot for each plate
    plot_list <- list()
    
    for (plate_id in plate_ids) {
        # Initialize matrices for this plate
        plate_matrix <- matrix(NA, nrow = 16, ncol = 24)
        rownames(plate_matrix) <- rows
        colnames(plate_matrix) <- as.character(cols)
        
        significance_matrix <- matrix(NA_character_, nrow = 16, ncol = 24)
        rownames(significance_matrix) <- rows
        colnames(significance_matrix) <- as.character(cols)
        
        # Filter sample_to_wells for this plate if multiple plates exist
        if (plate_id != "default") {
            # Only include samples that belong to this plate
            plate_sample_to_wells <- lapply(names(sample_to_wells), function(sample_name) {
                if (startsWith(sample_name, paste0(plate_id, "_"))) {
                    # Remove plate prefix from sample name for matching with comparisons
                    base_name <- sub(paste0("^", plate_id, "_"), "", sample_name)
                    list(base_name = base_name, wells = sample_to_wells[[sample_name]])
                } else {
                    NULL
                }
            })
            plate_sample_to_wells <- Filter(Negate(is.null), plate_sample_to_wells)
            
            # Create lookup: base_name -> wells for this plate
            current_sample_to_wells <- setNames(
                lapply(plate_sample_to_wells, function(x) x$wells),
                sapply(plate_sample_to_wells, function(x) x$base_name)
            )
        } else {
            # Single plate or no prefix - use all samples
            current_sample_to_wells <- sample_to_wells
        }
        
        # For each result, map the contrast condition to its well(s)
        # For each result, map the contrast condition to its well(s)
        for (comp_idx in 1:nrow(tf_results)) {
            result_row <- tf_results[comp_idx, ]
            comparison_id <- result_row$comparison_id
            logFC <- as.numeric(result_row$logFC)
            sig <- if ("sig" %in% colnames(result_row)) result_row$sig else "NotSig"
            
            # Extract reference and contrast from comparison_id (format: "ref_vs_contrast")
            comp_parts <- strsplit(comparison_id, "_vs_")[[1]]
            if (length(comp_parts) != 2) next
            
            reference_cond <- comp_parts[1]
            contrast_cond <- comp_parts[2]
            
            # Find wells that correspond to the contrast condition
            # The logFC is for this TF in contrast vs reference
            contrast_wells <- current_sample_to_wells[[contrast_cond]]
            
            if (!is.null(contrast_wells) && length(contrast_wells) > 0) {
                for (well in contrast_wells) {
                    coords <- well_to_coords(well)
                    if (!is.na(coords$row) && !is.na(coords$col)) {
                        # Only overwrite when current result is significant; keep strongest abs(logFC)
                        if (sig %in% c("Upregulated", "Downregulated")) {
                            existing_sig <- significance_matrix[coords$row, coords$col]
                            existing_val <- plate_matrix[coords$row, coords$col]
                            if (is.na(existing_sig) || abs(logFC) > abs(existing_val)) {
                                plate_matrix[coords$row, coords$col] <- logFC
                                significance_matrix[coords$row, coords$col] <- sig
                            }
                        }
                    }
                }
            }
        }
        
        # Check if we have any data to plot for this plate
        n_populated <- sum(!is.na(plate_matrix))
        
        if (n_populated == 0) {
            next  # Skip this plate if no data
        }
        
        # Create the heatmap visualization for this plate
        plate_df <- as.data.frame(as.table(plate_matrix))
        colnames(plate_df) <- c("Row", "Col", "logFC")
        plate_df$Row <- factor(plate_df$Row, levels = rows)
        plate_df$Col <- as.numeric(as.character(plate_df$Col))
        plate_df$logFC <- as.numeric(as.character(plate_df$logFC))  # Ensure numeric
        
        # Add significance information
        sig_df <- as.data.frame(as.table(significance_matrix))
        colnames(sig_df) <- c("Row", "Col", "sig")
        plate_df$sig <- sig_df$sig
        
        # Color significant values with logFC; show non-significant as 0 (white), missing as NA (gray)
        plate_df$logFC_display <- ifelse(is.na(plate_df$logFC), NA,
                                          ifelse(plate_df$sig %in% c("Upregulated", "Downregulated"),
                                                 plate_df$logFC, 0))
        
        plate_title <- if (plate_id == "default") {
            paste0("384-Well Plate Overview: ", tf_name)
        } else {
            paste0("384-Well Plate Overview: ", tf_name, " (Plate: ", plate_id, ")")
        }
        
        p <- ggplot(plate_df, aes(x = Col, y = Row, fill = logFC_display)) +
            geom_tile(color = "white", size = 0.2) +
            scale_fill_gradient2(
                low = "#377EB8", mid = "white", high = "#E41A1C",
                limits = c(-2, 2),
                oob = scales::squish,
                na.value = "#F0F0F0", breaks = scales::pretty_breaks(n = 5)
            ) +
            scale_x_continuous(breaks = cols, expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            coord_fixed() +
            labs(
                title = plate_title,
                x = "Column", y = "Row", fill = "logFC\n(sig only)"
            ) +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 7),
                plot.title = element_text(size = 11, face = "bold"),
                legend.position = "bottom"
            )
        
        plot_list[[plate_id]] <- p
    }
    
    # Return list of plots (one per plate) or NULL if no plots
    if (length(plot_list) == 0) {
        return(NULL)
    } else {
        return(plot_list)
    }
}

# Function to create lollipop plot for a TF
create_tf_lollipop <- function(tf_name, all_results_df) {
    # Extract data for this TF
    tf_results <- all_results_df[all_results_df$tf == tf_name, ]
    
    if (nrow(tf_results) == 0) {
        return(NULL)
    }
    
    # Create a formatted label for each comparison
    tf_results$comparison_label <- sapply(tf_results$comparison_id, function(comp_id) {
        parts <- strsplit(comp_id, "_vs_")[[1]]
        if (length(parts) == 2) {
            parts[2]  # Just the contrast condition name
        } else {
            comp_id
        }
    })
    
    # Get reference and contrast values if they exist in the results
    tf_results$reference_activity <- NA
    tf_results$contrast_activity <- NA
    
    # Extract reference and contrast condition names
    for (i in 1:nrow(tf_results)) {
        comp_parts <- strsplit(tf_results$comparison_id[i], "_vs_")[[1]]
        if (length(comp_parts) == 2) {
            ref_col <- comp_parts[1]
            contrast_col <- comp_parts[2]
            
            # Look for columns matching these names (the actual activity values)
            if (ref_col %in% colnames(tf_results)) {
                tf_results$reference_activity[i] <- as.numeric(tf_results[i, ref_col])
            }
            if (contrast_col %in% colnames(tf_results)) {
                tf_results$contrast_activity[i] <- as.numeric(tf_results[i, contrast_col])
            }
        }
    }
    
    # Ensure numeric
    tf_results$logFC <- as.numeric(tf_results$logFC)
    
    # Handle significance column
    if (!("sig" %in% colnames(tf_results))) {
        tf_results$sig <- "NS"
    }
    tf_results$sig[is.na(tf_results$sig)] <- "NS"
    
    # Set color for axis text
    tf_results$color_axis <- ifelse(tf_results$sig == "Upregulated", "#f37f80",
                                     ifelse(tf_results$sig == "NS", "gray30", "#6495ed"))
    
    # Order by logFC
    tf_results <- tf_results[order(tf_results$logFC, decreasing = TRUE), ]
    tf_results$comparison_label <- factor(tf_results$comparison_label, levels = unique(tf_results$comparison_label))
    
    # Create lollipop plot matching the style of per-comparison plots
    if (any(!is.na(tf_results$reference_activity)) && any(!is.na(tf_results$contrast_activity))) {
        # Plot with both reference and contrast points
        p <- ggplot(tf_results) +
            geom_segment(aes(
                y = reference_activity,
                yend = contrast_activity,
                x = comparison_label, xend = comparison_label,
                color = sig
            ), size = 1) +
            scale_color_manual(
                values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80"),
                name = "Significance"
            ) +
            geom_point(aes(x = comparison_label, y = reference_activity), color = "black", size = 3) +
            geom_point(aes(x = comparison_label, y = contrast_activity, color = sig), size = 3) +
            theme_bw() +
            theme(
                axis.text.x = element_text(
                    angle = 90, hjust = 1, vjust = 0.5,
                    color = tf_results$color_axis,
                    size = 7
                ),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                text = element_text(size = 12),
                legend.position = "right"
            ) +
            labs(
                y = "Activity (log2(RPM+1))",
                x = "",
                title = paste0("TF Activity: ", tf_name)
            )
    } else {
        # Simplified plot with just logFC
        p <- ggplot(tf_results, aes(x = comparison_label, y = logFC, color = sig)) +
            geom_segment(aes(xend = comparison_label, y = 0, yend = logFC), size = 1) +
            geom_point(size = 3) +
            scale_color_manual(
                values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80"),
                name = "Significance"
            ) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
            theme_bw() +
            theme(
                axis.text.x = element_text(
                    angle = 90, hjust = 1, vjust = 0.5,
                    color = tf_results$color_axis,
                    size = 8
                ),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                text = element_text(size = 12),
                legend.position = "right"
            ) +
            labs(
                y = "logFC",
                x = "",
                title = paste0("TF Activity: ", tf_name)
            )
    }
    
    return(p)
}

# Process each TF
message("Generating plots for each TF...")
plot_count <- 0

for (tf in all_tfs) {
    tryCatch({
        # Create lollipop plot
        lollipop_plot <- create_tf_lollipop(tf, all_results)
        
        if (!is.null(lollipop_plot)) {
            output_file <- file.path(opt$output_dir, paste0(tf, "_lollipop.pdf"))
            # Calculate width based on number of comparisons (0.35 inches per comparison, min 10, max 30)
            tf_data <- all_results[all_results$tf == tf, ]
            n_comparisons <- nrow(tf_data)
            plot_width <- max(10, min(30, n_comparisons * 0.35))
            ggsave(output_file, lollipop_plot, width = plot_width, height = 6, useDingbats = FALSE)
            message("  ✓ ", tf, " - lollipop plot")
            plot_count <- plot_count + 1
        }
        
        # Create plate overview
        plate_plots <- create_plate_overview(tf, all_results, all_barcodes, well_map, design)
        
        if (!is.null(plate_plots) && length(plate_plots) > 0) {
            output_file <- file.path(opt$output_dir, paste0(tf, "_plate_overview.pdf"))
            
            # Save as multi-page PDF if multiple plates, single page otherwise
            pdf(output_file, width = 12, height = 8, useDingbats = FALSE)
            for (plate_id in names(plate_plots)) {
                print(plate_plots[[plate_id]])
            }
            dev.off()
            
            n_plates <- length(plate_plots)
            message("  ✓ ", tf, " - plate overview (", n_plates, " plate", 
                    ifelse(n_plates > 1, "s", ""), ")")
            plot_count <- plot_count + n_plates
        } else {
            message("  - ", tf, " - no plate data")
        }
    }, error = function(e) {
        message("  ✗ ", tf, " - Error: ", e$message)
    })
}

message("")
message("Generated ", plot_count, " plots in: ", opt$output_dir)


