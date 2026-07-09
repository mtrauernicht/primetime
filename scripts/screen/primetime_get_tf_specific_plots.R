suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(dplyr)
    library(stringr)
    library(tidyr)
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
    reference_condition = get_arg("--reference-condition"),
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

find_existing_path <- function(paths) {
    for (candidate in paths) {
        if (file.exists(candidate)) return(candidate)
    }
    return(NA_character_)
}

# Auto-detect design file (for sample -> condition mapping)
design <- NULL
project_root <- dirname(dirname(opt$output_dir))
design_path <- find_existing_path(c(
    file.path(project_root, "tmp_primetime", "design.txt"),
    file.path(dirname(opt$output_dir), "tmp_primetime", "design.txt")
))
if (!is.na(design_path) && file.exists(design_path)) {
    message("Loading design file from: ", design_path)
    design <- tryCatch({
        read.delim(design_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
        warning("Failed to read design file: ", e$message)
        NULL
    })
} else {
    message("Design file not found at ", design_path, " - inferring treatment from replicate names")
}

# Fallback when design is unavailable: infer treatment from replicate names
infer_treatment <- function(x) {
    # Common suffixes: _rep1, _R1, .1
    sub("(_rep[0-9]+|_R[0-9]+|\\.[0-9]+)$", "", x, ignore.case = TRUE)
}

infer_condition_from_replicate <- function(x) {
    sample_name <- trimws(as.character(x))
    sample_name <- sub("_[0-9]+$", "", sample_name)
    sample_name <- sub("(_rep[0-9]+|_R[0-9]+|\\.[0-9]+)$", "", sample_name, ignore.case = TRUE)
    sample_name
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

# Try to auto-detect cDNA counts file (used for replicate comparison plots)
cdna_df <- NULL
pDNA_df <- NULL
sample_well_scores <- NULL
sample_level_activity_df <- NULL
cdna_path <- find_existing_path(c(
    file.path(project_root, "tmp_primetime", "activity", "cDNA_counts.txt"),
    file.path(dirname(opt$output_dir), "tmp_primetime", "activity", "cDNA_counts.txt")
))
if (!is.na(cdna_path) && file.exists(cdna_path)) {
    message("Loading cDNA counts file from: ", cdna_path)
    cdna_df <- tryCatch({
        read.table(cdna_path, header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
        warning("Failed to read cDNA counts: ", e$message)
        NULL
    })
} else {
    message("cDNA counts file not found at ", cdna_path, " - skipping replicate plots")
}

pDNA_path <- find_existing_path(c(
    file.path(project_root, "tmp_primetime", "activity", "pDNA_counts.txt"),
    file.path(dirname(opt$output_dir), "tmp_primetime", "activity", "pDNA_counts.txt")
))
if (!is.na(pDNA_path) && file.exists(pDNA_path)) {
    message("Loading pDNA counts file from: ", pDNA_path)
    pDNA_df <- tryCatch({
        read.table(pDNA_path, header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
        warning("Failed to read pDNA counts: ", e$message)
        NULL
    })
} else {
    message("pDNA counts file not found at ", pDNA_path, " - replicate plots will fall back to RPM")
}

barcode_activity_path <- find_existing_path(c(
    file.path(project_root, "tmp_primetime", "activity", "barcode_activity.txt"),
    file.path(dirname(opt$output_dir), "tmp_primetime", "activity", "barcode_activity.txt")
))
if (!is.na(barcode_activity_path) && file.exists(barcode_activity_path)) {
    message("Loading barcode activity file from: ", barcode_activity_path)
    barcode_activity_df <- tryCatch({
        read.table(barcode_activity_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
        warning("Failed to read barcode activity file: ", e$message)
        NULL
    })

    if (!is.null(barcode_activity_df) && nrow(barcode_activity_df) > 0) {
        if (!"log2_mean_RPM" %in% colnames(barcode_activity_df) && "mean_RPM" %in% colnames(barcode_activity_df)) {
            barcode_activity_df$log2_mean_RPM <- log2(as.numeric(barcode_activity_df$mean_RPM) + 1)
        }

        sample_level_activity <- barcode_activity_df
        if ("negative_control" %in% colnames(sample_level_activity)) {
            sample_level_activity <- sample_level_activity %>% filter(!negative_control)
        }

        sample_level_activity <- sample_level_activity %>%
            group_by(cDNA_sample, tf) %>%
            summarise(sample_log2_mean_RPM = mean(log2_mean_RPM, na.rm = TRUE), .groups = "drop")
        sample_level_activity_df <- sample_level_activity

        control_condition <- if (!is.na(opt$reference_condition) && nzchar(trimws(opt$reference_condition))) {
            opt$reference_condition
        } else {
            inferred_control <- sample_level_activity$cDNA_sample[str_detect(tolower(sample_level_activity$cDNA_sample), "control|ctrl|dmso")][1]
            inferred_control
        }
        if (is.na(control_condition) || control_condition == "") {
            control_condition <- sample_level_activity$cDNA_sample[1]
        }

        control_samples <- sample_level_activity$cDNA_sample[startsWith(tolower(sample_level_activity$cDNA_sample), tolower(control_condition))]
        if (length(control_samples) == 0) {
            control_samples <- sample_level_activity$cDNA_sample[str_detect(tolower(sample_level_activity$cDNA_sample), "control|ctrl|dmso")]
        }
        if (length(control_samples) == 0) {
            control_samples <- sample_level_activity$cDNA_sample[1]
        }

        control_activity <- sample_level_activity %>%
            filter(cDNA_sample %in% control_samples) %>%
            group_by(tf) %>%
            summarise(control_log2_mean_RPM = mean(sample_log2_mean_RPM, na.rm = TRUE), .groups = "drop")

        sample_well_scores <- sample_level_activity %>%
            left_join(control_activity, by = "tf") %>%
            mutate(
                score = sample_log2_mean_RPM - control_log2_mean_RPM,
                score = ifelse(is.na(score), sample_log2_mean_RPM, score)
            )
    }
} else {
    message("barcode_activity.txt not found - plate overview will fall back to comparison logFC values")
}


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

normalize_label <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("(_rep[0-9]+|_r[0-9]+|\\.[0-9]+|-[0-9]+)$", "", x, ignore.case = TRUE)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

infer_plate_id <- function(sample_name) {
    sample_name <- trimws(as.character(sample_name))
    gsub("(_rep[0-9]+|_r[0-9]+|\\.[0-9]+|-[0-9]+)$", "", sample_name, ignore.case = TRUE)
}

build_well_metadata <- function(well_mapping, design_data) {
    well_meta <- data.frame(
        well = names(well_mapping),
        sample = unname(as.character(well_mapping)),
        stringsAsFactors = FALSE
    )
    well_meta$plate_id <- vapply(well_meta$sample, infer_plate_id, character(1))
    well_meta$sample_key <- vapply(well_meta$sample, normalize_label, character(1))
    well_meta$condition <- well_meta$sample

    if (!is.null(design_data) && "replicate" %in% colnames(design_data)) {
        condition_col <- intersect(c("treatment", "condition", "sig"), colnames(design_data))
        if (length(condition_col) > 0) {
            design_lookup <- design_data[, c("replicate", condition_col[1]), drop = FALSE]
            colnames(design_lookup)[2] <- "design_condition"
            match_idx <- match(well_meta$sample, design_lookup$replicate)
            well_meta$condition[!is.na(match_idx)] <- as.character(design_lookup$design_condition[match_idx[!is.na(match_idx)]])
        }
    }

    well_meta$condition_key <- vapply(well_meta$condition, normalize_label, character(1))
    well_meta$plate_key <- vapply(well_meta$plate_id, normalize_label, character(1))
    well_meta
}

match_condition_to_wells <- function(condition, plate_meta) {
    condition_key <- normalize_label(condition)
    if (is.na(condition_key) || condition_key == "") {
        return(character(0))
    }

    wells <- plate_meta$well[plate_meta$condition_key == condition_key]
    if (length(wells) == 0) {
        wells <- plate_meta$well[plate_meta$sample_key == condition_key]
    }
    if (length(wells) == 0) {
        wells <- plate_meta$well[grepl(condition_key, plate_meta$sample_key, fixed = TRUE)]
    }

    unique(wells)
}

resolve_target_condition <- function(comp_parts, reference_condition) {
    if (length(comp_parts) != 2) {
        return(NA_character_)
    }

    ref_key <- normalize_label(reference_condition)
    first_key <- normalize_label(comp_parts[1])
    second_key <- normalize_label(comp_parts[2])

    if (!is.na(ref_key) && ref_key != "") {
        if (first_key == ref_key && second_key != ref_key) {
            return(comp_parts[2])
        }
        if (second_key == ref_key && first_key != ref_key) {
            return(comp_parts[1])
        }
    }

    comp_parts[1]
}

# Function to create 384-well plate visualization for a TF
create_plate_overview <- function(tf_name, all_results, barcode_list, well_mapping, design_data) {
    rows <- LETTERS[1:16]
    cols <- 1:24

    tf_results <- all_results[all_results$tf == tf_name, ]
    if (nrow(tf_results) == 0) {
        return(NULL)
    }

    well_meta <- build_well_metadata(well_mapping, design_data)
    if (nrow(well_meta) == 0) {
        return(NULL)
    }

    use_qc_scores <- FALSE
    qc_score_hits <- 0

    plate_matrix <- matrix(NA_real_, nrow = length(rows), ncol = length(cols), dimnames = list(rows, as.character(cols)))
    tooltip_matrix <- matrix("", nrow = length(rows), ncol = length(cols), dimnames = list(rows, as.character(cols)))
    sample_matrix <- matrix(NA_character_, nrow = length(rows), ncol = length(cols), dimnames = list(rows, as.character(cols)))
    condition_matrix <- matrix(NA_character_, nrow = length(rows), ncol = length(cols), dimnames = list(rows, as.character(cols)))
    sig_matrix <- matrix(NA_character_, nrow = length(rows), ncol = length(cols), dimnames = list(rows, as.character(cols)))

    if (!is.null(sample_well_scores) && nrow(sample_well_scores) > 0) {
        tf_scores <- sample_well_scores %>% filter(tf == tf_name)
        if (nrow(tf_scores) > 0) {
            use_qc_scores <- TRUE
            score_lookup <- setNames(tf_scores$score, normalize_label(tf_scores$cDNA_sample))

            for (well_idx in seq_len(nrow(well_meta))) {
                well <- well_meta$well[well_idx]
                coords <- well_to_coords(well)
                if (is.na(coords$row) || is.na(coords$col)) {
                    next
                }
                sample_name <- well_meta$sample[well_idx]
                score_value <- unname(score_lookup[normalize_label(sample_name)])
                if (length(score_value) == 0) {
                    score_value <- NA_real_
                }
                if (!is.na(score_value)) {
                    qc_score_hits <- qc_score_hits + 1
                }
                plate_matrix[coords$row, coords$col] <- score_value
                tooltip_matrix[coords$row, coords$col] <- paste0(
                    "Well: ", well,
                    "\nSample: ", sample_name,
                    "\nCondition: ", well_meta$condition[well_idx],
                    "\nScore: ", format(round(score_value, 3), nsmall = 3)
                )
            }

            if (qc_score_hits == 0) {
                use_qc_scores <- FALSE
                plate_matrix[,] <- NA_real_
                tooltip_matrix[,] <- ""
            }
        }
    }

    for (well_idx in seq_len(nrow(well_meta))) {
        well <- well_meta$well[well_idx]
        coords <- well_to_coords(well)
        if (is.na(coords$row) || is.na(coords$col)) {
            next
        }
        sample_matrix[coords$row, coords$col] <- well_meta$sample[well_idx]
        condition_matrix[coords$row, coords$col] <- well_meta$condition[well_idx]
        tooltip_matrix[coords$row, coords$col] <- paste0(
            "Well: ", well,
            "<br>Sample: ", well_meta$sample[well_idx],
            "<br>Condition: ", well_meta$condition[well_idx]
        )
    }

    for (comp_idx in seq_len(nrow(tf_results))) {
        result_row <- tf_results[comp_idx, ]
        comparison_id <- result_row$comparison_id
        logFC <- as.numeric(result_row$logFC)
        sig <- if ("sig" %in% colnames(result_row)) result_row$sig else "NotSig"
        comp_parts <- strsplit(comparison_id, "_vs_")[[1]]
        target_cond <- resolve_target_condition(comp_parts, opt$reference_condition)
        if (is.na(target_cond) || !sig %in% c("Upregulated", "Downregulated")) {
            next
        }

        target_wells <- match_condition_to_wells(target_cond, well_meta)
        if (length(target_wells) == 0) {
            next
        }

        for (well in target_wells) {
            coords <- well_to_coords(well)
            if (is.na(coords$row) || is.na(coords$col)) {
                next
            }
            if (use_qc_scores) {
                existing_sig <- sig_matrix[coords$row, coords$col]
                if (is.na(existing_sig) || existing_sig == "" || existing_sig == "NS") {
                    sig_matrix[coords$row, coords$col] <- sig
                    tooltip_matrix[coords$row, coords$col] <- paste0(
                        tooltip_matrix[coords$row, coords$col],
                        "<br>Comparison: ", comparison_id,
                        "<br>Primetime significance: ", sig
                    )
                }
            } else {
                existing_val <- plate_matrix[coords$row, coords$col]
                if (is.na(existing_val) || abs(logFC) > abs(existing_val)) {
                    plate_matrix[coords$row, coords$col] <- logFC
                    sig_matrix[coords$row, coords$col] <- sig
                    tooltip_matrix[coords$row, coords$col] <- paste0(
                        tooltip_matrix[coords$row, coords$col],
                        "<br>Comparison: ", comparison_id,
                        "<br>logFC: ", format(round(logFC, 3), nsmall = 3),
                        "<br>Significance: ", sig
                    )
                }
            }
        }
    }

    n_populated <- sum(!is.na(plate_matrix))
    if (n_populated == 0) {
        return(NULL)
    }

    plate_df <- expand.grid(Row = rows, Col = cols, stringsAsFactors = FALSE)
    plate_df$Row <- factor(plate_df$Row, levels = rows)
    plate_df$Col <- as.numeric(plate_df$Col)
    plate_df$well <- paste0(as.character(plate_df$Row), plate_df$Col)
    plate_df$logFC <- mapply(function(r, c) plate_matrix[r, as.character(c)], as.character(plate_df$Row), plate_df$Col)
    plate_df$sig <- mapply(function(r, c) sig_matrix[r, as.character(c)], as.character(plate_df$Row), plate_df$Col)
    plate_df$sample <- mapply(function(r, c) sample_matrix[r, as.character(c)], as.character(plate_df$Row), plate_df$Col)
    plate_df$condition <- mapply(function(r, c) condition_matrix[r, as.character(c)], as.character(plate_df$Row), plate_df$Col)
    plate_df$tooltip <- mapply(function(r, c) tooltip_matrix[r, as.character(c)], as.character(plate_df$Row), plate_df$Col)
    plate_df$logFC_display <- ifelse(is.na(plate_df$logFC), NA, ifelse(plate_df$sig %in% c("Upregulated", "Downregulated"), plate_df$logFC, 0))

    plate_title <- paste0("384-Well Plate Overview: ", tf_name)

    static_plot <- ggplot(plate_df, aes(x = Col, y = Row, fill = logFC_display)) +
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
    static_plot
}

# Precompute replicate plotting dataset once (major speedup vs rebuilding per TF)
replicate_by_tf <- NULL
if (!is.null(cdna_df)) {
    ref_cond <- ifelse(is.na(opt$reference_condition) || opt$reference_condition == "", "Control", opt$reference_condition)

    # Prefer QC-derived barcode activity, which is already cDNA_RPM/pDNA_RPM based.
    if (!is.null(barcode_activity_df) && nrow(barcode_activity_df) > 0 && all(c("cDNA_sample", "tf") %in% colnames(barcode_activity_df))) {
        merged <- barcode_activity_df
        merged$replicate <- as.character(merged$cDNA_sample)
        if ("log2_mean_RPM" %in% colnames(merged)) {
            merged$activity_log2 <- as.numeric(merged$log2_mean_RPM)
        } else if ("mean_RPM" %in% colnames(merged)) {
            merged$activity_log2 <- log2(as.numeric(merged$mean_RPM) + 1)
        } else {
            merged$activity_log2 <- NA_real_
        }

        # Merge with design if available; otherwise infer treatment from replicate names
        if (!is.null(design) && "replicate" %in% colnames(design)) {
            merged <- merge(merged, design, by.x = "replicate", by.y = "replicate", all.x = TRUE)
        } else {
            merged$treatment <- infer_treatment(merged$replicate)
        }

        merged$treatment_col <- if ("treatment" %in% colnames(merged)) {
            ifelse(merged$treatment == ref_cond, "Control condition", "Unknown")
        } else {
            "Unknown"
        }

        merged$condition <- infer_condition_from_replicate(merged$replicate)
        merged$condition_key <- toupper(merged$condition)
        replicate_by_tf <- split(as.data.frame(merged), merged$tf)
    } else {

        # Prepare pDNA values per barcode for log2(cDNA_RPM/pDNA_RPM)
        pDNA_mean <- NULL
        if (!is.null(pDNA_df) && "barcode" %in% colnames(pDNA_df)) {
            pDNA_counts <- pDNA_df
            pDNA_num_cols <- setdiff(colnames(pDNA_counts), "barcode")
            for (col in pDNA_num_cols) {
                pDNA_counts[[col]] <- as.numeric(pDNA_counts[[col]])
                total_col <- sum(pDNA_counts[[col]], na.rm = TRUE)
                pDNA_counts[[col]] <- if (is.finite(total_col) && total_col > 0) {
                    pDNA_counts[[col]] / total_col * 1e6
                } else {
                    NA_real_
                }
            }
            if (length(pDNA_num_cols) > 0) {
                pDNA_counts$pDNA_mean <- rowMeans(pDNA_counts[, pDNA_num_cols, drop = FALSE], na.rm = TRUE)
            } else {
                pDNA_counts$pDNA_mean <- NA_real_
            }
            pDNA_mean <- pDNA_counts[, c("barcode", "pDNA_mean")]
        }

        # Prepare counts: remove optional columns if present
        counts <- cdna_df
        drop_cols <- intersect(c("negative_control", "promoter"), colnames(counts))
        if (length(drop_cols) > 0) counts <- counts[, setdiff(colnames(counts), drop_cols), drop = FALSE]

        # Convert replicate columns to RPM (depth-normalized)
        num_cols <- setdiff(colnames(counts), c("tf", "negative_control", "promoter", "barcode"))
        for (col in num_cols) {
            counts[[col]] <- as.numeric(counts[[col]])
            total_col <- sum(counts[[col]], na.rm = TRUE)
            counts[[col]] <- if (is.finite(total_col) && total_col > 0) {
                counts[[col]] / total_col * 1e6
            } else {
                NA_real_
            }
        }

        # Melt to long format using data.table
        counts_dt <- data.table::as.data.table(counts)
        long_dt <- data.table::melt(counts_dt, id.vars = c("barcode", "tf", intersect(c("negative_control", "promoter"), colnames(counts_dt))), variable.name = "replicate", value.name = "cDNA_rpm")
        long_df <- as.data.frame(long_dt)

        if (!is.null(pDNA_mean)) {
            long_df <- merge(long_df, pDNA_mean, by = "barcode", all.x = TRUE)
            long_df$activity_log2 <- log2((as.numeric(long_df$cDNA_rpm) + 1) / (as.numeric(long_df$pDNA_mean) + 1))
        } else {
            long_df$activity_log2 <- log2(as.numeric(long_df$cDNA_rpm) + 1)
        }

        # Merge with design if available; otherwise infer treatment from replicate names
        if (!is.null(design) && "replicate" %in% colnames(design)) {
            merged <- merge(long_df, design, by.x = "replicate", by.y = "replicate", all.x = TRUE)
        } else {
            merged <- long_df
            merged$treatment <- infer_treatment(merged$replicate)
        }
        if ("pDNA" %in% colnames(merged)) {
            merged <- merged[merged$pDNA != "True", ]
        }

        if ("treatment" %in% colnames(merged) && "sig" %in% colnames(merged)) {
            merged$treatment_col <- ifelse(merged$treatment == ref_cond, "Control condition", merged$sig)
        } else if ("treatment" %in% colnames(merged)) {
            merged$treatment_col <- merged$treatment
        } else {
            merged$treatment_col <- "Unknown"
        }

        merged$condition <- infer_condition_from_replicate(merged$replicate)
        merged$condition_key <- toupper(merged$condition)
        replicate_by_tf <- split(merged, merged$tf)
    }
}

# Process each TF
message("Generating plots for each TF...")
plot_count <- 0

for (tf in all_tfs) {
    tryCatch({
        tf_results <- all_results[all_results$tf == tf, , drop = FALSE]

        # Create plate overview
        plate_plots <- create_plate_overview(tf, all_results, all_barcodes, well_map, design)

        if (!is.null(plate_plots)) {
            output_file <- file.path(opt$output_dir, paste0(tf, "_plate_overview.pdf"))

            # Save a static PDF only
            pdf(output_file, width = 12, height = 8, useDingbats = FALSE)
            print(plate_plots)
            dev.off()

            message("  ✓ ", tf, " - plate overview")
            plot_count <- plot_count + 1
        } else {
            message("  - ", tf, " - no plate data")
        }

        # Generate per-TF replicate comparison plot if precomputed replicate data is available
        if (!is.null(replicate_by_tf)) {
            tf_df <- replicate_by_tf[[tf]]
            if (!is.null(tf_df) && nrow(tf_df) > 0) {
                ref_key <- "dmso"
                ref_key_upper <- "DMSO"

                other_conditions <- unique(tf_df$condition[tf_df$condition_key != ref_key_upper & !is.na(tf_df$condition_key)])
                if (length(other_conditions) > 0) {
                    out_file <- file.path(opt$output_dir, paste0(tf, "_replicates.pdf"))
                    pdf(out_file, width = 8, height = 6, useDingbats = FALSE)

                    for (condition_name in other_conditions) {
                        plot_df <- tf_df[tf_df$condition_key %in% c(ref_key_upper, toupper(condition_name)), ]
                        if (nrow(plot_df) == 0) {
                            next
                        }

                        cond_key <- normalize_label(condition_name)

                        comp_matches <- vapply(tf_results$comparison_id, function(cid) {
                            parts <- strsplit(as.character(cid), "_vs_", fixed = TRUE)[[1]]
                            if (length(parts) != 2) {
                                return(FALSE)
                            }
                            part_keys <- vapply(parts, normalize_label, character(1))
                            has_condition <- cond_key %in% part_keys
                            has_reference <- any(part_keys %in% c(ref_key, "dmso", "control", "ctrl"))
                            has_condition && has_reference
                        }, logical(1))

                        sign <- "NS"
                        if (any(comp_matches, na.rm = TRUE)) {
                            candidate_rows <- tf_results[comp_matches, , drop = FALSE]
                            if ("logFC" %in% colnames(candidate_rows)) {
                                candidate_rows$logFC_abs <- abs(as.numeric(candidate_rows$logFC))
                                candidate_rows <- candidate_rows[order(candidate_rows$logFC_abs, decreasing = TRUE), , drop = FALSE]
                            }
                            candidate_sig <- as.character(candidate_rows$sig)
                            candidate_sig <- candidate_sig[!is.na(candidate_sig)]
                            if (length(candidate_sig) > 0) {
                                sign <- candidate_sig[1]
                            }
                        }

                        if (!sign %in% c("NS", "Downregulated", "Upregulated", "Control condition")) {
                            sign <- "Unknown"
                        }

                        plot_df$treatment_col <- ifelse(plot_df$condition_key == ref_key_upper, "Control condition", sign)


                        p_rep <- ggplot(plot_df, aes(x = replicate, y = activity_log2)) +
                            geom_point(aes(color = treatment_col)) +
                            scale_color_manual(values = c("Control condition" = "#000000", "NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80", "Unknown" = "grey")) +
                            labs(title = paste0("Replicate activity: ", tf, " - ", condition_name), x = "Replicate", y = "Log2(cDNA/pDNA)") +
                            theme_bw() +
                            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

                        # Add violin layer if treatment_col exists
                        p_rep <- p_rep + ggplot2::geom_violin(aes(fill = treatment_col, color = treatment_col), alpha = 0.4, width = 1) +
                            scale_fill_manual(values = c("Control condition" = "#00000033", "NS" = "#e4e4e455", "Downregulated" = "#6495ed33", "Upregulated" = "#f37f8033", "Unknown" = "#e4e4e455"))

                        print(p_rep)
                    }

                    dev.off()
                    message("  ✓ ", tf, " - replicate plot")
                    plot_count <- plot_count + 1
                }
            }
        }
    }, error = function(e) {
        message("  ✗ ", tf, " - Error: ", e$message)
    })
}

message("")
message("Generated ", plot_count, " plots in: ", opt$output_dir)


