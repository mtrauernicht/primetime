# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script runs the BCalm analysis for the TF reporter data.
# It reads the barcode counts, processes the data, and performs differential
# TF activity analysis comparing different conditions of the samples.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================

suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(reshape2)
        library(BiocParallel)
        library(optparse)
        library(tidyr)
        library(ggrepel)
        library(ggnewscale)
})

options(dplyr.width = Inf)

if (!require(BCalm, quietly = TRUE)){
        cat("--------------------------- Installing BCalm\n\n")
        remotes::install_github("kircherlab/BCalm")
        suppressPackageStartupMessages(library(BCalm))
} else {
        suppressPackageStartupMessages(library(BCalm))
}

##########################################################################################
# Read arguments #########################################################################
##########################################################################################

option_list <- list(
        make_option(c("-p", "--pdna"), type = "character", default = NULL, help = "Path to the pDNA counts file"),
        make_option(c("-c", "--cdna"), type = "character", default = NULL, help = "Path to the cDNA counts file"),
        make_option(c("--cdna_baseline"), type = "character", default = NULL, help = "Path to uncorrected cDNA counts used to estimate banana correction independently of other corrections"),
        make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to the output directory"),
        make_option(c("-t", "--threads"), type = "integer", default = 1, help = "Number of threads to use"),
        make_option(c("-q", "--pval_threshold"), type = "numeric", default = 0.05, help = "P-value threshold for significance"),
        make_option(c("--contrast_condition"), type = "character", default = NULL, help = "Condition to contrast"),
        make_option(c("--reference_condition"), type = "character", default = NULL, help = "Condition to contrast"),
        make_option(c("--contrast_replicates"), type = "character", default = NULL, help = "Comma-separated replicate column names for the contrast condition"),
        make_option(c("--reference_replicates"), type = "character", default = NULL, help = "Comma-separated replicate column names for the reference condition"),
        make_option(c("--plot_output"), type = "character", default = NULL, help = "Path to the output directory for the plots"),
        make_option(c("--num_replicates_contrast"), type = "integer", default = 1, help = "Number of replicates for the contrast condition"),
        make_option(c("--num_replicates_reference"), type = "integer", default = 1, help = "Number of replicates for the reference condition"),
        make_option(c("--single_model"), type = "logical", default = FALSE, help = "Whether to make a single model for all the TFs or not"),
        make_option(c("--normalize"), type = "logical", default = TRUE, help = "Whether to normalize counts per promoter and replicate using the RANDOM reporters (count-level)"),
        make_option(c("--shrinkage"), type = "logical", default = FALSE, help = "Add empirical-Bayes shrunken logFC (computed after BCalm from logFC and its SE) to the output and show shrunk and unshrunk values in the plots"),
        make_option(c("--plate_ids"), type = "character", default = "", help = "Comma-separated sample plate IDs (PLATES_FILE) used to assign replicates to plates for the per-plate banana anchors"),
        make_option(c("--banana_correction"), type = "logical", default = TRUE, help = "Whether to apply the banana-shaped bias correction to contrast vs reference barcode counts"),
        make_option(c("--banana_reference_file"), type = "character", default = NULL, help = "Optional TSV (column 'replicate') listing the low-bleedthrough control replicates used as the banana-correction anchor. If missing, the reference replicates are used")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

find_existing_path <- function(paths) {
        for (candidate in paths) {
                if (!is.na(candidate) && file.exists(candidate)) {
                        return(candidate)
                }
        }
        return(NA_character_)
}

find_project_root <- function(start_path) {
        current <- normalizePath(start_path)
        if (file.exists(current) && !dir.exists(current)) {
                current <- dirname(current)
        }

        repeat {
                if (file.exists(file.path(current, "misc", "tf_functions.tsv"))) {
                        return(current)
                }
                parent <- dirname(current)
                if (identical(parent, current)) {
                        break
                }
                current <- parent
        }

        return(NA_character_)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else getwd()
project_root <- find_project_root(script_path)
if (is.na(project_root)) {
        project_root <- find_project_root(getwd())
}


# Extract "replicate" and "barcode suffix" (bcN) from a this_cdna_matrix/this_pdna_matrix
# column name of the form "<replicate>_bc<N>".
split_replicate_bc_colname <- function(col_names) {
        bc_suffix <- sub("^.*_(bc[0-9]+)$", "\\1", col_names)
        replicate_name <- sub("_bc[0-9]+$", "", col_names)
        data.frame(col_name = col_names, replicate = replicate_name, bc = bc_suffix, stringsAsFactors = FALSE)
}

# Normalize counts within each promoter and replicate using that promoter's
# negative-control (RANDOM) reporters as the reference.
normalize_counts_by_random <- function(cdna_matrix, pdna_matrix, negative_control_flag, promoter_vector, target = 1e7) {
        if (!identical(dim(cdna_matrix), dim(pdna_matrix)) ||
            !identical(rownames(cdna_matrix), rownames(pdna_matrix))) {
                stop("cDNA and pDNA matrices must have matching dimensions and row order")
        }
        if (length(negative_control_flag) != nrow(cdna_matrix) ||
            length(promoter_vector) != nrow(cdna_matrix)) {
                stop("Promoter and negative-control metadata must match the count matrix rows")
        }
        if (length(target) != 1 || !is.finite(target) || target <= 0) {
                stop("Normalization target must be a single finite positive number")
        }

        promoter_values <- as.character(promoter_vector)
        if (anyNA(promoter_values) || any(!nzchar(promoter_values))) {
                stop("Promoter metadata contains missing or empty values")
        }
        if (anyNA(negative_control_flag)) {
                stop("Negative-control metadata contains missing values")
        }
        if (is.logical(negative_control_flag)) {
                is_random <- negative_control_flag
        } else {
                is_random <- toupper(trimws(as.character(negative_control_flag))) %in% c("TRUE", "T", "1", "YES", "Y")
        }

        normalize_matrix <- function(count_matrix, matrix_name) {
                normalized <- count_matrix * 1.0
                warnings <- character(0)

                for (promoter in unique(promoter_values)) {
                        promoter_rows <- which(promoter_values == promoter)
                        random_rows <- promoter_rows[is_random[promoter_rows]]
                        if (length(random_rows) == 0) {
                                warnings <- c(warnings, paste0("no RANDOM reporters for promoter '", promoter, "'"))
                                next
                        }

                        random_sums <- colSums(count_matrix[random_rows, , drop = FALSE], na.rm = TRUE)
                        valid_sums <- is.finite(random_sums) & random_sums > 0
                        scale_factors <- rep(1, ncol(count_matrix))
                        scale_factors[valid_sums] <- target / random_sums[valid_sums]
                        normalized[promoter_rows, ] <- sweep(
                                count_matrix[promoter_rows, , drop = FALSE],
                                2,
                                scale_factors,
                                FUN = "*"
                        )

                        if (any(!valid_sums)) {
                                warnings <- c(
                                        warnings,
                                        paste0(
                                                "zero or non-finite RANDOM count sum for promoter '", promoter,
                                                "' in column(s): ",
                                                paste(colnames(count_matrix)[!valid_sums], collapse = ", ")
                                        )
                                )
                        }
                }

                if (length(warnings) > 0) {
                        warning(
                                paste0("Count-level normalization for ", matrix_name, " skipped scaling where needed: ",
                                       paste(warnings, collapse = "; "))
                        )
                }
                normalized
        }

        list(
                cdna = normalize_matrix(cdna_matrix, "cDNA"),
                pdna = normalize_matrix(pdna_matrix, "pDNA")
        )
}

# Loess-corrects contrast-condition barcode counts toward the reference condition,
# fitting on log2(cDNA_RPM / pDNA_RPM) so barcode representation in the input library
# does not masquerade as a contrast-specific abundance trend.
correct_banana_counts_activity <- function(cdna_count_matrix, pdna_count_matrix, ref_replicates, contrast_replicates, pseudocount = 1, span = 0.75, degree = 2, baseline_cdna_matrix = cdna_count_matrix) {
        col_info <- split_replicate_bc_colname(colnames(cdna_count_matrix))
        ref_cols <- col_info$col_name[col_info$replicate %in% ref_replicates]
        contrast_cols <- col_info$col_name[col_info$replicate %in% contrast_replicates]

        if (length(ref_cols) == 0 || length(contrast_cols) == 0) {
                message("---- Banana correction: no matching reference/contrast columns found, skipping")
                return(cdna_count_matrix)
        }

        if (!identical(dim(cdna_count_matrix), dim(pdna_count_matrix)) ||
            !identical(colnames(cdna_count_matrix), colnames(pdna_count_matrix)) ||
            !identical(rownames(cdna_count_matrix), rownames(pdna_count_matrix)) ||
            !identical(dim(cdna_count_matrix), dim(baseline_cdna_matrix)) ||
            !identical(colnames(cdna_count_matrix), colnames(baseline_cdna_matrix)) ||
            !identical(rownames(cdna_count_matrix), rownames(baseline_cdna_matrix))) {
                stop("Banana correction requires row- and column-aligned cDNA, baseline cDNA, and pDNA matrices")
        }

        cdna_rpm <- sweep(baseline_cdna_matrix, 2, colSums(baseline_cdna_matrix, na.rm = TRUE), FUN = function(x, s) x / s * 1e6)
        pdna_rpm <- sweep(pdna_count_matrix, 2, colSums(pdna_count_matrix), FUN = function(x, s) x / s * 1e6)
        log2_activity <- log2((cdna_rpm + pseudocount) / (pdna_rpm + pseudocount))

        long_df <- as.data.frame(log2_activity) %>%
                mutate(tf = rownames(cdna_count_matrix)) %>%
                pivot_longer(-tf, names_to = "col_name", values_to = "log2_activity") %>%
                left_join(col_info, by = "col_name")

        reference_means <- long_df %>%
                filter(col_name %in% ref_cols) %>%
                group_by(tf, bc) %>%
                summarise(reference_mean = mean(log2_activity, na.rm = TRUE), .groups = "drop")

        contrast_means <- long_df %>%
                filter(col_name %in% contrast_cols) %>%
                group_by(tf, bc) %>%
                summarise(contrast_mean = mean(log2_activity, na.rm = TRUE), .groups = "drop")

        banana_df <- reference_means %>%
                inner_join(contrast_means, by = c("tf", "bc")) %>%
                mutate(
                        A = reference_mean,
                        M = contrast_mean - reference_mean
                ) %>%
                filter(is.finite(A), is.finite(M))

        # Fit and evaluate the loess on TF-level averages, not per-barcode values: a
        # single barcode's own reference activity is a noisy estimate of that TF's true
        # activity, and fitting/predicting per-barcode causes regression dilution that
        # leaves a residual per-TF bias uncorrected - invisible in the noisy per-barcode
        # scatter but dominant once barcodes are averaged downstream (e.g. by mpralm).
        tf_level_df <- banana_df %>%
                group_by(tf) %>%
                summarise(A_tf = mean(A, na.rm = TRUE), M_tf = mean(M, na.rm = TRUE), .groups = "drop") %>%
                filter(is.finite(A_tf), is.finite(M_tf))

        n_unique_A <- n_distinct(tf_level_df$A_tf)
        fit <- if (nrow(tf_level_df) >= 10 && n_unique_A >= 5) {
                tryCatch(loess(M_tf ~ A_tf, data = tf_level_df, span = span, degree = degree), error = function(e) NULL)
        } else {
                NULL
        }

        if (is.null(fit)) {
                message("---- Banana correction: loess fit failed or insufficient TFs, skipping correction")
                return(cdna_count_matrix)
        }

        tf_level_df$M_fitted <- predict(fit, newdata = tf_level_df)
        tf_level_df$M_fitted[!is.finite(tf_level_df$M_fitted)] <- 0
        tf_level_df$correction_factor <- 2^(-tf_level_df$M_fitted)

        # Same per-TF correction factor applied to every barcode of that TF's contrast
        # columns - barcode-level M is treated as noise around the TF-level trend.
        cols_to_correct <- intersect(col_info$col_name[col_info$replicate %in% contrast_replicates], contrast_cols)
        corrected_matrix <- cdna_count_matrix
        for (i in seq_len(nrow(tf_level_df))) {
                # rownames(corrected_matrix) has one entry per barcode, so a TF's name is
                # duplicated across rows - name-based indexing (`matrix[tf_name, ]`) would
                # only touch the first match; use which() to hit every row for that TF.
                row_idx <- which(rownames(corrected_matrix) == tf_level_df$tf[i])
                corrected_matrix[row_idx, cols_to_correct] <- cdna_count_matrix[row_idx, cols_to_correct] * tf_level_df$correction_factor[i]
        }

        corrected_matrix
}

# Replicates are named <sample>_<plate_id>_<id>; the plate is the longest listed plate ID
# embedded in the name, or "all" when none matches.
banana_plate_ids <- trimws(strsplit(if (is.null(opt$plate_ids)) "" else opt$plate_ids, ",")[[1]])
banana_plate_ids <- banana_plate_ids[banana_plate_ids != ""]
banana_plate_ids <- banana_plate_ids[order(-nchar(banana_plate_ids))]
banana_plate_key <- function(replicates) {
        vapply(replicates, function(name) {
                hit <- banana_plate_ids[vapply(banana_plate_ids, function(p) grepl(paste0("_", p, "_[^_]*$"), name), logical(1))]
                if (length(hit) == 0) "all" else hit[1]
        }, character(1), USE.NAMES = FALSE)
}

# Corrects every group of replicates (e.g. non-anchor reference replicates and the
# contrast replicates) separately against the mean of the anchor replicates.
# With per_plate = TRUE each group is further split by plate and corrected against the
# anchors of its own plate (all anchors if that plate has none).
# Anchor replicates themselves are left untouched.
correct_banana_against_anchor <- function(cdna_count_matrix, pdna_count_matrix, anchor_replicates, replicate_groups, baseline_cdna_matrix = cdna_count_matrix, per_plate = FALSE) {
        anchor_keys <- banana_plate_key(anchor_replicates)
        for (group in replicate_groups) {
                group <- setdiff(group, anchor_replicates)
                if (length(group) == 0) next
                subgroups <- if (per_plate) split(group, banana_plate_key(group)) else list(group)
                for (key in names(subgroups)) {
                        anchors <- if (per_plate && any(anchor_keys == key)) anchor_replicates[anchor_keys == key] else anchor_replicates
                        cdna_count_matrix <- correct_banana_counts_activity(
                                cdna_count_matrix,
                                pdna_count_matrix,
                                anchors,
                                subgroups[[key]],
                                baseline_cdna_matrix = baseline_cdna_matrix
                        )
                }
        }
        cdna_count_matrix
}

# Per-replicate MA data (activity level) before vs after banana correction, written to a TSV that
# primetime_get_banana_qc.R turns into one QC page per replicate.
# A = mean log2 activity of the anchor/reference replicates, M = replicate - A.
banana_replicate_ma_df <- function(cdna_matrix_before, cdna_matrix_after, pdna_matrix, anchor_replicates, replicate_condition, pseudocount = 1) {
        col_info <- split_replicate_bc_colname(colnames(cdna_matrix_before))
        anchor_cols <- col_info$col_name[col_info$replicate %in% anchor_replicates]
        replicates <- names(replicate_condition)
        if (length(anchor_cols) == 0 || length(replicates) == 0) return(data.frame())

        pdna_rpm <- sweep(pdna_matrix, 2, colSums(pdna_matrix), FUN = function(x, s) x / s * 1e6)

        make_ma_df <- function(cdna_matrix, stage) {
                cdna_rpm <- sweep(cdna_matrix, 2, colSums(cdna_matrix, na.rm = TRUE), FUN = function(x, s) x / s * 1e6)
                log2_activity <- log2((cdna_rpm + pseudocount) / (pdna_rpm + pseudocount))
                long_df <- as.data.frame(log2_activity) %>%
                        mutate(tf = rownames(cdna_matrix)) %>%
                        pivot_longer(-tf, names_to = "col_name", values_to = "log2_activity") %>%
                        left_join(col_info, by = "col_name")
                anchor_means <- long_df %>%
                        filter(col_name %in% anchor_cols) %>%
                        group_by(tf, bc) %>%
                        summarise(A = mean(log2_activity, na.rm = TRUE), .groups = "drop")
                long_df %>%
                        filter(replicate %in% replicates) %>%
                        group_by(replicate, tf, bc) %>%
                        summarise(activity = mean(log2_activity, na.rm = TRUE), .groups = "drop") %>%
                        inner_join(anchor_means, by = c("tf", "bc")) %>%
                        mutate(M = activity - A, stage = stage) %>%
                        select(replicate, tf, bc, A, M, stage)
        }

        bind_rows(
                make_ma_df(cdna_matrix_before, "Before correction"),
                make_ma_df(cdna_matrix_after, "After correction")
        ) %>%
                filter(is.finite(A), is.finite(M)) %>%
                mutate(condition = unname(replicate_condition[replicate]))
}

plot_bcalm_normalization_diagnostic <- function(results_before, results_after, reference_condition, contrast_condition, output_file) {
        diagnostic_df <- bind_rows(
                results_before %>% mutate(stage = "Before normalization"),
                results_after %>% mutate(stage = "After normalization")
        ) %>%
                pivot_longer(cols = all_of(c(reference_condition, contrast_condition)), names_to = "condition", values_to = "activity") %>%
                mutate(
                        reporter_type = ifelse(grepl("RANDOM", tf, ignore.case = TRUE), "RANDOM reporter", "TF reporter"),
                        stage = factor(stage, levels = c("Before normalization", "After normalization")),
                        condition = factor(condition, levels = c(reference_condition, contrast_condition))
                ) %>%
                filter(is.finite(activity))

        if (nrow(diagnostic_df) == 0) return(invisible(NULL))

        pdf(output_file, width = 12, height = 7)
        print(
                ggplot(diagnostic_df, aes(x = activity, y = condition)) +
                        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                        geom_boxplot(aes(fill = reporter_type), alpha = 0.25, outlier.shape = NA, width = 0.5) +
                        geom_jitter(aes(color = reporter_type), height = 0.12, alpha = 0.3, size = 0.9) +
                        stat_summary(data = function(data) filter(data, reporter_type == "RANDOM reporter"), fun = median, geom = "point", shape = 23, fill = "white", color = "black", size = 3) +
                        facet_wrap(~stage) +
                        scale_color_manual(values = c("RANDOM reporter" = "#e15759", "TF reporter" = "#4c78a8")) +
                        scale_fill_manual(values = c("RANDOM reporter" = "#e15759", "TF reporter" = "#4c78a8")) +
                        ggpubr::theme_pubr(border = TRUE) +
                        labs(title = "BCalm activity normalization by RANDOM reporter median", x = "BCalm activity (log2 scale)", y = "Condition", color = NULL, fill = NULL) +
                        theme(text = element_text(size = 14), legend.position = "bottom")
        )
        invisible(dev.off())
}

plot_bcalm_normalization_ma_diagnostic <- function(results_before, results_after, reference_condition, contrast_condition, output_file) {
        diagnostic_df <- bind_rows(
                results_before %>% mutate(stage = "Before normalization"),
                results_after %>% mutate(stage = "After normalization")
        ) %>%
                transmute(
                        tf,
                        stage = factor(stage, levels = c("Before normalization", "After normalization")),
                        A = (!!sym(reference_condition) + !!sym(contrast_condition)) / 2,
                        M = !!sym(contrast_condition) - !!sym(reference_condition),
                        reporter_type = ifelse(grepl("RANDOM", tf, ignore.case = TRUE), "RANDOM reporter", "TF reporter")
                ) %>%
                filter(is.finite(A), is.finite(M))

        if (nrow(diagnostic_df) == 0) return(invisible(NULL))

        pdf(output_file, width = 12, height = 6)
        print(
                ggplot(diagnostic_df, aes(x = A, y = M, color = reporter_type)) +
                        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                        geom_point(alpha = 0.45, size = 1.2) +
                        geom_smooth(data = function(data) filter(data, reporter_type == "TF reporter"), method = "loess", se = FALSE, color = "#2c7fb8", span = 0.75) +
                        facet_wrap(~stage) +
                        scale_color_manual(values = c("RANDOM reporter" = "#e15759", "TF reporter" = "#4c78a8")) +
                        ggpubr::theme_pubr(border = TRUE) +
                        labs(title = "MA diagnostic of BCalm activity normalization", x = "A (mean reference and contrast BCalm activity)", y = "M (contrast - reference BCalm activity)", color = NULL) +
                        theme(text = element_text(size = 14), legend.position = "bottom")
        )
        invisible(dev.off())
}

# Removes promoter-level baseline offset WITHIN each replicate independently, using the
# negative-control (RANDOM) barcodes as that replicate's own "core promoter, no TF"
# baseline. For each replicate and promoter, computes that replicate's own median
# negative-control log2 RPM and rescales ALL barcodes on that promoter (real TFs
# included) in that replicate so its own core-promoter baseline sits at log2(RPM+1) = 0.
# No cross-replicate reference is used - each replicate is normalized purely to its own
# core promoter, independently of every other replicate.
correct_core_promoter_counts <- function(count_matrix, negative_control_flag, promoter_vector, pseudocount = 1) {
        rpm_matrix <- sweep(count_matrix, 2, colSums(count_matrix, na.rm = TRUE), FUN = function(x, s) x / s * 1e6)
        log2_rpm <- log2(rpm_matrix + pseudocount)

        long_df <- as.data.frame(log2_rpm) %>%
                mutate(
                        tf_row = rownames(count_matrix),
                        promoter = promoter_vector,
                        negative_control = negative_control_flag
                ) %>%
                pivot_longer(-c(tf_row, promoter, negative_control), names_to = "col_name", values_to = "log2_rpm")

        # Each replicate's own core-promoter (negative-control) baseline, corrected
        # independently - no grand/common reference across replicates.
        replicate_promoter_baseline <- long_df %>%
                filter(negative_control) %>%
                group_by(col_name, promoter) %>%
                summarise(baseline = median(log2_rpm, na.rm = TRUE), .groups = "drop") %>%
                filter(is.finite(baseline)) %>%
                mutate(correction_factor = 2^(-baseline))

        if (nrow(replicate_promoter_baseline) == 0) {
                message("---- Core promoter correction: no valid negative-control baselines found, skipping")
                return(count_matrix)
        }

        corrected_matrix <- count_matrix
        for (i in seq_len(nrow(replicate_promoter_baseline))) {
                this_col <- replicate_promoter_baseline$col_name[i]
                this_promoter <- replicate_promoter_baseline$promoter[i]
                rows_to_correct <- which(promoter_vector == this_promoter)
                corrected_matrix[rows_to_correct, this_col] <- count_matrix[rows_to_correct, this_col] * replicate_promoter_baseline$correction_factor[i]
        }

        corrected_matrix
}

# Diagnostic plot (before vs after core-promoter correction): each replicate's
# negative-control median log2 RPM per promoter, so you can visually confirm the
# per-replicate promoter baselines converge toward the common reference level after
# correction instead of staying spread out.
plot_core_promoter_diagnostic <- function(count_matrix_before, count_matrix_after, negative_control_flag, promoter_vector, title_suffix, pseudocount = 1) {
        make_baseline_df <- function(cm) {
                rpm <- sweep(cm, 2, colSums(cm, na.rm = TRUE), FUN = function(x, s) x / s * 1e6)
                log2_rpm <- log2(rpm + pseudocount)
                as.data.frame(log2_rpm) %>%
                        mutate(promoter = promoter_vector, negative_control = negative_control_flag) %>%
                        filter(negative_control) %>%
                        select(-negative_control) %>%
                        pivot_longer(-promoter, names_to = "col_name", values_to = "log2_rpm") %>%
                        group_by(col_name, promoter) %>%
                        summarise(baseline = median(log2_rpm, na.rm = TRUE), .groups = "drop")
        }

        before_df <- make_baseline_df(count_matrix_before) %>% mutate(stage = "Before correction")
        after_df <- make_baseline_df(count_matrix_after) %>% mutate(stage = "After correction")
        combined_df <- bind_rows(before_df, after_df) %>% filter(is.finite(baseline))

        if (nrow(combined_df) == 0) return(invisible(NULL))

        combined_df$stage <- factor(combined_df$stage, levels = c("Before correction", "After correction"))

        print(
                ggplot(combined_df, aes(x = col_name, y = baseline, color = promoter)) +
                        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                        geom_point(size = 2) +
                        facet_wrap(~stage) +
                        ggpubr::theme_pubr(border = TRUE) +
                        labs(
                                title = paste("Core promoter correction diagnostic -", title_suffix),
                                x = "Replicate x barcode column",
                                y = "Negative-control median log2(RPM+1) (per promoter)"
                        ) +
                        theme(
                                text = element_text(size = 12),
                                axis.text.x = element_blank(),
                                axis.ticks.x = element_blank()
                        )
        )
}

##########################################################################################
# Preparing the data #####################################################################
##########################################################################################

pdna <- read.table(opt$pdna, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
name_of_the_pdna_replicate <- colnames(pdna) %>% setdiff(c("barcode", "negative_control"))

contrast_condition <- opt$contrast_condition
reference_condition <- opt$reference_condition
normalize_counts <- isTRUE(opt$normalize)

# Adaptive-shrinkage-style empirical Bayes (ashr is not installed): zero-centred mixture of normals
# prior learned from the observed logFC and their standard errors; returns posterior mean/sd.
shrink_logfc <- function(beta, se, grid_mult = c(0, 0.25, 0.5, 1, 2, 4, 8, 16), n_iter = 500) {
        ok <- is.finite(beta) & is.finite(se) & se > 0
        out <- list(mean = rep(NA_real_, length(beta)), weight = rep(NA_real_, length(beta)))
        if (sum(ok) < 10) {
                return(out)
        }
        b <- beta[ok]
        s <- se[ok]
        sds <- sort(unique(c(0, min(s) * grid_mult[-1], median(s) * grid_mult[-1])))
        marg_sd <- sqrt(outer(s^2, sds^2, "+"))
        lik <- dnorm(b, 0, marg_sd)
        pi_k <- rep(1 / length(sds), length(sds))
        for (i in seq_len(n_iter)) {
                resp <- sweep(lik, 2, pi_k, "*")
                resp <- resp / pmax(rowSums(resp), 1e-300)
                pi_k <- colMeans(resp)
        }
        resp <- sweep(lik, 2, pi_k, "*")
        resp <- resp / pmax(rowSums(resp), 1e-300)
        comp_w <- outer(rep(1, length(s)), sds^2) / marg_sd^2
        w <- rowSums(resp * comp_w)
        out$mean[ok] <- w * b
        out$weight[ok] <- w
        out
}
banana_correction <- isTRUE(opt$banana_correction)
shrinkage <- isTRUE(opt$shrinkage)
logfc_threshold <- if (normalize_counts) 0.263 else 0

parse_replicate_arg <- function(value) {
        if (is.null(value) || is.na(value)) {
                return(character(0))
        }
        pieces <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
        pieces[pieces != ""]
}

# Setup the replicates that are gonna be used here
ref_replicates <- parse_replicate_arg(opt$reference_replicates)
contrast_replicates <- parse_replicate_arg(opt$contrast_replicates)

if (length(ref_replicates) == 0) {
        ref_replicates <- paste0(reference_condition, "_", seq_len(opt$num_replicates_reference))
}
if (length(contrast_replicates) == 0) {
        contrast_replicates <- paste0(contrast_condition, "_", seq_len(opt$num_replicates_contrast))
}


# Low-bleedthrough control replicates used as the banana-correction anchor. The anchors
# are loaded as extra columns so they can serve as reference, then dropped before modelling.
banana_anchor_replicates <- character(0)
if (banana_correction && !is.null(opt$banana_reference_file) && file.exists(opt$banana_reference_file)) {
        # Anchors are chosen per comparison: only the flagged replicates of this comparison's reference condition
        banana_reference_df <- read.table(opt$banana_reference_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
        banana_anchor_replicates <- as.character(banana_reference_df$replicate[banana_reference_df$used_as_reference & banana_reference_df$replicate %in% ref_replicates])
}
extra_anchor_replicates <- character(0)
extra_anchor_replicates <- setdiff(banana_anchor_replicates, c(ref_replicates, contrast_replicates))
cdna_all <- read.table(opt$cdna, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
if (is.null(opt$cdna_baseline) || !file.exists(opt$cdna_baseline)) {
        stop("The uncorrected cDNA baseline file is required for order-independent banana correction: ", opt$cdna_baseline)
}
cdna_baseline_all <- read.table(opt$cdna_baseline, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
cdna_key_columns <- c("tf", "negative_control", "promoter", "barcode")
if (!all(cdna_key_columns %in% colnames(cdna_all)) ||
    !all(cdna_key_columns %in% colnames(cdna_baseline_all))) {
        stop("Corrected and uncorrected cDNA files must include TF, negative_control, promoter, and barcode columns")
}
cdna_keys <- paste(cdna_all$tf, cdna_all$negative_control, cdna_all$promoter, cdna_all$barcode, sep = "\r")
baseline_keys <- paste(cdna_baseline_all$tf, cdna_baseline_all$negative_control, cdna_baseline_all$promoter, cdna_baseline_all$barcode, sep = "\r")
if (anyDuplicated(cdna_keys) ||
    anyDuplicated(baseline_keys) ||
    length(cdna_keys) != length(baseline_keys) ||
    anyNA(match(cdna_keys, baseline_keys))) {
        stop("Corrected and uncorrected cDNA files must contain the same unique TF/promoter/barcode rows")
}
cdna_baseline_all <- cdna_baseline_all[match(cdna_keys, baseline_keys), , drop = FALSE]
extra_anchor_replicates <- intersect(extra_anchor_replicates, colnames(cdna_all))

# Replicates removed upstream (e.g. MIN_READS_PER_REPLICATE) are no longer in the count tables.
missing_replicates <- setdiff(c(ref_replicates, contrast_replicates), colnames(cdna_all))
if (length(missing_replicates) > 0) {
        message("Dropping replicates absent from the count tables (excluded upstream): ", paste(missing_replicates, collapse = ", "))
        ref_replicates <- intersect(ref_replicates, colnames(cdna_all))
        contrast_replicates <- intersect(contrast_replicates, colnames(cdna_all))
}
if (length(ref_replicates) == 0 || length(contrast_replicates) == 0) {
        message("No replicates left for ", opt$contrast_condition, " vs ", opt$reference_condition, "; writing empty outputs")
        for (d in c("volcano_plots", "lollipop_plots", "circular_lollipop_plots", "output_data")) {
                dir.create(file.path(opt$plot_output, d), showWarnings = FALSE, recursive = TRUE)
        }
        nm <- paste0(opt$contrast_condition, "_vs_", opt$reference_condition)
        file.create(file.path(opt$plot_output, "output_data", paste0(nm, ".txt")))
        for (d in c("volcano_plots", "lollipop_plots", "circular_lollipop_plots")) {
                pdf(file.path(opt$plot_output, d, paste0(nm, ".pdf")), width = 6, height = 3)
                plot.new()
                text(0.5, 0.5, paste(nm, "\nno replicates passed the replicate QC filters"))
                invisible(dev.off())
        }
        banana_qc_empty <- file.path(opt$output, "banana_qc", paste0(nm, ".tsv"))
        dir.create(dirname(banana_qc_empty), showWarnings = FALSE, recursive = TRUE)
        file.create(banana_qc_empty)
        quit(status = 0)
}
banana_anchor_replicates <- intersect(banana_anchor_replicates, c(ref_replicates, contrast_replicates, extra_anchor_replicates))
cdna <- cdna_all %>%
        select(
                tf,
                negative_control,
                promoter,
                barcode,
                all_of(c(ref_replicates, contrast_replicates, extra_anchor_replicates))
        )
cdna_baseline <- cdna_baseline_all %>%
        select(
                tf,
                negative_control,
                promoter,
                barcode,
                all_of(c(ref_replicates, contrast_replicates, extra_anchor_replicates))
        )

# Now, biuld the design df, with columns: sample (conditions and pDNA), replicate (names ), pDNA (whether is pDNA or not) and treatment (condition name; None for pDNA)
design_df <- data.frame(
        replicate = c(ref_replicates, contrast_replicates, extra_anchor_replicates),
        treatment = c(
                rep(reference_condition, length(ref_replicates)),
                rep(contrast_condition, length(contrast_replicates)),
                rep("banana_anchor_only", length(extra_anchor_replicates))
        ),
        sample = c(
                rep(reference_condition, length(ref_replicates)),
                rep(contrast_condition, length(contrast_replicates)),
                rep("banana_anchor_only", length(extra_anchor_replicates))
        )
)

output_pdf_name <- paste0(contrast_condition, "_vs_", reference_condition, ".pdf")
output_txt_name <- paste0(contrast_condition, "_vs_", reference_condition, ".txt")
p_threshold <- opt$pval_threshold

dir.create(file.path(opt$plot_output, "volcano_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$plot_output, "lollipop_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$plot_output, "circular_lollipop_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$plot_output, "output_data"), showWarnings = FALSE, recursive = TRUE)
banana_qc_file <- file.path(opt$output, "banana_qc", sub("\\.pdf$", ".tsv", output_pdf_name))
dir.create(dirname(banana_qc_file), showWarnings = FALSE, recursive = TRUE)
file.create(banana_qc_file)

tf_function_path <- find_existing_path(c(
        if (!is.na(project_root)) file.path(project_root, "misc", "tf_functions.tsv") else NA_character_,
        file.path(dirname(opt$plot_output), "misc", "tf_functions.tsv"),
        file.path(dirname(dirname(opt$plot_output)), "misc", "tf_functions.tsv")
))
tf_function_map <- NULL
tf_function_order <- character(0)
tf_order_in_map <- character(0)
if (!is.na(tf_function_path) && file.exists(tf_function_path)) {
        tf_function_map <- read.table(tf_function_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        if ("tf" %in% colnames(tf_function_map)) {
                tf_function_map$tf <- trimws(as.character(tf_function_map$tf))
        }
        if ("biological_function" %in% colnames(tf_function_map)) {
                tf_function_map$biological_function <- trimws(as.character(tf_function_map$biological_function))
        }
        if ("biological_function" %in% colnames(tf_function_map)) {
                tf_function_order <- unique(tf_function_map$biological_function)
        }
        if ("tf" %in% colnames(tf_function_map)) {
                tf_order_in_map <- unique(tf_function_map$tf)
        }
} else {
        warning("Could not find tf_functions.tsv - circular plot will fall back to alphabetical TF ordering")
}

df <-
        merge(pdna, cdna %>% select(-negative_control), by = "barcode") %>%
        group_by(tf) %>%
        # reframe and keep all the other columns
        mutate(
                barcode_id = paste0("bc", row_number())
        ) %>%
        # melt the data to have the columns: tf, promomer, barcode_id, replicate, count
        melt(
                id.vars = c("tf", "barcode_id", "promoter", "barcode", "negative_control"),
                variable.name = "replicate",
                value.name = "count"
        )

parsed_cdna <-
        df %>%
        filter(replicate != name_of_the_pdna_replicate) %>%
        mutate(replicate_id = paste(replicate, barcode_id, sep = "_")) %>%
        select(-barcode_id, -barcode, -replicate) %>%
        # Transform back again to every column being a replicate, but now we have the Sample_RepID_bcID
        pivot_wider(names_from = replicate_id, values_from = count) %>%
        as.data.frame()

baseline_barcode_ids <- df %>%
        distinct(tf, negative_control, promoter, barcode, barcode_id)
if (anyDuplicated(baseline_barcode_ids[c("tf", "negative_control", "promoter", "barcode")])) {
        stop("Could not uniquely map cDNA baseline rows to barcode identifiers")
}
cdna_baseline_with_ids <- cdna_baseline %>%
        left_join(
                baseline_barcode_ids,
                by = c("tf", "negative_control", "promoter", "barcode")
        )
if (anyNA(cdna_baseline_with_ids$barcode_id)) {
        stop("Some uncorrected cDNA rows could not be matched to barcode identifiers")
}
parsed_cdna_baseline <-
        cdna_baseline_with_ids %>%
        pivot_longer(
                cols = all_of(c(ref_replicates, contrast_replicates, extra_anchor_replicates)),
                names_to = "replicate",
                values_to = "count"
        ) %>%
        mutate(replicate_id = paste(replicate, barcode_id, sep = "_")) %>%
        select(tf, negative_control, promoter, replicate_id, count) %>%
        pivot_wider(names_from = replicate_id, values_from = count) %>%
        as.data.frame()

parsed_cdna_keys <- paste(parsed_cdna$tf, parsed_cdna$negative_control, parsed_cdna$promoter, sep = "\r")
parsed_baseline_keys <- paste(parsed_cdna_baseline$tf, parsed_cdna_baseline$negative_control, parsed_cdna_baseline$promoter, sep = "\r")
if (anyDuplicated(parsed_cdna_keys) ||
    anyDuplicated(parsed_baseline_keys) ||
    length(parsed_cdna_keys) != length(parsed_baseline_keys) ||
    anyNA(match(parsed_cdna_keys, parsed_baseline_keys))) {
        stop("Parsed corrected and uncorrected cDNA matrices do not have matching TF/promoter rows")
}
parsed_cdna_baseline <- parsed_cdna_baseline[match(parsed_cdna_keys, parsed_baseline_keys), , drop = FALSE]

# Keep a tf -> promoter lookup, since this_cDNA/this_cdna_matrix below drop the
# promoter column (needed for the core-promoter correction pre-mpralm).
tf_promoter_map <- parsed_cdna %>% distinct(tf, promoter)

annotation_df <-
        data.frame(
                obs = colnames(parsed_cdna)
        ) %>%
        as_tibble() %>%
        mutate(
                barcode = gsub(".*_", "", obs),
                replicate = gsub("_bc.*", "", obs),
                # condition for U2OS_Calcitriol_12W_1 is U2OS_Calcitriol_12W. So just remove the last number and underscore
                bio_replicate = gsub("_\\d$", "", gsub(".*_", "", gsub("_bc.*", "", replicate)))
        ) %>%
        merge(design_df, by = "replicate")

rownames(annotation_df) <- annotation_df$obs
annotation_df %>% select(-sample, -replicate, -obs) -> annotation_df

# 1. Make the pDNA df
df %>%
        filter(replicate == name_of_the_pdna_replicate) %>%
        mutate(replicate_id = paste(replicate, barcode_id, sep = "_")) %>%
        select(-barcode_id, -barcode, -replicate) %>%
        pivot_wider(names_from = replicate_id, values_from = count) -> tmp_parsed_pdna


# now get the values from tmp_parsed_pdna but with the names of the columns of
# parsed_cdna, matching the barcode id
# =================================
# the idea is that I will iterate over all the columns of the tmp_parsed_pdna (which are
# like pDNA_bc1, pDNA_bc2, etc) then: see which barcode it is (bc1, bc2, etc) and then get
# the column of parsed_cdna that has the same barcode. The output is a df with the same
# columns as parsed_cdna, but with the values of the pDNA
parsed_pdna <- list(negative_controls = tmp_parsed_pdna$negative_control, promoter = tmp_parsed_pdna$promoter)
columns_of_parsed_cdna <- rownames(annotation_df)
for (col in colnames(tmp_parsed_pdna) %>% setdiff(c("tf", "promoter", "negative_control"))) {
        barcode <- gsub(".*_", "", col)
        # See which column of parsed_cdna has the same barcode and return them
        columns_w_same_bc <- columns_of_parsed_cdna[endsWith(columns_of_parsed_cdna, paste0("_", barcode))]
        for (col_w_same_bc in columns_w_same_bc) {
                parsed_pdna[[col_w_same_bc]] <- as.numeric(tmp_parsed_pdna[[col]])
        }
}
parsed_pdna <- do.call(cbind, parsed_pdna) %>% as.data.frame()

cdna_row_keys <- paste(parsed_cdna$tf, parsed_cdna$negative_control, parsed_cdna$promoter, sep = "\r")
pdna_row_keys <- paste(tmp_parsed_pdna$tf, tmp_parsed_pdna$negative_control, tmp_parsed_pdna$promoter, sep = "\r")
if (anyDuplicated(cdna_row_keys) || anyDuplicated(pdna_row_keys)) {
        stop("cDNA/pDNA TF-promoter row keys are not unique")
}
pdna_row_order <- match(cdna_row_keys, pdna_row_keys)
if (anyNA(pdna_row_order)) {
        stop("Some cDNA TF-promoter rows have no matching pDNA rows")
}
parsed_pdna <- parsed_pdna[pdna_row_order, , drop = FALSE]
tmp_parsed_pdna <- tmp_parsed_pdna[pdna_row_order, , drop = FALSE]

# ============ SINGLE MODEL - ALL PROMOTERS TOGETHER =====================================
# ==========================================================================================
message("==== Fitting a single BCalm model for all promoters together\n")

this_cdna <- parsed_cdna %>%
        select(-promoter)
this_pdna <- parsed_pdna %>%
        select(-promoter)
this_tmp_parsed_pdna <- tmp_parsed_pdna

# Updating rownames
rownames(this_pdna) <- this_tmp_parsed_pdna$tf
rownames(this_cdna) <- this_cdna$tf
message("this_pdna")
print(head(this_pdna))
message("this_cdna")
print(head(this_cdna))

# Convert to matrix
this_pdna %>%
        select(-negative_controls) %>%
        mutate_all(as.numeric) %>%
        # arrange columns alphabetically
        select(order(colnames(this_pdna %>% select(-negative_controls)))) %>%
        as.matrix() -> this_pdna_matrix


this_cdna %>%
        select(-tf, -negative_control) %>%
        mutate_all(as.numeric) %>%
        # arrange columns alphabetically
        select(colnames(this_pdna_matrix)) %>%
        as.matrix() -> this_cdna_matrix

parsed_cdna_baseline %>%
        select(-tf, -negative_control, -promoter) %>%
        mutate_all(as.numeric) %>%
        select(colnames(this_pdna_matrix)) %>%
        as.matrix() -> this_cdna_baseline_matrix

rownames(this_pdna_matrix) = this_tmp_parsed_pdna$tf
rownames(this_cdna_matrix) <- this_cdna$tf
rownames(this_cdna_baseline_matrix) <- this_cdna$tf

print("this pDNA matrix")
print(head(this_pdna_matrix))
print("This cDNA matrix (before any correction)")
print(head(this_cdna_matrix))

this_cdna_matrix_precorrection <- this_cdna_matrix
# Banana fitted on the depth/offset/shrinkage-corrected counts, so it only removes the
# bias those corrections have not already explained.
this_cdna_baseline_matrix <- this_cdna_matrix_precorrection
if (banana_correction) {
        message("---- Correcting banana-shaped bias (activity level: log2(cDNA_RPM/pDNA_RPM)) in contrast vs reference barcode counts")
        if (length(banana_anchor_replicates) > 0) {
                message("---- Banana anchor (low-bleedthrough controls): ", paste(banana_anchor_replicates, collapse = ", "))
                diag_ref_replicates <- banana_anchor_replicates
                this_cdna_matrix <- correct_banana_against_anchor(
                        this_cdna_matrix, this_pdna_matrix, banana_anchor_replicates,
                        list(ref_replicates, contrast_replicates),
                        baseline_cdna_matrix = this_cdna_baseline_matrix,
                        per_plate = TRUE
                )
        } else {
                diag_ref_replicates <- ref_replicates
                this_cdna_matrix <- correct_banana_counts_activity(
                        this_cdna_matrix,
                        this_pdna_matrix,
                        ref_replicates,
                        contrast_replicates,
                        baseline_cdna_matrix = this_cdna_baseline_matrix
                )
        }

        # Replicates that were corrected (everything except the anchors) with their condition
        banana_qc_replicates <- setdiff(c(ref_replicates, contrast_replicates), diag_ref_replicates)
        banana_qc_condition <- setNames(
                ifelse(banana_qc_replicates %in% contrast_replicates, contrast_condition, reference_condition),
                banana_qc_replicates
        )
        write.table(
                banana_replicate_ma_df(this_cdna_matrix_precorrection, this_cdna_matrix, this_pdna_matrix, diag_ref_replicates, banana_qc_condition),
                banana_qc_file, sep = "\t", row.names = FALSE, quote = FALSE
        )

        print("This cDNA matrix (after banana correction)")
        print(head(this_cdna_matrix))
} else {
        message("---- Skipping banana-shaped bias correction (banana_correction = FALSE)")
}

if (length(extra_anchor_replicates) > 0) {
        keep_cols <- !split_replicate_bc_colname(colnames(this_cdna_matrix))$replicate %in% extra_anchor_replicates
        this_cdna_matrix <- this_cdna_matrix[, keep_cols, drop = FALSE]
        this_pdna_matrix <- this_pdna_matrix[, keep_cols, drop = FALSE]
        this_cdna_baseline_matrix <- this_cdna_baseline_matrix[, keep_cols, drop = FALSE]
        this_cdna_matrix_precorrection <- this_cdna_matrix_precorrection[, keep_cols, drop = FALSE]
}

# Count-level normalization: scale each promoter/column using its RANDOM reporter sum.
# This is applied BEFORE the model, so variance estimates and p-values see the corrected data.
if (normalize_counts) {
        message("---- Normalizing counts at count level by per-promoter RANDOM reporter sums")
        norm_result <- normalize_counts_by_random(
                this_cdna_matrix,
                this_pdna_matrix,
                this_tmp_parsed_pdna$negative_control,
                this_tmp_parsed_pdna$promoter,
                target = 1e7
        )
        this_cdna_matrix <- norm_result$cdna
        this_pdna_matrix <- norm_result$pdna
        
        print("This cDNA matrix (after count-level normalization)")
        print(head(this_cdna_matrix))
}

corrected_matrix_long <- as.data.frame(this_cdna_matrix, check.names = FALSE) %>%
        mutate(tf = rownames(this_cdna_matrix)) %>%
        pivot_longer(-tf, names_to = "replicate_id", values_to = "count") %>%
        mutate(
                replicate = sub("_bc[0-9]+$", "", replicate_id),
                barcode_id = sub("^.*_(bc[0-9]+)$", "\\1", replicate_id)
        ) %>%
        select(tf, barcode_id, replicate, count) %>%
        pivot_wider(names_from = replicate, values_from = count)

corrected_cdna_output <- cdna_all %>%
        left_join(
                df %>% distinct(tf, barcode, barcode_id),
                by = c("tf", "barcode")
        ) %>%
        distinct(tf, negative_control, promoter, barcode, barcode_id) %>%
        left_join(
                corrected_matrix_long,
                by = c("tf", "barcode_id"),
                suffix = c("", ".corrected")
        )
for (replicate_name in c(ref_replicates, contrast_replicates)) {
        corrected_cdna_output[[replicate_name]] <- corrected_cdna_output[[paste0(replicate_name, ".corrected")]]
}
corrected_cdna_output %>%
        select(-barcode_id, -ends_with(".corrected")) -> corrected_cdna_output
write.table(
        corrected_cdna_output,
        file.path(dirname(opt$cdna), "cDNA_counts_corrected.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
)

# ========================================================================================
BcVariantMPRASet <- MPRASet(
        DNA = this_pdna_matrix,
        RNA = this_cdna_matrix,
        eid = rownames(this_pdna_matrix),
        barcode = NULL
)

# Preparing the design_bcalm data, where each row is a sample (same order as the matrix)
# and we have whether it is the contrast condition
ordered_annotation = annotation_df[colnames(this_pdna_matrix), ]
if (is.na(contrast_condition) || is.na(reference_condition) ||
    contrast_condition == reference_condition ||
    anyNA(ordered_annotation$treatment) ||
    !all(ordered_annotation$treatment %in% c(reference_condition, contrast_condition)) ||
    !all(c(reference_condition, contrast_condition) %in% ordered_annotation$treatment)) {
        stop("Invalid BCalm design: both distinct reference and contrast conditions must have annotated observations")
}

# Now create the design data
design_bcalm =
        data.frame(
                intcpt = 1,
                ordered_annotation$treatment == contrast_condition
        )

colnames(design_bcalm) = c("intcpt", contrast_condition)

print("Design BCalm")
print(design_bcalm)

bcalm_normalize <- !normalize_counts
if (!bcalm_normalize) {
        message("---- Disabling BCalm global library-size normalization; per-promoter count-level normalization is already applied")
}

mpralm_fit_var <- mpralm(
        object = BcVariantMPRASet,
        design = design_bcalm,
        aggregate = "none",
        normalize = bcalm_normalize,
        model_type = "indep_groups",
        plot = FALSE
)
# For the real comparison, retrieve the second coefficient of the model
results_bcalm = topTable(mpralm_fit_var, coef = 2, number = Inf)

logfc_se <- tryCatch(
        mpralm_fit_var$stdev.unscaled[, 2] * sqrt(mpralm_fit_var$s2.post),
        error = function(e) NULL
)
reporter_se <- if (is.null(logfc_se)) setNames(abs(results_bcalm$logFC / results_bcalm$t), rownames(results_bcalm)) else logfc_se

# For the activity of the reference condition, retrieve the first coefficient of the model
# (later, we sum the LogFC column of results bcalm to this to get the activity of
# the contrast condition)
reference_results =
        topTable(mpralm_fit_var, coef = 1, number = Inf) %>%
        mutate(
                tf = rownames(.),
                !!reference_condition := logFC, .keep = "none"
        )
# merge both to get the activity of contrast condition
results_bcalm =
        results_bcalm %>%
        mutate(tf = rownames(.)) %>%
        left_join(reference_results, by = "tf") %>%
        mutate(!!contrast_condition := logFC + !!sym(reference_condition))

all_results = results_bcalm

all_results <- all_results %>% mutate(old_logFC = logFC)

# Correcting the p-values for multiple testing
message("==== Correcting p-values")
all_results <- all_results %>%
        mutate(p_adjusted = p.adjust(P.Value, method = "BH"))

if (normalize_counts) {
        # For count-level normalization: p-values already reflect the corrected data,
        # so we do not need a sign-agreement guard. Use the logFC directly.
        message("---- Using count-level normalized results (p-values see corrected data)")
        all_results <- all_results %>%
                mutate(
                        sig = ifelse(p_adjusted <= p_threshold,
                                ifelse(logFC >= logfc_threshold, "Upregulated", 
                                       ifelse(logFC <= -logfc_threshold, "Downregulated", "NS")),
                                "NS"
                        )
                )
} else {
        all_results <- all_results %>%
                mutate(
                        sig = ifelse(p_adjusted <= p_threshold,
                                ifelse(logFC >= 0, "Upregulated", "Downregulated"),
                                "NS"
                        )
                )
}

# Post-BCalm shrinkage: inference (p-values, sig) stays on the unshrunk logFC.
all_results$logFC_se <- unname(reporter_se[all_results$tf])
if (shrinkage) {
        message("==== Shrinking logFC after BCalm (empirical Bayes, prior learned from non-RANDOM TFs)")
        real_mask <- !grepl("RANDOM", all_results$tf, ignore.case = TRUE)
        # The BCalm SE is not calibrated (barcodes are shared between conditions but modelled as
        # independent groups), so z = logFC / SE of the RANDOM reporters (true nulls) is rescaled to unit spread.
        random_z <- (all_results$logFC / all_results$logFC_se)[!real_mask]
        random_z <- random_z[is.finite(random_z)]
        se_scale <- if (length(random_z) >= 10) min(max(1.4826 * mad(random_z, constant = 1), 0.3), 1) else 1
        all_results$se_scale <- se_scale
        message("SE calibration factor from RANDOM reporters: ", round(se_scale, 3))
        shrunk <- shrink_logfc(all_results$logFC[real_mask], all_results$logFC_se[real_mask] * se_scale)
        all_results$logFC_shrunk <- NA_real_
        all_results$shrink_weight <- NA_real_
        all_results$logFC_shrunk[real_mask] <- shrunk$mean
        all_results$shrink_weight[real_mask] <- shrunk$weight
        all_results$logFC_shrunk[!real_mask] <- all_results$logFC[!real_mask]
        all_results[[paste0(contrast_condition, "_shrunk")]] <- all_results[[reference_condition]] + all_results$logFC_shrunk
        message("Mean shrink weight: ", round(mean(shrunk$weight, na.rm = TRUE), 3),
                "; fraction of TFs with w < 0.5: ", round(mean(shrunk$weight < 0.5, na.rm = TRUE), 3))
}

# Each plot is drawn for the unshrunk values and, if SHRINKAGE, once more for the shrunk values.
all_results_full <- all_results
plot_variants <- list(unshrunk = all_results_full)
if (shrinkage) {
        plot_variants$shrunk <- all_results_full %>%
                mutate(
                        logFC = ifelse(is.finite(logFC_shrunk), logFC_shrunk, logFC),
                        !!contrast_condition := !!sym(reference_condition) + logFC
                )
}
variant_label <- function(name) if (name == "shrunk") " [shrunk]" else if (shrinkage) " [unshrunk]" else ""

##########################################################################################
# Plot results ###########################################################################
##########################################################################################
p_threshold <- opt$pval_threshold
message(paste0("==== Applying signifficance thresholds: fdr <= ", p_threshold, "\n"))
pdf(file.path(opt$plot_output, "volcano_plots/", output_pdf_name), width = 10, height = 10)

message("==== Writing BCalm results")
message("Reference condition:", opt$reference_condition)
message("Contrast condition:", opt$contrast_condition)

plot_title <- paste(opt$contrast_condition, "vs.", opt$reference_condition, "(p.adjusted <=", p_threshold, ")")
message("Plot title:", plot_title)


for (variant in names(plot_variants)) {
all_results <- plot_variants[[variant]]
print(all_results %>%
        ggplot(aes(x = logFC, y = -log10(p_adjusted), color = sig, alpha = sig, tf=tf)) +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        geom_point() +
        # sig will be tab-blue
        scale_color_manual(values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
        scale_alpha_manual(values = c("NS" = 0.4, "Upregulated" = 1, "Downregulated" = 1)) +
        ggpubr::theme_pubr(border = T) + 
        # annotate the significant points
        geom_label_repel(data = . %>% filter(sig != "NS"), size=5, aes(label = tf), box.padding = 0.5) +
        # TeX on labs
        labs(
                x = expression(log[2]("Fold Change")),
                y = expression(-log[10]("Adjusted p-value")),
                title = paste0(plot_title, variant_label(variant)),
        ) +
        guides(color = "none", alpha = "none") +
        theme(text = element_text(size = 18)))
}
invisible(dev.off())

# Plotting lollipop plot
pdf(file.path(opt$plot_output, "lollipop_plots/", output_pdf_name), width = 21, height = 7)

for (variant in names(plot_variants)) {
all_results <- plot_variants[[variant]]
all_results %>%
        filter(!grepl("RANDOM", tf)) %>%
        mutate(color_axis = ifelse(sig == "Upregulated", "#f37f80",
               ifelse(sig == "NS", "gray30",
                   "#6495ed"
               )
           )) %>%
        arrange(desc(!!sym(reference_condition))) %>%
        mutate(tf = factor(tf, levels = unique(tf))) %>%
           distinct(tf, .keep_all = T) %>%
           arrange(desc(!!sym(reference_condition))) %>%
           mutate(tf = factor(tf, levels = unique(tf))) -> plot_df


print(plot_df %>%
            ggplot() +
            geom_segment(aes(
                y = !!sym(reference_condition), 
                yend = !!sym(contrast_condition),
                x = tf, xend = tf,
                color = sig
            ), size = 1) +
            scale_color_manual(
                values = c("NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80"),
                name = "Result of\nComparative\nAnalysis"
            ) +
            new_scale_color() +
            geom_point(aes(x = tf, y = !!sym(reference_condition), color = "C"), size = 3) +
            geom_point(aes(x = tf, y = !!sym(contrast_condition), color = sig), size = 3) +
            scale_color_manual(values = c('C'= "black", "NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
            guides(color = "none") +
            theme_bw() +
            theme(
                axis.text.x = element_text(
                    angle = 90, hjust = 1, vjust = 0.5,
                    color = plot_df$color_axis
                ),
                # Remove inner lines of y-axis
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                text = element_text(size = 14)
            ) +
            labs(
                y = "Activity (log2(RPM+1))",
                x = "",
                                title= paste0(contrast_condition, " vs. ", reference_condition, variant_label(variant))
            ))
}
invisible(dev.off())

pdf(file.path(opt$plot_output, "circular_lollipop_plots/", output_pdf_name), width = 14, height = 14)

for (variant in names(plot_variants)) {
all_results <- plot_variants[[variant]]
circular_plot_df <-
        all_results %>%
        filter(!grepl("RANDOM", tf)) %>%
        mutate(tf_lookup = trimws(as.character(tf))) %>%
        left_join(
                if (!is.null(tf_function_map)) {
                        tf_function_map %>%
                                mutate(tf_lookup = trimws(as.character(tf))) %>%
                                select(tf_lookup, biological_function)
                } else {
                        data.frame(tf_lookup = character(0), biological_function = character(0))
                },
                by = "tf_lookup"
        ) %>%
        mutate(
                biological_function = ifelse(is.na(biological_function), "Unknown", biological_function),
                biological_function_order = match(biological_function, c(tf_function_order, "Unknown")),
                biological_function_order = ifelse(is.na(biological_function_order), length(c(tf_function_order, "Unknown")) + 1, biological_function_order),
                biological_function = factor(
                        biological_function,
                        levels = c(tf_function_order, setdiff("Unknown", tf_function_order))
                ),
                tf_map_order = match(tf_lookup, tf_order_in_map),
                tf_map_order = ifelse(is.na(tf_map_order), Inf, tf_map_order)
        ) %>%
        arrange(biological_function_order, tf_map_order, tf_lookup) %>%
        distinct(tf, .keep_all = TRUE) %>%
        mutate(color_axis = ifelse(sig == "Upregulated", "#f37f80",
                                    ifelse(sig == "NS", "gray30", "#6495ed"))
        ) %>%
        mutate(x_index = row_number())

if (nrow(circular_plot_df) > 0) {
                activity_limits <- range(
                        circular_plot_df$logFC,
                        na.rm = TRUE
                )
                base_radius <- abs(min(activity_limits, na.rm = TRUE)) + 1
                lower_fold_limit <- if (activity_limits[1] < -3) activity_limits[1] - 1 else -3
                upper_fold_limit <- if (activity_limits[2] > 3) activity_limits[2] + 1 else 3
                y_breaks <- seq(lower_fold_limit, upper_fold_limit, by = 1)
                circular_plot_df <- circular_plot_df %>%
                  mutate(
                    fold_change_radius = base_radius + logFC,
                    x_index = row_number(),
                    y_axis_label = "Fold change (log2)",
                    sig_group = cumsum(c(TRUE, sig[-1] != sig[-n()])),
                    function_group = 
                      cumsum(c(TRUE, as.character(biological_function)[-1] != as.character(biological_function)[-n()]))
                        )

                n_tf <- nrow(circular_plot_df)
                circular_plot_df <- circular_plot_df %>%
                        mutate(
                                bar_xmin = x_index - 0.45,
                                bar_xmax = x_index + 0.45,
                                bar_ymin = pmin(base_radius, fold_change_radius),
                                bar_ymax = pmax(base_radius, fold_change_radius)
                        )

                leader_line_df <- circular_plot_df %>%
                        mutate(label_radius = base_radius + upper_fold_limit + 1.15)
                
                label_df <- circular_plot_df %>%
                  mutate(
                    label_radius = base_radius + upper_fold_limit + 1.15,
                    text_radius = label_radius + 0.35,
                    label_angle_raw = 90 - 360 * (x_index - 0.5) / nrow(circular_plot_df),
                    label_angle = ifelse(label_angle_raw < -90, label_angle_raw + 180, label_angle_raw),
                    label_hjust = ifelse(label_angle_raw < -90, 1, 0)
                  )
                
                # small stub segments for functions with only one TF, since geom_path
                # needs at least 2 points per group to draw anything
                group_sizes <- label_df %>% dplyr::count(function_group)
                singleton_groups <- group_sizes$function_group[group_sizes$n == 1]
                
                singleton_df <- label_df %>%
                  dplyr::filter(function_group %in% singleton_groups) %>%
                  mutate(x_start = x_index - 0.3, x_end = x_index + 0.3)

                function_levels <- unique(as.character(circular_plot_df$biological_function))
                function_levels <- function_levels[!is.na(function_levels)]
                function_colors <- setNames(grDevices::hcl.colors(length(function_levels), "Dark 3"), function_levels)
                if ("Unknown" %in% names(function_colors)) {
                        function_colors["Unknown"] <- "grey50"
                }

                # Save a standalone legend for TF annotation group colors.
                circular_legend_pdf <- file.path(opt$plot_output, "circular_lollipop_plots", "legend.pdf")
                legend_df <- data.frame(
                        biological_function = factor(function_levels, levels = function_levels),
                        x = 1,
                        y = 1
                )
                pdf(circular_legend_pdf, width = 4.5, height = 3)
                print(
                        ggplot(legend_df, aes(x = x, y = y, color = biological_function)) +
                                geom_point(size = 4) +
                                scale_color_manual(values = function_colors) +
                                guides(color = guide_legend(title = "TF Annotation Group")) +
                                theme_void() +
                                theme(legend.position = "center")
                )
                invisible(dev.off())

        p_circular <- circular_plot_df %>%
                ggplot() +
                                geom_segment(data = leader_line_df, aes(x = x_index, xend = x_index, y = fold_change_radius, yend = label_radius), color = "grey85", linewidth = 0.3) +
                                geom_rect(aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bar_ymin, ymax = bar_ymax, fill = sig), color = "grey35", linewidth = 0.15) +
          geom_path(data = label_df, aes(x = x_index, y = label_radius, group = function_group, color = biological_function),
                    linewidth = 2.2, lineend = "round") +
          geom_segment(data = singleton_df, aes(x = x_start, xend = x_end, y = label_radius, yend = label_radius, color = biological_function),
                       linewidth = 2.2, lineend = "round") +
          geom_text(data = label_df, aes(x = x_index, y = text_radius, label = tf,
                                         angle = label_angle, hjust = label_hjust),
                    size = 6, vjust = 0.5, color = "grey15") +
                scale_fill_manual(values = c("NS" = "grey80", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
                                scale_color_manual(values = function_colors, guide = "none") +
                                guides(fill = "none") +
                                scale_x_continuous(
                                        breaks = circular_plot_df$x_index,
                                        limits = c(0.5, max(circular_plot_df$x_index) + 0.5),
                                        expand = expansion(mult = c(0, 0))
                                ) +
          scale_y_continuous(
            limits = c(base_radius + lower_fold_limit, base_radius + upper_fold_limit + 1.65),
            breaks = base_radius + y_breaks,
            labels = y_breaks,
            expand = expansion(mult = c(0, 0))
          ) +
                coord_polar(theta = "x") +
                theme_bw() +
                theme(
                        axis.text.x = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        panel.border = element_blank(),
                        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
                        text = element_text(size = 14)
                ) +
                labs(
                                        y = "Fold change (log2)",
                        x = "",
                        title = paste0(contrast_condition, " vs. ", reference_condition, variant_label(variant))
                )

        print(p_circular)
}
}

invisible(dev.off())

all_results <- all_results_full
all_results %>%
        select(-AveExpr, -t, -adj.P.Val, -B) %>%
        write.table(file.path(opt$plot_output, 'output_data', output_txt_name), sep = "\t", quote = FALSE, row.names = FALSE)