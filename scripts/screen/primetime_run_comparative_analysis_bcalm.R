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
        make_option(c("--split_by_promoter"), type = "logical", default = TRUE, help = "Whether to split the analysis by promoter or not"),
        make_option(c("--normalize"), type = "logical", default = TRUE, help = "Whether to normalize activities by promoter-specific negative controls")
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


##########################################################################################
# Preparing the data #####################################################################
##########################################################################################

pdna <- read.table(opt$pdna, header = TRUE, sep = "\t", check.names = FALSE)
name_of_the_pdna_replicate <- colnames(pdna) %>% setdiff(c("barcode", "negative_control"))

contrast_condition <- opt$contrast_condition
reference_condition <- opt$reference_condition
normalize_counts <- isTRUE(opt$normalize)
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

cdna <- read.table(opt$cdna, header = TRUE, sep = "\t", check.names = FALSE) %>%
        select(
                tf,
                negative_control,
                promoter,
                barcode,
                all_of(c(ref_replicates, contrast_replicates))
        )

# Now, biuld the design df, with columns: sample (conditions and pDNA), replicate (names ), pDNA (whether is pDNA or not) and treatment (condition name; None for pDNA)
design_df <- data.frame(
        replicate = c(ref_replicates, contrast_replicates),
        treatment = c(
                rep(reference_condition, length(ref_replicates)),
                rep(contrast_condition, length(contrast_replicates))
        ),
        sample = c(
                rep(reference_condition, length(ref_replicates)),
                rep(contrast_condition, length(contrast_replicates))
        )
)

output_pdf_name <- paste0(contrast_condition, "_vs_", reference_condition, ".pdf")
output_txt_name <- paste0(contrast_condition, "_vs_", reference_condition, ".txt")
p_threshold <- opt$pval_threshold

dir.create(file.path(opt$plot_output, "volcano_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$plot_output, "lollipop_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$plot_output, "circular_lollipop_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$plot_output, "output_data"), showWarnings = FALSE, recursive = TRUE)

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
        # print(barcode)
        # print('===========')
        # See which column of parsed_cdna has the same barcode and return them
        columns_w_same_bc <- columns_of_parsed_cdna[columns_of_parsed_cdna %>% grep(barcode, .)]
        # print(columns_w_same_bc)
        for (col_w_same_bc in columns_w_same_bc) {
                parsed_pdna[[col_w_same_bc]] <- as.numeric(tmp_parsed_pdna[[col]])
        }
}
parsed_pdna <- do.call(cbind, parsed_pdna) %>% as.data.frame()

# ============ SPLIT BY PROMOTERS ========================================================
# ========================================================================================
all_promoters = unique(parsed_cdna$promoter)
all_results = data.frame()
if (opt$split_by_promoter) {
        message("==== Splitting analysis by promoter\n")
        for (this_promoter in all_promoters) {
                message(paste0("---- Analyzing promoter: ", this_promoter, "\n"))
                this_cdna <- parsed_cdna %>%
                        filter(promoter == this_promoter) %>%
                        select(-promoter)
                this_pdna <- parsed_pdna %>%
                        filter(promoter == this_promoter) %>%
                        select(-promoter)
                this_tmp_parsed_pdna <- tmp_parsed_pdna %>% filter(promoter == this_promoter)

                # Updating rownames
                rownames(this_pdna) <- this_tmp_parsed_pdna$tf
                rownames(this_cdna) <- this_cdna$tf

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

                rownames(this_pdna_matrix) = this_tmp_parsed_pdna$tf
                rownames(this_cdna_matrix) <- this_cdna$tf


                print("this pDNA matrix")
                print(head(this_pdna_matrix))
                print("This cDNA matrix")
                print(head(this_cdna_matrix))

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

                # Now create the design data
                design_bcalm =
                        data.frame(
                                intcpt = 1,
                                grepl(contrast_condition, ordered_annotation$treatment)
                        )

                colnames(design_bcalm) = c("intcpt", contrast_condition)

                block_vector = sapply(colnames(this_pdna_matrix),
                        FUN = function(x) {
                                y = strsplit(x, "_")[[1]]
                                as.numeric(y[length(y) - 1])
                        }
                )

                print("Design BCalm")
                print(design_bcalm)

                mpralm_fit_var <- mpralm(
                        object = BcVariantMPRASet,
                        design = design_bcalm,
                        aggregate = "none",
                        normalize = TRUE,
                        model_type = "indep_groups",
                        plot = FALSE,
                        block = block_vector
                )
                # For the real comparison, retrieve the second coefficient of the model
                results_bcalm = topTable(mpralm_fit_var, coef = 2, number = Inf)

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

                all_results = rbind(all_results, results_bcalm)
        }
} else {
        #### NOT SPLITTING BY PROMOTER - SINGLE MODEL FOR ALL THE PROMOTERS ####
        message("==== NOT Splitting analysis by promoter - Single model for all the promoters\n")
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

        rownames(this_pdna_matrix) = this_tmp_parsed_pdna$tf
        rownames(this_cdna_matrix) <- this_cdna$tf


        print("this pDNA matrix")
        print(head(this_pdna_matrix))
        print("This cDNA matrix")
        print(head(this_cdna_matrix))

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

        # Now create the design data
        design_bcalm =
                data.frame(
                        intcpt = 1,
                        grepl(contrast_condition, ordered_annotation$treatment)
                )

        colnames(design_bcalm) = c("intcpt", contrast_condition)

        block_vector = sapply(colnames(this_pdna_matrix),
                FUN = function(x) {
                        y = strsplit(x, "_")[[1]]
                        as.numeric(y[length(y) - 1])
                }
        )

        print("Design BCalm")
        print(design_bcalm)

        mpralm_fit_var <- mpralm(
                object = BcVariantMPRASet,
                design = design_bcalm,
                aggregate = "none",
                normalize = TRUE,
                model_type = "indep_groups",
                plot = FALSE,
                block = block_vector
        )
        # For the real comparison, retrieve the second coefficient of the model
        results_bcalm = topTable(mpralm_fit_var, coef = 2, number = Inf)

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

        all_results = rbind(all_results, results_bcalm)
}


if (normalize_counts) {
                                # Correcting the activities by dividing by the median of the negative controls per promoter
                                message("==== Correcting activities by negative controls median")
                                # Get the median of the negative controls per condition and promoter
                                tf_promoter_map <- df %>% distinct(tf, promoter)

                                negative_control_medians <- all_results %>%
                                        filter(grepl("RANDOM", tf)) %>%
                                        left_join(tf_promoter_map, by = "tf") %>%
                                        group_by(promoter) %>%
                                        summarise(
                                                median_reference = median(!!sym(reference_condition), na.rm = TRUE),
                                                median_contrast  = median(!!sym(contrast_condition),  na.rm = TRUE),
                                                .groups = "drop"
                                        )

                                all_results <- all_results %>%
                                        left_join(tf_promoter_map, by = "tf") %>%                # add promoter column
                                        left_join(negative_control_medians, by = "promoter") %>%
                                        mutate(
                                                !!reference_condition := !!sym(reference_condition) - coalesce(median_reference, 0),
                                                !!contrast_condition := !!sym(contrast_condition) - coalesce(median_contrast, 0)
                                        ) %>%
                                        select(-median_reference, -median_contrast, -promoter)


                                # Compute logFC again after correction
                                all_results <- all_results %>%
                                                                mutate(
                                                                                                # first save the old logFC
                                                                                                old_logFC = logFC,
                                                                                                # then compute the new logFC
                                                                                                logFC = !!sym(contrast_condition) - !!sym(reference_condition)
                                                                )
} else {
                                all_results <- all_results %>% mutate(old_logFC = logFC)
}

# Correcting the p-values for multiple testing
message("==== Correcting p-values")
all_results <- all_results %>%
        mutate(
                # adjust p-values
                p_adjusted = p.adjust(P.Value, method = "BH"),
                # define significance: original BCalm logFC should be in the correct direction, new logFC should have a certain magnitude
                sig = ifelse(logFC > logfc_threshold & old_logFC > 0 & p_adjusted <= p_threshold, "Upregulated",
                        ifelse(logFC < -logfc_threshold & old_logFC < 0 & p_adjusted <= p_threshold, "Downregulated",
                                "NS"
                        )
                )
        )

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


all_results %>%
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
                title = plot_title,
        ) +
        guides(color = "none", alpha = "none") +
        theme(text = element_text(size = 18))
invisible(dev.off())

# Plotting lollipop plot
pdf(file.path(opt$plot_output, "lollipop_plots/", output_pdf_name), width = 21, height = 7)

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


plot_df %>%
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
                                title= paste(contrast_condition, "vs.", reference_condition)
            )
invisible(dev.off())

pdf(file.path(opt$plot_output, "circular_lollipop_plots/", output_pdf_name), width = 14, height = 14)

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
                        title = paste(contrast_condition, "vs.", reference_condition)
                )

        print(p_circular)
}

invisible(dev.off())

all_results %>%
        select(-AveExpr, -t, -adj.P.Val, -B) %>%
        write.table(file.path(opt$plot_output, 'output_data', output_txt_name), sep = "\t", quote = FALSE, row.names = FALSE)
