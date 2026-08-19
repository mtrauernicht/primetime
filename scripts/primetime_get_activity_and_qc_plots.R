# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script generates activity and quality control (QC) plots for the TF reporter data.
# It reads the barcode counts, processes the data, and generates various plots
# including correlation plots, density plots, and scatter plots.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(GGally)
    library(tidyr)
    library(ggridges)
    library(ggbeeswarm)
    library(stringr)
    library(patchwork)
    library(ggrepel)
    library(rlang)
    library(purrr)
    library(ggrastr)
})
options(
    dplyr.width = Inf,
    error = rlang::entrace,
    show.error.locations = TRUE
)
# Read the arguments: list of dfs for this sample, this sample name, and output directory
option_list <- list(
    make_option(c("--list_of_annotated_files"), help = "Path to all annotated files separated by commas", type = "character"),
    make_option(c("--plots_basedir"), type = "character", help = "Basedir for the plots"),
    make_option(c("--activity_basedir"), type = "character", help = "Basedir for the MPRAnalyze files"),
    make_option(c("--design"), type = "character", help = "Design DF with sample names"),
    make_option(c("--expected_pdna"), type = "character", help = "Path to expected pDNA counts"),
    make_option(c("--cdna_output"), type = "character", help = "Path to save the cDNA counts for MPRAnalyze"),
    make_option(c("--barcode_activity_output"), type = "character", help = "Path to save barcode-level activity used for barcode correlations"),
    make_option(c("--viability_file"), type = "character", default = "", help = "Optional comma-separated viability matrix files")
)

# Functions for the plots

upper_diag_plot <- function(data, mapping, color = I("black"), sizeRange = c(1, 3), ...) {
    boundaries <- seq(from = 0.8, by = 0.05, length.out = 4)

    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
    if (!is.finite(r)) {
        r <- NA_real_
    }
    rt <- if (is.na(r)) "NA" else format(r, digits = 3)
    tt <- as.character(rt)
    cex <- max(sizeRange)

    # helper function to calculate a useable size
    percent_of_range <- function(percent, range) {
        percent * diff(range) + min(range, na.rm = TRUE)
    }

    # plot correlation coefficient
    p <- ggally_text(
        label = tt, mapping = aes(), xP = 0.5, yP = 0.5,
        size = I(percent_of_range(cex * ifelse(is.na(r), 0, abs(r)), sizeRange)) + 5, color = color, ...
    ) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
        )

    corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")[2:6]

    if (is.na(r)) {
        corCol <- "grey90"
    } else if (r <= boundaries[1]) {
        corCol <- corColors[1]
    } else if (r <= boundaries[2]) {
        corCol <- corColors[2]
    } else if (r < boundaries[3]) {
        corCol <- corColors[3]
    } else if (r < boundaries[4]) {
        corCol <- corColors[4]
    } else {
        corCol <- corColors[5]
    }

    p <- p +
        theme(panel.background = element_rect(fill = corCol))

    return(p)
}


lower_diag_plot <- function(data, mapping, ...) {
    ggally_points(data = data, mapping = mapping, alpha = 0.5, size = 0.7) +
        geom_abline(slope = 1, lty = "dashed", col = "red") +
        theme_pubr(border = T)
}

diag_plot <- function(data, mapping, ...) {
    ggally_densityDiag(data = data, mapping = mapping, alpha = 0.3, fill = "red") +
        theme_pubr(border = T)
}

corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")[2:6]

counts_df <- data.frame()

normalize_well_id <- function(value) {
    value <- toupper(trimws(as.character(value)))
    value <- gsub("[^A-Z0-9]", "", value)
    value
}

read_viability_matrix <- function(path) {
    if (is.null(path) || !nzchar(trimws(path)) || !file.exists(path) || file.info(path)$size == 0) {
        return(NULL)
    }

    viability_df <- tryCatch(
        read.delim(path, sep = ";", header = TRUE, check.names = FALSE, quote = "", comment.char = "", na.strings = c("", "NA")),
        error = function(e) data.frame()
    )
    if (is.null(viability_df) || ncol(viability_df) < 2) {
        return(data.frame(well = character(0), viability = numeric(0)))
    }

    colnames(viability_df)[1] <- "well_row"
    viability_long <- viability_df %>%
        pivot_longer(-well_row, names_to = "well_column", values_to = "viability") %>%
        mutate(
            well_row = normalize_well_id(well_row),
            well_column = normalize_well_id(well_column),
            well = paste0(well_row, well_column),
            viability = suppressWarnings(as.numeric(viability))
        ) %>%
        filter(well != "", !is.na(viability), is.finite(viability)) %>%
        distinct(well, .keep_all = TRUE)

    if (nrow(viability_long) == 0) {
        return(data.frame(well = character(0), viability = numeric(0)))
    }

    viability_long
}

read_design_map <- function(path) {
    if (is.null(path) || !nzchar(trimws(path)) || !file.exists(path) || file.info(path)$size == 0) {
        return(NULL)
    }

    first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character(0))
    if (length(first_line) == 0) {
        return(NULL)
    }

    has_header <- grepl("sample|well|row|column|plate", tolower(first_line[1]))
    design_df <- tryCatch(
        read.csv(path, header = has_header, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )
    if (is.null(design_df) || ncol(design_df) < 2) {
        return(NULL)
    }

    lower_names <- tolower(names(design_df))
    well_col <- which(lower_names %in% c("well", "wells", "plate_well"))[1]
    sample_col <- which(lower_names %in% c("sample", "condition", "name", "replicate"))[1]
    row_col <- which(lower_names %in% c("row", "well_row"))[1]
    column_col <- which(lower_names %in% c("column", "col", "well_column"))[1]

    if (!is.na(well_col) && !is.na(sample_col)) {
        return(
            design_df %>%
                transmute(
                    sample = trimws(as.character(.data[[names(design_df)[sample_col]]])),
                    well = normalize_well_id(.data[[names(design_df)[well_col]]])
                ) %>%
                filter(sample != "", well != "") %>%
                distinct(sample, well, .keep_all = TRUE)
        )
    }

    if (!is.na(row_col) && !is.na(column_col) && !is.na(sample_col)) {
        return(
            design_df %>%
                transmute(
                    sample = trimws(as.character(.data[[names(design_df)[sample_col]]])),
                    well = normalize_well_id(paste0(.data[[names(design_df)[row_col]]], .data[[names(design_df)[column_col]]]))
                ) %>%
                filter(sample != "", well != "") %>%
                distinct(sample, well, .keep_all = TRUE)
        )
    }

    design_df <- design_df[, 1:2]
    colnames(design_df) <- c("well", "sample")
    design_df %>%
        mutate(
            sample = trimws(as.character(sample)),
            well = normalize_well_id(well)
        ) %>%
        filter(sample != "", well != "") %>%
        distinct(sample, well, .keep_all = TRUE)
}

##########################################################################################
## Start of the script ###################################################################
##########################################################################################

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$plots_basedir) && nzchar(opt$plots_basedir)) {
    dir.create(opt$plots_basedir, recursive = TRUE, showWarnings = FALSE)
}
if (!is.null(opt$activity_basedir) && nzchar(opt$activity_basedir)) {
    dir.create(opt$activity_basedir, recursive = TRUE, showWarnings = FALSE)
}

list_of_annotated_files <- if (file.exists(opt$list_of_annotated_files)) {
    readLines(opt$list_of_annotated_files, warn = FALSE)
} else {
    strsplit(opt$list_of_annotated_files, "[[:space:],]+")[[1]]
}
list_of_annotated_files <- trimws(list_of_annotated_files)
list_of_annotated_files <- list_of_annotated_files[list_of_annotated_files != ""]
unique_files <- unique(list_of_annotated_files)

counts_df = data.frame()
for(l in unique_files){
    # message("- Processing file ", l)
    basename = tools::file_path_sans_ext(basename(l))
    # Remove the .cluster.annotated part from the end
    sample_plus_replicate = sub("\\.cluster\\.annotated$", "", basename)
    # Split at last underscore to get sample and replicate
    underscore_pos = max(unlist(gregexpr("_", sample_plus_replicate)))
    this_sample = substr(sample_plus_replicate, 1, underscore_pos - 1)
    replicate = substr(sample_plus_replicate, underscore_pos + 1, nchar(sample_plus_replicate))
    this_replicate = paste(this_sample, replicate, sep = "_")
    this_replicate_df <- read.table(l, header = TRUE, sep = "\t")
    count_col <- if ("count" %in% colnames(this_replicate_df)) {
        "count"
    } else if ("raw_count" %in% colnames(this_replicate_df)) {
        "raw_count"
    } else {
        stop("Annotated file is missing a count column: ", l)
    }
    neg_ctrl_col <- if ("neg_ctrls" %in% colnames(this_replicate_df)) {
        "neg_ctrls"
    } else if ("negative_control" %in% colnames(this_replicate_df)) {
        "negative_control"
    } else {
        stop("Annotated file is missing a negative-control column: ", l)
    }
    # ADD PSEUDOCOUNT OF 1!!!!
    this_replicate_df[[count_col]] <- as.numeric(this_replicate_df[[count_col]]) + 1
    is_pDNA = this_sample == "pDNA"
    total_read_count = sum(this_replicate_df[[count_col]])
    this_df = data.frame(
        sample = this_sample,
        replicate = this_replicate,
        pDNA = is_pDNA,
        tf = this_replicate_df$tf,
        negative_control = ifelse(this_replicate_df[[neg_ctrl_col]] == "Yes", T, F),
        promoter = this_replicate_df$promoter,
        barcode = this_replicate_df$barcode,
        raw_count = this_replicate_df[[count_col]],
        log2_count = log2(this_replicate_df[[count_col]]),
        total_read_count=total_read_count,
        RPM = this_replicate_df[[count_col]] / total_read_count * 1e6
    )
    # ---- Update counts ---------------------------------------------------
    counts_df <- rbind(
        counts_df,
        this_df
    )

}

# ==============================================================================

# Get number of obvservations: barcodes x TF x promoters
counts_df %>%
    group_by(replicate, sample) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    head(1) %>%
    pull(n_obs) -> n_observations
# Ideal number of reads is 10 times the number of barcodes
ideal_reads <- 100 * n_observations[1]

# Get the pDNA rows and make it wider
pDNA <- counts_df %>%
    # mutate(RPM = log2(RPM)) %>%
    filter(pDNA == T) %>%
    select(-pDNA, -tf, -promoter, -total_read_count, -log2_count) %>%
    # Take the average of the RPMs in case of multiple replicates
    dplyr::group_by(replicate, barcode) %>%
    dplyr::reframe(mean_RPM = mean(RPM), sample = sample, barcode = barcode) %>%
    pivot_wider(names_from = sample, values_from = mean_RPM) %>%
    # Assign 1 to NAs
    mutate_all(~ ifelse(is.na(.), 1, .)) %>%
    select(-replicate) %>%
    ##### VF250404: TAKE THE MEAN OF THE pDNA REPLICATES ## NOW DEALS WITH MORE THAN 1 pDNA rep
    group_by(barcode) %>%
    summarise(pDNA = mean(pDNA))

# print('Check if any NA in pDNA')
# pDNA %>%
#     summarise(across(everything(), ~ sum(is.na(.)))) %>%
#     pivot_longer(everything(), names_to = "sample", values_to = "n_NA") %>%
#     filter(n_NA > 0) %>% print()

cDNA <- counts_df %>%
    filter(pDNA == F) %>%
    select(-sample, -pDNA, -log2_count, -total_read_count, -raw_count) %>%
    pivot_wider(
        names_from = replicate,
        values_from = RPM,
        id_cols = c('barcode', 'tf', 'negative_control', 'promoter')
    ) %>%
    # Assign 1 to NAs
    mutate_all(~ ifelse(is.na(.), 1, .))

# Save the cDNA df for MPRAnalyze (raw_counts)
# message ('============= savign CDNA ==============')
message("==== Saving cDNA and pDNA raw counts (for comparative analysis)")
counts_df %>%
    filter(pDNA == F) %>%
    select(-sample, -pDNA, -total_read_count, -log2_count, -RPM) %>%
    pivot_wider(names_from = replicate, values_from = raw_count) %>%
    # Assign 1 to NAs
    mutate_all(~ ifelse(is.na(.), 1, .)) %>%
    write.table(
        file = opt$cdna_output,
        row.names = FALSE, quote = F, sep = "\t"
    )
# Save the pDNA df for MPRAnalyze (raw_counts)
# message ('============= savign PDNA ==============')
counts_df %>%
    filter(pDNA == T) %>%
    select(-pDNA, -tf, -promoter, -total_read_count, -log2_count, -RPM) %>%
    pivot_wider(names_from = sample, values_from = raw_count) %>%
    # Assign 1 to NAs
    mutate_all(~ ifelse(is.na(.), 1, .)) %>%
    select(-replicate) %>%
    write.table(
        file = file.path(opt$activity_basedir, "pDNA_counts.txt"),
        row.names = FALSE, quote = F, sep = "\t"
    )

# Merge the two
merged <-
    merge(pDNA,
        cDNA,
        by = "barcode"
    )


pdna_sample <- colnames(pDNA) %>%
    unique() %>%
    setdiff("barcode")
cdna_sample <- colnames(cDNA) %>%
    unique() %>%
    setdiff(c("barcode", "tf", "promoter", "negative_control"))
# message ("cdna_sample: ", cdna_sample)
# message ("pdna_sample: ", pdna_sample)

all_cdna_conditions =
    cdna_sample %>%
    # split by _ remove the last one and concatenate again
    str_split("_") %>%
    map_chr(~ paste(.x[-length(.x)], collapse = "_")) %>%
    unique()

# message('All cDNA conditions: ', all_cdna_conditions)

print("COUNT DF")
print(counts_df%>% head)
print("Check if any log2 values are NA or inf. If so, here they are:")
print(counts_df %>% filter(is.na(log2_count) | !is.finite(log2_count)))


# ==============================================================================
# Estimate bleed-through using the random tfs
message("==== Estimating bleed-through")
equation_df <- data.frame()
for (this_cdna_sample in cdna_sample) {
    # message ("Processing ", this_cdna_sample)
    for (this_pdna_sample in pdna_sample) {
        # message ("    Processing ", this_pdna_sample)
        # Get the equation of the line
        this_df <- merged %>% filter(negative_control)
        x <- this_df[[this_pdna_sample]]
        y <- this_df[[this_cdna_sample]]
        # Calculate the slope and intercept: y = slope * x + intercept
        fit <- lm(y ~ x)
        slope <- coef(fit)[2]
        intercept <- coef(fit)[1]
        # Add to the equation df -----------------------------------------------
        equation_df <- rbind(
            equation_df,
            data.frame(
                cDNA_sample = this_cdna_sample,
                pDNA_sample = this_pdna_sample,
                slope = slope,
                intercept = intercept
            )
        )
    }
}

# Merge the cDNA and pDNA dfs with the equation df. For this, needs to melt the cDNA df first

bleed_through_slope_df <-
    cDNA %>%
    # filter(negative_control) %>%
    pivot_longer(
        cols = -c(barcode, tf, negative_control, promoter),
        names_to = "sample",
        values_to = "cDNA_RPM"
    ) %>%
    merge(.,
        pDNA,
        by = "barcode"
    ) %>%
    left_join(equation_df %>% select(-pDNA_sample),
        by = c("sample" = "cDNA_sample")
    ) %>%
    filter(negative_control)

bleedthrough_output <- file.path(opt$activity_basedir, "bleedthrough_per_condition.txt")
write.table(
    bleed_through_slope_df %>%
        select(sample, slope, intercept) %>%
        distinct(),
    file = bleedthrough_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)



message("==== Plotting distribution of BC counts")

# Plot the BC counts: as ggridges density plots ==============================================
pdf(file.path(opt$plots_basedir, "distribution_of_BC_counts.pdf"),
    width = 10,
    height = 0.2*length(cdna_sample)
)
order_of_replicates =
    counts_df %>%
    group_by(replicate) %>%
    summarise(median_log2 = median(log2_count)) %>%
    arrange(median_log2) %>%
    pull(replicate)

# Alarm colouring: will check how many barcodes with low counts (10) each replicate has.
#                  if more than 10, the sample will be highlighted in red
counts_df %>%
    group_by(replicate) %>%
    summarise(low_count_barcodes = sum(raw_count < 5)) %>%
    mutate(alarm = ifelse(low_count_barcodes > 10, "T", "F")) %>%
    ungroup() -> alarm_df



counts_df %>%
    # Add the alarming information
    left_join(alarm_df, by = "replicate") %>%
    # Refactor to follow the order_of_replicates
    mutate(replicate = factor(replicate, levels = order_of_replicates)) %>%
        filter(!is.na(log2_count), is.finite(log2_count)) %>%
        ggplot(aes(x = log2_count, y = replicate, fill = alarm)) +
    geom_density_ridges(alpha = 0.7, scale = 3,
        rel_min_height = 0.01,
        quantile_lines = TRUE, 
        quantiles = 2
    ) +
    scale_fill_manual(values = c(
        "F" = "#7f7f7f",
        "T" = "#e15759"
    ),
    name = "More than 10\nlow count barcodes?") +
    theme_minimal() +
    xlab("BC count (log2)") +
    ylab("Replicate") +
    ggtitle("Distribution of BC counts per replicate") +
    theme(text = element_text(size = 18))
invisible(dev.off())


# Plot the read counts =========================================================
pdf(file.path(opt$plots_basedir, "read_counts.pdf"),
    width = 10, 
    height = 0.5*length(cdna_sample)
)


# Sort by total_read_count
order_of_replicates =
    counts_df %>%
    select(replicate, total_read_count) %>%
    distinct() %>%
    arrange(total_read_count) %>%
    pull(replicate)

counts_df %>%
    select(replicate, sample, total_read_count) %>%
    # Sort by total_read_count
    mutate(replicate = factor(replicate, levels = order_of_replicates)) %>%
    distinct() %>%
    ggplot(aes(
        x = replicate,
        y = total_read_count,
        fill = ifelse(total_read_count < ideal_reads, "below", "above")
    )) +
    geom_bar(stat = "identity", alpha = 0.5) +
    geom_text(
        aes(
            label = total_read_count,
            color = ifelse(total_read_count < ideal_reads, "below", "above")
        ),
        size = 5, fontface = "bold"
    ) +
    coord_flip() +
    theme_pubr() +
    # Add line for ideal number of reads
    geom_hline(yintercept = ideal_reads, linetype = "dashed", color = "black") +
    xlab("") +
    scale_fill_manual(
        values = c(
            "above" = "#7f7f7f",
            "below" = "#e15759"
        ),
        aesthetics = c("fill", "color"),
    )+
    ylab("Read count") +
    ggtitle("Read counts per sample") +
    theme(
        text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none"
    )

invisible(dev.off())

read_count_summary_df <- counts_df %>%
    filter(!pDNA) %>%
    select(replicate, sample, total_read_count) %>%
    distinct() %>%
    mutate(
        condition = sub("_[0-9]+$", "", replicate),
        read_count_lt_25000 = total_read_count < 25000
    )

write.table(
    read_count_summary_df,
    file = file.path(opt$activity_basedir, "read_count_per_sample.tsv"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

viability_paths <- if (is.null(opt$viability_file) || !nzchar(trimws(opt$viability_file))) {
    character(0)
} else {
    trimws(unlist(strsplit(opt$viability_file, "[,;\n]+")))
}
viability_paths <- viability_paths[viability_paths != ""]
viability_matrices <- lapply(seq_along(viability_paths), function(plate_index) {
    viability_df <- read_viability_matrix(viability_paths[plate_index])
    if (is.null(viability_df)) {
        message("==== Skipping missing or empty viability file: ", viability_paths[plate_index])
        return(NULL)
    }
    viability_df %>%
        mutate(
            plate_index = .env$plate_index,
            plate_name = paste0("Plate ", .env$plate_index)
        )
})
viability_matrices <- Filter(Negate(is.null), viability_matrices)
viability_matrices <- lapply(seq_along(viability_matrices), function(plate_index) {
    viability_matrices[[plate_index]] %>%
        mutate(
            plate_index = .env$plate_index,
            plate_name = paste0("Plate ", .env$plate_index)
        )
})

if (length(viability_matrices) > 0) {
    message("==== Plotting viability heatmaps")
    pdf(file.path(opt$plots_basedir, "viability_well_heatmap.pdf"), width = 12, height = 7)
    for (viability_df in viability_matrices) {
        viability_heatmap_df <- viability_df %>%
            mutate(
                row = factor(substr(well, 1, 1), levels = LETTERS[1:16]),
                column = suppressWarnings(as.integer(sub("^[A-Za-z]+", "", well)))
            ) %>%
            filter(!is.na(column), !is.na(row)) %>%
            complete(row, column = 1:24)

        print(
            ggplot(viability_heatmap_df, aes(x = column, y = row, fill = viability)) +
                geom_tile(color = "white", linewidth = 0.25) +
                scale_x_continuous(breaks = 1:24, expand = c(0, 0)) +
                scale_y_discrete(limits = rev(levels(viability_heatmap_df$row)), expand = c(0, 0)) +
                scale_fill_gradientn(colors = viridisLite::viridis(256), na.value = "grey90") +
                coord_fixed() +
                theme_pubr(border = T) +
                labs(
                    title = paste("Viability per well -", unique(viability_df$plate_name)),
                    x = "Column",
                    y = "Row",
                    fill = "Viability"
                ) +
                theme(
                    panel.grid = element_blank(),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.text = element_text(size = 10),
                    legend.key.width = grid::unit(1.2, "cm")
                )
        )
    }
    invisible(dev.off())

    design_map <- read_design_map(opt$design)
    read_count_by_replicate_df <- counts_df %>%
        filter(!pDNA) %>%
        group_by(sample, replicate) %>%
        summarise(total_read_count = first(total_read_count), .groups = "drop") %>%
        mutate(plate_index = suppressWarnings(as.integer(sub(".*_([0-9]+)$", "\\1", replicate)))) %>%
        filter(total_read_count > 0)

    if (length(viability_matrices) == 1) {
        read_count_wells_df <- read_count_by_replicate_df %>%
            group_by(sample) %>%
            summarise(total_read_count = first(total_read_count), .groups = "drop") %>%
            mutate(plate_index = 1L)
    } else {
        read_count_wells_df <- read_count_by_replicate_df %>%
            filter(plate_index %in% seq_along(viability_matrices))
    }

    if (!is.null(design_map)) {
        read_count_wells_df <- read_count_wells_df %>%
            inner_join(design_map, by = "sample") %>%
            distinct(plate_index, well, .keep_all = TRUE)
        message("==== Viability design map rows: ", nrow(design_map))
    } else {
        read_count_wells_df <- read_count_wells_df %>%
            mutate(well = normalize_well_id(sample)) %>%
            filter(well != "") %>%
            distinct(plate_index, well, .keep_all = TRUE)
        message("==== Viability design map unavailable; falling back to sample-derived wells")
    }

    viability_read_count_df <- map_dfr(viability_matrices, function(viability_df) {
        this_plate_index <- unique(viability_df$plate_index)
        read_count_wells_df %>%
            filter(plate_index == this_plate_index) %>%
            inner_join(viability_df %>% select(well, viability), by = "well") %>%
            mutate(plate_name = unique(viability_df$plate_name))
    }) %>%
        group_by(plate_index) %>%
        mutate(
            dmso_mean_viability = mean(viability[grepl("dmso", sample, ignore.case = TRUE)], na.rm = TRUE),
            dmso_mean_viability = ifelse(!is.finite(dmso_mean_viability) | dmso_mean_viability <= 0, 1, dmso_mean_viability),
            viability_rel_dmso = viability / dmso_mean_viability
        ) %>%
        ungroup()

    message("==== Plotting read count versus viability")
    pdf(file.path(opt$plots_basedir, "read_counts_vs_viability.pdf"), width = 8, height = 6)
    for (this_plate_name in unique(viability_read_count_df$plate_name)) {
        this_plate_df <- viability_read_count_df %>% filter(.data$plate_name == this_plate_name)
        max_viability_rel <- max(this_plate_df$viability_rel_dmso, na.rm = TRUE)
        if (!is.finite(max_viability_rel) || max_viability_rel <= 0) {
            max_viability_rel <- 1
        }
        if (nrow(this_plate_df) >= 2) {
            outlier_fit <- tryCatch(
                lm(log1p(total_read_count) ~ viability_rel_dmso, data = this_plate_df),
                error = function(e) NULL
            )
            if (!is.null(outlier_fit)) {
                y_outlier_score <- abs(rstandard(outlier_fit))
                x_outlier_score <- abs(as.numeric(scale(this_plate_df$viability_rel_dmso)))
                x_outlier_score[!is.finite(x_outlier_score)] <- 0
                this_plate_df$outlier_score <- pmax(y_outlier_score, x_outlier_score)
                label_df <- this_plate_df %>%
                    filter(is.finite(outlier_score)) %>%
                    slice_max(order_by = outlier_score, n = 10, with_ties = FALSE)
            } else {
                label_df <- this_plate_df[0, , drop = FALSE]
            }
            print(
                ggplot(this_plate_df, aes(x = viability_rel_dmso, y = total_read_count)) +
                    geom_point(alpha = 0.35, size = 1.4) +
                    geom_text_repel(
                        data = label_df,
                        aes(label = sample),
                        max.overlaps = 10,
                        box.padding = 0.35,
                        point.padding = 0.2
                    ) +
                    geom_smooth(method = "lm", se = FALSE, color = "#4c78a8") +
                    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
                    scale_x_continuous(limits = c(0, max_viability_rel), expand = c(0, 0)) +
                    theme_pubr(border = T) +
                    labs(
                        title = paste("Read counts versus viability -", this_plate_name),
                        x = "Viability / mean DMSO viability",
                        y = "Total read count"
                    )
            )
        } else {
            plot.new()
            text(0.5, 0.5, paste("No overlapping wells found for", this_plate_name))
        }
    }
    invisible(dev.off())

    message("==== Plotting viability correlation between plates")
    pdf(file.path(opt$plots_basedir, "viability_plate_correlations.pdf"), width = 8, height = 6)
    if (length(viability_matrices) >= 2) {
        plate_pairs <- combn(seq_along(viability_matrices), 2, simplify = FALSE)
        for (plate_pair in plate_pairs) {
            first_plate <- viability_matrices[[plate_pair[1]]]
            second_plate <- viability_matrices[[plate_pair[2]]]
            correlation_df <- first_plate %>%
                select(well, viability_first = viability) %>%
                inner_join(second_plate %>% select(well, viability_second = viability), by = "well")
            if (nrow(correlation_df) >= 2) {
                if (!is.null(design_map)) {
                    correlation_df <- correlation_df %>%
                        left_join(design_map %>% select(well, sample), by = "well")
                } else {
                    correlation_df <- correlation_df %>%
                        mutate(sample = well)
                }
                correlation_fit <- tryCatch(
                    lm(viability_second ~ viability_first, data = correlation_df),
                    error = function(e) NULL
                )
                if (!is.null(correlation_fit)) {
                    y_outlier_score <- abs(rstandard(correlation_fit))
                    x_outlier_score <- abs(as.numeric(scale(correlation_df$viability_first)))
                    x_outlier_score[!is.finite(x_outlier_score)] <- 0
                    correlation_df$outlier_score <- pmax(y_outlier_score, x_outlier_score)
                    correlation_labels <- correlation_df %>%
                        filter(is.finite(outlier_score)) %>%
                        slice_max(order_by = outlier_score, n = 10, with_ties = FALSE)
                } else {
                    correlation_labels <- correlation_df[0, , drop = FALSE]
                }
                print(
                    ggplot(correlation_df, aes(x = viability_first, y = viability_second)) +
                        geom_point(alpha = 0.5, size = 1.4) +
                        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
                        geom_text_repel(
                            data = correlation_labels,
                            aes(label = sample),
                            max.overlaps = 10,
                            box.padding = 0.35,
                            point.padding = 0.2
                        ) +
                        stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
                        theme_pubr(border = T) +
                        labs(
                            title = paste("Viability correlation:", first_plate$plate_name[1], "vs", second_plate$plate_name[1]),
                            x = paste(first_plate$plate_name[1], "viability"),
                            y = paste(second_plate$plate_name[1], "viability")
                        )
                )
            } else {
                plot.new()
                text(0.5, 0.5, "No overlapping wells found between plates")
            }
        }
    } else {
        plot.new()
        text(0.5, 0.5, "At least two viability files are required")
    }
    invisible(dev.off())
} else {
    message("==== No viability file provided or file missing; writing placeholder viability QC PDFs")
    for (output_file in c("viability_well_heatmap.pdf", "read_counts_vs_viability.pdf", "viability_plate_correlations.pdf")) {
        pdf(file.path(opt$plots_basedir, output_file), width = 8, height = 6)
        plot.new()
        text(0.5, 0.5, "No viability file provided")
        invisible(dev.off())
    }
}

# Plot correlation between cDNA and pDNA =======================================


n_conditions = length(all_cdna_conditions)
# Now, plot cDNA_RPM vs pDNA facetting by sample
pdf(file.path(opt$plots_basedir, "bleedthrough_estimation.pdf"),
    width = 5, height = 5
)
bleed_through_slope_df <-
    bleed_through_slope_df %>%
    mutate(slope_col = ifelse(slope <= 0.1, "positive", ifelse(slope <= 0.2, "medium", "negative")))
bleed_through_labels <-
    bleed_through_slope_df %>%
    select(slope, slope_col, sample) %>%
    mutate(slope = round(slope * 100, 2)) %>%
    mutate(color_axis = ifelse(slope_col == "positive", corColors[5],
        ifelse(slope_col == "medium",
            corColors[2],
            corColors[1]
        )
    )) %>%
    # order by sample
    arrange(sample) %>%
    distinct()

# Bleed-through plot will be ordered by slope (higher first)
order_of_samples =
    bleed_through_slope_df %>%
    select(sample, slope) %>%
    distinct() %>%
    arrange(desc(slope)) %>%
    pull(sample)

for (smp in order_of_samples) {
    message("----- ploting ", smp)
    df <- bleed_through_slope_df %>%
        filter(sample == smp) %>%
        distinct()
    df_lbl <- bleed_through_labels %>%
        filter(sample == smp) %>%
        distinct()
    x_median <- max(df$pDNA)
    y_max <- min(df$cDNA_RPM)
    p = ggplot(
        data = df %>% filter(sample == smp)
    ) +
        aes(x = !!sym(pdna_sample), y = cDNA_RPM, color = slope_col, fill = slope_col) +
        geom_point(alpha = 0.2, color = "black") +
        geom_abline(aes(intercept = intercept, slope = slope, color = slope_col),
            linewidth = 1.5
        ) +
        facet_wrap(~sample) +
        theme_pubr(border = T) +
        # Add the slope value, very big, above the plot
        geom_label_repel(
            data = df_lbl,
            aes(label = paste0("bleedt. = ", slope, "%")),
            x = x_median,
            color = "black",
            segment.color = "white",
            y = y_max,
            hjust = 0,
            vjust = 1,
            size = 5
        ) +
        scale_color_manual(values = c(
            "positive" = corColors[5], # "#26a74a",
            "medium" = corColors[2], # "#fdc010",
            "negative" = corColors[1] # "#dc3644"
        ), aesthetics = c("color", "fill")) +
        xlab("pDNA count") +
        ylab("cDNA count") +
        guides(color = "none", fill = "none") +
        theme(
            strip.background = element_rect(colour = "black", fill = NA),
            strip.text = element_text(face = "bold", size = 12),
            strip.background.x = element_rect(fill = df_lbl$color_axis)
        )
    # Remove the y-label for all but the first plot
    if (smp != unique(bleed_through_slope_df$sample)[1]) {
        p <- p + theme(axis.title.y = element_blank())
    }
    print(p)
    }   

# n_conditions = length(unique(bleed_through_slope_df$sample))
# # use patchwork to arrange the plots
# plots %>% wrap_plots(nrow = 1) + plot_annotation(title = "Bleedthrough estimation (Percentage of cDNA counts coming from pDNA)")

invisible(dev.off())



# ==============================================================================
# Get the activity of the TFs by deviding the values of the cDNA by the values of the pDNA
# message("==== Calculating corrected activity")
activity_df <- data.frame()
# corrected_activity_df <- data.frame() # corrected for bleed-through
for (this_cdna_sample in cdna_sample) {
    # Get the values
    cDNA_values <- merged[[this_cdna_sample]]
    for (this_pdna_sample in pdna_sample) {
        # message ("Processing ", this_cdna_sample, " and ", this_pdna_sample)
        pDNA_values <- merged[[this_pdna_sample]]
        # Calculate the activity
        activity <- cDNA_values / pDNA_values

        # Add to the activity df -----------------------------------------------
        activity_df <- rbind(
            activity_df,
            data.frame(
                barcode = merged$barcode,
                tf = merged$tf,
                negative_control = merged$negative_control,
                promoter = merged$promoter,
                cDNA_sample = this_cdna_sample,
                pDNA_sample = this_pdna_sample,
                activity_RPM = activity
            )
        )
    }
}

summarise_barcode_activity <- function(df) {
    df %>%
        filter(!negative_control) %>%
        group_by(cDNA_sample, barcode, tf, promoter) %>%
        summarise(
            mean_RPM = mean(activity_RPM),
            log2_mean_RPM = log2(mean_RPM),
            .groups = "drop"
        )
}

control_condition <-
    all_cdna_conditions[str_detect(str_to_lower(all_cdna_conditions), "control|ctrl|dmso")][1]
if (is.na(control_condition) || control_condition == "") {
    control_condition <- all_cdna_conditions[1]
}

control_replicates <-
    cdna_sample[str_detect(cdna_sample, paste0("^", control_condition, "_[0-9]+$"))]
if (length(control_replicates) == 0) {
    control_replicates <- cdna_sample[str_detect(str_to_lower(cdna_sample), "control|ctrl|dmso")]
}
if (length(control_replicates) == 0) {
    control_replicates <- cdna_sample[1]
}

build_barcode_control_comparison <- function(barcode_activity_df, control_replicates) {
    barcode_activity_df %>%
        filter(!cDNA_sample %in% control_replicates) %>%
        inner_join(
            barcode_activity_df %>%
                filter(cDNA_sample %in% control_replicates) %>%
                group_by(barcode, tf, promoter) %>%
                summarise(
                    control_mean_RPM = mean(mean_RPM),
                    control_log2_mean_RPM = log2(control_mean_RPM),
                    .groups = "drop"
                ),
            by = c("barcode", "tf", "promoter")
        ) %>%
        mutate(
            comparison = cDNA_sample,
            sample_log2_mean_RPM = log2_mean_RPM,
            deviation = sample_log2_mean_RPM - control_log2_mean_RPM,
            residual = sample_log2_mean_RPM - control_log2_mean_RPM
        )
}

barcode_activity_df <- summarise_barcode_activity(activity_df)

barcode_control_comparison_df <- build_barcode_control_comparison(barcode_activity_df, control_replicates)

if (!is.null(opt$barcode_activity_output)) {
    write.table(
        barcode_activity_df,
        file = opt$barcode_activity_output,
        row.names = FALSE, quote = F, sep = "\t"
    )
}

activity_summary_df <- activity_df %>%
    group_by(cDNA_sample, barcode, tf, negative_control, promoter) %>%
    summarise(
        mean_RPM = mean(activity_RPM),
        log2_mean_RPM = log2(mean_RPM),
        .groups = "drop"
    )

barcode_control_labels <-
    barcode_control_comparison_df %>%
    group_by(comparison, tf) %>%
    slice_max(order_by = abs(deviation), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(comparison) %>%
    slice_max(order_by = abs(deviation), n = 5, with_ties = FALSE) %>%
    ungroup()

barcode_control_pages <- split(
    sort(unique(barcode_control_comparison_df$comparison)),
    ceiling(seq_along(sort(unique(barcode_control_comparison_df$comparison))) / 25)
)

message("==== Plotting barcode activities vs control")
pdf(file.path(opt$plots_basedir, "control_barcode_correlations.pdf"), width = 16, height = 16)
for (page_idx in seq_along(barcode_control_pages)) {
    page_comparisons <- barcode_control_pages[[page_idx]]
    page_df <- barcode_control_comparison_df %>%
        filter(comparison %in% page_comparisons) %>%
        mutate(comparison = factor(comparison, levels = page_comparisons))
    page_labels <- barcode_control_labels %>%
        filter(comparison %in% page_comparisons) %>%
        filter(control_log2_mean_RPM > -2 | sample_log2_mean_RPM > -2) %>%
        mutate(comparison = factor(comparison, levels = page_comparisons))

    if (nrow(page_df) == 0) {
        next
    }

    p <- ggplot(page_df, aes(x = control_log2_mean_RPM, y = sample_log2_mean_RPM)) +
        geom_point(alpha = 0.25, size = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        geom_text_repel(
            data = page_labels,
            aes(label = tf),
            size = 4,
            max.overlaps = Inf,
            min.segment.length = 0,
            box.padding = 0.05,
            point.padding = 0.05,
            segment.color = "grey50"
        ) +
        facet_wrap(~comparison, ncol = 5) +
        coord_equal() +
        theme_pubr(border = T) +
        ggtitle(paste0("Barcode activities vs ", control_condition, " control (page ", page_idx, "/", length(barcode_control_pages), ")")) +
        xlab(paste0(control_condition, " control barcode activity (log2 mean RPM)")) +
        ylab("Sample barcode activity (log2 mean RPM)") +
        theme(
            text = element_text(size = 14),
            strip.text = element_text(size = 8, face = "bold"),
            axis.text = element_text(size = 7)
        )
    print(p)
}
invisible(dev.off())

control_residual_correlation_stats <-
    barcode_control_comparison_df %>%
    filter(is.finite(control_log2_mean_RPM), is.finite(residual)) %>%
    group_by(comparison) %>%
    group_modify(~ {
        this_df <- .x
        n_points <- nrow(this_df)
        n_unique_x <- n_distinct(this_df$control_log2_mean_RPM)

        if (n_points < 10 || n_unique_x < 5) {
            return(data.frame(
                n_points = n_points,
                n_unique_x = n_unique_x,
                slope = NA_real_,
                corr = NA_real_,
                corr_p_two_sided = NA_real_,
                test_status = "insufficient_points"
            ))
        }

        fit <- tryCatch(
            lm(residual ~ control_log2_mean_RPM, data = this_df),
            error = function(e) NULL
        )
        corr_test <- tryCatch(
            cor.test(this_df$control_log2_mean_RPM, this_df$residual, method = "pearson", alternative = "two.sided"),
            error = function(e) NULL
        )

        if (is.null(fit) || is.null(corr_test)) {
            return(data.frame(
                n_points = n_points,
                n_unique_x = n_unique_x,
                slope = NA_real_,
                corr = NA_real_,
                corr_p_two_sided = NA_real_,
                test_status = "fit_failed"
            ))
        }

        fit_coef <- summary(fit)$coefficients
        slope_term <- "control_log2_mean_RPM"
        if (!slope_term %in% rownames(fit_coef)) {
            return(data.frame(
                n_points = n_points,
                n_unique_x = n_unique_x,
                slope = NA_real_,
                corr = NA_real_,
                corr_p_two_sided = NA_real_,
                test_status = "slope_term_missing"
            ))
        }

        data.frame(
            n_points = n_points,
            n_unique_x = n_unique_x,
            slope = unname(fit_coef[slope_term, "Estimate"]),
            corr = unname(corr_test$estimate),
            corr_p_two_sided = unname(corr_test$p.value),
            test_status = "ok"
        )
    }) %>%
    ungroup() %>%
    mutate(
        corr_p_two_sided_adj = ifelse(is.na(corr_p_two_sided), NA_real_, p.adjust(corr_p_two_sided, method = "BH")),
        has_correlation = test_status == "ok" & !is.na(corr_p_two_sided_adj) & corr_p_two_sided_adj < 0.0001,
        has_negative_correlation = has_correlation & corr < 0,
        has_positive_correlation = has_correlation & !is.na(slope) & !is.na(corr) & slope > 0.125 & corr > 0,
        correlation_label = case_when(
            has_positive_correlation ~ "Correlation: POSITIVE",
            has_negative_correlation ~ "Correlation: NEGATIVE",
            has_correlation ~ "Correlation: YES",
            test_status != "ok" ~ "Correlation: test unavailable",
            TRUE ~ "Correlation: NO"
        ),
        positive_slope_label = case_when(
            test_status != "ok" ~ "Positive slope > 0.125: test unavailable",
            is.na(slope) ~ "Positive slope > 0.125: NA",
            slope > 0.125 ~ "Positive slope > 0.125: YES",
            TRUE ~ "Positive slope > 0.125: NO"
        ),
        p_label = ifelse(
            is.na(corr_p_two_sided_adj),
            "adj p(two-sided)=NA",
            paste0("adj p(two-sided)=", formatC(corr_p_two_sided_adj, format = "e", digits = 2))
        ),
        corr_label = ifelse(
            is.na(corr),
            "r=NA",
            paste0("r=", formatC(corr, format = "f", digits = 3))
        ),
        slope_label = ifelse(
            is.na(slope),
            "slope=NA",
            paste0("slope=", formatC(slope, format = "f", digits = 4))
        ),
        facet_label = paste(correlation_label, positive_slope_label, p_label, corr_label, slope_label, sep = "\n")
    )

write.table(
    control_residual_correlation_stats %>% select(-facet_label),
    file = file.path(opt$plots_basedir, "control_barcode_residual_correlation_stats.tsv"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

# Keep legacy output path for downstream compatibility.
write.table(
    control_residual_correlation_stats %>% select(-facet_label),
    file = file.path(opt$plots_basedir, "control_barcode_residual_negative_correlation_stats.tsv"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

message("==== Plotting residuals vs control with correlation test")
pdf(file.path(opt$plots_basedir, "control_barcode_residuals_with_correlation_test.pdf"), width = 16, height = 16)
for (page_idx in seq_along(barcode_control_pages)) {
    page_comparisons <- barcode_control_pages[[page_idx]]
    page_df <- barcode_control_comparison_df %>%
        filter(comparison %in% page_comparisons) %>%
        mutate(comparison = factor(comparison, levels = page_comparisons))
    page_labels <- barcode_control_labels %>%
        filter(comparison %in% page_comparisons) %>%
        mutate(comparison = factor(comparison, levels = page_comparisons))
    page_stats <- control_residual_correlation_stats %>%
        filter(comparison %in% page_comparisons) %>%
        mutate(comparison = factor(comparison, levels = page_comparisons))

    if (nrow(page_df) == 0) {
        next
    }

    p <- ggplot(page_df, aes(x = control_log2_mean_RPM, y = residual)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_point(alpha = 0.25, size = 0.7) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#2c7fb8", linewidth = 0.7) +
        geom_text(
            data = page_stats,
            aes(x = -Inf, y = Inf, label = facet_label),
            inherit.aes = FALSE,
            hjust = -0.02,
            vjust = 1.05,
            size = 2.9,
            lineheight = 0.95
        ) +
        facet_wrap(~comparison, ncol = 5) +
        theme_pubr(border = T) +
        ggtitle(paste0("Residual barcode activity vs ", control_condition, " control (page ", page_idx, "/", length(barcode_control_pages), ")")) +
        xlab(paste0(control_condition, " control barcode activity (log2 mean RPM)")) +
        ylab("Residual activity (sample - control, log2 mean RPM)") +
        theme(
            text = element_text(size = 14),
            strip.text = element_text(size = 8, face = "bold"),
            axis.text = element_text(size = 7)
        )
    print(p)
}
invisible(dev.off())

# Keep legacy plot filename for downstream compatibility.
file.copy(
    from = file.path(opt$plots_basedir, "control_barcode_residuals_with_correlation_test.pdf"),
    to = file.path(opt$plots_basedir, "control_barcode_residuals_with_negative_correlation_test.pdf"),
    overwrite = TRUE
)

################################ PLOT REPLICATE CORRELATIONS #############################

message("==== Plotting replicate correlations")
# # Plot replicate correlation of the activity
pdf(file.path(opt$plots_basedir, "replicate_correlations.pdf"),
    width = 15, height = 15
)
for (this_condition in all_cdna_conditions) {
    message("----- ploting, ", this_condition)
    replicates_of_this_sample <- cdna_sample[str_detect(cdna_sample, paste0("^", this_condition, "_[0-9]+$"))]
    message("    Replicates: ", paste(replicates_of_this_sample, collapse = ", "))

    # Choose the y-label based on the number of replicates
    if(length(replicates_of_this_sample) > 1){
        y_lbl = "Log2(cDNA/pDNA)"
    } else {
        y_lbl = "Density"
    }

    this_pdna_sample <- pdna_sample[1]
    # message ("Processing ", this_cdna_sample)
    # Get the acitivty of the samples and take the average
    activity_df %>%
        filter(
            cDNA_sample %in% replicates_of_this_sample,
            !negative_control
        ) %>%
        dplyr::group_by(barcode, tf, promoter, cDNA_sample) %>%
        dplyr::reframe(
            activity_RPM = log2(mean(activity_RPM)),
            cDNA_sample = cDNA_sample
        ) %>%
        pivot_wider(names_from = cDNA_sample, values_from = activity_RPM) %>%
        select(-barcode, -tf, -promoter) %>%
        # print()
        # Plot correlations
        ggpairs(
            upper = list(continuous = upper_diag_plot),
            lower = list(continuous = lower_diag_plot),
            diag = list(continuous = diag_plot)
        ) +
        ggtitle(paste("Replicate correlation of ", this_condition)) +
        xlab("Log2(cDNA/pDNA)") +
        ylab(y_lbl) +
        theme(text = element_text(size = 18)) -> p
    print(p)
}
invisible(dev.off())

######################## PLOT BC CORRELATION WITH ACTIVITY ###############################

message("==== Plotting barcode correlations")
pdf(file.path(opt$plots_basedir, "barcode_correlations.pdf"), width = 17, height = 17)
# Start the plot
for (this_cdna_sample in cdna_sample) {
    #print(this_cdna_sample)
    message("----- ploting, ", this_cdna_sample)
    activity_df %>%
        filter(cDNA_sample == this_cdna_sample) %>%
        # Remove the random barcodes
        filter(!negative_control) %>%
        # Group by barcode and tf, calculate the mean count
        group_by(barcode, tf) %>%
        # Assign 1 to NAs
        summarise(mean_RPM = mean(activity_RPM), .groups = "drop") %>%
        # Group by tf, to make each barcode a column
        group_by(tf) %>%
        reframe(
            barcode_id = paste0("barcode", row_number()),
            log2_RPM = log2(mean_RPM)
        ) %>%
        pivot_wider(names_from = barcode_id, values_from = log2_RPM) %>%
        select(-tf) %>%
        # Plot correlations
        ggpairs(
            upper = list(continuous = upper_diag_plot),
            lower = list(continuous = lower_diag_plot),
            diag = list(continuous = diag_plot)
        ) +
        ggtitle(paste("Barcode correlation", this_cdna_sample)) +
        xlab("Log2(cDNA/pDNA)") +
        ylab("Log2(cDNA/pDNA)") +
        theme(text = element_text(size = 20)) -> p
    print(p)
}
invisible(dev.off())

##### VF240404: Compare the pDNA counts with the expected (the counts from our lab)
expected_pdna_counts <-
    read.table(opt$expected_pdna, header = T, sep = "\t") %>%
    # Normalize to RPM
    mutate(expected_RPM = pDNA / sum(pDNA) * 1e6) %>%
    select(expected_RPM, barcode)

pdf(file.path(opt$plots_basedir, "expected_vs_observed_pDNA_counts.pdf"),
    width = 8, height = 8
)
counts_df %>%
    filter(pDNA == T) %>%
    select(-pDNA, -tf, -promoter, -total_read_count, -log2_count, -RPM) %>%
    pivot_wider(names_from = sample, values_from = raw_count) %>%
    # Assign 1 to NAs
    mutate_all(~ ifelse(is.na(.), 1, .)) %>%
    mutate(observed_RPM = pDNA / sum(pDNA) * 1e6) %>%
    select(observed_RPM, barcode) %>%
    inner_join(expected_pdna_counts, by = "barcode") %>%
    ggplot(aes(x = observed_RPM, y = expected_RPM)) +
    geom_point(alpha = 0.2) +
    geom_abline(slope = 1, lty = "dashed", col = "red") +
    theme_pubr(border = T) +
    stat_cor(method = "pearson", size = 6) +
    ggtitle("Expected vs. observed pDNA counts per barcode") +
    xlab("Observed pDNA read count (RPM)") +
    ylab("Expected pDNA read count (RPM)") +
    theme(text = element_text(size = 18))
invisible(dev.off())