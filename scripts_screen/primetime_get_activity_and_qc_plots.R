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
    make_option(c("--cdna_output"), type = "character", help = "Path to save the cDNA counts for MPRAnalyze")
)

# Functions for the plots

upper_diag_plot <- function(data, mapping, color = I("black"), sizeRange = c(1, 3), ...) {
    boundaries <- seq(from = 0.8, by = 0.05, length.out = 4)

    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    r <- cor(x, y, "pairwise.complete.obs")
    rt <- format(r, digits = 3)
    tt <- as.character(rt)
    cex <- max(sizeRange)

    # helper function to calculate a useable size
    percent_of_range <- function(percent, range) {
        percent * diff(range) + min(range, na.rm = TRUE)
    }

    # plot correlation coefficient
    p <- ggally_text(
        label = tt, mapping = aes(), xP = 0.5, yP = 0.5,
        size = I(percent_of_range(cex * abs(r), sizeRange)) + 5, color = color, ...
    ) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
        )

    corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")[2:6]

    if (r <= boundaries[1]) {
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

##########################################################################################
## Start of the script ###################################################################
##########################################################################################

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
list_of_annotated_files = strsplit(opt$list_of_annotated_files, " ")[[1]]
unique_files = unique(list_of_annotated_files)

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
    # ADD PSEUDOCOUNT OF 1!!!!
    this_replicate_df$count <- this_replicate_df$count + 1
    is_pDNA = this_sample == "pDNA"
    total_read_count = sum(this_replicate_df$count)
    this_df = data.frame(
        sample = this_sample,
        replicate = this_replicate,
        pDNA = is_pDNA,
        tf = this_replicate_df$tf,
        negative_control = ifelse(this_replicate_df$neg_ctrls == "Yes", T, F),
        promoter = this_replicate_df$promoter,
        barcode = this_replicate_df$barcode,
        raw_count = this_replicate_df$count,
        log2_count = log2(this_replicate_df$count),
        total_read_count=total_read_count,
        RPM = this_replicate_df$count / total_read_count * 1e6
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
    rowwise() %>%
    mutate(pDNA = mean(c_across(-barcode))) %>%
    ungroup() %>%
    select(barcode, pDNA)

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



message("==== Saving bleedthrough data")
# Save bleedthrough slope data per condition for later use
# Remove replicate suffix (_1, _2, etc.) from sample names to match comparison format
bleedthrough_summary <- bleed_through_slope_df %>%
    select(sample, slope) %>%
    distinct() %>%
    mutate(condition = sub("_\\d+$", "", sample)) %>%
    select(condition, bleedthrough = slope)

write.table(
    bleedthrough_summary,
    file.path(opt$activity_basedir, "bleedthrough_per_condition.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
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