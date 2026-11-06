#!/usr/bin/env Rscript
# QC for count-level normalization: log2 shift applied to each condition and promoter,
# relative to plain library-size normalization.
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

option_list <- list(
    make_option("--cdna", type = "character", help = "cDNA_counts.txt"),
    make_option("--pdna", type = "character", help = "pDNA_counts.txt"),
    make_option("--read_count_summary", type = "character", help = "read_count_per_sample.tsv (replicate -> condition)"),
    make_option("--output", type = "character", help = "Output PDF")
)
opt <- parse_args(OptionParser(option_list = option_list))

cdna <- read.delim(opt$cdna, check.names = FALSE, stringsAsFactors = FALSE)
pdna <- read.delim(opt$pdna, check.names = FALSE, stringsAsFactors = FALSE)
condition_map <- read.delim(opt$read_count_summary, check.names = FALSE, stringsAsFactors = FALSE) %>%
    select(replicate, condition) %>%
    distinct()

meta_cols <- c("tf", "negative_control", "promoter", "barcode")
replicates <- setdiff(colnames(cdna), meta_cols)
pdna_col <- setdiff(colnames(pdna), c("negative_control", "barcode"))[1]
data <- inner_join(cdna, pdna[, c("barcode", pdna_col)] %>% rename(pDNA = all_of(pdna_col)), by = "barcode")
is_random <- toupper(trimws(as.character(data$negative_control))) %in% c("TRUE", "T", "1", "YES", "Y")

pdna_total <- sum(data$pDNA, na.rm = TRUE)
shift_df <- bind_rows(lapply(replicates, function(rep) {
    cdna_total <- sum(data[[rep]], na.rm = TRUE)
    bind_rows(lapply(unique(data$promoter[is_random]), function(promoter) {
        rows <- is_random & data$promoter == promoter
        cdna_random <- sum(data[[rep]][rows], na.rm = TRUE)
        pdna_random <- sum(data$pDNA[rows], na.rm = TRUE)
        data.frame(
            replicate = rep,
            promoter = promoter,
            # Count-level scaling multiplies cDNA/pDNA by 1/RANDOM sum (the target constant cancels);
            # library-size normalization would use 1/total. Their ratio is the extra correction.
            log2_shift = log2((pdna_random / pdna_total) / (cdna_random / cdna_total))
        )
    }))
})) %>%
    left_join(condition_map, by = "replicate") %>%
    mutate(condition = ifelse(is.na(condition), replicate, condition)) %>%
    filter(is.finite(log2_shift))

condition_order <- shift_df %>%
    group_by(condition) %>%
    summarise(m = median(log2_shift), .groups = "drop") %>%
    arrange(m) %>%
    pull(condition)
shift_df$condition <- factor(shift_df$condition, levels = condition_order)

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
n_cond <- length(condition_order)
pdf(opt$output, width = max(8, 0.35 * n_cond + 3), height = 4 + 2.2 * n_distinct(shift_df$promoter))
print(
    ggplot(shift_df, aes(x = condition, y = log2_shift)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(aes(color = promoter), position = position_jitter(width = 0.15, height = 0), alpha = 0.8) +
        stat_summary(fun = median, geom = "crossbar", width = 0.6, linewidth = 0.3, color = "black") +
        facet_grid(promoter ~ .) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
        labs(
            title = "Count-level normalization: log2 shift per condition and promoter",
            subtitle = "Relative to library-size normalization. > 0: activities corrected up (RANDOM reporters take a smaller share of cDNA than of pDNA); < 0: corrected down.\nOne point per replicate, black bar = median.",
            x = "Condition (ordered by median shift)",
            y = "log2 shift"
        )
)
invisible(dev.off())
