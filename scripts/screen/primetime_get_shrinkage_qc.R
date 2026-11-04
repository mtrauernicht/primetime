#!/usr/bin/env Rscript
# QC of the post-BCalm logFC shrinkage (SHRINKAGE: True). Reads all output_data tables.
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
    library(tidyr)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option(c("--results-file"), type = "character", help = "File with newline-separated results files"),
    make_option(c("--output"), type = "character", help = "Output PDF"),
    make_option(c("--logfc-threshold"), type = "numeric", default = 0.263, help = "logFC cutoff used for the hit calls")
)))

files <- readLines(opt$`results-file`)
files <- files[nzchar(files) & file.exists(files)]
df <- bind_rows(lapply(files, function(f) {
    d <- tryCatch(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !all(c("tf", "logFC", "logFC_shrunk", "shrink_weight", "logFC_se") %in% colnames(d))) {
        return(NULL)
    }
    d <- d[!grepl("RANDOM", d$tf, ignore.case = TRUE), ]
    d$comparison <- sub("\\.txt$", "", basename(f))
    d$se_scale <- if ("se_scale" %in% colnames(d)) d$se_scale else NA_real_
    d[, c("comparison", "tf", "logFC", "logFC_shrunk", "shrink_weight", "logFC_se", "se_scale", "sig")]
}))

pdf(opt$output, width = 11, height = 8)
if (is.null(df) || nrow(df) == 0) {
    plot.new()
    text(0.5, 0.5, "No shrunk results found (SHRINKAGE is probably off)")
    invisible(dev.off())
    quit(status = 0)
}

df <- df %>% mutate(se_cal = logFC_se * ifelse(is.finite(se_scale), se_scale, 1))
set.seed(1)
sub_df <- df %>% filter(is.finite(logFC), is.finite(logFC_shrunk))
sub_df <- sub_df[sample(nrow(sub_df), min(nrow(sub_df), 30000)), ]

# 1. unshrunk vs shrunk
print(ggplot(sub_df, aes(x = logFC, y = logFC_shrunk, color = shrink_weight)) +
    geom_abline(color = "grey50") +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_viridis_c(limits = c(0, 1)) +
    theme_bw() +
    labs(title = "Shrunk vs unshrunk logFC (all contrasts, subsample)", x = "logFC (unshrunk)", y = "logFC (shrunk)", color = "weight"))

# 2. weight vs calibrated SE
print(ggplot(sub_df, aes(x = se_cal, y = shrink_weight)) +
    geom_point(size = 0.6, alpha = 0.4) +
    theme_bw() +
    labs(title = "Shrinkage weight vs standard error", subtitle = "weight 1 = unchanged; large SE and small |logFC| are shrunk hardest", x = "Calibrated SE of logFC", y = "Shrinkage weight"))

# 3. weight vs |logFC|
print(ggplot(sub_df, aes(x = abs(logFC), y = shrink_weight, color = sig)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_color_manual(values = c("NS" = "grey60", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
    theme_bw() +
    labs(title = "Shrinkage weight vs |logFC|", x = "|logFC| (unshrunk)", y = "Shrinkage weight", color = "Call"))

# 4. per-contrast summary
comp_df <- df %>%
    group_by(comparison) %>%
    summarise(
        mean_weight = mean(shrink_weight, na.rm = TRUE),
        frac_w_lt_half = mean(shrink_weight < 0.5, na.rm = TRUE),
        se_scale = first(se_scale),
        sd_unshrunk = sd(logFC, na.rm = TRUE),
        sd_shrunk = sd(logFC_shrunk, na.rm = TRUE),
        n_hits = sum(sig != "NS"),
        n_down_hits = sum(sig == "Downregulated"),
        n_down_hits_lost = sum(sig == "Downregulated" & abs(logFC_shrunk) < opt$`logfc-threshold`),
        n_up_hits = sum(sig == "Upregulated"),
        n_up_hits_lost = sum(sig == "Upregulated" & abs(logFC_shrunk) < opt$`logfc-threshold`),
        .groups = "drop"
    )

print(ggplot(comp_df, aes(x = mean_weight)) +
    geom_histogram(bins = 30, fill = "#4e79a7", color = "white") +
    theme_bw() +
    labs(title = "Mean shrinkage weight per contrast", subtitle = paste0("Overall fraction of TF values with w < 0.5: ", round(mean(df$shrink_weight < 0.5, na.rm = TRUE), 3)), x = "Mean weight", y = "Contrasts"))

print(ggplot(comp_df, aes(x = sd_unshrunk, y = sd_shrunk, color = n_hits)) +
    geom_abline(color = "grey50") +
    geom_point(size = 1.5) +
    scale_color_viridis_c(trans = "sqrt") +
    theme_bw() +
    labs(title = "Spread of TF logFC per contrast before and after shrinkage", x = "SD of logFC (unshrunk)", y = "SD of logFC (shrunk)", color = "# hits"))

if (any(is.finite(comp_df$se_scale))) {
    print(ggplot(comp_df, aes(x = se_scale)) +
        geom_histogram(bins = 30, fill = "#59a14f", color = "white") +
        theme_bw() +
        labs(title = "SE calibration factor per contrast", subtitle = "Spread of RANDOM-reporter z-scores (capped at 1); < 1 means the BCalm SE is conservative", x = "SE scale", y = "Contrasts"))
}

# 5. hits that would not pass the logFC cutoff after shrinkage (unshrunk calls are the reported ones)
hit_df <- comp_df %>%
    transmute(comparison, `Down kept` = n_down_hits - n_down_hits_lost, `Down lost` = n_down_hits_lost, `Up kept` = n_up_hits - n_up_hits_lost, `Up lost` = n_up_hits_lost) %>%
    pivot_longer(-comparison, names_to = "class", values_to = "n") %>%
    group_by(class) %>%
    summarise(n = sum(n), .groups = "drop")
print(ggplot(hit_df, aes(x = class, y = n, fill = grepl("lost", class))) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "#e15759", "FALSE" = "grey60"), guide = "none") +
    theme_bw() +
    labs(title = "Significant hits with |shrunk logFC| below the cutoff", subtitle = "Hit-list sensitivity: 'lost' hits would not pass the logFC threshold after shrinkage", x = "", y = "Hits (all contrasts)"))

invisible(dev.off())
write.table(comp_df, sub("\\.pdf$", ".tsv", opt$output), sep = "\t", quote = FALSE, row.names = FALSE)
