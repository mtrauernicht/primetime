#!/usr/bin/env Rscript
# QC for the banana-shaped bias correction: one page per corrected replicate showing
# M (replicate - anchor activity) vs A (anchor activity), before and after correction.
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
})

option_list <- list(
    make_option("--input-file", type = "character", dest = "input_file", help = "Text file listing the per-comparison banana_qc TSVs (one path per line)"),
    make_option("--output", type = "character", help = "Output PDF")
)
opt <- parse_args(OptionParser(option_list = option_list))

paths <- readLines(opt$input_file, warn = FALSE)
paths <- paths[nzchar(paths) & file.exists(paths) & file.info(paths)$size > 0]

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

data <- if (length(paths) == 0) {
    data.frame()
} else {
    bind_rows(lapply(seq_along(paths), function(i) read.delim(paths[i], check.names = FALSE) %>% mutate(source_idx = i)))
}

pdf(opt$output, width = 12, height = 6)
if (nrow(data) == 0) {
    plot.new()
    text(0.5, 0.5, "No banana correction data available")
} else {
    # A replicate is corrected identically in every comparison it appears in; keep the first copy.
    data <- data %>%
        group_by(replicate) %>%
        filter(source_idx == min(source_idx)) %>%
        ungroup() %>%
        mutate(stage = factor(stage, levels = c("Before correction", "After correction"))) %>%
        arrange(condition, replicate)

    for (rep_name in unique(data$replicate)) {
        rep_df <- data %>% filter(replicate == rep_name)
        print(
            ggplot(rep_df, aes(x = A, y = M)) +
                geom_point(alpha = 0.15, size = 0.5) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#2c7fb8", span = 0.75) +
                facet_wrap(~stage) +
                ggpubr::theme_pubr(border = TRUE) +
                labs(
                    title = paste0("Banana correction: ", unique(rep_df$condition), " | ", rep_name),
                    x = "Anchor mean log2 activity",
                    y = "M (replicate - anchor, log2 activity)"
                ) +
                theme(text = element_text(size = 14))
        )
    }
}
invisible(dev.off())
