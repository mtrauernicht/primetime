# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script normalizes the raw cDNA counts from a TF reporter by dividing
# the median count of negative control constructs with the same minimal promoter
# as the reporter constructs.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(reshape2)
    library(tidyr)
})
options(
    dplyr.width = Inf,
    error = rlang::entrace,
    show.error.locations = TRUE
)

option_list <- list(
    make_option(c("--input_cDNA"),
        help = "Path to the input cDNA file (with raw counts)",
        type = "character"
    ),
    make_option(c("--output_cDNA"),
        help = "Path to the output normalized cDNA file",
        type = "character"
    )
)

##########################################################################################
## Start of the script ###################################################################
##########################################################################################

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cdna <- read.table(opt$input_cDNA,
    header = TRUE,
    sep = "\t"
) %>% as_tibble()

# Calculate median count for every minimal promoter in the negative controls
cdna %>%
    filter(negative_control) %>%
    select(-barcode, -tf) %>%
    group_by(promoter) %>%
    summarise_all(median) -> ctrl_medians

# Now divide the cDNA counts by the median control counts
cdna %>%
    left_join(ctrl_medians, by = "promoter", suffix = c("", "_ctrl_median")) %>%
    mutate(
        across(
            matches("_\\d+$"),
            ~ . / get(paste0(cur_column(), "_ctrl_median"))
        )
    ) %>%
    select(-ends_with("_ctrl_median")) -> cdna_normalized

# Now multiply by the median read count per sequencing sample and round numbers to reconstitute proper read counts
cdna %>%
    select(-promoter, -barcode, -tf, -negative_control) %>%
    summarise_all(median) %>%
    pivot_longer(
        cols = everything(),
        names_to = "sample",
        values_to = "median_count"
    ) -> sample_medians

## Export sample medians to file
write.table(sample_medians,
    file = sub("\\.txt$", "_sample_medians.txt", opt$output_cDNA),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cdna_normalized %>%
    mutate(across(
        matches("_\\d+$"),
        ~ round(. * sample_medians$median_count[match(cur_column(), sample_medians$sample)])
    )) -> cdna_normalized

# Plot beeswarm to show the negative values. Later we will set a floor at 1

# cdna_normalized %>%
#     filter(!negative_control) %>%
#     select(-promoter, -barcode, -tf, -negative_control) %>%
#     melt() %>%
#     ggplot(aes(
#         x = variable,
#         y = value,
#         group = variable,
#         color = ifelse(value >= 0, "b", "a")
#     )) +
#     geom_beeswarm(cex = 0.12, alpha = 0.5) +
#     geom_hline(yintercept = 0) +
#     labs(
#         title = "cDNA counts normalized by subtracting control median",
#         x = "",
#         y = "Normalized cDNA counts"
#     ) +
#     scale_color_manual(values = c("a" = "#e41a1c", "b" = "gray40")) +
#     ggpubr::theme_pubr() +
#     theme(
#         legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1)
#     )

# Now, we set a minimum value of 1 for all negative values
# cdna_normalized %>%
#     mutate(
#         across(
#             matches("_\\d+$"),
#             ~ ifelse(. < 1, 1, .)
#         )
#     ) -> cdna_normalized_floored


# cdna_normalized_floored %>%
#     filter(!negative_control) %>%
#     select(-promoter, -barcode, -tf, -negative_control) %>%
#     melt() %>%
#     ggplot(aes(
#         x = variable,
#         y = value,
#         group = variable,
#         color = ifelse(value >= 0, "b", "a")
#     )) +
#     geom_beeswarm(cex = 0.05, alpha = 0.5) +
#     labs(
#         title = "cDNA counts normalized by subtracting control median (floored)",
#         x = "",
#         y = "Normalized cDNA counts"
#     ) +
#     scale_color_manual(values = c("a" = "#e41a1c", "b" = "gray40")) +
#     ggpubr::theme_pubr() +
#     theme(
#         legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1)
#     )

# Write the normalized cDNA counts to file
write.table(cdna_normalized,
    file = opt$output_cDNA,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)