suppressPackageStartupMessages({
    library(tidyverse)
    library(ggthemes)
    library(ggnewscale)
    library(optparse)
    library(forcats)
    library(ggplot2)
    library(reshape2)
    library(ggtext)
})

# Read the arguments:
option_list <- list(
    make_option(c("--results"), type = "character", help = "Results file from PrimeTime"),
    make_option(c("--output"), type = "character", help = "Output file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("DEBUG: opt$results =", opt$results, "\n")

# Construct the df
all_comparisons_df = data.frame()
all_comparisons_randoms = data.frame()
for(this_comparison in opt$results %>% strsplit(split=",") %>% unlist) {
    cat("DEBUG: Processing file:", this_comparison, "\n")
    tmp = basename(this_comparison) %>% str_remove(".txt") %>% strsplit(split="_vs_") %>% unlist
    comparison_name = sprintf("**%s** vs. %s", tmp[1], tmp[2])
    cat("DEBUG: Comparison name:", comparison_name, "\n")
    read_table(this_comparison) %>%
        select(tf, sig, logFC) %>%
        mutate(comparison = comparison_name) -> tmp_df
    cat("DEBUG: tmp_df nrow =", nrow(tmp_df), "\n")
    cat("DEBUG: sig values unique:", paste(unique(tmp_df$sig), collapse=", "), "\n")
    randoms = tmp_df %>% filter(grepl("RANDOM", tf))
    not_randoms = tmp_df %>% filter(!grepl("RANDOM", tf))
    cat("DEBUG: randoms nrow =", nrow(randoms), ", not_randoms nrow =", nrow(not_randoms), "\n")
    all_comparisons_df = rbind(all_comparisons_df, not_randoms)
    all_comparisons_randoms = rbind(all_comparisons_randoms, randoms)
}

cat("DEBUG: all_comparisons_df total nrow =", nrow(all_comparisons_df), "\n")
cat("DEBUG: all_comparisons_df head:\n")
print(head(all_comparisons_df))

# Cluster the comparisons based on the logFC values
filtered_for_clustering <- all_comparisons_df %>%
    filter(sig != 'NS') %>%
    select(-sig)
cat("DEBUG: Rows after filtering sig != 'NS':", nrow(filtered_for_clustering), "\n")
cat("DEBUG: Unique comparisons:", paste(unique(filtered_for_clustering$comparison), collapse=", "), "\n")

filtered_for_clustering %>%
    pivot_wider(names_from = comparison, values_from = logFC, values_fill = 0) %>%
    column_to_rownames(var = "tf") %>%
    t() %>%
    scale() %>%
    dist() %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram() %>%
    labels() -> comparison_order
cat("DEBUG: comparison_order length:", length(comparison_order), "\n")

# Cluster the TFs based on the logFC values
all_comparisons_df %>%
    filter(sig != 'NS') %>%
    select(-sig) %>%
    pivot_wider(names_from = tf, values_from = logFC, values_fill = 0) %>%
    column_to_rownames(var = "comparison") %>%
    t() %>%
    dist() %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram() %>%
    labels() -> tf_order


# Reorder the comparisons based on the clustering
all_comparisons_df %>%
    mutate(comparison = factor(comparison, levels = comparison_order)) %>%
    mutate(tf = factor(tf, levels = tf_order)) -> plot_df

# Plot the heatmap
up = plot_df %>% filter(sig=='Upregulated')
down = plot_df %>% filter(sig=='Downregulated')

cat("DEBUG: up nrow =", nrow(up), "\n")
cat("DEBUG: down nrow =", nrow(down), "\n")
cat("DEBUG: plot_df before NA filter nrow =", nrow(plot_df), "\n")

# Filter NA plot_Df
plot_df %>%
    filter(!is.na(comparison)) %>%
    filter(!is.na(tf)) -> plot_df

cat("DEBUG: plot_df after NA filter nrow =", nrow(plot_df), "\n")
cat("DEBUG: plot_df after NA filter head:\n")
print(head(plot_df))

pdf(opt$output, width=16, height=8)
ggplot(plot_df, aes(x = tf, y = comparison)) +
    geom_tile(fill = "white", color = NA) +
    geom_tile(data = up, aes(fill = logFC), size = 0.1, color = NA) +
    scale_fill_distiller(palette = "Reds", name = "LogFC\n(Upregulated)", direction = 1) +
    ggnewscale::new_scale_fill() +
    theme_few() +
    geom_tile(data = down, aes(fill = logFC), size = 0.1, color = NA) +
    scale_fill_distiller(palette = "Blues", name = "LogFC\n(Downregulated)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
        title = "TF Activity Changes Across Comparisons",
        x = "",
        y = ""
    ) +
    theme(
        axis.text.x = element_text(size = 6),
        axis.text.y = element_markdown(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    coord_fixed(ratio = 1)
invisible(dev.off())