# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script generates a plot of barcode counts for the top responding TFs
# based on the highest LogFC * p_adjusted values.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggthemes)
  library(ggnewscale)
    library(optparse)
    library(forcats)
    library(ggplot2)
    library(reshape2)
})
# Read the arguments:
option_list <- list(
    make_option(c("--result"), type = "character", help = "Results file from PrimeTime"),
    make_option(c("--cdna"), type = "character", help = "cDNA counts file"),
  make_option(c("--pdna"), type = "character", default = NA_character_, help = "pDNA counts file (optional; auto-detected if omitted)"),
    make_option(c("--output"), type = "character", help = "Output directory for the plots"),
    make_option(c("--design"), type = "character", help = "Design DF with sample names"),
    make_option(c("--reference_condition"), type = "character", default = "Control", help = "Reference condition for the experiment (default: Control)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

find_existing_path <- function(paths) {
  for (candidate in paths) {
    if (!is.na(candidate) && file.exists(candidate)) return(candidate)
  }
  return(NA_character_)
}

pdna_path <- if (!is.na(opt$pdna) && nzchar(opt$pdna)) {
  opt$pdna
} else {
  find_existing_path(c(
    file.path(dirname(dirname(dirname(opt$output))), "tmp_primetime", "activity", "pDNA_counts.txt"),
    file.path(dirname(dirname(opt$output)), "tmp_primetime", "activity", "pDNA_counts.txt"),
    file.path(dirname(opt$output), "tmp_primetime", "activity", "pDNA_counts.txt")
  ))
}

if (is.na(pdna_path)) {
  stop("Could not find pDNA counts file. Pass --pdna explicitly or ensure tmp_primetime/activity/pDNA_counts.txt exists next to the output directory.")
}


read.table(opt$result, header=T) %>% 
  as_tibble() %>%
  # remove rows that tf contains 'RANDOM'
  filter(!grepl("RANDOM", tf)) %>%
  mutate(FC_x_Padj = -abs(logFC) * p_adjusted) %>%
  arrange(desc(FC_x_Padj)) -> top10_tfs

# Read design DF 
read.table(opt$design, header=T) %>% 
  as_tibble() %>% filter(pDNA != "True") -> design_df


# Read counts df
cdna_df = 
  read.table(opt$cdna, 
             header=T) %>% 
  as_tibble() %>%
  select(-any_of(c("negative_control", "promoter")))

pdna_df = 
  read.table(pdna_path, 
             header=T) %>% 
  as_tibble() %>%
  select(any_of(c("barcode")), where(is.numeric))

if (!("barcode" %in% colnames(cdna_df))) {
  stop("cDNA counts file must contain a 'barcode' column to compute log2(cDNA/pDNA).")
}

if (!("barcode" %in% colnames(pdna_df))) {
  stop("pDNA counts file must contain a 'barcode' column to compute log2(cDNA/pDNA).")
}

pdna_num_cols <- setdiff(colnames(pdna_df), "barcode")
if (length(pdna_num_cols) == 0) {
  stop("pDNA counts file has no numeric pDNA columns.")
}

# Match QC activity scaling: normalize pDNA counts to RPM per pDNA replicate, then average.
pdna_mean <- pdna_df %>%
  mutate(across(all_of(pdna_num_cols), ~ as.numeric(.) / sum(as.numeric(.), na.rm = TRUE) * 1e6)) %>%
  mutate(pDNA_mean_RPM = rowMeans(across(all_of(pdna_num_cols)), na.rm = TRUE)) %>%
  select(barcode, pDNA_mean_RPM) %>%
  group_by(barcode) %>%
  summarize(pDNA_mean_RPM = mean(pDNA_mean_RPM, na.rm = TRUE), .groups = 'drop')

df = 
  melt(cdna_df, id.vars=c('tf', 'barcode'), variable.name = 'replicate', value.name='cDNA_count') %>% as_tibble() %>%
  mutate(cDNA_count = as.numeric(cDNA_count)) %>%
  group_by(replicate) %>%
  mutate(cDNA_RPM = cDNA_count / sum(cDNA_count, na.rm = TRUE) * 1e6) %>%
  ungroup() %>%
  left_join(pdna_mean, by = 'barcode') %>%
  mutate(read_count_log2 = log2(cDNA_RPM / pDNA_mean_RPM))

inner_join(df, design_df, by = 'replicate') %>%
  inner_join(top10_tfs, by = 'tf') %>%
  mutate(treatment_col=ifelse(treatment==opt$reference_condition, "Control condition", sig)) %>%
  # Make sure sig has the desired factor level order
  mutate(sig = factor(sig, levels = c('Upregulated', 'Downregulated', 'NS'))) %>%
  arrange(sig) %>%
  # Create tf_lbl and set factor levels based on sig ordering
  mutate(tf_lbl = paste0(tf, ' (', sig, ')\nlogFC: ', round(logFC, 2))) %>%
  mutate(tf_lbl = fct_inorder(tf_lbl)) -> plot_df
  

pdf(opt$output, width=20, height=20)
plot_df %>%
# plot
  ggplot(aes(x=replicate, y=read_count_log2)) +
  geom_point(aes(color=treatment_col)) +
  scale_color_manual(values = c("Control condition" = "#000000", "NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
  facet_wrap(~factor(tf_lbl, levels(plot_df$tf_lbl)), scales='free_y') +
  labs(title='Barcode counts for all TFs (sorted by LogFC*p_adj)', x='Replicate', y='Log2(cDNA/pDNA)') +
  ggthemes::theme_few() +
  new_scale_color() +
  geom_violin(aes(fill=treatment_col, color=treatment_col), alpha=0.4, width=1) +
    scale_fill_manual(values = c("Control condition" = "#00000033", "NS" = "#e4e4e455", "Downregulated" = "#6495ed33", "Upregulated" = "#f37f8033")) +
    scale_color_manual(values = c("Control condition" = "#00000022", "NS" = "#e4e4e455", "Downregulated" = "#6495ed33", "Upregulated" = "#f37f8033")) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    strip.text = element_text(size=8),
    aspect.ratio = 1,
    legend.position = "none"
  ) 

invisible(dev.off())