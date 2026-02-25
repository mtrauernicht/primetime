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
    make_option(c("--output"), type = "character", help = "Output directory for the plots"),
    make_option(c("--design"), type = "character", help = "Design DF with sample names"),
    make_option(c("--reference_condition"), type = "character", default = "Control", help = "Reference condition for the experiment (default: Control)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


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
df = 
  read.table(opt$cdna, 
             header=T) %>% 
  as_tibble() %>%
  select(-c(negative_control, promoter, barcode)) %>% 
  # Convert to RPM
  mutate_if(is.numeric, function(x) x / sum(x) * 1e6) %>%
  melt(id.vars='tf', variable.name = 'replicate', value.name='read_count_RPM') %>% as_tibble()

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
  ggplot(aes(x=replicate, y=read_count_RPM)) +
  geom_point(aes(color=treatment_col)) +
  scale_color_manual(values = c("Control condition" = "#000000", "NS" = "grey", "Downregulated" = "#6495ed", "Upregulated" = "#f37f80")) +
  facet_wrap(~factor(tf_lbl, levels(plot_df$tf_lbl)), scales='free_y') +
  labs(title='Barcode counts for all TFs (sorted by LogFC*p_adj)', x='Replicate', y='Read Count (RPM)') +
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