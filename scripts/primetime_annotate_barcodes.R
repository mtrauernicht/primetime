# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script overlaps the barcode counts with the barcode annotation.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================

suppressPackageStartupMessages({
    library("optparse")
    library(dplyr)
})
options(warn = -1)

option_list <- list(
    make_option(c("--bc_counts"),
        type = "character", default = NULL,
        help = "path to the barcode counts, after clustering",
        metavar = "path.cluster.txt"
    ),
    make_option(c("--bc_annotation"),
        type = "character", default = NULL,
        help = "Path to the barcode annotation file",
        metavar = "path.csv"
    ),
    make_option(c("--sample"),
        type = "character", default = NULL,
        help = "Sample name, for metadata",
        metavar = "sample_name"
    ),
    make_option(c("--replicate_number"),
        type = "numeric", default = NULL,
        help = "Replicate number, for metadata", metavar = "replicate_no"
    ),
    make_option(c("--treatment"),
        type = "numeric", default = NULL,
        help = "Treatment name, for metadata", metavar = "Treatment"
    ),
    make_option(c("--is_pdna"),
        type = "logical", default = NULL,
        help = "Is this sample pDNA?", metavar = "T/F"
    ),
    make_option(c("--output"),
        type = "character", default = NULL,
        help = "Path to the output file", metavar = "path.csv"
    ),
    make_option(c("--amount_spikes"),
        type = "numeric", default = 100000,
        help = "Number of spikes to add", metavar = "amount_spikes"
    ),
    make_option(c('--stats'),
        type = "character", default = NULL,
        help = "Path to the stats file", metavar = "path.txt"
    )
)

# Parse the arguments
opt_parser <- OptionParser(
    usage = "\tRscript %prog [options]",
    option_list = option_list,
    description = "Annotate the barcode counts with the barcode annotation"
)
opt <- parse_args(opt_parser)

# Read input
readLines(opt$bc_counts) %>%
    strsplit("\t") %>%
    lapply(function(x) if(length(x) == 3) x else NULL) %>%
    Filter(Negate(is.null), .) %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("barcode", "count", "unclustered_bc")) %>%
    mutate(count = as.numeric(count)) %>%
    as_tibble() -> bc_counts

# Read annotation
read.table(opt$bc_annotation, sep = ",", header = TRUE) %>%
    as_tibble() -> bc_annotation

# Merge
inner_join(bc_counts, bc_annotation, by = join_by("barcode")) %>%
    select(-unclustered_bc) %>%
    # Remove the amount of spikes inserted for the clustering
    mutate(count = count - opt$amount_spikes) -> merged_df

# Stats file starts with "Total reads: 180957, of those: ". use this to get the read number

read_stats <- readLines(opt$stats)
total_reads <- as.numeric(gsub("Total reads: ([0-9]+), of those: .*", "\\1", read_stats[1]))

# Now, get the number of reads after merged (sum of counts)
merged_reads <- sum(merged_df$count)
pct_match <- merged_reads / total_reads * 100
lost_reads <- total_reads - merged_reads
pct_unmatch <- 100 - pct_match

# Add this to the stats file
rows_to_write = c(
    "------------------ After merging Barcode annotation ------------------",
    sprintf("Merged reads: %d (%.2f%%)", merged_reads, pct_match),
    sprintf("Lost reads: %d (%.2f%%)", lost_reads, pct_unmatch)
)
cat(rows_to_write, file = opt$stats, append = TRUE, sep = "\n")

# Write output table
write.table(merged_df, file = opt$output, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
