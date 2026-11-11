test_dir <- tempfile("heatmap-transpose-")
dir.create(test_dir)
on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

script <- normalizePath("scripts/primetime_get_heatmap_of_conditions.R")
wrapper <- file.path(test_dir, "wrapper.R")
writeLines(c(
    "suppressPackageStartupMessages(library(pheatmap))",
    "heatmap_calls <- 0L",
    "trace('pheatmap', where = as.environment('package:pheatmap'), print = FALSE, tracer = quote({",
    "    heatmap_calls <<- heatmap_calls + 1L",
    "    saveRDS(list(mat = mat, annotation_col = annotation_col, annotation_row = annotation_row,",
    "                 cluster_rows = cluster_rows, cluster_cols = cluster_cols, angle_col = angle_col,",
    "                 cellwidth = cellwidth, fontsize_col = fontsize_col),",
    "            file.path(Sys.getenv('HEATMAP_CAPTURE_DIR'), paste0(heatmap_calls, '.rds')))",
    "}))",
    "source(Sys.getenv('HEATMAP_SCRIPT'))"
), wrapper)

samples <- c("Drug'A", "DrugB", "DrugC")
shrunk <- rbind(c(2, -1, 0.2), c(-1.5, 0.5, 1.3), c(0.1, 0.15, 0.25))
pvalues <- rbind(c(0.01, 0.2, 0.3), c(0.001, 0.2, 0.02), c(0.9, 0.9, 0.9))
significant <- rbind(c("UP", "NS", "NS"), c("DOWN", "NS", "UP"), c("NS", "NS", "NS"))
result_files <- vapply(seq_along(samples), function(i) {
    path <- file.path(test_dir, paste0(samples[i], "_vs_Control.txt"))
    result <- data.frame(tf = c("TF1", "TF2", "TF3"), logFC = 2 * shrunk[i, ],
                         logFC_shrunk = shrunk[i, ], p_adjusted = pvalues[i, ], sig = significant[i, ])
    result[[samples[i]]] <- shrunk[i, ]
    write.table(result, path, sep = "\t", quote = FALSE, row.names = FALSE)
    path
}, character(1))

metadata <- file.path(test_dir, "metadata.csv")
write.csv(data.frame(sample = samples, target = c("TargetA", "TargetB", "TargetC"),
                     class = c("classA", "classB", "classC")), metadata, row.names = FALSE)
read_counts <- file.path(test_dir, "reads.tsv")
write.table(data.frame(condition = samples, total_read_count = c(10000, 20000, 30000)),
            read_counts, sep = "\t", quote = FALSE, row.names = FALSE)
viability <- file.path(test_dir, "viability.tsv")
write.table(data.frame(sample = samples, mean_viability_pct = c(50, 75, 100)),
            viability, sep = "\t", quote = FALSE, row.names = FALSE)

capture <- file.path(test_dir, "capture")
dir.create(capture)
output <- file.path(test_dir, "heatmap.pdf")
args <- c(wrapper, "--results", paste(result_files, collapse = ","), "--sample_info", metadata,
          "--read-count-summary", read_counts, "--sample-viability-summary", viability, "--output", output)
status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args),
                  env = c(paste0("HEATMAP_CAPTURE_DIR=", shQuote(capture)),
                          paste0("HEATMAP_SCRIPT=", shQuote(script))))
stopifnot(status == 0, file.info(output)$size > 0)

plots <- lapply(file.path(capture, paste0(seq_len(7), ".rds")), readRDS)
labels <- paste(samples, "->", c("TargetA", "TargetB", "TargetC"))
for (i in seq_along(plots)) {
    plot <- plots[[i]]
    expected_labels <- if (i %in% c(4, 5, 7)) labels[1:2] else labels
    stopifnot(
        identical(rownames(plot$mat), c("TF1", "TF2", "TF3")),
        identical(colnames(plot$mat), expected_labels),
        identical(plot$angle_col, "270"),
        identical(plot$cellwidth, 6),
        identical(plot$fontsize_col, 6),
        all(expected_labels %in% rownames(plot$annotation_col)),
        all(is.na(plot$annotation_row)),
        identical(plot$annotation_col[expected_labels, "sample_class"], c("classA", "classB", "classC")[seq_along(expected_labels)]),
        identical(plot$annotation_col[expected_labels, "read_count"], c(10000, 20000, 30000)[seq_along(expected_labels)]),
        identical(plot$annotation_col[expected_labels, "mean_viability"], c(50, 75, 100)[seq_along(expected_labels)]),
        identical(plot$cluster_rows$labels, rownames(plot$mat)),
        identical(plot$cluster_cols$labels, colnames(plot$mat))
    )
}
stopifnot(
    isTRUE(all.equal(unname(plots[[1]]$mat), t(2 * shrunk))),
    isTRUE(all.equal(unname(plots[[2]]$mat), t(shrunk))),
    isTRUE(all.equal(unname(plots[[3]]$mat), t(sweep(shrunk, 2, apply(shrunk, 2, median))))),
    is.na(plots[[5]]$mat["TF2", labels[1]]),
    plots[[5]]$mat["TF1", labels[1]] == shrunk[1, 1],
    isTRUE(all.equal(unname(plots[[6]]$mat), t(shrunk * -log10(pvalues))))
)

single_result <- file.path(test_dir, "Single_vs_Control.txt")
write.table(data.frame(tf = "TF1", logFC = 0.1, p_adjusted = 0.9, sig = "NS"),
            single_result, sep = "\t", quote = FALSE, row.names = FALSE)
single_capture <- file.path(test_dir, "single-capture")
dir.create(single_capture)
single_output <- file.path(test_dir, "single.pdf")
args <- c(wrapper, "--results", single_result, "--output", single_output)
status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args),
                  env = c(paste0("HEATMAP_CAPTURE_DIR=", shQuote(single_capture)),
                          paste0("HEATMAP_SCRIPT=", shQuote(script))))
stopifnot(status == 0, file.info(single_output)$size > 0)
single_plots <- lapply(file.path(single_capture, paste0(seq_len(3), ".rds")), readRDS)
stopifnot(
    length(list.files(single_capture, pattern = "\\.rds$")) == 3,
    all(vapply(single_plots, function(plot) identical(colnames(plot$mat), "Single"), logical(1))),
    all(vapply(single_plots, function(plot) identical(plot$angle_col, "270"), logical(1))),
    all(vapply(single_plots, function(plot) identical(plot$cellwidth, 6), logical(1))),
    all(vapply(single_plots, function(plot) identical(plot$cluster_cols, FALSE), logical(1))),
    identical(single_plots[[1]]$annotation_col["Single", "sample_class"], "unknown")
)

pdf(NULL)
example <- pheatmap::pheatmap(
    matrix(c(1, 2, 3, 5), nrow = 2, dimnames = list(c("TF1", "TF2"), c("A -> X", "B -> Y"))),
    annotation_col = data.frame(group = c("one", "two"), row.names = c("A -> X", "B -> Y")),
    main = "Example", angle_col = "270", silent = TRUE
)
layout <- example$gtable
location <- function(name) layout$layout$t[match(name, layout$layout$name)]
label_grob <- layout$grobs[[match("col_names", layout$layout$name)]]
stopifnot(
    location("main") < location("col_tree"),
    location("col_tree") < location("col_annotation"),
    location("col_annotation") < location("matrix"),
    location("matrix") < location("col_names"),
    location("col_annotation_names") == location("col_annotation"),
    label_grob$rot == 270, label_grob$hjust == 0,
    identical(label_grob$label, c("A -> X", "B -> Y")[example$tree_col$order])
)
invisible(dev.off())
cat("Heatmap transpose tests passed\n")
