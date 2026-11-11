suppressPackageStartupMessages(library(dplyr))

expressions <- parse("scripts/primetime_get_activity_and_qc_plots.R")
for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% c("estimate_early_bleedthrough", "apply_bleedthrough_cutoff")) {
        eval(expr)
    }
}

counts <- data.frame(
    pDNA = c(rep(TRUE, 3), rep(FALSE, 9)),
    sample = c(rep("pDNA", 3), rep("A", 3), rep("B", 3), rep("C", 3)),
    replicate = c(rep("pDNA_1", 3), rep("A_1", 3), rep("B_1", 3), rep("C_1", 3)),
    barcode = rep(letters[1:3], 4),
    negative_control = TRUE,
    RPM = c(1, 2, 3, 1.5, 2, 2.5, 1.6, 2.2, 2.8, 1.4, 1.8, 2.2)
)
estimates <- apply_bleedthrough_cutoff(counts, 0.5)
stopifnot(
    isTRUE(all.equal(estimates$slope, c(0.5, 0.6, 0.4))),
    identical(estimates$replicate[estimates$excluded], "B_1"),
    nrow(apply_bleedthrough_cutoff(counts, -1)) == 0
)
expect_error <- function(expr, expected) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(expected, conditionMessage(error), fixed = TRUE))
}
expect_error(apply_bleedthrough_cutoff(counts, -0.5), "--max_bleedthrough")
expect_error(apply_bleedthrough_cutoff(counts[counts$barcode == "a", ], 0.5), "A_1")
expect_error(apply_bleedthrough_cutoff(counts[!counts$pDNA, ], 0.5), "no pDNA")
cat("Bleed-through filter tests passed\n")
