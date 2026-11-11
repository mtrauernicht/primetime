test_dir <- tempfile("comparison-counts-")
dir.create(test_dir)
on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

cdna <- data.frame(
    tf = c("TF1", "TF2"), negative_control = c(FALSE, TRUE),
    promoter = c("mCMV", "hBGm"), barcode = c("AAAA", "CCCC"),
    check.names = FALSE
)
cdna[["GSK'963_1"]] <- c(5, NA)
cdna[["sample#2_1"]] <- c(3, 10)
expected_pdna <- data.frame(barcode = c("AAAA", "CCCC"), negative_control = c(FALSE, TRUE), pDNA = c(4, 9))
opt <- list(
    cdna = file.path(test_dir, "cdna.tsv"),
    cdna_baseline = file.path(test_dir, "baseline.tsv"),
    pdna = file.path(test_dir, "pdna.tsv")
)
write.table(cdna, opt$cdna, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cdna, opt$cdna_baseline, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expected_pdna, opt$pdna, sep = "\t", quote = FALSE, row.names = FALSE)

scripts <- c("scripts/primetime_run_comparative_analysis_bcalm.R",
             "scripts/screen/primetime_run_comparative_analysis_bcalm.R")
for (script in scripts) {
    assignments <- list()
    for (expr in parse(script)) {
        if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
            is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% c("pdna", "cdna_all", "cdna_baseline_all") &&
            is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("read.table"))) {
            assignments[[as.character(expr[[2]])]] <- expr
        }
    }
    expected <- c("pdna", "cdna_all", if (grepl("/screen/", script, fixed = TRUE)) "cdna_baseline_all")
    stopifnot(setequal(names(assignments), expected))
    for (assignment in assignments) eval(assignment)
    stopifnot(isTRUE(all.equal(cdna_all, cdna)), isTRUE(all.equal(pdna, expected_pdna)))
    if ("cdna_baseline_all" %in% expected) stopifnot(isTRUE(all.equal(cdna_baseline_all, cdna)))
}
cat("Comparison count-reading tests passed\n")
