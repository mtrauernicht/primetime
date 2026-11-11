scripts <- c("scripts/primetime_run_comparative_analysis_bcalm.R",
             "scripts/screen/primetime_run_comparative_analysis_bcalm.R")

for (script in scripts) {
    expressions <- parse(script)
    position <- which(vapply(expressions, function(expr) {
        is.call(expr) && identical(expr[[1]], as.name("=")) &&
            identical(expr[[2]], as.name("ordered_annotation"))
    }, logical(1)))
    stopifnot(length(position) == 1)
    design_expressions <- expressions[position + 0:2]
    make_design <- function(reference, contrast, treatments) {
        env <- new.env(parent = baseenv())
        env$reference_condition <- reference
        env$contrast_condition <- contrast
        env$annotation_df <- data.frame(treatment = treatments, extra = 1, row.names = paste0("obs", seq_along(treatments)))
        env$this_pdna_matrix <- matrix(0, nrow = 1, ncol = length(treatments),
                                       dimnames = list(NULL, rev(rownames(env$annotation_df))))
        for (expr in design_expressions) eval(expr, env)
        env$design_bcalm
    }

    for (contrast in c("Lorlatinib?PF-6463922", "ENMD-2076_L-+-Tartaric_acid_", "Drug.A", "Drug")) {
        design <- make_design(paste0(contrast, "_control"), contrast,
                              c(paste0(contrast, "_control"), contrast, contrast))
        stopifnot(identical(design[[2]], c(TRUE, TRUE, FALSE)), all(design[[1]] == 1))
    }
    for (treatments in list(c("DMSO", "DMSO"), c("DMSO", NA), c("DMSO", "Other"))) {
        err <- tryCatch(make_design("DMSO", "Drug", treatments), error = identity)
        stopifnot(inherits(err, "error"), grepl("Invalid BCalm design", conditionMessage(err), fixed = TRUE))
    }
    err <- tryCatch(make_design("Drug", "Drug", c("Drug", "Drug")), error = identity)
    stopifnot(inherits(err, "error"), grepl("Invalid BCalm design", conditionMessage(err), fixed = TRUE))
}
cat("Comparison design tests passed\n")
