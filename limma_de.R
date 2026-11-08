suppressPackageStartupMessages({
  library(limma)
})

# ---------------------------------------------------------------------
# Differential expression with limma (control vs. treatment contrast).
#
# Arguments:
#   expr_mat     : genes x samples matrix; sample names <batch>_<sampleId>_<group>
#   control, treat : group labels used in the contrast
#   logfc_filter : absolute logFC threshold
#   p_filter     : raw p-value threshold (adjust = "none", as in original analysis)
#
# Returns: list(all, diff, expr)
# ---------------------------------------------------------------------
run_limma_de <- function(expr_mat, control = "Control", treat = "PMOP",
                         logfc_filter = 0, p_filter = 0.05) {
  data <- avereps(as.matrix(expr_mat))
  
  group <- sub(".*_", "", colnames(data))
  data  <- data[, order(group), drop = FALSE]
  group <- sub(".*_", "", colnames(data))
  colnames(data) <- sub("^[^_]+_([^_]+)_.*$", "\\1", colnames(data))
  
  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- levels(factor(group))
  
  fit <- lmFit(data, design)
  contrast <- makeContrasts(paste0(treat, "-", control), levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  all_diff <- topTable(fit2, adjust = "none", number = Inf)
  diff_sig <- all_diff[abs(all_diff$logFC) > logfc_filter &
                         all_diff$P.Value < p_filter, , drop = FALSE]
  diff_exp <- data[rownames(diff_sig), , drop = FALSE]
  
  list(all = all_diff, diff = diff_sig, expr = diff_exp)
}