suppressPackageStartupMessages({
  library(limma)
  library(sva)
})

# ---------------------------------------------------------------------
# Cross-study batch correction via ComBat.
#
# Arguments:
#   input_files : character vector of expression-matrix file paths
#                 (first column = gene ID, remaining columns = samples)
#
# Returns:
#   matrix : batch-corrected expression (genes x samples),
#            sample names prefixed with the study tag
# ---------------------------------------------------------------------
batch_correct <- function(input_files) {
  gene_list <- lapply(input_files, function(f) {
    rt <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
    unique(as.vector(rt[, 1]))
  })
  inter_genes <- Reduce(intersect, gene_list)
  
  all_tab    <- NULL
  batch_type <- c()
  
  for (i in seq_along(input_files)) {
    rt  <- read.table(input_files[i], header = TRUE, sep = "\t",
                      check.names = FALSE)
    rt  <- as.matrix(rt)
    rownames(rt) <- rt[, 1]
    exp <- rt[, 2:ncol(rt), drop = FALSE]
    
    data <- matrix(as.numeric(as.matrix(exp)),
                   nrow = nrow(exp),
                   dimnames = list(rownames(exp), colnames(exp)))
    data <- avereps(data)
    
    prefix <- sub("\\..*$", "", basename(input_files[i]))
    colnames(data) <- paste0(prefix, "_", colnames(data))
    
    if (is.null(all_tab)) {
      all_tab <- data[inter_genes, , drop = FALSE]
    } else {
      all_tab <- cbind(all_tab, data[inter_genes, , drop = FALSE])
    }
    batch_type <- c(batch_type, rep(i, ncol(data)))
  }
  
  ComBat(all_tab, batch_type, par.prior = TRUE)
}