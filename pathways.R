suppressPackageStartupMessages({
  library(limma)
  library(GSEABase)
  library(GSVA)
})

# ---------------------------------------------------------------------
# ssGSEA per sample, then t-test between high/low target-gene groups.
#
# Arguments:
#   expr_mat     : genes x samples matrix
#   gmt_file     : gene-set file
#   target_gene  : gene used to split samples at the median
#   sample_group : restrict to a subset of samples (e.g. disease only)
#
# Returns:
#   data.frame(Pathway, t, pvalue, Sig)
# ---------------------------------------------------------------------
run_gsva <- function(expr_mat, gmt_file, target_gene,
                     sample_group = "PMOP") {
  data <- avereps(as.matrix(expr_mat))
  group <- sub(".*_", "", colnames(data))
  data  <- data[, group == sample_group, drop = FALSE]
  
  gene_sets <- getGmt(gmt_file, geneIdType = SymbolIdentifier())
  par   <- ssgseaParam(data, gene_sets, minSize = 5, maxSize = 500)
  score <- gsva(par)
  
  normalize <- function(x) (x - min(x)) / (max(x) - min(x))
  score <- normalize(score)
  score <- score[apply(score, 1, sd) > 0.01, , drop = FALSE]
  
  low  <- colnames(data)[data[target_gene, ] <  median(data[target_gene, ])]
  high <- colnames(data)[data[target_gene, ] >= median(data[target_gene, ])]
  
  score_group <- cbind(score[, low, drop = FALSE],
                       score[, high, drop = FALSE])
  group <- c(rep("Low", length(low)), rep("High", length(high)))
  
  out <- data.frame()
  for (pathway in rownames(score_group)) {
    test  <- t.test(score_group[pathway, ] ~ group)
    t_val <- test$statistic
    if (test$estimate[2] > test$estimate[1]) t_val <- abs(t_val) else t_val <- -abs(t_val)
    sig <- ifelse(test$p.value > 0.05, "Not",
                  ifelse(t_val > 0, "Up", "Down"))
    out <- rbind(out, data.frame(Pathway = pathway, t = t_val,
                                 pvalue = test$p.value, Sig = sig))
  }
  
  out$Pathway <- vapply(out$Pathway, function(x) {
    x <- sub("^KEGG_", "", x)
    x <- gsub("_", " ", x)
    x <- tolower(x)
    words <- strsplit(x, " ")[[1]]
    paste0(toupper(substr(words, 1, 1)),
           substr(words, 2, nchar(words)),
           collapse = " ")
  }, character(1))
  
  out
}




suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
})

# ---------------------------------------------------------------------
# Gene-set enrichment analysis between high/low target-gene groups.
#
# Returns:
#   data.frame with GSEA columns (ID, NES, pvalue, p.adjust,
#   core_enrichment, …), filtered to pvalue < 0.05.
# ---------------------------------------------------------------------
run_gsea <- function(expr_mat, gene_list, gmt_file, target_gene,
                     sample_group = "PMOP") {
  data <- avereps(as.matrix(expr_mat))
  
  keep <- rowMeans(data) > 0
  keep[intersect(gene_list, rownames(data))] <- TRUE
  data <- data[keep, , drop = FALSE]
  
  group <- sub(".*_", "", colnames(data))
  data  <- data[, group == sample_group, drop = FALSE]
  
  low  <- data[, data[target_gene, ] <  median(data[target_gene, ]), drop = FALSE]
  high <- data[, data[target_gene, ] >= median(data[target_gene, ]), drop = FALSE]
  
  mean_low  <- rowMeans(low);  mean_low[mean_low   < 1e-5] <- 1e-5
  mean_high <- rowMeans(high); mean_high[mean_high < 1e-5] <- 1e-5
  logfc <- sort(mean_high - mean_low, decreasing = TRUE)
  
  gmt <- read.gmt(gmt_file)
  kk  <- GSEA(logfc, TERM2GENE = gmt, pvalueCutoff = 1)
  kk_tab <- as.data.frame(kk)
  kk_tab[kk_tab$pvalue < 0.05, , drop = FALSE]
}