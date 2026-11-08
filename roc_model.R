suppressPackageStartupMessages({
  library(pROC)
  library(glmnet)
})

# ---------------------------------------------------------------------
# Per-gene ROC AUC and multivariate logistic-regression ROC.
#
# Arguments:
#   expr_mat      : genes x samples matrix
#   genes         : gene vector for single-gene AUC + model
#   negative_group: group encoded as 0 (default "Control")
#
# Returns: list(gene_auc, model_auc, model_ci, pred)
# ---------------------------------------------------------------------
run_roc_model <- function(expr_mat, genes, negative_group = "Control") {
  group <- sub(".*_", "", colnames(expr_mat))
  y <- as.integer(group != negative_group)
  genes <- intersect(genes, rownames(expr_mat))
  
  gene_auc <- do.call(rbind, lapply(genes, function(g) {
    r <- roc(y, as.numeric(expr_mat[g, ]), quiet = TRUE)
    data.frame(gene = g, AUC = as.numeric(r$auc))
  }))
  
  dat   <- as.data.frame(t(expr_mat[genes, , drop = FALSE]))
  logit <- glm(y ~ ., family = binomial(link = "logit"), data = dat)
  pred  <- predict(logit, newx = dat, type = "response")
  
  roc_model <- roc(y, as.numeric(pred), quiet = TRUE)
  ci <- ci.auc(roc_model, method = "bootstrap")
  
  list(gene_auc = gene_auc,
       model_auc = as.numeric(roc_model$auc),
       model_ci  = c(lower = as.numeric(ci[1]),
                     upper = as.numeric(ci[3])),
       pred      = pred)
}