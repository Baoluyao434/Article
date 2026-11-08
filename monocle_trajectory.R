suppressPackageStartupMessages({
  library(monocle)
  library(Seurat)
  library(dplyr)
})

# ---------------------------------------------------------------------
# Monocle 2 pseudotime ordering and cell-type differential test.
#
# Arguments:
#   seurat_obj     : processed Seurat object
#   group_col      : cell-type metadata column
#   ordering_genes : genes used to build the trajectory
#   root_type      : group label used to pick the root State
#
# Returns: list(cds, diff_test)
# ---------------------------------------------------------------------
run_monocle <- function(seurat_obj, group_col = "celltype",
                        ordering_genes, root_type = "AP") {
  monocle_matrix <- as.matrix(seurat_obj@assays$RNA$data)
  sample_df <- seurat_obj@meta.data
  gene_ann  <- data.frame(gene_short_name = rownames(monocle_matrix),
                          row.names = rownames(monocle_matrix))
  
  data <- as(as.matrix(monocle_matrix), "sparseMatrix")
  pd <- new("AnnotatedDataFrame", data = sample_df)
  fd <- new("AnnotatedDataFrame", data = gene_ann)
  cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
  
  pData(cds)$cell_type <- pData(cds)[[group_col]]
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- setOrderingFilter(cds, ordering_genes)
  cds <- reduceDimension(cds, max_components = 2,
                         reduction_method = "DDRTree")
  cds <- orderCells(cds)
  
  state_n <- pData(cds) %>%
    dplyr::filter(Type == root_type) %>%
    dplyr::group_by(State) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n))
  if (nrow(state_n) > 0) {
    cds <- orderCells(cds, root_state = as.character(state_n$State[1]))
  }
  
  diff_test <- differentialGeneTest(
    cds,
    fullModelFormulaStr    = " ~ cell_type + orig.ident",
    reducedModelFormulaStr = " ~ orig.ident",
    cores = 8
  )
  
  list(cds = cds, diff_test = diff_test)
}