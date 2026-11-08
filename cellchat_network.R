suppressPackageStartupMessages({
  library(CellChat)
})

# ---------------------------------------------------------------------
# CellChat network inference.
#
# Arguments:
#   seurat_obj  : processed Seurat object
#   group_col   : metadata column used as cell groups
#   db_category : CellChatDB category (default "Secreted Signaling")
#   min_cells   : filter interactions supported by < min_cells cells
#
# Returns: aggregated CellChat object (cc@net, cc@netP available)
# ---------------------------------------------------------------------
run_cellchat <- function(seurat_obj, group_col = "celltype",
                         db_category = "Secreted Signaling",
                         min_cells = 10) {
  exp_matrix <- as.matrix(seurat_obj@assays$RNA$data)
  meta <- data.frame(labels = seurat_obj[[group_col, drop = TRUE]])
  rownames(meta) <- colnames(seurat_obj)
  
  cellchat <- createCellChat(object = exp_matrix, meta = meta,
                             group.by = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  
  cellchat@DB <- subsetDB(CellChatDB.human, search = db_category)
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  cellchat
}