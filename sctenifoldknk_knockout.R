suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(Seurat)
})

# ---------------------------------------------------------------------
# In-silico gene knockout with scTenifoldKnk on the disease subgroup.
#
# Arguments:
#   seurat_obj    : processed Seurat object
#   target_gene   : gene to knock out
#   disease_type  : value of the Type metadata used to subset cells
#   nfeatures     : number of HVGs to include
#   n_net, n_cells: scTenifoldKnk subnetwork parameters
#   seed          : random seed
#
# Returns: list(result, sig)
#   result$diffRegulation (full table)
#   sig  (p.adj < 0.05, target gene excluded)
# ---------------------------------------------------------------------
run_sctenifoldknk <- function(seurat_obj, target_gene,
                              disease_type = "Disease",
                              nfeatures = 2000,
                              n_net = 10, n_cells = 500,
                              seed = 1234) {
  set.seed(seed)
  
  pbmc <- subset(seurat_obj, subset = Type == disease_type)
  count_mat <- GetAssayData(pbmc, layer = "counts")
  
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                               nfeatures = nfeatures)
  hvgs <- VariableFeatures(pbmc)
  data <- as.data.frame(count_mat[unique(c(target_gene, hvgs)), ])
  
  result <- scTenifoldKnk(
    countMatrix    = data,
    gKO            = target_gene,
    qc_mtThreshold = 0.1,
    qc_minLSize    = 500,
    nc_nNet        = n_net,
    nc_nCells      = n_cells
  )
  
  df  <- result$diffRegulation
  df  <- df[df$gene != target_gene, ]
  sig <- df[df$p.adj < 0.05, ]
  
  list(result = result, sig = sig)
}