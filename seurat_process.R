suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
})

# ---------------------------------------------------------------------
# Standard Seurat pipeline with Harmony batch correction and
# manual annotation via a provided cluster -> cell-type map.
#
# Arguments:
#   data_dir          : directory of 10x subdirectories
#   resolution        : Louvain resolution for FindClusters
#   pc_select         : number of PCs used for neighbors / UMAP / tSNE
#   logfc_threshold   : FindAllMarkers logFC threshold
#   adjp_filter       : FindAllMarkers adjusted-p threshold
#   cluster_id        : named character vector (cluster -> celltype)
#
# Returns: Seurat object with metadata columns celltype, Type, group
# ---------------------------------------------------------------------
run_seurat <- function(data_dir, resolution = 1, pc_select = 10,
                       logfc_threshold = 0.5, adjp_filter = 0.05,
                       cluster_id = NULL) {
  sample_dirs <- list.dirs(data_dir, recursive = FALSE)
  names(sample_dirs) <- basename(sample_dirs)
  
  counts <- Read10X(data.dir = sample_dirs)
  pbmc   <- CreateSeuratObject(counts, min.cells = 3, min.features = 100)
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc[["percent.hb"]] <- PercentageFeatureSet(pbmc, pattern = "^HB[AB]")
  pbmc <- subset(pbmc, nFeature_RNA > 200 & percent.hb < 3 & percent.mt < 15)
  
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                        scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                               nfeatures = 2000)
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc, npcs = 20, features = VariableFeatures(pbmc))
  pbmc <- RunHarmony(pbmc, "orig.ident", project.dim = FALSE, max_iter = 50)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:pc_select)
  pbmc <- FindClusters(pbmc, resolution = resolution)
  pbmc <- RunUMAP(pbmc, dims = 1:pc_select)
  pbmc <- RunTSNE(pbmc, dims = 1:pc_select)
  
  if (!is.null(cluster_id)) {
    pbmc <- RenameIdents(pbmc, cluster_id)
    pbmc$celltype <- pbmc@active.ident
  }
  
  type   <- sub("\\..*$", "", colnames(pbmc))
  groups <- paste0(type, "_", as.character(Idents(pbmc)))
  names(type)   <- colnames(pbmc)
  names(groups) <- colnames(pbmc)
  pbmc <- AddMetaData(pbmc, metadata = type,   col.name = "Type")
  pbmc <- AddMetaData(pbmc, metadata = groups, col.name = "group")
  
  pbmc
}