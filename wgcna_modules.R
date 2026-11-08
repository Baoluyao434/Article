suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
})

# ---------------------------------------------------------------------
# WGCNA module detection and module-trait association.
#
# Arguments:
#   expr_mat          : genes x samples matrix
#   soft_power        : if NULL, chosen automatically by pickSoftThreshold
#   min_module_size   : minimum genes per module (dynamic tree cut)
#   merge_cut_height  : threshold for merging correlated modules
#
# Returns: list with elements
#   soft_power, module_colors, MEs,
#   module_trait_cor, module_trait_p,
#   gene_info (moduleColor + MM/GS + p-values),
#   module_genes (named list: module -> gene IDs)
# ---------------------------------------------------------------------
run_wgcna <- function(expr_mat, soft_power = NULL,
                      min_module_size = 50, merge_cut_height = 0.3) {
  enableWGCNAThreads()
  
  data <- avereps(as.matrix(expr_mat))
  data <- data[apply(data, 1, sd) > 0.2, , drop = FALSE]
  
  group <- sub(".*_", "", colnames(data))
  data  <- data[, order(group), drop = FALSE]
  group <- sub(".*_", "", colnames(data))
  
  dat_expr0 <- t(data)
  gsg <- goodSamplesGenes(dat_expr0, verbose = 0)
  dat_expr0 <- dat_expr0[gsg$goodSamples, gsg$goodGenes]
  
  groups <- unique(group)
  traits <- sapply(groups, function(g) as.integer(group == g))
  if (is.null(dim(traits))) traits <- matrix(traits, ncol = 1)
  colnames(traits) <- groups
  rownames(traits) <- colnames(data)
  
  same <- intersect(rownames(dat_expr0), rownames(traits))
  dat_expr0  <- dat_expr0[same, , drop = FALSE]
  dat_traits <- traits[same, , drop = FALSE]
  
  if (is.null(soft_power)) {
    sft <- pickSoftThreshold(dat_expr0, powerVector = 1:20, verbose = 0)
    soft_power <- sft$powerEstimate
    if (is.na(soft_power)) soft_power <- 9
  }
  
  adjacency <- adjacency(dat_expr0, power = soft_power)
  diss_tom  <- 1 - TOMsimilarity(adjacency)
  
  gene_tree    <- hclust(as.dist(diss_tom), method = "average")
  dynamic_mods <- cutreeDynamic(dendro = gene_tree, distM = diss_tom,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = min_module_size)
  dynamic_colors <- labels2colors(dynamic_mods)
  
  me_list <- moduleEigengenes(dat_expr0, colors = dynamic_colors)
  mes     <- me_list$eigengenes
  
  merged <- mergeCloseModules(dat_expr0, dynamic_colors,
                              cutHeight = merge_cut_height, verbose = 0)
  module_colors <- merged$colors
  mes           <- merged$newMEs
  
  module_trait_cor <- cor(mes, dat_traits, use = "p")
  module_trait_p   <- corPvalueStudent(module_trait_cor, nrow(dat_expr0))
  
  mod_names <- substring(names(mes), 3)
  gmm  <- as.data.frame(cor(dat_expr0, mes, use = "p"))
  mm_p <- as.data.frame(corPvalueStudent(as.matrix(gmm), nrow(dat_expr0)))
  names(gmm)  <- paste0("MM", mod_names)
  names(mm_p) <- paste0("p.MM", mod_names)
  
  gs   <- as.data.frame(cor(dat_expr0, dat_traits, use = "p"))
  gs_p <- as.data.frame(corPvalueStudent(as.matrix(gs), nrow(dat_expr0)))
  names(gs)   <- paste0("GS.", colnames(dat_traits))
  names(gs_p) <- paste0("p.GS.", colnames(dat_traits))
  
  gene_info <- data.frame(probes = colnames(dat_expr0),
                          moduleColor = module_colors)
  gene_info <- cbind(gene_info, gs, gs_p, gmm, mm_p)
  gene_info <- gene_info[order(gene_info$moduleColor), ]
  
  module_genes <- split(colnames(dat_expr0), module_colors)
  
  list(soft_power = soft_power, module_colors = module_colors,
       MEs = mes, module_trait_cor = module_trait_cor,
       module_trait_p = module_trait_p, gene_info = gene_info,
       module_genes = module_genes)
}