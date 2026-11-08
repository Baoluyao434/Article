# Reproducible multi-omics analysis pipeline

This repository contains the computational code accompanying our study.
It is split into three independent modules:

| Module | Folder | Language |
|---|---|---|
| Bulk transcriptomics | `transcriptome/` | R, Python |
| Mendelian randomization | `mendelian_randomization/` | R |
| Single-cell RNA-seq | `single_cell/` | R |

Every function is a **pure analysis function** that takes inputs as arguments
and **returns results as R/Python objects**. No function writes files or
creates plots. Users are expected to chain functions together and decide how
to persist the returned objects (e.g. `saveRDS()` for R objects,
`adata.write_h5ad()` for AnnData objects). This design keeps the code
reusable, testable and free of hard-coded paths.

---

## 1. Requirements

### 1.1 R (≥ 4.2)

Install packages per module as needed. `environment.R` lists them all:

```r
install.packages(c(
  "BiocManager",
  # transcriptomics
  "limma", "sva", "WGCNA", "GSEABase", "GSVA",
  "clusterProfiler", "pROC", "glmnet", "dplyr",
  # mendelian randomization
  "TwoSampleMR", "coloc", "MRPRESSO", "data.table",
  # single-cell
  "Seurat", "harmony", "monocle", "CellChat",
  "scTenifoldKnk"
))

BiocManager::install(c(
  "limma", "sva", "GSEABase", "GSVA",
  "clusterProfiler", "monocle", "VariantAnnotation"
))

# MR and single-cell packages from GitHub
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("rondolab/MR-PRESSO")
# devtools::install_github("cailab-tamu/scTenifoldKnk")
```
---

## 2. Expected input formats

### 2.1 Expression matrices (transcriptomics, GSVA, GSEA, ROC, WGCNA)

Tab-separated `.txt` file:

- **First column**: gene symbol or probe ID (no name required).
- **Header row**: sample IDs.
- **Sample ID convention**: `<batch>_<sampleId>_<group>`, e.g.
  `GSE1_S01_Control`, `GSE1_S02_PMOP`. The **last underscore-delimited
  field** must be the group label used by downstream functions
  (`Control`, `PMOP`, …).

### 2.2 Gene set files (GSVA / GSEA)

Standard GMT format, e.g. `c2.cp.kegg.Hs.symbols.gmt`,
`c5.go.Hs.symbols.gmt`.

### 2.3 Gene lists (GSVA target / ROC / GSEA filtering)

Plain text, one gene symbol per line.

### 2.4 MR exposure file (CSV)

Must be readable by `TwoSampleMR::read_exposure_data()`, with columns
`SNP, beta.exposure, se.exposure, pval.exposure,
effect_allele.exposure, other_allele.exposure, eaf.exposure,
id.exposure, exposure, samplesize.exposure, chr.exposure, pos.exposure`.

### 2.5 MR outcome file (FinnGen-style TSV)

Columns: `rsids, #chrom, pos, beta, sebeta, alt, ref, af_alt, pval`
(compressed `.gz` supported).

### 2.6 Colocalization exposure file (TSV)

Columns: `SNP, Chr, BP, b, SE, A1, A2, Freq, p, Samplesize, Probe, Gene`.

### 2.7 10x Genomics directory

A directory of subdirectories, each a 10x `filtered_feature_bc_matrix/`
produced by `cellranger`; subdirectory names become `orig.ident`.

---

## 3. Usage

### 3.1 Transcriptomics

```r
source("transcriptome/batch_correct.R")
source("transcriptome/limma_de.R")
source("transcriptome/wgcna_modules.R")
source("transcriptome/gsva_pathways.R")
source("transcriptome/gsea_pathways.R")
source("transcriptome/roc_model.R")

# Step 1: batch correction across multiple studies
expr_files <- list.files("data/transcriptome", pattern = "\\.txt$",
                         full.names = TRUE)
expr_norm  <- batch_correct(expr_files)         # matrix, genes x samples

# Step 2: differential expression
de <- run_limma_de(expr_norm,
                   control = "Control", treat = "PMOP",
                   logfc_filter = 0, p_filter = 0.05)
# de$all  : all genes
# de$diff : filtered DE genes
# de$expr : DE gene expression matrix

# Step 3: WGCNA (returns module assignments + MM/GS table)
wgcna_res <- run_wgcna(expr_norm, soft_power = 9)
# wgcna_res$gene_info, $module_colors, $MEs, $module_genes, ...

# Step 4: GSVA on a single target gene
gsva_tab <- run_gsva(expr_norm,
                     gmt_file    = "data/gmt/c2.cp.kegg.Hs.symbols.gmt",
                     target_gene = "CAMLG",
                     sample_group = "PMOP")

# Step 5: GSEA (requires a gene list, e.g. DE gene IDs)
gsea_tab <- run_gsea(expr_norm,
                     gene_list    = de$diff$ID,
                     gmt_file     = "data/gmt/c5.go.Hs.symbols.gmt",
                     target_gene  = "CAMLG",
                     sample_group = "PMOP")

# Step 6: ROC / logistic model on selected genes
roc_res <- run_roc_model(expr_norm, genes = de$diff$ID)
# roc_res$gene_auc, $model_auc, $model_ci, $pred
```

### 3.2 Mendelian randomization

```r
source("mendelian_randomization/mr_analysis.R")
source("mendelian_randomization/coloc_analysis.R")

mr_res <- run_mr(
  exposure_file = "data/mr/exposure.eqtl.csv",
  outcome_file  = "data/mr/finngen_R12_OSTPOPATFRCTURE_POSTEMENO.gz",
  outcome_name  = "PMOP",
  ncase = 2001, ncontrol = 279785
)
# mr_res$mr, $dat, $heterogeneity, $pleiotropy, $steiger

coloc_res <- run_coloc(
  exposure_file = "data/mr/exposure_gene.tsv",
  outcome_file  = "data/mr/finngen_R12_OSTPOPATFRCTURE_POSTEMENO.gz",
  window = 5e5,
  outcome_type = "cc",
  N_total = 281786,
  s_case_prop = 2001 / 281786
)
# coloc_res is a data.table: one row per probe, columns PP.H0..H4, Ratio, n_snps
```

### 3.3 Single-cell RNA-seq

```r
source("single_cell/seurat_process.R")
source("single_cell/cellchat_network.R")
source("single_cell/monocle_trajectory.R")
source("single_cell/sctenifoldknk_knockout.R")

# Step 1: Seurat processing + manual annotation
cluster_map <- c("0" = "OC", "1" = "OC", "2" = "TC", "3" = "OC",
                 "4" = "MC", "5" = "MC", "6" = "MON", "7" = "pre-BC",
                 "8" = "BC", "9" = "PC", "10" = "MC", "11" = "MC",
                 "12" = "TC", "13" = "MAC", "14" = "EC", "15" = "OB",
                 "16" = "OB", "17" = "PC", "18" = "MAC", "19" = "pre-BC",
                 "20" = "pDC", "21" = "SMC", "22" = "pre-BC")

seurat_obj <- run_seurat("data/scRNA/", cluster_id = cluster_map)

# Step 2: CellChat
cc <- run_cellchat(seurat_obj, group_col = "celltype")
# cc@net, cc@netP, cc@netP$centr

# Step 3: Monocle pseudotime
mono <- run_monocle(
  seurat_obj     = seurat_obj,
  group_col      = "celltype",
  ordering_genes = rownames(mono_result$diff_test)[1:200]
)
# mono$cds, mono$diff_test

# Step 4: scTenifoldKnk in silico knockout
knk <- run_sctenifoldknk(seurat_obj,
                         target_gene = "CAMLG",
                         disease_type = "Disease")
# knk$sig : significant downstream genes
```

### 3.4 V-score (Python)

```python
import pandas as pd
from transcriptome.vscore import run_vscore

df = pd.read_csv("data/merge.normalize.txt", sep="\t", index_col=0)
adata, vscores_df = run_vscore(df)
# adata       : AnnData, one observation "Treat_vs_Control"
# vscores_df  : DataFrame with columns gene, v_score, mean_control, mean_treat
# Persist yourself, e.g. adata.write_h5ad("input.h5ad")
```

---

## 4. Design conventions

1. **No file I/O inside functions** other than reading inputs.
   All results are returned as R lists, `data.frame`/`data.table`, S4
   objects (Seurat, CellChat, Monocle, `coloc`) or `AnnData`. Saving is
   the caller's responsibility.
2. **No plotting code** is included. Users can build plots from the
   returned objects with any toolkit (ggplot2, pheatmap, matplotlib,
   scanpy, …).
3. **No hard-coded paths** or per-study sample prefixes. Sample grouping
   is inferred from the `<batch>_<sampleId>_<group>` convention; any
   other scheme requires adjusting the `sub(".*_", "", colnames(...))`
   calls.
4. **Group labels** (`Control`, `PMOP`, `Disease`, …) are function
   arguments, not globals.
5. **Random seed**: functions that require it (`scTenifoldKnk`,
   bootstrap CIs in `pROC`) accept a seed argument or call `set.seed()`
   locally.

---

## 5. Reproducibility notes

- Record R and Python versions (`sessionInfo()`, `pip freeze`) alongside
  any results.
- `WGCNA` soft power selection, `Seurat` clustering resolution and
  `CellChat` `min.cells` are the main user-tunable parameters; all are
  exposed as arguments and should be reported in the Methods section.
- Colocalization results depend heavily on the `N_total` and
  `s_case_prop` supplied. Always report these with the results.
- `scTenifoldKnk` is stochastic; set `set.seed()` before running and
  report the seed.
- `v-score` uses population variance (`ddof=0`), matching
  `drugreflector`'s `compute_vscore_two_groups`.

---

## 6. License

Copyright (c) 2025 The Authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.