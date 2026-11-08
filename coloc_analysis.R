suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
})

# ---------------------------------------------------------------------
# Colocalization (coloc.abf) around the lead exposure SNP in a
# ±window-bp region. One coloc.abf call per probe/gene.
#
# Arguments:
#   exposure_file, outcome_file : file paths
#   window                      : half-window size in bp
#   outcome_type                : "cc" or "quant"
#   N_total, s_case_prop        : outcome sample size / case proportion
#   remove_MHC                  : drop the MHC region before coloc
#   build                       : "37" or "38" (applied to both datasets)
#
# Returns: data.table, one row per probe
#   Probe, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, Ratio, n_snps
# ---------------------------------------------------------------------
.remove_mhc <- function(data, chr_col = "chr", pos_col = "pos",
                        build = c("37", "38")) {
  build <- match.arg(build)
  data[[chr_col]] <- as.character(data[[chr_col]])
  if (build == "37") {
    d1 <- data[data[[chr_col]] != "6"]
    d2 <- data[data[[chr_col]] == "6" & data[[pos_col]] > 33448354]
    d3 <- data[data[[chr_col]] == "6" & data[[pos_col]] < 28477797]
  } else {
    d1 <- data[data[[chr_col]] != "6"]
    d2 <- data[data[[chr_col]] == "6" & data[[pos_col]] > 33480577]
    d3 <- data[data[[chr_col]] == "6" & data[[pos_col]] < 28510120]
  }
  rbind(d1, d2, d3)
}

run_coloc <- function(exposure_file, outcome_file,
                      window = 500000, outcome_type = "cc",
                      N_total = 281786, s_case_prop = 2001 / 281786,
                      remove_MHC = FALSE, build = "38") {
  
  exposure <- fread(exposure_file)
  setnames(exposure,
           old = c("SNP","Chr","BP","b","SE","A1","A2",
                   "Freq","p","Samplesize","Probe","Gene"),
           new = c("snp","chr","pos","beta","se","ALT","REF",
                   "eaf","pval","N","cell","gene"),
           skip_absent = TRUE)
  exposure[, chr := gsub("^chr", "", as.character(chr))]
  exposure[, pos := as.integer(pos)]
  exposure[, effect_allele := toupper(ALT)]
  exposure[, other_allele  := toupper(REF)]
  exposure[, Probe := gene]
  
  outcome <- fread(outcome_file)
  setnames(outcome,
           old = c("rsids","#chrom","pos","beta","sebeta",
                   "alt","ref","af_alt","pval"),
           new = c("snp","chr","pos","beta","se",
                   "effect_allele","other_allele","eaf","pval"),
           skip_absent = TRUE)
  outcome[, chr := gsub("^chr", "", as.character(chr))]
  outcome[, pos := as.integer(pos)]
  outcome[, effect_allele := toupper(effect_allele)]
  outcome[, other_allele  := toupper(other_allele)]
  
  results <- list()
  for (i in seq_along(unique(exposure$Probe))) {
    probe   <- unique(exposure$Probe)[i]
    exp_sub <- exposure[Probe == probe]
    if (nrow(exp_sub) == 0) next
    
    if (remove_MHC) {
      exp_sub <- .remove_mhc(exp_sub, "chr", "pos", build = build)
      if (nrow(exp_sub) == 0) next
    }
    
    lead_pos <- exp_sub$pos[which.min(exp_sub$pval)]
    this_chr <- unique(exp_sub$chr)
    if (length(this_chr) != 1) next
    
    start <- max(1, lead_pos - window)
    end   <- lead_pos + window
    outc_sub <- outcome[chr == this_chr & pos >= start & pos <= end]
    if (nrow(outc_sub) == 0) next
    
    if (remove_MHC) {
      outc_sub <- .remove_mhc(outc_sub, "chr", "pos", build = build)
      if (nrow(outc_sub) == 0) next
    }
    
    merged <- merge(exp_sub, outc_sub, by = "snp",
                    suffixes = c(".exposure", ".outcome"))
    if (nrow(merged) == 0) next
    
    merged[, flip := effect_allele.exposure != effect_allele.outcome &
             effect_allele.exposure == other_allele.outcome]
    merged[flip == TRUE, beta.outcome := -beta.outcome]
    if ("eaf.outcome" %in% names(merged)) {
      merged[flip == TRUE, eaf.outcome := 1 - eaf.outcome]
    }
    merged <- merged[
      effect_allele.exposure == effect_allele.outcome |
        (effect_allele.exposure == other_allele.outcome &
           !is.na(other_allele.outcome))
    ]
    merged <- merged[!is.na(snp) & snp != ""]
    merged[, snp := toupper(snp)]
    merged <- merged[!duplicated(snp)]
    setorder(merged, snp)
    if (nrow(merged) == 0) next
    
    eaf_exp <- if ("eaf.exposure" %in% names(merged)) merged$eaf.exposure else rep(NA_real_, nrow(merged))
    eaf_out <- if ("eaf.outcome"  %in% names(merged)) merged$eaf.outcome  else rep(NA_real_, nrow(merged))
    
    maf_exp <- eaf_exp
    maf_exp[!is.na(maf_exp) & maf_exp > 0.5] <-
      1 - maf_exp[!is.na(maf_exp) & maf_exp > 0.5]
    maf_exp[is.na(maf_exp)] <- 0.5
    maf_out <- eaf_out
    maf_out[!is.na(maf_out) & maf_out > 0.5] <-
      1 - maf_out[!is.na(maf_out) & maf_out > 0.5]
    maf_out[is.na(maf_out)] <- maf_exp[is.na(maf_out)]
    maf_out[is.na(maf_out)] <- 0.5
    
    D1 <- list(snp = as.character(merged$snp),
               position = as.integer(merged$pos.exposure),
               beta = merged$beta.exposure,
               varbeta = merged$se.exposure^2,
               type = "quant",
               N = unique(exp_sub$N)[1],
               MAF = maf_exp)
    D2 <- list(snp = as.character(merged$snp),
               position = as.integer(merged$pos.exposure),
               beta = merged$beta.outcome,
               varbeta = merged$se.outcome^2,
               type = outcome_type,
               N = N_total,
               MAF = maf_out)
    if (outcome_type == "cc") D2$s <- s_case_prop
    
    res <- tryCatch(coloc.abf(D1, D2), error = function(e) NULL)
    if (is.null(res)) next
    
    pp_h3 <- res$summary["PP.H3.abf"]
    pp_h4 <- res$summary["PP.H4.abf"]
    results[[i]] <- data.table(
      Probe = probe,
      PP.H0 = res$summary["PP.H0.abf"],
      PP.H1 = res$summary["PP.H1.abf"],
      PP.H2 = res$summary["PP.H2.abf"],
      PP.H3 = pp_h3,
      PP.H4 = pp_h4,
      Ratio = pp_h4 / (pp_h3 + pp_h4),
      n_snps = nrow(merged)
    )
  }
  
  rbindlist(results, fill = TRUE)
}