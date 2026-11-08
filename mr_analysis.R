suppressPackageStartupMessages({
  library(TwoSampleMR)
})

# ---------------------------------------------------------------------
# Two-sample MR + sensitivity analyses.
#
# Arguments:
#   exposure_file, outcome_file : file paths (CSV / TSV.gz)
#   outcome_name                : label for the outcome trait
#   ncase, ncontrol             : optional; if provided, Steiger test runs
#
# Returns: list(dat, mr, heterogeneity, pleiotropy, steiger)
# ---------------------------------------------------------------------
run_mr <- function(exposure_file, outcome_file, outcome_name,
                   ncase = NULL, ncontrol = NULL) {
  exposure_dat <- read_exposure_data(
    filename = exposure_file, sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    pval_col = "pval.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col  = "other_allele.exposure",
    eaf_col           = "eaf.exposure",
    id_col            = "id.exposure",
    phenotype_col     = "exposure",
    samplesize_col    = "samplesize.exposure",
    chr_col           = "chr.exposure",
    pos_col           = "pos.exposure",
    clump             = FALSE
  )
  
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP, filename = outcome_file, sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col  = "ref",
    pval_col = "pval",
    eaf_col  = "af_alt"
  )
  outcome_dat$outcome <- outcome_name
  
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  mr_result <- mr(dat)
  mr_or     <- generate_odds_ratios(mr_result)
  
  heter <- mr_heterogeneity(dat)
  pleio <- mr_pleiotropy_test(dat)
  
  steiger <- NULL
  if (!is.null(ncase) && !is.null(ncontrol)) {
    dat$ncase.outcome    <- ncase
    dat$ncontrol.outcome <- ncontrol
    dat$samplesize.outcome <- ncase + ncontrol
    dat$r.outcome <- get_r_from_lor(
      lor = dat$beta.outcome, af = dat$eaf.outcome,
      ncase = ncase, ncontrol = ncontrol,
      prevalence = ncase / (ncase + ncontrol)
    )
    dat$r.exposure <- get_r_from_pn(p = dat$pval.exposure,
                                    n = dat$samplesize.exposure)
    steiger <- directionality_test(dat)
  }
  
  list(dat = dat, mr = mr_or, heterogeneity = heter,
       pleiotropy = pleio, steiger = steiger)
}