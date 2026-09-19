#!/usr/bin/env Rscript

library(data.table)

#----------------------------
# Read BCAC subtype summary statistics
#----------------------------
dat <- fread("MR/incidence/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")

#----------------------------
# Define intrinsic subtypes
#----------------------------
subtypes <- c(
  "Luminal_A",
  "Luminal_B",
  "Luminal_B_HER2Neg",
  "HER2_Enriched",
  "Triple_Neg"
)

#----------------------------
# Generate one outcome file for each subtype
#----------------------------
for (st in subtypes) {

  beta_col <- paste0(st, "_log_or_meta")
  se_col   <- paste0(st, "_se_meta")

  out <- data.table(

    SNP = dat$var_name,

    beta = dat[[beta_col]],

    se = dat[[se_col]],

    effect_allele = dat$Effect.Meta,

    other_allele = dat$Baseline.Meta,

    chr = dat$chr.iCOGs,

    pos = dat$Position.iCOGs
  )

  # Calculate P value
  out[, pval := 2 * pnorm(-abs(beta / se))]

  # Remove missing values
  out <- out[
    !is.na(beta) &
    !is.na(se) &
    !is.na(effect_allele) &
    !is.na(other_allele)
  ]

  # Remove duplicated SNPs
  out <- unique(out, by = "SNP")

  # Write output
  fwrite(
    out,
    file = paste0("MR/incidence/", st, ".txt"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat(st, ":", nrow(out), "SNPs exported\n")
}

cat("\nFinished!\n")