#!/usr/bin/env Rscript

library(SlopeHunter)
library(data.table)
library(tools)

# ------------------------------------------------------------------
# 1. Read and standardise each GWAS
# ------------------------------------------------------------------

incidence <- read_incidence(
  filename = "MR/incidence/for_slopehunter/Luminal_B_rsid.txt",
  snp_col = "rsid",
  beta_col = "beta", 
  se_col = "se",
  pval_col = "pval",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  chr_col = "chr",
  pos_col = "pos"
)

prognosis <- read_prognosis(
  filename = "MR/outcome/survival_2021_gwas_results_files_with_rsid_no_chr23/h_ERPRHER2posi_15_years_with_rsid.txt",
  snp_col = "rsid",
  beta_col = "Beta",
  se_col = "SE",
  pval_col = "P",
  eaf_col = "eaf",
  effect_allele_col = "a1",
  other_allele_col = "a0",
  chr_col = "Chromosome",
  pos_col = "Position"
)

# Remove rows already found to lack essential data.
incidence <- incidence[
  incidence$SH_keep.incidence %in% TRUE,
]

prognosis <- prognosis[
  prognosis$SH_keep.prognosis %in% TRUE,
]

# ------------------------------------------------------------------
# 2. Match variants and harmonise alleles
# ------------------------------------------------------------------

harmonised <- harmonise_effects(
  incidence_dat = incidence,
  prognosis_dat = prognosis,

  # Prefer rsID matching. See caveat about by.pos below.
  by.pos = FALSE
)

harmonisation_info <- attr(harmonised, "info")
print(harmonisation_info)

harmonised <- as.data.frame(harmonised)

# harmonise_effects() marks exclusions but does not remove them.
analysis_data <- harmonised[
  harmonised$remove %in% FALSE &
    is.finite(harmonised$BETA.incidence) &
    is.finite(harmonised$SE.incidence) &
    harmonised$SE.incidence > 0 &
    is.finite(harmonised$PVAL.incidence) &
    harmonised$PVAL.incidence >= 0 &
    harmonised$PVAL.incidence <= 1 &
    is.finite(harmonised$BETA.prognosis) &
    is.finite(harmonised$SE.prognosis) &
    harmonised$SE.prognosis > 0 &
    is.finite(harmonised$PVAL.prognosis),
]

analysis_data <- analysis_data[!duplicated(analysis_data$SNP), ]

# Setup output directory and logging
out_dir <- file.path("MR", "slopehunter", "slopehunter_out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "error_log.txt")
if (file.exists(log_file)) file.remove(log_file)
cat("time\tfile\terror\n", file = log_file)

# Initialize summary stats file
summary_file <- file.path(out_dir, "slopehunter_summary_stats.tsv")
if (file.exists(summary_file)) file.remove(summary_file)
cat(paste(c("phenotype", "slope", "slope_se", "ci_lower", "ci_upper", "mixture_pi", "iterations"), collapse = "\t"), "\n", file = summary_file)

# ------------------------------------------------------------------
# 3. LD prune before fitting SlopeHunter
# ------------------------------------------------------------------

# Perform clumping locally using PLINK (avoids OpenGWAS API payload limits).
# Requires a local PLINK executable and an EUR reference at:
# MR/reference/EUR.bed, MR/reference/EUR.bim, MR/reference/EUR.fam

# Write PLINK clump input (SNP + P)
clump_input_file <- file.path(out_dir, "plink_clump_input.txt")
write.table(
  analysis_data[, c("SNP", "PVAL.incidence")],
  file = clump_input_file,
  row.names = FALSE,
  col.names = c("SNP", "P"),
  quote = FALSE,
  sep = "\t"
)

# Prepare PLINK command
plink_exec <- Sys.which("plink")
if (is.na(plink_exec) || plink_exec == "") {
  plink_exec <- "plink"
}
plink_bfile <- file.path("MR", "reference", "EUR")
plink_out_prefix <- file.path(out_dir, "plink_clump")

plink_args <- c(
  "--bfile", plink_bfile,
  "--clump", clump_input_file,
  "--clump-field", "P",
  "--clump-p1", "1",
  "--clump-p2", "1",
  "--clump-r2", "0.1",
  "--clump-kb", "250",
  "--out", plink_out_prefix
)

message("Running PLINK clumping (this may take a while)...")
plink_status <- system2(plink_exec, args = plink_args)
if (!is.null(plink_status) && plink_status != 0) {
  stop("PLINK clumping failed with exit code: ", plink_status)
}

# Read PLINK .clumped output and extract lead SNPs
clumped_file <- paste0(plink_out_prefix, ".clumped")
if (!file.exists(clumped_file)) {
  stop("Expected PLINK output not found: ", clumped_file)
}
clumped <- data.table::fread(clumped_file, data.table = FALSE)
if ("SNP" %in% colnames(clumped)) {
  lead_snps <- as.character(clumped$SNP)
} else if (ncol(clumped) >= 3) {
  lead_snps <- as.character(clumped[[3]])
} else {
  stop("Cannot locate SNP column in PLINK .clumped output")
}

# Subset the analysis data to the clumped lead SNPs
pruned <- analysis_data[analysis_data$SNP %in% lead_snps, ]


# ------------------------------------------------------------------
# 4. Fit the correction factor
# ------------------------------------------------------------------

lambda <- 1e-3

n_for_mixture <- sum(
  pruned$PVAL.incidence <= lambda,
  na.rm = TRUE
)

message(
  "Independent variants entering the mixture model: ",
  n_for_mixture
)

fit <- hunt(
  dat = pruned,

  # Always specify these explicitly.
  snp_col = "SNP",
  xbeta_col = "BETA.incidence",
  xse_col = "SE.incidence",
  xp_col = "PVAL.incidence",
  ybeta_col = "BETA.prognosis",
  yse_col = "SE.prognosis",
  yp_col = "PVAL.prognosis",

  xp_thresh = lambda,
  init_pi = 0.6,
  init_sigmaIP = 1e-5,

  Bootstrapping = TRUE,
  M = 100,                # Consider more for final analysis
  seed = 20260813,

  Plot = TRUE,
  show_adjustments = FALSE
)

name <- file_path_sans_ext(basename("MR/outcome/survival_2021_gwas_results_files_with_rsid_no_chr23/h_ERPRHER2posi_15_years_with_rsid.txt"))

out_stats <- data.frame(
  phenotype = name,
  slope = fit$b,
  slope_se = fit$bse,
  ci_lower = fit$b_CI[1],
  ci_upper = fit$b_CI[2],
  mixture_pi = fit$pi,
  iterations = fit$iter,
  stringsAsFactors = FALSE
)
data.table::fwrite(out_stats, file = summary_file, append = TRUE, sep = "\t", col.names = FALSE)
message("Saved stats for ", name)

# Save model object and plot
fit_file <- file.path(out_dir, paste0(name, "_fit.rds"))
saveRDS(fit, file = fit_file)
pdf_file <- file.path(out_dir, paste0(name, "_plot.pdf"))
pdf(pdf_file, width = 8, height = 6)
if (!is.null(fit$plot)) {
  print(fit$plot)
} else {
  tryCatch(
    plot(fit),
    error = function(e) {
      msg <- paste0(Sys.time(), "\t", name, "\tplot_error\t", conditionMessage(e), "\n")
      cat(msg, file = log_file, append = TRUE)
      message("Plotting failed for ", name, ": ", conditionMessage(e))
    }
  )
}
dev.off()

print(table(fit$Fit$clusters))
print(summary(fit$Fit$pt))

# ------------------------------------------------------------------
# 5. Apply the fitted factor to the complete harmonised dataset
# ------------------------------------------------------------------

adjusted <- SHadj(
  x = fit,
  dat = analysis_data,
  snp_col = "SNP",
  xbeta_col = "BETA.incidence",
  xse_col = "SE.incidence",
  ybeta_col = "BETA.prognosis",
  yse_col = "SE.prognosis"
)

data.table::fwrite(
  adjusted,
  file = "conditional_outcome.slopehunter.tsv.gz",
  sep = "\t",
  quote = FALSE
)

# Plot produced during fitting
if (!is.null(fit$plot)) {
  tryCatch({ print(fit$plot) }, error = function(e) {
    msg <- paste0(Sys.time(), "\t", name, "\tplot_error\t", conditionMessage(e), "\n")
    cat(msg, file = log_file, append = TRUE)
    message("Plot display failed for ", name, ": ", conditionMessage(e))
  })
} else {
  tryCatch({ plot(fit) }, error = function(e) {
    msg <- paste0(Sys.time(), "\t", name, "\tplot_error\t", conditionMessage(e), "\n")
    cat(msg, file = log_file, append = TRUE)
    message("Plot display failed for ", name, ": ", conditionMessage(e))
  })
}