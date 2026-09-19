#!/usr/bin/env Rscript

library(SlopeHunter)
library(data.table)
library(tools)

# ------------------------------------------------------------------
# 1. Read and standardise each GWAS
# ------------------------------------------------------------------

incidence <- read_incidence(
  filename = "MR/incidence/for_slopehunter/overall_incidence_rsid.txt",
  snp_col = "rsid",
  beta_col = "Beta.meta", 
  se_col = "sdE.meta",
  pval_col = "p.meta",
  eaf_col = "Freq.Gwas",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  chr_col = "chr.Onco",
  pos_col = "Position.Onco"
)

## Batch mode: prognosis files are processed in `process_one()` below.
## Only `incidence` is read here; per-file `prognosis` is read inside the loop.

# Setup output directories and logging
out_dir <- file.path("MR", "slopehunter", "slopehunter_out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir3, "error_log.txt")

# Initialize/clear the error log for this run
if (file.exists(log_file)) file.remove(log_file)
cat("time\tfile\terror\n", file = log_file)

# Initialize summary stats file
summary_file <- file.path(out_dir, "slopehunter_summary_stats.tsv")
if (file.exists(summary_file)) file.remove(summary_file)
cat(paste(c("phenotype", "slope", "slope_se", "ci_lower", "ci_upper", "mixture_pi", "iterations"), collapse = "\t"), "\n", file = summary_file)

# Quick checks: PLINK and EUR reference panel must exist
plink_exec_check <- Sys.which("plink")
if (is.na(plink_exec_check) || plink_exec_check == "") {
  stop("plink not found on PATH. Please install PLINK or add it to PATH.")
}
ref_files <- file.path("MR", "reference", c("EUR.bed", "EUR.bim", "EUR.fam"))
missing_ref <- ref_files[!file.exists(ref_files)]
if (length(missing_ref) > 0) {
  stop("Missing reference files: ", paste(missing_ref, collapse = ", "),
       "; place EUR.bed/EUR.bim/EUR.fam in MR/reference/")
}

# Directory containing prognosis GWAS files
prognosis_dir <- file.path("MR", "outcome", "survival_2021_gwas_results_files_with_rsid_no_chr23")
# Files to exclude
exclude_files <- c(
  "h_ERPRHER2posi_15_years_with_rsid.txt",
  "i_ERPRnega_HER2posi_15_years_with_rsid.txt",
  "j_TNBC_15_years_with_rsid.txt"
)

prognosis_files <- list.files(prognosis_dir, pattern = "\\.txt$", full.names = TRUE)
prognosis_files <- prognosis_files[!basename(prognosis_files) %in% exclude_files]

message("Found ", length(prognosis_files), " prognosis files to process.")

process_one <- function(prog_file) {
  name <- file_path_sans_ext(basename(prog_file))

  tryCatch({
    message("Processing: ", basename(prog_file))

    prognosis <- read_prognosis(
      filename = prog_file,
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

    # Remove rows lacking essential data
    incidence_sub <- incidence[incidence$SH_keep.incidence %in% TRUE, ]
    prognosis_sub <- prognosis[prognosis$SH_keep.prognosis %in% TRUE, ]

    harmonised <- harmonise_effects(
      incidence_dat = incidence_sub,
      prognosis_dat = prognosis_sub,
      by.pos = FALSE
    )
    harmonisation_info <- attr(harmonised, "info")
    print(harmonisation_info)
    harmonised <- as.data.frame(harmonised)

    analysis_data_local <- harmonised[
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
    analysis_data_local <- analysis_data_local[!duplicated(analysis_data_local$SNP), ]

    # PLINK clumping files unique per run
    clump_input_file <- file.path("MR", "slopehunter", paste0("plink_clump_input_", name, ".txt"))
    write.table(
      analysis_data_local[, c("SNP", "PVAL.incidence")],
      file = clump_input_file,
      row.names = FALSE,
      col.names = c("SNP", "P"),
      quote = FALSE,
      sep = "\t"
    )

    plink_exec <- Sys.which("plink")
    if (is.na(plink_exec) || plink_exec == "") {
      plink_exec <- "plink"
    }
    plink_bfile <- file.path("MR", "reference", "EUR")
    plink_out_prefix <- file.path("MR", "slopehunter", paste0("plink_clump_", name))

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

    message("Running PLINK clumping for ", name)
    plink_status <- system2(plink_exec, args = plink_args)
    if (!is.null(plink_status) && plink_status != 0) {
      stop("PLINK clumping failed with exit code: ", plink_status)
    }

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
      stop("Cannot locate SNP column in PLINK .clumped output for: ", name)
    }

    pruned <- analysis_data_local[analysis_data_local$SNP %in% lead_snps, ]

    lambda <- 1e-3
    n_for_mixture <- sum(pruned$PVAL.incidence <= lambda, na.rm = TRUE)
    message("Independent variants entering the mixture model: ", n_for_mixture)

    fit <- hunt(
      dat = pruned,
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
      M = 100,
      seed = 20260813,
      Plot = FALSE,
      show_adjustments = FALSE
    )

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
      tryCatch({ print(fit$plot) }, error = function(e) {
        msg <- paste0(Sys.time(), "\t", name, "\tplot_error\t", conditionMessage(e), "\n")
        cat(msg, file = log_file, append = TRUE)
        message("Plotting failed for ", name, ": ", conditionMessage(e))
      })
    } else {
      tryCatch({ plot(fit) }, error = function(e) {
        msg <- paste0(Sys.time(), "\t", name, "\tplot_error\t", conditionMessage(e), "\n")
        cat(msg, file = log_file, append = TRUE)
        message("Plotting failed for ", name, ": ", conditionMessage(e))
      })
    }
    dev.off()

    adjusted <- SHadj(
      x = fit,
      dat = analysis_data_local,
      snp_col = "SNP",
      xbeta_col = "BETA.incidence",
      xse_col = "SE.incidence",
      ybeta_col = "BETA.prognosis",
      yse_col = "SE.prognosis"
    )

    out_file <- file.path(out_dir, paste0(name, ".slopehunter.tsv.gz"))
    data.table::fwrite(adjusted, file = out_file, sep = "\t", quote = FALSE)

    TRUE
  }, error = function(e) {
    msg <- paste0(Sys.time(), "\t", basename(prog_file), "\t", conditionMessage(e), "\n")
    cat(msg, file = log_file, append = TRUE)
    message("Error processing ", basename(prog_file), ": ", conditionMessage(e))
    FALSE
  })
}

# Iterate over files
results <- vapply(prognosis_files, process_one, logical(1))
message(sum(results), " succeeded, ", sum(!results), " failed. See ", log_file, " for errors.")