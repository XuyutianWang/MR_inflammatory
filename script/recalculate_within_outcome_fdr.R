#!/usr/bin/env Rscript

# Recalculate BH-FDR within each outcome from completed MR primary results.
# This does NOT re-run MR, harmonisation, or sensitivity analyses.
# Usage: Rscript MR/script/recalculate_within_outcome_fdr.R <mr_output_dir>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript MR/script/recalculate_within_outcome_fdr.R <mr_output_dir>", call. = FALSE)
}

results_dir <- args[1]
input_file <- file.path(results_dir, "MR_primary_results.csv")
if (!file.exists(input_file)) stop("Cannot find: ", input_file, call. = FALSE)

dt <- fread(input_file)
required <- c("pval", "outcome_id_file")
if (!all(required %in% names(dt))) {
  stop("MR_primary_results.csv must contain pval and outcome_id_file.", call. = FALSE)
}

# Preserve prior global q-values if present; otherwise calculate them.
if (!"qval_bh_global" %in% names(dt)) dt[, qval_bh_global := p.adjust(pval, method = "BH")]
dt[, fdr_significant_global := qval_bh_global < 0.05]

# Correct only across the successfully tested proteins for each outcome/subtype.
dt[, qval_bh_within_outcome := p.adjust(pval, method = "BH"), by = outcome_id_file]
dt[, fdr_significant_within_outcome := qval_bh_within_outcome < 0.05]

# Compatibility fields used by plot_mr_results.R now represent within-outcome FDR.
dt[, `:=`(
  qval_bh = qval_bh_within_outcome,
  fdr_significant = fdr_significant_within_outcome
)]

fwrite(dt, file.path(results_dir, "MR_primary_results_within_outcome_FDR.csv"))
fwrite(dt[fdr_significant_within_outcome == TRUE],
       file.path(results_dir, "MR_primary_within_outcome_FDR_significant.csv"))

summary_dt <- dt[, .(
  n_tested = .N,
  n_within_outcome_FDR_significant = sum(fdr_significant_within_outcome),
  min_p = min(pval),
  min_q_within_outcome = min(qval_bh_within_outcome)
), by = outcome_id_file][order(outcome_id_file)]
fwrite(summary_dt, file.path(results_dir, "MR_within_outcome_FDR_summary.csv"))

message("Completed FDR recalculation without rerunning MR. Within-outcome FDR hits: ",
        sum(dt$fdr_significant_within_outcome))
