#!/usr/bin/env Rscript

# Batch univariable MR for one harmonised-result directory.
# Usage: Rscript MR/script/mr.R <harmonised_input_dir> <output_dir> [--overwrite] [--run-steiger]

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript MR/script/mr.R <harmonised_input_dir> <output_dir> [--overwrite] [--run-steiger]", call. = FALSE)
}
input_dir <- args[1]
output_dir <- args[2]
flags <- args[-c(1, 2)]
if (!all(flags %in% c("--overwrite", "--run-steiger")) || anyDuplicated(flags)) {
  stop("Optional flags are --overwrite and --run-steiger.", call. = FALSE)
}
overwrite <- "--overwrite" %in% flags
run_steiger <- "--run-steiger" %in% flags

if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir, call. = FALSE)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE, no.. = TRUE)) > 0 && !overwrite) {
  stop("Output directory is not empty: ", output_dir, ". Choose a new directory or pass --overwrite.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(input_dir, pattern = "_harmonised\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(csv_files)) stop("No *_harmonised.csv files found in: ", input_dir, call. = FALSE)

pair_from_filename <- function(path) {
  tag <- sub("_harmonised\\.csv$", "", basename(path))
  parts <- strsplit(tag, "__VS__", fixed = TRUE)[[1]]
  list(exposure_id = parts[1], outcome_id = if (length(parts) > 1) parts[2] else NA_character_)
}

add_pair_metadata <- function(result, file) {
  if (is.null(result) || !nrow(result)) return(NULL)
  pair <- pair_from_filename(file)
  result <- as.data.table(result)
  result[, `:=`(filename = basename(file), exposure_id_file = pair$exposure_id,
                outcome_id_file = pair$outcome_id)]
  result
}

main_results <- list(); primary_results <- list(); heterogeneity_results <- list()
pleiotropy_results <- list(); steiger_results <- list(); loo_results <- list(); failures <- list()

bind_or_empty <- function(items) {
  items <- Filter(Negate(is.null), items)
  if (!length(items)) return(data.table())
  rbindlist(items, fill = TRUE)
}

write_checkpoint <- function(i) {
  # These files are partial: BH-FDR correction is valid only after completion.
  fwrite(bind_or_empty(main_results), file.path(output_dir, "MR_all_method_results.partial.csv"))
  fwrite(bind_or_empty(primary_results), file.path(output_dir, "MR_primary_results.partial.csv"))
  fwrite(bind_or_empty(failures), file.path(output_dir, "MR_failures.partial.csv"))
  writeLines(sprintf("Processed %d of %d files; run is still in progress.", i, length(csv_files)),
             file.path(output_dir, "MR_run_status.txt"))
}

writeLines(sprintf("Started: 0 of %d files processed.", length(csv_files)), file.path(output_dir, "MR_run_status.txt"))
message("Found ", length(csv_files), " harmonised files.")

for (i in seq_along(csv_files)) {
  file <- csv_files[i]
  message(sprintf("[%d/%d] %s", i, length(csv_files), basename(file)))
  tryCatch({
    dat <- fread(file)
    if (!"mr_keep" %in% names(dat)) stop("Missing mr_keep column")
    dat_clean <- dat[mr_keep %in% TRUE]
    nsnp <- nrow(dat_clean)
    if (!nsnp) {
      failures[[length(failures) + 1]] <- data.table(filename = basename(file), reason = "no_valid_instruments", nsnp_retained = 0L)
      next
    }

    # Package defaults are used; primary result is Wald ratio for one SNP and
    # IVW for multiple SNPs, selected below with subset_on_method().
    res <- mr(dat_clean)
    if (is.null(res) || !nrow(res)) {
      failures[[length(failures) + 1]] <- data.table(filename = basename(file), reason = "mr_returned_no_result", nsnp_retained = nsnp)
      next
    }
    res <- add_pair_metadata(res, file)
    res[, `:=`(exp_b = exp(b), exp_lci95 = exp(b - 1.96 * se), exp_uci95 = exp(b + 1.96 * se))]
    main_results[[length(main_results) + 1]] <- res
    primary <- subset_on_method(res)
    if (nrow(primary)) primary_results[[length(primary_results) + 1]] <- as.data.table(primary)

    if (nsnp >= 2) heterogeneity_results[[length(heterogeneity_results) + 1]] <- add_pair_metadata(tryCatch(mr_heterogeneity(dat_clean), error = function(e) NULL), file)
    if (nsnp >= 3) {
      pleiotropy_results[[length(pleiotropy_results) + 1]] <- add_pair_metadata(tryCatch(mr_pleio_test(dat_clean), error = function(e) NULL), file)
      loo_results[[length(loo_results) + 1]] <- add_pair_metadata(tryCatch(mr_leaveoneout(dat_clean), error = function(e) NULL), file)
    }
    # This is opt-in: survival GWAS summary data often do not permit a valid
    # Steiger correlation calculation.
    if (run_steiger) steiger_results[[length(steiger_results) + 1]] <- add_pair_metadata(tryCatch(directionality_test(dat_clean), error = function(e) NULL), file)
  }, error = function(e) {
    failures[[length(failures) + 1]] <<- data.table(filename = basename(file), reason = conditionMessage(e), nsnp_retained = NA_integer_)
  })
  if (i %% 100 == 0 || i == length(csv_files)) write_checkpoint(i)
}

main_dt <- bind_or_empty(main_results)
primary_dt <- bind_or_empty(primary_results)
if (!nrow(primary_dt)) stop("No primary MR results were generated; inspect MR_failures.csv.", call. = FALSE)

# BH correction applies only to this one outcome family. Run incidence, raw
# survival and Slope-Hunter survival into separate output directories.
primary_dt[, qval_bh := p.adjust(pval, method = "BH")]
primary_dt[, fdr_significant := qval_bh < 0.05]

fwrite(main_dt, file.path(output_dir, "MR_all_method_results.csv"))
fwrite(primary_dt, file.path(output_dir, "MR_primary_results.csv"))
fwrite(primary_dt[fdr_significant == TRUE], file.path(output_dir, "MR_primary_FDR_significant.csv"))
fwrite(bind_or_empty(heterogeneity_results), file.path(output_dir, "MR_heterogeneity.csv"))
fwrite(bind_or_empty(pleiotropy_results), file.path(output_dir, "MR_Egger_pleiotropy.csv"))
fwrite(bind_or_empty(steiger_results), file.path(output_dir, "MR_Steiger.csv"))
fwrite(bind_or_empty(loo_results), file.path(output_dir, "MR_leave_one_out.csv"))
fwrite(bind_or_empty(failures), file.path(output_dir, "MR_failures.csv"))

message("Completed ", nrow(primary_dt), " primary tests; ", sum(primary_dt$fdr_significant),
        " pass BH FDR < 0.05. Results: ", output_dir)
