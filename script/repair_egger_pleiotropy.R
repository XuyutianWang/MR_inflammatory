#!/usr/bin/env Rscript

# Regenerate MR-Egger intercept results only, from existing harmonised files.
# Usage:
# Rscript script/repair_egger_pleiotropy.R <harmonised_dir> <mr_output_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript repair_egger_pleiotropy.R <harmonised_dir> <mr_output_dir>", call. = FALSE)
}
input_dir <- args[1]
output_dir <- args[2]
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pair_from_filename <- function(path) {
  tag <- sub("_harmonised\\.csv$", "", basename(path))
  parts <- strsplit(tag, "__VS__", fixed = TRUE)[[1]]
  list(exposure_id = parts[1], outcome_id = if (length(parts) > 1) parts[2] else NA_character_)
}

files <- list.files(input_dir, pattern = "_harmonised\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No harmonised files found in: ", input_dir, call. = FALSE)

results <- list()
status <- vector("list", length(files))
for (i in seq_along(files)) {
  file <- files[i]
  pair <- pair_from_filename(file)
  base_status <- data.table(filename = basename(file), exposure_id_file = pair$exposure_id,
                            outcome_id_file = pair$outcome_id, nsnp_retained = NA_integer_,
                            status = NA_character_, error = NA_character_)
  result <- tryCatch({
    dat <- fread(file)
    if (!"mr_keep" %in% names(dat)) stop("missing mr_keep column")
    clean <- dat[mr_keep %in% TRUE]
    base_status$nsnp_retained <- nrow(clean)
    if (nrow(clean) < 3) {
      base_status$status <- "not_estimable_fewer_than_3_instruments"
      NULL
    } else {
      pleio <- as.data.table(TwoSampleMR::mr_pleiotropy_test(clean))
      pleio[, `:=`(filename = basename(file), exposure_id_file = pair$exposure_id,
                   outcome_id_file = pair$outcome_id, nsnp_retained = nrow(clean))]
      base_status$status <- "tested"
      pleio
    }
  }, error = function(e) {
    base_status$status <<- "failed"
    base_status$error <<- conditionMessage(e)
    NULL
  })
  status[[i]] <- base_status
  if (!is.null(result) && nrow(result)) results[[length(results) + 1]] <- result
}

result_dt <- if (length(results)) rbindlist(results, fill = TRUE) else data.table()
status_dt <- rbindlist(status, fill = TRUE)
fwrite(result_dt, file.path(output_dir, "MR_Egger_pleiotropy.csv"))
fwrite(status_dt, file.path(output_dir, "MR_Egger_pleiotropy_status.csv"))
message("MR-Egger: ", sum(status_dt$status == "tested"), " tested; ",
        sum(grepl("not_estimable", status_dt$status)), " not estimable; ",
        sum(status_dt$status == "failed"), " failed.")
