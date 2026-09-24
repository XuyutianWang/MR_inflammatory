args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop(
    "Usage: Rscript script/filter_cis_lead_snps.R [labelled_dir] [cis_output_dir]",
    call. = FALSE
  )
}

source_dir <- if (length(args) >= 1L) args[[1L]] else "exposure_all/UKB_1_22_all_rsid_cis_trans"
target_dir <- if (length(args) >= 2L) args[[2L]] else "exposure_all/UKB_1_22_all_rsid_cis"

if (!dir.exists(source_dir)) {
  stop("Labelled input directory does not exist: ", source_dir, call. = FALSE)
}

dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- list.files(source_dir, pattern = "\\.txt$", full.names = TRUE)
for (source_file in source_files) {
  dat <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"qtl_type" %in% names(dat)) {
    stop("Missing qtl_type column: ", source_file)
  }

  cis_dat <- dat[tolower(trimws(dat$qtl_type)) == "cis", , drop = FALSE]
  write.table(
    cis_dat,
    file = file.path(target_dir, basename(source_file)),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
