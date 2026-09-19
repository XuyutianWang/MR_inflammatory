source_dir <- "exposure_all/lead_snps_cis_trans"
target_dir <- "exposure_all/lead_snps_cis"

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
