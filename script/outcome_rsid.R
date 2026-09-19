#!/usr/bin/env Rscript

# Add dbSNP144/GRCh37 rsIDs to every survival GWAS text file.
#
# Usage (paths are optional):
#   Rscript MR/script/outcome_rsid.R [input_dir] [output_dir] [--overwrite]

required_packages <- c(
  "data.table",
  "GenomicRanges",
  "IRanges",
  "SNPlocs.Hsapiens.dbSNP144.GRCh37"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
})

get_project_dir <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) == 1L) {
    script_file <- sub("^--file=", "", file_argument)
    return(dirname(dirname(dirname(normalizePath(script_file)))))
  }
  normalizePath(getwd())
}

arguments <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% arguments
arguments <- arguments[arguments != "--overwrite"]
if (length(arguments) > 2L) {
  stop(
    "Usage: Rscript MR/script/outcome_rsid.R ",
    "[input_dir] [output_dir] [--overwrite]",
    call. = FALSE
  )
}

project_dir <- get_project_dir()
input_dir <- if (length(arguments) >= 1L) {
  arguments[[1L]]
} else {
  file.path(project_dir, "MR", "outcome", "survival_2021_gwas_results_files")
}
output_dir <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  file.path(project_dir, "MR", "outcome", "survival_2021_gwas_results_files_with_rsid")
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, call. = FALSE)
}
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}

input_files <- sort(list.files(
  input_dir,
  pattern = "\\.txt$",
  full.names = TRUE,
  ignore.case = TRUE
))
if (length(input_files) == 0L) {
  stop("No .txt files found in: ", input_dir, call. = FALSE)
}

snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
required_columns <- c("Chromosome", "Position")

make_snp_map <- function(outcome) {
  chromosome <- suppressWarnings(as.integer(outcome[["Chromosome"]]))
  position <- suppressWarnings(as.integer(outcome[["Position"]]))
  valid <- !is.na(chromosome) & chromosome %between% c(1L, 22L) &
    !is.na(position) & position > 0L

  coordinates <- unique(data.table(
    Chromosome = chromosome[valid],
    Position = position[valid]
  ))
  if (nrow(coordinates) == 0L) {
    return(data.table(
      Chromosome = integer(),
      Position = integer(),
      rsid = character()
    ))
  }

  # SNPlocs.Hsapiens.dbSNP144.GRCh37 uses "1", ..., "22", not
  # "ch1", ..., "ch22" or UCSC-style "chr1", ..., "chr22".
  query_ranges <- GRanges(
    seqnames = as.character(coordinates[["Chromosome"]]),
    ranges = IRanges(start = coordinates[["Position"]], width = 1L)
  )
  known_snps <- snpsByOverlaps(snps, query_ranges)
  if (length(known_snps) == 0L) {
    return(data.table(
      Chromosome = integer(),
      Position = integer(),
      rsid = character()
    ))
  }

  snp_map <- as.data.table(known_snps)[, .(
    Chromosome = as.integer(as.character(seqnames)),
    Position = as.integer(pos),
    RefSNP_id = as.character(RefSNP_id)
  )]
  snp_map[, rsid := ifelse(
    startsWith(RefSNP_id, "rs"),
    RefSNP_id,
    paste0("rs", RefSNP_id)
  )]

  duplicate_positions <- snp_map[, .N, by = .(Chromosome, Position)][N > 1L]
  if (nrow(duplicate_positions) > 0L) {
    message(
      "  Note: ", nrow(duplicate_positions),
      " position(s) have multiple dbSNP IDs; retaining the first ID."
    )
  }

  # dbSNP can assign more than one RefSNP ID to the same position. Keep the
  # first database result so that joining never duplicates outcome rows.
  unique(snp_map[, .(Chromosome, Position, rsid)], by = c("Chromosome", "Position"))
}

process_file <- function(input_file) {
  output_name <- sub("\\.txt$", "_with_rsid.txt", basename(input_file), ignore.case = TRUE)
  output_file <- file.path(output_dir, output_name)

  if (file.exists(output_file) && !overwrite) {
    message("Skipping existing file: ", output_file)
    return(invisible(FALSE))
  }

  message("Reading: ", input_file)
  outcome <- fread(input_file, nThread = getDTthreads(), showProgress = TRUE)
  missing_columns <- setdiff(required_columns, names(outcome))
  if (length(missing_columns) > 0L) {
    stop("missing column(s): ", paste(missing_columns, collapse = ", "))
  }

  # Normalize join columns once. Integer positions are sufficient for GRCh37.
  outcome[, Chromosome := suppressWarnings(as.integer(Chromosome))]
  outcome[, Position := suppressWarnings(as.integer(Position))]
  outcome[, rsid := NA_character_]

  message("Matching rsIDs against dbSNP144/GRCh37...")
  snp_map <- make_snp_map(outcome)

  # An update join is memory-efficient and preserves the original row order.
  outcome[snp_map, on = .(Chromosome, Position), rsid := i.rsid]
  setcolorder(outcome, c("rsid", setdiff(names(outcome), "rsid")))

  matched_count <- sum(!is.na(outcome[["rsid"]]))
  matched_percent <- if (nrow(outcome) == 0L) {
    0
  } else {
    100 * matched_count / nrow(outcome)
  }
  message(sprintf(
    "Matched %s of %s rows (%.2f%%).",
    format(matched_count, big.mark = ",", scientific = FALSE),
    format(nrow(outcome), big.mark = ",", scientific = FALSE),
    matched_percent
  ))

  temporary_file <- tempfile(
    pattern = paste0(".", output_name, "."),
    tmpdir = output_dir
  )
  on.exit(unlink(temporary_file), add = TRUE)
  fwrite(
    outcome,
    temporary_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    nThread = getDTthreads()
  )

  # The temporary file is in the destination directory, so this rename is an
  # atomic replacement on the Linux filesystem used for these analyses.
  if (!file.rename(temporary_file, output_file)) {
    stop("could not move temporary output to: ", output_file)
  }
  message("Saved: ", output_file)
  invisible(TRUE)
}

message("=== Matching outcome rsIDs ===")
message("Input files: ", length(input_files))
message("Output directory: ", output_dir)

failed_files <- character()
for (input_file in input_files) {
  message("\n--- ", basename(input_file), " ---")
  tryCatch(
    process_file(input_file),
    error = function(error) {
      failed_files <<- c(failed_files, input_file)
      message("ERROR: ", conditionMessage(error))
    }
  )
  gc(verbose = FALSE)
}

if (length(failed_files) > 0L) {
  message(
    "\nFinished with ", length(failed_files), " failed file(s):\n  ",
    paste(basename(failed_files), collapse = "\n  ")
  )
  quit(save = "no", status = 1L)
}

message("\n=== Completed all ", length(input_files), " outcome files ===")
