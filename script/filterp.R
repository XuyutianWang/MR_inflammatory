#!/usr/bin/env Rscript

# Filter UKB-PPP exposure files at the genome-wide significance threshold.
# LOG10P is -log10(P), so P < 5e-8 is equivalent to LOG10P > -log10(5e-8).
#
# Usage (paths are optional):
#   Rscript MR/script/filterp.R [input_dir] [output_dir] [--overwrite]

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Missing required R package: data.table", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))

gzip_executable <- Sys.which("gzip")
if (!nzchar(gzip_executable)) {
  stop("The gzip executable was not found on PATH.", call. = FALSE)
}

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
    "Usage: Rscript MR/script/filterp.R ",
    "[input_dir] [output_dir] [--overwrite]",
    call. = FALSE
  )
}

project_dir <- get_project_dir()
input_dir <- if (length(arguments) >= 1L) {
  arguments[[1L]]
} else {
  file.path(project_dir, "MR", "exposure_all", "UKB_1_22_all")
}
output_dir <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  file.path(project_dir, "MR", "exposure_all", "UKB_1_22_all_pfiltered")
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, call. = FALSE)
}
input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)

if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

if (identical(input_dir, output_dir)) {
  stop("Input and output directories must be different.", call. = FALSE)
}
if (startsWith(output_dir, paste0(input_dir, "/"))) {
  stop("Output directory cannot be inside the recursively scanned input directory.", call. = FALSE)
}

gz_files <- sort(list.files(
  input_dir,
  pattern = "\\.gz$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
))
if (length(gz_files) == 0L) {
  stop("No .gz files found in: ", input_dir, call. = FALSE)
}

p_threshold <- -log10(5e-8)

relative_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (!startsWith(path, prefix)) {
    stop("File is not under the input directory: ", path)
  }
  substring(path, nchar(prefix) + 1L)
}

process_file <- function(input_file, index, total) {
  rel_path <- relative_path(input_file, input_dir)
  output_file <- file.path(output_dir, rel_path)

  if (file.exists(output_file) && !overwrite) {
    message(sprintf("[%d/%d] Skipping existing file: %s", index, total, rel_path))
    return(invisible(FALSE))
  }

  output_parent <- dirname(output_file)
  if (!dir.exists(output_parent) && !dir.create(output_parent, recursive = TRUE)) {
    stop("could not create output directory: ", output_parent)
  }

  message(sprintf("[%d/%d] Reading: %s", index, total, rel_path))
  # Reading a .gz path directly makes data.table depend on the optional
  # R.utils package. Streaming through gzip keeps the script self-contained.
  read_command <- paste(
    shQuote(gzip_executable),
    "-cd --",
    shQuote(input_file)
  )
  exposure <- fread(
    cmd = read_command,
    nThread = getDTthreads(),
    showProgress = TRUE
  )
  row_before <- nrow(exposure)

  p_columns <- intersect(c("LOG10P", "LOG10_P"), names(exposure))
  if (length(p_columns) == 0L) {
    stop("LOG10P or LOG10_P column not found")
  }
  p_column <- p_columns[[1L]]
  log10p <- suppressWarnings(as.numeric(exposure[[p_column]]))
  invalid_count <- sum(is.na(log10p) & !is.na(exposure[[p_column]]))
  if (invalid_count > 0L) {
    warning(
      sprintf(
        "%s contains %s non-numeric value(s); these rows will be excluded.",
        p_column,
        format(invalid_count, big.mark = ",", scientific = FALSE)
      ),
      call. = FALSE
    )
  }

  exposure_filtered <- exposure[!is.na(log10p) & log10p > p_threshold]
  row_after <- nrow(exposure_filtered)

  # A temporary gzip in the destination directory prevents a failed run from
  # leaving a truncated file that would be mistaken for completed output.
  temporary_file <- tempfile(
    pattern = paste0(".", basename(output_file), "."),
    tmpdir = output_parent,
    fileext = ".gz"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  fwrite(
    exposure_filtered,
    temporary_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    nThread = getDTthreads()
  )

  # The temporary and final files share a filesystem, making replacement
  # atomic on the Linux environment used for this analysis.
  if (!file.rename(temporary_file, output_file)) {
    stop("could not move temporary output to: ", output_file)
  }

  retention_rate <- if (row_before == 0L) 0 else 100 * row_after / row_before
  message(sprintf(
    "  Saved %s: %s of %s rows retained (%.4f%%).",
    rel_path,
    format(row_after, big.mark = ",", scientific = FALSE),
    format(row_before, big.mark = ",", scientific = FALSE),
    retention_rate
  ))
  invisible(TRUE)
}

message("=== Filtering exposure data ===")
message("Files: ", length(gz_files))
message(sprintf("Threshold: LOG10P > %.10f (P < 5e-8)", p_threshold))
message("Output directory: ", output_dir)

failed_files <- character()
for (index in seq_along(gz_files)) {
  input_file <- gz_files[[index]]
  tryCatch(
    process_file(input_file, index, length(gz_files)),
    error = function(error) {
      failed_files <<- c(failed_files, input_file)
      message("  ERROR [", basename(input_file), "]: ", conditionMessage(error))
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

message("\n=== Filtering completed for all ", length(gz_files), " files ===")
