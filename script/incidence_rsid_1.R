#!/usr/bin/env Rscript

library(data.table)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

input_files <- c(
  "MR/incidence/overall_incidence.txt"
)

chunk_size <- 250000L

normalize_chr <- function(x) {
  gsub("^chr", "", as.character(x), ignore.case = TRUE)
}

pick_first <- function(candidates, available) {
  hit <- intersect(candidates, available)
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  hit[[1L]]
}

message("Loading GRCh37 dbSNP144 reference ...")
snps_ref <- SNPlocs.Hsapiens.dbSNP144.GRCh37
ref_style <- seqlevelsStyle(snps_ref)[1]

message(sprintf("Processing %d subtype files for rsID mapping...\n", length(input_files)))

process_file <- function(file, output_file, snps_ref, ref_style, chunk_size) {
  header_dt <- fread(file, nrows = 0, showProgress = FALSE)
  all_cols <- names(header_dt)

  chr_col <- pick_first(grep("^chr(\\.|$)", all_cols, value = TRUE), all_cols)
  pos_col <- pick_first(grep("^Position(\\.|$)", all_cols, value = TRUE), all_cols)
  id_col <- pick_first(c("SNP", "var_name", grep("^SNP(\\.|$)", all_cols, value = TRUE)), all_cols)

  if (is.na(chr_col) || is.na(pos_col)) {
    stop(sprintf("Could not detect chr/pos columns in %s", file))
  }
  if (is.na(id_col)) {
    stop(sprintf("Could not detect an ID column in %s", file))
  }

  if (file.exists(output_file)) {
    file.remove(output_file)
  }

  con <- file(file, open = "r")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) == 0L) {
    message("   -> File is empty, skipping.")
    return(invisible(NULL))
  }

  first_write <- TRUE
  total_rows <- 0L
  total_matched <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    chunk <- fread(
      text = c(header_line, lines),
      showProgress = FALSE,
      nThread = 4
    )

    setDT(chunk)
    chunk[, rsid := as.character(get(id_col))]
    chunk[, `:=`(
      chr_tmp = normalize_chr(get(chr_col)),
      pos_tmp = as.integer(get(pos_col))
    )]

    gr <- GRanges(
      seqnames = chunk$chr_tmp,
      ranges = IRanges(start = chunk$pos_tmp, width = 1)
    )
    seqlevelsStyle(gr) <- ref_style

    locs <- snpsByOverlaps(snps_ref, gr)
    if (length(locs) > 0L) {
      snp_map <- as.data.table(locs)[, .(
        chr = normalize_chr(seqnames),
        pos = as.integer(pos),
        rsid_ref = as.character(RefSNP_id)
      )]
      snp_map <- unique(snp_map, by = c("chr", "pos"))

      chunk[snp_map, rsid := fifelse(!is.na(i.rsid_ref) & i.rsid_ref != "", i.rsid_ref, rsid), on = .(chr_tmp = chr, pos_tmp = pos)]
    }

    chunk[, c("chr_tmp", "pos_tmp") := NULL]
    setcolorder(chunk, c("rsid", setdiff(names(chunk), "rsid")))

    fwrite(
      chunk,
      output_file,
      sep = "\t",
      quote = FALSE,
      append = !first_write,
      col.names = first_write,
      nThread = 4
    )

    matched_num <- sum(grepl("^rs", chunk$rsid))
    total_rows <- total_rows + nrow(chunk)
    total_matched <- total_matched + matched_num
    first_write <- FALSE

    message(sprintf("   -> Processed %d rows so far, matched %d rsIDs (%.2f%%).", 
                    total_rows, total_matched, (total_matched / total_rows) * 100))
  }

  message(sprintf(
    "   -> Complete! Total %d sites, %d successfully matched to rsid (match rate %.2f%%).\n",
    total_rows,
    total_matched,
    ifelse(total_rows > 0L, (total_matched / total_rows) * 100, 0)
  ))
}

for (i in seq_along(input_files)) {
  file <- input_files[i]
  output_file <- sub("\\.txt$", "_rsid.txt", file)

  message(sprintf("[%d/%d] Processing: %s ...", i, length(input_files), file))
  process_file(file, output_file, snps_ref, ref_style, chunk_size)
}

message("=== All subtype data rsID mapping completed! ===")