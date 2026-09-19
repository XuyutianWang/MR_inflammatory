#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

user_dir <- Sys.getenv("R_LIBS_USER")
if (dir.exists(user_dir)) .libPaths(c(user_dir, .libPaths()))

input_dir  <- "MR/exposure_all/UKB_1_22_all_pfiltered_5108"
output_dir <- "MR/exposure_all/UKB_1_22_all_rsid"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste0("Setup dir: ", output_dir))
}

message("Loading GRCh38 dbSNP155 ...")
snps_38 <- SNPlocs.Hsapiens.dbSNP155.GRCh38

gz_files <- list.files(input_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)

message(sprintf("Processing %d files for rsid mapping...\n", length(gz_files)))

for (i in seq_along(gz_files)) {
  file <- gz_files[i]
  
  rel_path <- sub(paste0(input_dir, "/"), "", file, fixed = TRUE)
  rel_path <- sub(input_dir, "", rel_path, fixed = TRUE)
  rel_path <- gsub("^/", "", rel_path)
  
  output_file <- file.path(output_dir, rel_path)
  
  out_file_dir <- dirname(output_file)
  if (!dir.exists(out_file_dir)) {
    dir.create(out_file_dir, recursive = TRUE)
  }
  
  message(sprintf("[%d/%d] Processing: %s ...", i, length(gz_files), rel_path))
  
  df <- fread(file, nThread = 4)
  
  if (nrow(df) == 0) {
    message("   -> File is empty, skipping.")

    next
  }
  
gr <- GRanges(
  seqnames = as.character(df$CHROM),
  ranges = IRanges(start = df$GENPOS, width = 1)
)
  
  locs <- snpsByOverlaps(snps_38, gr)
  
  if (length(locs) > 0) {
    snp_map <- as.data.table(locs) %>%
      mutate(
        CHROM = as.integer(as.character(seqnames)),
        GENPOS = as.integer(pos),
        rsid = RefSNP_id
      ) %>%
      select(CHROM, GENPOS, rsid) %>%
      unique(by = c("CHROM", "GENPOS"))
    
    df_matched <- df %>%
      left_join(snp_map, by = c("CHROM", "GENPOS")) %>%
      mutate(
        rsid = ifelse(!is.na(rsid) & rsid != "", rsid, ID)
      )
  } else {
    df_matched <- df %>% mutate(rsid = ID)
  }
  
  setcolorder(df_matched, c("rsid", names(df)))
  
  fwrite(df_matched, output_file, sep = "\t", quote = FALSE, compress = "gzip", nThread = 4)
  
  matched_num <- sum(grepl("^rs", df_matched$rsid))
  message(sprintf("   -> Complete! Total %d sites, %d successfully matched to rsid (match rate %.2f%%).\n", 
                  nrow(df_matched), matched_num, (matched_num / nrow(df_matched)) * 100))
}

message("=== All exposure data rsid mapping completed! Data saved to: ", output_dir, " ===")