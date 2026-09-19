library(data.table)
library(dplyr)

input_dir  <- "MR/outcome/survival_2021_gwas_results_files_with_rsid"
output_dir <- "MR/outcome/survival_2021_gwas_results_files_with_rsid_no_chr23"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste0("made dir: ", output_dir))
}

txt_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

message(sprintf("found %d files...\n", length(txt_files)))

for (file in txt_files) {
  file_name <- basename(file)
  message(paste0("processing: ", file_name, " ..."))
  
  df <- fread(file, nThread = 4)
  
  row_before <- nrow(df)
  
  df_filtered <- df %>% 
    filter(!(Chromosome %in% c(23, "23", "X", "chr23", "chrX")))
  
  row_after <- nrow(df_filtered)
  removed_rows <- row_before - row_after
  
  output_file <- file.path(output_dir, file_name)
  
  fwrite(df_filtered, output_file, sep = "\t", quote = FALSE, nThread = 4)
  
  message(sprintf(" -> removed %d rows, kept %d rows\n", 
                  removed_rows, row_after))
}

message("=== done: ", output_dir, " ===")