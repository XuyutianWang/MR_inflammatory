library(data.table)
library(dplyr)

input_dir  <- "MR/incidence/rsid"
output_dir <- "MR/incidence/rsid_no_chr23"

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
  
  # 检查是否存在 chr 列，如果不存在则使用 chr.iCOGs
  if ("chr" %in% colnames(df)) {
    message(" -> using 'chr' column")
    df_filtered <- df %>% 
      filter(!(chr %in% c(23, "23", "X", "chr23", "chrX")))
  } else if ("chr.iCOGs" %in% colnames(df)) {
    message(" -> using 'chr.iCOGs' column")
    df_filtered <- df %>% 
      filter(!(chr.iCOGs %in% c(23, "23", "X", "chr23", "chrX")))
  } else {
    # 如果两列都不存在，发出警告并跳过该文件
    warning(paste0("No chromosome column found in ", file_name, ", skipping..."))
    next
  }
  
  row_after <- nrow(df_filtered)
  removed_rows <- row_before - row_after
  
  output_file <- file.path(output_dir, file_name)
  
  fwrite(df_filtered, output_file, sep = "\t", quote = FALSE, nThread = 4)
  
  message(sprintf(" -> removed %d rows, kept %d rows\n", 
                  removed_rows, row_after))
}

message("=== done: ", output_dir, " ===")