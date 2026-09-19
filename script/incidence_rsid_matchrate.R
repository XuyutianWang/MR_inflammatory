library(data.table)
library(dplyr)

input_dir <- "MR/incidence/rsid_no_chr23"
output_file <- "incidence_rsid_matchrate.txt"

# 获取所有txt文件
txt_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

message(sprintf("found %d files...\n", length(txt_files)))

# 创建结果数据框
results <- data.frame(
  file_name = character(),
  total_rows = integer(),
  rsid_not_empty = integer(),
  rsid_start_with_rs = integer(),
  real_rsid_percentage = numeric(),
  stringsAsFactors = FALSE
)

for (file in txt_files) {
  file_name <- basename(file)
  message(paste0("processing: ", file_name, " ..."))
  
  # 读取文件
  df <- fread(file, nThread = 4)
  
  total_rows <- nrow(df)
  
  # 检查是否存在rsid列（假设第一列是rsid）
  if (!"rsid" %in% colnames(df)) {
    warning(paste0("No 'rsid' column found in ", file_name, ", skipping..."))
    next
  }
  
  # 计算rsid不为空的行数（排除NA、空字符串）
  rsid_not_empty <- sum(!is.na(df$rsid) & df$rsid != "" & !grepl("^\\s*$", df$rsid))
  
  # 计算以"rs"开头的真正rsid数量
  rsid_start_with_rs <- sum(!is.na(df$rsid) & grepl("^rs", df$rsid))
  
  # 计算真正rsid的比例（以总行数为分母）
  real_rsid_pct <- (rsid_start_with_rs / total_rows) * 100
  
  # 保存结果
  results <- rbind(results, data.frame(
    file_name = file_name,
    total_rows = total_rows,
    rsid_not_empty = rsid_not_empty,
    rsid_start_with_rs = rsid_start_with_rs,
    real_rsid_percentage = round(real_rsid_pct, 2)
  ))
  
  message(sprintf(" -> total: %d, non-empty: %d, real rsid (starting with 'rs'): %d, real rsid percentage: %.2f%%\n", 
                  total_rows, rsid_not_empty, rsid_start_with_rs, real_rsid_pct))
}

# 按文件名字母顺序排序
results <- results[order(results$file_name), ]

# 输出汇总结果
message("\n=== Summary ===")
print(results)

# 保存到文件
fwrite(results, output_file, sep = "\t", quote = FALSE)
message(paste0("\nResults saved to: ", output_file))

# 额外显示统计信息
message("\n=== Statistics ===")
message(sprintf("Total files processed: %d", nrow(results)))
message(sprintf("Mean real rsid percentage: %.2f%%", mean(results$real_rsid_percentage)))
message(sprintf("Median real rsid percentage: %.2f%%", median(results$real_rsid_percentage)))
message(sprintf("Min real rsid percentage: %.2f%%", min(results$real_rsid_percentage)))
message(sqrt("Max real rsid percentage: %.2f%%", max(results$real_rsid_percentage)))