#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(TwoSampleMR)

input_dir  <- "MR/harmonised_results"
output_dir <- "MR/mr_results_summary"
pleio_file <- file.path(output_dir, "ALL_MR_Egger_Pleiotropy.csv")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
message(sprintf("找到 %d 个文件，开始补跑 MR-Egger 多效性检验...\n", length(csv_files)))

# 1. 显式初始化创建带有表头的空 CSV 文件
init_pleio <- data.table(
  id.exposure = character(),
  id.outcome = character(),
  outcome = character(),
  exposure = character(),
  egger_intercept = numeric(),
  se = numeric(),
  pval = numeric(),
  filename = character()
)
fwrite(init_pleio, pleio_file)

success_count <- 0
fail_count    <- 0

# 2. 循环专门运行多效性检验
for (i in seq_along(csv_files)) {
  file <- csv_files[i]
  file_name <- basename(file)
  
  dat <- fread(file)
  dat_clean <- subset(dat, mr_keep == TRUE)
  
  # 必须 >= 3 个 SNP 才能算多效性截距
  if (nrow(dat_clean) < 3) {
    fail_count <- fail_count + 1
    next
  }
  
  # 显式使用 TwoSampleMR::mr_pleio_test
  res_pleio <- tryCatch({
    TwoSampleMR::mr_pleio_test(dat_clean)
  }, error = function(e) {
    # 如果 TwoSampleMR::mr_pleio_test 仍找不到，降级调用 mr_egger 拟合提取截距
    tryCatch({
      fit <- lm(b_out ~ b_exp, weights = 1/se_out^2, data = data.frame(b_out = dat_clean$beta.outcome, b_exp = dat_clean$beta.exposure, se_out = dat_clean$se.outcome))
      coef_fit <- summary(fit)$coefficients
      data.frame(
        id.exposure = dat_clean$id.exposure[1],
        id.outcome = dat_clean$id.outcome[1],
        outcome = dat_clean$outcome[1],
        exposure = dat_clean$exposure[1],
        egger_intercept = coef_fit[1, 1],
        se = coef_fit[1, 2],
        pval = coef_fit[1, 4]
      )
    }, error = function(e2) NULL)
  })
  
  if (!is.null(res_pleio) && nrow(res_pleio) > 0) {
    res_pleio$filename <- file_name
    fwrite(res_pleio, pleio_file, append = TRUE)
    success_count <- success_count + 1
  } else {
    fail_count <- fail_count + 1
  }
}


df <- fread("MR/mr_results_summary/ALL_MR_Egger_Pleiotropy.csv")

df_filtered <- df[df$pval < 0.05, ]

output_file <- "MR/sensitivity/Pleiotropy.csv"
fwrite(df_filtered, output_file)


message(sprintf("\n 补跑完成！成功提取 %d 组多效性结果，跳过/不满足条件 %d 组。", success_count, fail_count))
message("多效性结果已保存至: ", pleio_file)