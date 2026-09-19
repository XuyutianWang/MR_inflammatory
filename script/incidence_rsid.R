#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# 配置文件路径 / Configure file paths
input_files <- c(
  "MR/incidence/Luminal_A.txt",
  "MR/incidence/HER2_Enriched.txt",
  "MR/incidence/Luminal_B_HER2Neg.txt",
  "MR/incidence/Luminal_B.txt",
  "MR/incidence/TNBC.txt",
  "MR/incidence/overall_incidence.txt"
)

message("Loading GRCh37 dbSNP144 reference ...")
snps_ref <- SNPlocs.Hsapiens.dbSNP144.GRCh37
ref_style <- seqlevelsStyle(snps_ref)[1] # 获取参考库的染色体命名风格 (如 "dbSNP" 风格) / Get reference chromosome naming style

message(sprintf("Processing %d subtype files for rsID mapping...\n", length(input_files)))

for (i in seq_along(input_files)) {
  file <- input_files[i]
  
  # 设置输出文件名 (如在原文件名后加 _rsid)
  # Set output filename (e.g., append _rsid to original filename)
  output_file <- sub("\\.txt$", "_rsid.txt", file)
  
  message(sprintf("[%d/%d] Processing: %s ...", i, length(input_files), file))
  
  df <- fread(file, nThread = 4)
  
  if (nrow(df) == 0) {
    message("   -> File is empty, skipping.")
    next
  }
  
  # 构建 GRanges 对象并统一染色体命名风格
  # Build GRanges object and harmonize chromosome naming styles
  gr <- GRanges(
    seqnames = as.character(df$chr),
    ranges = IRanges(start = df$pos, width = 1)
  )
  # 自动匹配参考基因组的 chr 格式 (防止 "10" 和 "ch10" 匹配不上的问题)
  # Automatically match chr formats to prevent mismatch between "10" and "ch10"
  seqlevelsStyle(gr) <- ref_style 
  
  locs <- snpsByOverlaps(snps_ref, gr)
  
  if (length(locs) > 0) {
    # 提取匹配结果 / Extract match results
    snp_map <- as.data.table(locs) %>%
      mutate(
        chr = as.integer(gsub("[^0-9]", "", as.character(seqnames))), # 提取纯数字 chr / Extract purely numeric chr
        pos = as.integer(pos),
        rsid = RefSNP_id
      ) %>%
      select(chr, pos, rsid) %>%
      unique(by = c("chr", "pos"))
    
    # 与原数据合并 / Merge with original data
    df_matched <- df %>%
      left_join(snp_map, by = c("chr", "pos")) %>%
      mutate(
        rsid = ifelse(!is.na(rsid) & rsid != "", rsid, SNP)
      )
  } else {
    # 如果全没匹配上，用原 SNP 标识符兜底 / Fallback to original SNP ID if no matches are found
    df_matched <- df %>% mutate(rsid = SNP) 
  }
  
  # 将 rsid 挪到第一列 / Move rsid to the first column
  setcolorder(df_matched, c("rsid", names(df)))
  
  # 导出文件 / Export file
  fwrite(df_matched, output_file, sep = "\t", quote = FALSE, nThread = 4)
  
  matched_num <- sum(grepl("^rs", df_matched$rsid))
  message(sprintf("   -> Complete! Total %d sites, %d successfully matched to rsid (match rate %.2f%%).\n", 
                  nrow(df_matched), matched_num, (matched_num / nrow(df_matched)) * 100))
}

message("=== All subtype data rsID mapping completed! ===")