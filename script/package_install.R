#!/usr/bin/env Rscript

# 1. 确保个人库路径正常
user_dir <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_dir)) dir.create(user_dir, recursive = TRUE)
.libPaths(c(user_dir, .libPaths()))

# 2. 一次性补全报错信息中缺失的所有底层依赖包
deps <- c(
  "BiocGenerics", "S4Vectors", "IRanges", "XVector", 
  "GenomeInfoDb", "GenomicRanges", "Biostrings", "BiocIO", 
  "XML", "restfulr", "Rsamtools", "GenomicAlignments", 
  "rtracklayer", "BSgenome"
)

# 强制安装依赖包 (不更新无关的旧包，加快速度)
BiocManager::install(deps, lib = user_dir, update = FALSE, ask = FALSE)