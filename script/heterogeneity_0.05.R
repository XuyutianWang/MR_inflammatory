#!/usr/bin/env Rscript

library(data.table)

df <- fread("MR/mr_results_summary/ALL_MR_Heterogeneity.csv")

df_filtered <- df[df$Q_pval < 0.05, ]

output_file <- "MR/sensitivity/Heterogeneity.csv"
fwrite(df_filtered, output_file)