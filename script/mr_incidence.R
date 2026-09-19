#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)

# 1. Define input and output directories
input_dir  <- "MR/harmonised_results"
output_dir <- "MR/mr_results_summary"
fig_dir    <- file.path(output_dir, "figures_significant")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(fig_dir))    dir.create(fig_dir, recursive = TRUE)

# define 5 output files for different MR results
main_file    <- file.path(output_dir, "ALL_MR_Main_Results.csv")
hetero_file  <- file.path(output_dir, "ALL_MR_Heterogeneity.csv")
pleio_file   <- file.path(output_dir, "ALL_MR_Egger_Pleiotropy.csv")
steiger_file <- file.path(output_dir, "ALL_MR_Steiger_Directionality.csv")
loo_file     <- file.path(output_dir, "ALL_MR_LeaveOneOut_Details.csv")

# Delete existing output files to avoid appending to old results
for (f in c(main_file, hetero_file, pleio_file, steiger_file, loo_file)) {
  if (file.exists(f)) file.remove(f)
}

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
message(sprintf("Found %d harmonised csv files, starting batch processing (real-time disk writing mode)...\n", length(csv_files)))

for (i in seq_along(csv_files)) {
  file <- csv_files[i]
  file_name <- basename(file)
  pair_tag  <- gsub("\\.csv$", "", file_name)
  
  message(sprintf("[%d/%d] Processing: %s ...", i, length(csv_files), file_name))
  
  # Used tryCatch to ensure that any error in processing a single file does not stop the entire batch process
  tryCatch({
    
    dat <- fread(file)
    dat_clean <- subset(dat, mr_keep == TRUE)
    
    if (nrow(dat_clean) < 1) {
      message("   -> [Skip] No valid instruments (mr_keep == TRUE)。")
      next
    }
    
    # --- A. Main MR analsis ---
    res <- mr(dat_clean)
    if (is.null(res) || nrow(res) == 0) next
    
    res_or <- generate_odds_ratios(res)
    res_or$filename <- file_name
    
    fwrite(res_or, main_file, append = file.exists(main_file))
    
    # --- B. Heterogeneity test (Cochran's Q) ---
    if (nrow(dat_clean) >= 2) {
      het <- tryCatch({ mr_heterogeneity(dat_clean) }, error = function(e) NULL)
      if (!is.null(het) && nrow(het) > 0) {
        het$filename <- file_name
        fwrite(het, hetero_file, append = file.exists(hetero_file))
      }
    }
    
    # --- C. Pleiotropy Test (MR Egger) & Steiger Directionality ---
    if (nrow(dat_clean) >= 3) {
      pleio <- tryCatch({ mr_pleio_test(dat_clean) }, error = function(e) NULL)
      if (!is.null(pleio) && nrow(pleio) > 0) {
        pleio$filename <- file_name
        fwrite(pleio, pleio_file, append = file.exists(pleio_file))
      }
      
      st <- tryCatch({ directionality_test(dat_clean) }, error = function(e) NULL)
      if (!is.null(st) && nrow(st) > 0) {
        st$filename <- file_name
        fwrite(st, steiger_file, append = file.exists(steiger_file))
      }
    }
    
    # --- D. Leave-one-out  ---
    loo <- tryCatch({ mr_leaveoneout(dat_clean) }, error = function(e) NULL)
    if (!is.null(loo) && nrow(loo) > 0) {
      loo$filename <- file_name
      fwrite(loo, loo_file, append = file.exists(loo_file))
    }
    
    # --- E. Significant (IVW P < 0.05) automatically export plots ---
    ivw_p <- res_or %>% filter(method == "Inverse variance weighted") %>% pull(pval)
    
    if (length(ivw_p) > 0 && !is.na(ivw_p) && ivw_p < 0.05) {
      message("   --> [significant]  Funnel plot 与 Leave-one-out plot ...")
      
      # 1. Funnel plot
      res_single <- tryCatch({ mr_singlesnp(dat_clean) }, error = function(e) NULL)
      if (!is.null(res_single)) {
        p_funnel <- mr_funnel_plot(res_single)
        ggsave(file.path(fig_dir, paste0("Funnel_", pair_tag, ".pdf")), 
               p_funnel[[1]], width = 6, height = 6)
      }
      
      # 2. Leave-one-out 
      if (!is.null(loo)) {
        p_loo <- mr_leaveoneout_plot(loo)
        fig_height <- max(5, nrow(dat_clean) * 0.3)
        ggsave(file.path(fig_dir, paste0("LOO_", pair_tag, ".pdf")), 
               p_loo[[1]], width = 7, height = fig_height)
      }
    }
    
  }, error = function(e) {
    warning(sprintf("Processing file %s caught an error and skipped: %s", file_name, e$message))
  })
}

# 3. Final summary of significant IVW results
if (file.exists(main_file)) {
  all_main <- fread(main_file)
  sig_main <- all_main %>% filter(method == "Inverse variance weighted" & pval < 0.05)
  fwrite(sig_main, file.path(output_dir, "SIGNIFICANT_IVW_Results.csv"))
  
  message(sprintf("\n🎉 Run completed! A total of %d IVW significant positive results (P < 0.05) were identified!", nrow(sig_main)))
  message("All data and plots have been saved to directory: ", output_dir)
} else {
  message("\nNo valid analysis results were generated. Please check the harmonised input data.")
}