library(TwoSampleMR)
library(dplyr)
library(readr)

exposure_dir <- "MR/exposure_all/lead_snps_all"
outcome_dir  <- "MR/slopehunter/slopehunter_out_gzfile_formr"
output_dir   <- "MR/harmonised_results_sh"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

exp_files <- list.files(exposure_dir, pattern = "\\.txt$", full.names = TRUE)
out_files <- list.files(outcome_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)

cat("found", length(exp_files), "exposure files\n")
cat("found", length(out_files), "outcome files\n\n")

for (exp_path in exp_files) {
  exp_name <- gsub("_lead_snps.txt|\\.txt$", "", basename(exp_path))
  cat("==================================================\n")
  cat(" Exposure:", exp_name, "\n")
  
  # A. read Exposure data
  exp_dat <- read_exposure_data(
    filename = exp_path,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",  
    other_allele_col = "ALLELE0",  
    eaf_col = "A1FREQ",              
    pval_col = "P",
    samplesize_col = "N"
  )
  
  if (nrow(exp_dat) == 0) {
    cat("  [Skip] Exposure is empty, skipping.\n")
    next
  }

  for (out_path in out_files) {
    out_name <- gsub("\\.tsv\\.gz$", "", basename(out_path))
    cat("  -> matching Outcome:", out_name, "\n")
    
    # B. read Outcome data (only SNPs present in Exposure, significantly speeding up execution)
    tryCatch({
      out_dat <- read_outcome_data(
        snps = exp_dat$SNP,
        filename = out_path,
        sep = "\t", 
        snp_col = "SNP",       
        beta_col = "ybeta",       
        se_col = "yse",           
        effect_allele_col = "EA.prognosis",  
        other_allele_col = "OA.prognosis",   
        eaf_col = "EAF.prognosis",           
        pval_col = "PVAL.prognosis"            
      )
      
      if (is.null(out_dat) || nrow(out_dat) == 0) {
        cat("     [Warning] No matching SNPs found in Outcome data, skipping.\n")
        next
      }

      # C. Doing Harmonisation 
      # action = 2 means flipping the strand and discarding palindromic or unalignable SNPs.
      dat <- harmonise_data(
        exposure_dat = exp_dat,
        outcome_dat = out_dat,
        action = 2
      )

      # D. Save the harmonised data to CSV
      if (nrow(dat) > 0) {
        save_file <- file.path(output_dir, paste0(exp_name, "__VS__", out_name, "_harmonised.csv"))
        write_csv(dat, save_file)
        cat("     [Success] successfully matched ", nrow(dat), " SNPs, results saved!\n")
      } else {
        cat("     [Warning] No available SNPs after Harmonisation.\n")
      }
      
    }, error = function(e) {
      cat("     [Error]  Harmonise failed:", conditionMessage(e), "\n")
    })
  }
}

cat("\n Harmonisation Done, results saved to:", output_dir, "\n")