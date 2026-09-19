library(dplyr)
library(readr)
library(stringr)

base_dir <- "MR/slopehunter/slopehunter_out"

# Find all SlopeHunter fit.rds files
fit_files <- list.files(
  base_dir,
  pattern = "_fit\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(fit_files), "fit.rds files\n")

results <- lapply(fit_files, function(f) {

  x <- readRDS(f)

  # Phenotype = name of the folder containing the rds file
  phenotype <- basename(dirname(f))

  data.frame(
    phenotype = phenotype,
    slope = as.numeric(x$b),
    slope_se = as.numeric(x$bse),
    ci_lower = as.numeric(x$b_CI[1]),
    ci_upper = as.numeric(x$b_CI[2]),
    mixture_pi = as.numeric(x$pi[1]),
    iterations = as.numeric(x$iter)
  )
})

summary_stats <- bind_rows(results)

# Save as a NEW file
output_file <- file.path(
  base_dir,
  "slopehunter_summary_stats_recovered.tsv"
)

write_tsv(summary_stats, output_file)

cat("\nSaved to:\n", output_file, "\n")
cat("Number of phenotypes:", nrow(summary_stats), "\n")
