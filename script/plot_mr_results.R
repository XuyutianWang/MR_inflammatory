#!/usr/bin/env Rscript

# Scalable overview plots for one MR analysis family.
# Usage: Rscript MR/script/plot_mr_results.R <mr_output_dir> <plot_dir> <analysis_label> <effect_label>
# Example effect labels: OR (incidence) or HR (survival).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript MR/script/plot_mr_results.R <mr_output_dir> <plot_dir> <analysis_label> <effect_label>", call. = FALSE)
}
results_dir <- args[1]; plot_dir <- args[2]; analysis_label <- args[3]; effect_label <- args[4]
input_file <- file.path(results_dir, "MR_primary_results.csv")
if (!file.exists(input_file)) stop("Cannot find: ", input_file, call. = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file)
required <- c("filename", "b", "se", "pval", "qval_bh", "exp_b", "exp_lci95", "exp_uci95")
if (!all(required %in% names(dt))) stop("MR_primary_results.csv has unexpected columns.", call. = FALSE)

parse_labels <- function(filename) {
  tag <- sub("_harmonised\\.csv$", "", filename)
  parts <- strsplit(tag, "__VS__", fixed = TRUE)[[1]]
  exposure_id <- parts[1]
  outcome_id <- if (length(parts) > 1) parts[2] else "Unknown outcome"
  list(protein = sub("_.*$", "", exposure_id), outcome = gsub("_", " ", outcome_id, fixed = TRUE))
}
labels <- lapply(dt$filename, parse_labels)
dt[, protein := vapply(labels, `[[`, character(1), "protein")]
dt[, outcome_label := vapply(labels, `[[`, character(1), "outcome")]
dt[, neglog10q := pmin(-log10(pmax(qval_bh, .Machine$double.xmin)), 12)]
dt[, fdr_label := ifelse(qval_bh < 0.05, "FDR < 0.05", "Not FDR significant")]

# 1. Faceted volcano overview.
volcano <- ggplot(dt, aes(x = b, y = -log10(pmax(pval, .Machine$double.xmin)))) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(aes(colour = fdr_label), alpha = 0.75, size = 1.4) +
  scale_colour_manual(values = c("FDR < 0.05" = "#c51b7d", "Not FDR significant" = "grey70")) +
  facet_wrap(~outcome_label, scales = "free_y") +
  labs(title = paste("MR volcano plots:", analysis_label), x = "MR beta", y = expression(-log[10](P)), colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "top", strip.text = element_text(size = 7))
ggsave(file.path(plot_dir, "01_volcano_by_outcome.pdf"), volcano, width = 14, height = 10)

# 2. Rain plot: fill is effect direction/magnitude; circle area is BH-FDR evidence.
protein_order <- dt[, .(best_q = min(qval_bh, na.rm = TRUE)), by = protein][order(best_q, protein)]
dt[, protein := factor(protein, levels = rev(protein_order$protein))]
rain <- ggplot(dt, aes(x = outcome_label, y = protein)) +
  geom_point(aes(size = neglog10q, fill = b, colour = fdr_label), shape = 21, stroke = 0.25) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "MR beta", oob = squish) +
  scale_colour_manual(values = c("FDR < 0.05" = "black", "Not FDR significant" = "grey80"), guide = "none") +
  scale_size_continuous(range = c(0.4, 5), name = expression(-log[10]("BH q"))) +
  labs(title = paste("MR rain plot:", analysis_label), x = NULL, y = NULL,
       subtitle = "Fill: direction/magnitude; circle size: BH-FDR evidence; black outline: q < 0.05") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6), axis.text.y = element_text(size = 5))
ggsave(file.path(plot_dir, "02_rain_plot.pdf"), rain, width = 16,
       height = max(10, length(levels(dt$protein)) * 0.12), limitsize = FALSE)

# 3. Forest plot of at most 30 FDR-significant results.
hits <- dt[qval_bh < 0.05][order(qval_bh)][1:min(.N, 30)]
if (nrow(hits)) {
  hits[, hit_label := paste(protein, "—", outcome_label)]
  hits[, hit_label := factor(hit_label, levels = rev(hit_label))]
  forest <- ggplot(hits, aes(y = hit_label, x = exp_b, xmin = exp_lci95, xmax = exp_uci95)) +
    geom_vline(xintercept = 1, colour = "grey60", linetype = 2) + geom_errorbarh(height = 0.2) +
    geom_point(size = 1.8, colour = "#2c7fb8") + scale_x_log10() +
    labs(title = paste("Top FDR-significant MR results:", analysis_label), x = paste(effect_label, "(95% CI; log scale)"), y = NULL) +
    theme_bw(base_size = 9)
  ggsave(file.path(plot_dir, "03_forest_top_FDR_hits.pdf"), forest, width = 10,
         height = max(5, nrow(hits) * 0.28), limitsize = FALSE)
}

# 4. Method-consistency heatmap for the 30 lowest-q primary associations.
all_methods_file <- file.path(results_dir, "MR_all_method_results.csv")
if (file.exists(all_methods_file)) {
  methods <- fread(all_methods_file)
  top_pairs <- dt[order(qval_bh), .(filename)][1:min(.N, 30)]
  methods <- methods[filename %in% top_pairs$filename]
  if (nrow(methods)) {
    methods[, signed_logp := sign(b) * pmin(-log10(pmax(pval, .Machine$double.xmin)), 12)]
    methods[, pair_label := vapply(lapply(filename, parse_labels), function(x) paste(x$protein, "—", x$outcome), character(1))]
    methods[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
    consistency <- ggplot(methods, aes(x = method, y = pair_label, fill = signed_logp)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                           name = "signed\n−log10(P)", limits = c(-12, 12), oob = squish) +
      labs(title = paste("Method-consistency heatmap:", analysis_label), x = NULL, y = NULL) +
      theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    ggsave(file.path(plot_dir, "04_method_consistency_heatmap.pdf"), consistency, width = 12,
           height = max(5, length(unique(methods$pair_label)) * 0.28), limitsize = FALSE)
  }
}
message("Plots written to: ", plot_dir)
