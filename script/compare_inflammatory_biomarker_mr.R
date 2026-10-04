#!/usr/bin/env Rscript

# Compare primary (IVW) inflammatory-biomarker MR results from UKB-PPP, deCODE,
# and Fenland.  The q-values are those supplied in each original MR run (BH
# correction across all successful tests in that run).  Shared hits are matched
# by UniProt accession and endpoint, not by assay identifier; when several
# assays have the same accession, the lowest-P nominal assay is retained for a
# database in the shared-hit summaries.
#
# Usage from the repository root:
#   Rscript script/compare_inflammatory_biomarker_mr.R

suppressPackageStartupMessages({
  library(data.table)
})

out_dir <- "comparison_results/inflammatory_biomarkers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sources <- data.table(
  database = rep(c("UKB-PPP", "deCODE", "Fenland"), 3),
  outcome_group = rep(c("incidence", "survival", "slopehunter_adjusted_survival"), each = 3),
  file = c(
    "cis/results/mr_results_cis_incidence/MR_primary_results.csv",
    "replicate/decode/mr_F10_incidence/MR_primary_results.csv",
    "replicate/fenland/mr_F10_incidence/MR_primary_results.csv",
    "cis/results/mr_results_cis_survival/MR_primary_results.csv",
    "replicate/decode/mr_F10_survival/MR_primary_results.csv",
    "replicate/fenland/mr_F10_survival/MR_primary_results.csv",
    "cis/results/mr_results_cis_slopehunter/MR_primary_results.csv",
    "replicate/decode/mr_F10_sh/MR_primary_results.csv",
    "replicate/fenland/mr_F10_sh/MR_primary_results.csv"
  )
)

accession_from_id <- function(x) {
  # Canonical six-character UniProt accessions present in all three assay IDs.
  m <- regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]", x, perl = TRUE)
  ans <- regmatches(x, m)
  ans[m < 0] <- NA_character_
  ans
}
biomarker_label_from_id <- function(x) {
  # deCODE/Fenland put a readable assay name before "__ACCESSION__"; UKB-PPP
  # puts a gene/protein symbol immediately before "_ACCESSION_OID".
  ans <- sub("__[OPQ][0-9][A-Z0-9]{3}[0-9].*$", "", x, perl = TRUE)
  needs_ukb_rule <- !grepl("__", x, fixed = TRUE)
  ans[needs_ukb_rule] <- sub("_[OPQ][0-9][A-Z0-9]{3}[0-9].*$", "", x[needs_ukb_rule], perl = TRUE)
  ans
}
endpoint_key <- function(x) {
  x <- sub("\\.slopehunter$", "", x, ignore.case = TRUE)
  x <- tolower(x)
  gsub("[^a-z0-9]+", "_", x)
}
pretty_endpoint <- function(x) {
  x <- sub("\\.slopehunter$", "", x, ignore.case = TRUE)
  x <- gsub("_with_rsid$", "", x, ignore.case = TRUE)
  gsub("_", " ", x, fixed = TRUE)
}

read_one <- function(database, outcome_group, file) {
  if (!file.exists(file)) stop("Missing input file: ", file)
  x <- fread(file)
  needed <- c("exposure_id_file", "outcome_id_file", "b", "se", "pval", "qval_bh",
              "exp_b", "exp_lci95", "exp_uci95", "nsnp")
  if (!all(needed %in% names(x))) stop("Unexpected columns in: ", file)
  x <- x[, ..needed]
  x[, `:=`(
    database = database,
    outcome_group = outcome_group,
    biomarker_accession = accession_from_id(exposure_id_file),
    biomarker_label = biomarker_label_from_id(exposure_id_file),
    endpoint_key = endpoint_key(outcome_id_file),
    endpoint = pretty_endpoint(outcome_id_file),
    fdr_significant = qval_bh < 0.05,
    evidence = fifelse(qval_bh < 0.05, "FDR significant (q < 0.05)",
                       fifelse(pval < 0.05, "Suggestive (P < 0.05)", "Not significant"))
  )]
  setcolorder(x, c("outcome_group", "database", "biomarker_accession", "biomarker_label", "exposure_id_file",
                   "endpoint", "endpoint_key", "nsnp", "b", "se", "exp_b", "exp_lci95",
                   "exp_uci95", "pval", "qval_bh", "fdr_significant", "evidence"))
  x
}

all_results <- rbindlist(Map(read_one, sources$database, sources$outcome_group, sources$file),
                         use.names = TRUE)
if (anyNA(all_results$biomarker_accession)) {
  warning("Some assay identifiers had no recognised UniProt accession.")
}

# A single table gives the requested lists and makes thresholding auditable.
fwrite(all_results, file.path(out_dir, "all_primary_results_with_evidence_tier.csv"))

db_cols <- c("UKB-PPP" = "#0072B2", "deCODE" = "#D55E00", "Fenland" = "#009E73")
group_titles <- c(incidence = "Breast cancer incidence",
                  survival = "Breast cancer survival",
                  slopehunter_adjusted_survival = "SlopeHunter-adjusted breast cancer survival")

make_forest <- function(x, title, filename, effect_label) {
  if (!nrow(x)) {
    pdf(filename, width = 10, height = 4)
    plot.new(); title(main = title)
    text(0.5, 0.5, "No results meet this evidence threshold", cex = 1.2)
    dev.off()
    return(invisible())
  }
  x <- copy(x)
  setorder(x, endpoint, biomarker_accession, pval, database)
  canonical_labels <- x[, .(canonical_biomarker_label = biomarker_label[1]), by = biomarker_accession]
  x <- merge(x, canonical_labels, by = "biomarker_accession", all.x = TRUE, sort = FALSE)
  x[, label := paste0(canonical_biomarker_label, " (", biomarker_accession, ") - ", endpoint)]
  labs_order <- unique(x$label)
  x[, y_base := match(label, labs_order)]
  db_offset <- c("UKB-PPP" = -0.22, "deCODE" = 0, "Fenland" = 0.22)
  x[, y := y_base + unname(db_offset[database])]
  finite_x <- c(x$exp_lci95, x$exp_uci95)
  finite_x <- finite_x[is.finite(finite_x) & finite_x > 0]
  xlim <- range(finite_x)
  xlim <- exp(log(xlim) + c(-0.15, 0.15) * diff(log(xlim) + c(-1e-5, 1e-5)))
  pdf(filename, width = 13, height = max(5, length(labs_order) * 0.25 + 1.5))
  par(mar = c(4.5, 14, 3.5, 1), xpd = FALSE)
  plot(NA, xlim = xlim, ylim = c(0.4, length(labs_order) + 0.6), log = "x", yaxt = "n",
       xlab = paste0(effect_label, " (95% CI; log scale)"), ylab = "", main = title)
  axis(2, at = seq_along(labs_order), labels = labs_order, las = 1, cex.axis = 0.62)
  abline(v = 1, lty = 2, col = "grey55")
  segments(x$exp_lci95, x$y, x$exp_uci95, x$y, col = db_cols[x$database], lwd = 1.1)
  points(x$exp_b, x$y, pch = ifelse(x$fdr_significant, 16, 1), col = db_cols[x$database], cex = 0.9)
  legend("top", inset = c(0, -0.03), legend = names(db_cols), col = db_cols, pch = 16,
         horiz = TRUE, bty = "n", title = "Biomarker GWAS")
  legend("bottomright", legend = c("FDR significant", "Suggestive only"), pch = c(16, 1),
         bty = "n", cex = 0.75)
  dev.off()
}

make_replication_dotplot <- function(x, title, filename) {
  if (!nrow(x)) return(invisible())
  x <- copy(x)
  # One label per protein-endpoint; x is already restricted to the best nominal
  # assay for each of the three databases.
  labs <- x[, .(biomarker_label = biomarker_label[1],
                endpoint = endpoint[1], mean_logp = mean(-log10(pval))),
            by = .(biomarker_accession, endpoint_key)]
  setorder(labs, -mean_logp, biomarker_accession, endpoint_key)
  labs[, result_label := paste0(biomarker_label, " (", biomarker_accession, ") - ", endpoint)]
  x <- merge(x, labs[, .(biomarker_accession, endpoint_key, result_label)],
             by = c("biomarker_accession", "endpoint_key"), all.x = TRUE, sort = FALSE)
  x[, y := match(result_label, labs$result_label)]
  x[, xpos := match(database, names(db_cols))]
  # Blue indicates a protective/lower-hazard association and red a harmful/
  # higher-hazard association.  Saturation is capped at the 95th percentile.
  cap <- max(quantile(abs(x$b), 0.95, na.rm = TRUE), 0.05)
  palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  idx <- pmax(1, pmin(101, round((pmax(-cap, pmin(cap, x$b)) + cap) / (2 * cap) * 100) + 1))
  point_col <- palette[idx]
  point_cex <- 0.7 + pmin(-log10(x$pval), 5) * 0.42
  pdf(filename, width = 12, height = max(5, nrow(labs) * 0.42 + 1.6))
  par(mar = c(4.3, 17, 3.5, 1), xpd = FALSE)
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.4, nrow(labs) + 0.6), xaxt = "n", yaxt = "n",
       xlab = "Biomarker GWAS database", ylab = "", main = title)
  abline(v = 1:3, col = "grey88")
  axis(1, at = 1:3, labels = names(db_cols))
  axis(2, at = seq_len(nrow(labs)), labels = labs$result_label, las = 1, cex.axis = 0.68)
  points(x$xpos, x$y, pch = ifelse(x$fdr_significant, 16, 21), bg = point_col,
         col = "grey25", cex = point_cex)
  legend("topright", legend = c("higher effect estimate", "lower effect estimate",
                                 "larger = smaller P", "filled = FDR significant"),
         pch = c(21, 21, 21, 16), pt.bg = c("#B2182B", "#2166AC", "grey70", "grey25"),
         col = "grey25", bty = "n", cex = 0.75)
  dev.off()
}

shared_rows <- list()
summary_rows <- list()
for (grp in names(group_titles)) {
  d <- all_results[outcome_group == grp]
  fdr <- d[fdr_significant == TRUE]
  suggestive <- d[pval < 0.05 & fdr_significant == FALSE]
  fwrite(fdr, file.path(out_dir, paste0(grp, "_fdr_significant_biomarkers.csv")))
  fwrite(suggestive, file.path(out_dir, paste0(grp, "_suggestive_biomarkers.csv")))
  fwrite(d[pval < 0.05], file.path(out_dir, paste0(grp, "_all_nominal_hits.csv")))

  effect_label <- if (grp == "incidence") "Odds ratio" else "Hazard ratio"
  make_forest(fdr, paste(group_titles[[grp]], "- FDR-significant results"),
              file.path(out_dir, paste0(grp, "_fdr_significant_forest.pdf")), effect_label)
  make_forest(suggestive, paste(group_titles[[grp]], "- suggestive results (P < 0.05, q >= 0.05)"),
              file.path(out_dir, paste0(grp, "_suggestive_forest.pdf")), effect_label)

  # A replicated nominal hit must be present at P < 0.05 in every database for
  # the same protein and endpoint.  This is deliberately distinct from FDR
  # replication, which is reported separately below.
  nominal_counts <- d[pval < 0.05 & !is.na(biomarker_accession),
                      .(n_databases = uniqueN(database)), by = .(biomarker_accession, endpoint_key)]
  fdr_counts <- d[fdr_significant == TRUE & !is.na(biomarker_accession),
                  .(n_fdr_databases = uniqueN(database)), by = .(biomarker_accession, endpoint_key)]
  replicated_keys <- nominal_counts[n_databases == 3]
  # Multiple related assays can share an accession (for example C3, C3a and
  # C3d).  Retain the strongest nominal assay for each database so each plotted
  # point represents one protein-level result per database and endpoint.
  replicated <- d[pval < 0.05 & !is.na(biomarker_accession)]
  replicated <- replicated[, .SD[which.min(pval)], by = .(biomarker_accession, endpoint_key, database)]
  replicated <- merge(replicated, replicated_keys, by = c("biomarker_accession", "endpoint_key"))
  replicated <- merge(replicated, fdr_counts, by = c("biomarker_accession", "endpoint_key"), all.x = TRUE)
  replicated[is.na(n_fdr_databases), n_fdr_databases := 0L]
  setorder(replicated, endpoint, biomarker_accession, database)
  fwrite(replicated, file.path(out_dir, paste0(grp, "_nominally_significant_in_all_three_databases.csv")))
  shared_rows[[grp]] <- replicated
  summary_rows[[grp]] <- d[, .(
    tests = .N,
    fdr_significant = sum(fdr_significant == TRUE),
    suggestive = sum(pval < 0.05 & fdr_significant == FALSE),
    nominally_significant = sum(pval < 0.05)
  ), by = database][, outcome_group := grp]

  # Summary: every displayed protein-endpoint combination is nominally
  # significant in all three datasets; point shape preserves its FDR status.
  make_forest(replicated, paste(group_titles[[grp]], "- nominal P < 0.05 in all three databases"),
              file.path(out_dir, paste0(grp, "_shared_nominal_hits_all_three_databases.pdf")), effect_label)
  make_replication_dotplot(replicated,
                           paste(group_titles[[grp]], "- replicated nominal evidence across databases"),
                           file.path(out_dir, paste0(grp, "_replicated_nominal_evidence_dotplot.pdf")))
}

summary_table <- rbindlist(summary_rows)
summary_table[, outcome := group_titles[outcome_group]]
setcolorder(summary_table, c("outcome_group", "outcome", "database", "tests", "fdr_significant",
                             "suggestive", "nominally_significant"))
fwrite(summary_table, file.path(out_dir, "evidence_count_summary.csv"))

shared_all <- rbindlist(shared_rows, idcol = "outcome_group_from_list", fill = TRUE)
if (nrow(shared_all)) {
  fwrite(shared_all, file.path(out_dir, "all_shared_nominal_hits_all_three_databases.csv"))
  priority <- shared_all[, .(
    biomarker_labels = paste(sort(unique(biomarker_label)), collapse = " / "),
    endpoint = endpoint[1],
    databases = uniqueN(database),
    direction_across_databases = if (all(b > 0)) "positive in all three" else if (all(b < 0)) "negative in all three" else "mixed",
    min_p = min(pval), max_p = max(pval), mean_neg_log10_p = mean(-log10(pval)),
    min_q = min(qval_bh),
    median_effect_ratio = median(exp_b),
    effect_ratio_range = paste0(formatC(min(exp_b), format = "f", digits = 2), " to ",
                                formatC(max(exp_b), format = "f", digits = 2))
  ), by = .(outcome_group, biomarker_accession, endpoint_key)]
  setorder(priority, outcome_group, -mean_neg_log10_p, max_p)
  priority[, replication_rank_within_outcome := seq_len(.N), by = outcome_group]
  fwrite(priority, file.path(out_dir, "replicated_nominal_priority_summary.csv"))

  # --- Focused cross-platform replication displays and exploratory meta-analysis ---
  # Pooling log OR/HRs is appropriate only if protein effects are on comparable
  # SD scales.  We therefore retain every source estimate and report I2 beside
  # the random-effects estimate rather than treating the pooled result as the
  # primary analysis.
  priority_accessions <- c("P01009", "Q9Y240", "P15692", "P02776", "P01024")
  pretty_targets <- c(
    P01009 = "SERPINA1 / alpha-1-antitrypsin", Q9Y240 = "CLEC11A / SCGF",
    P15692 = "VEGFA", P02776 = "PF4 / CXCL4", P01024 = "C3 / C3a"
  )
  focus <- shared_all[biomarker_accession %in% priority_accessions]
  focus[, target := unname(pretty_targets[biomarker_accession])]
  focus[is.na(target), target := biomarker_accession]
  fwrite(focus, file.path(out_dir, "priority_targets_source_estimates.csv"))

  make_replication_heatmap <- function(x, filename) {
    x <- copy(x)
    pairs <- x[, .(target = target[1], endpoint = endpoint[1], mean_logp = mean(-log10(pval))),
               by = .(outcome_group, biomarker_accession, endpoint_key)]
    setorder(pairs, outcome_group, -mean_logp, target, endpoint)
    short_outcome <- c(incidence = "Incidence", survival = "Survival",
                       slopehunter_adjusted_survival = "SlopeHunter survival")
    pairs[, label := paste0("[", unname(short_outcome[outcome_group]), "] ", target, " - ", endpoint)]
    x <- merge(x, pairs[, .(outcome_group, biomarker_accession, endpoint_key, label)],
               by = c("outcome_group", "biomarker_accession", "endpoint_key"), all.x = TRUE, sort = FALSE)
    x[, ypos := match(label, pairs$label)]
    x[, xpos := match(database, names(db_cols))]
    cap <- max(quantile(abs(x$b), 0.95, na.rm = TRUE), 0.05)
    pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
    idx <- pmax(1, pmin(101, round((pmax(-cap, pmin(cap, x$b)) + cap) / (2 * cap) * 100) + 1))
    pdf(filename, width = 13, height = max(7, nrow(pairs) * 0.36 + 2))
    par(mar = c(4.5, 20, 3.5, 1))
    plot(NA, xlim = c(0.5, 3.5), ylim = c(0.4, nrow(pairs) + 0.6), xaxt = "n", yaxt = "n",
         xlab = "Protein GWAS database", ylab = "",
         main = "Replication heatmap: prioritised protein-endpoint associations")
    for (i in seq_len(nrow(x))) {
      rect(x$xpos[i] - 0.43, x$ypos[i] - 0.37, x$xpos[i] + 0.43, x$ypos[i] + 0.37,
           col = pal[idx[i]], border = "white")
      text(x$xpos[i], x$ypos[i], sprintf("%.3g", x$pval[i]), cex = 0.6)
    }
    axis(1, at = 1:3, labels = names(db_cols))
    axis(2, at = seq_len(nrow(pairs)), labels = pairs$label, las = 1, cex.axis = 0.62)
    abline(h = seq(1.5, nrow(pairs) - 0.5, by = 1), col = "grey88")
    legend("topright", legend = c("protective / lower hazard", "higher risk / hazard", "cell text = P value"),
           fill = c("#2166AC", "#B2182B", NA), border = NA, bty = "n", cex = 0.75)
    mtext("Tile colour: source-specific log OR/HR (blue negative, red positive)", side = 3, line = 0.3, cex = 0.75)
    dev.off()
  }
  make_replication_heatmap(focus, file.path(out_dir, "priority_targets_replication_heatmap.pdf"))

  random_effects <- function(z) {
    w <- 1 / z$se^2
    fixed_b <- sum(w * z$b) / sum(w)
    Q <- sum(w * (z$b - fixed_b)^2)
    df <- nrow(z) - 1
    C <- sum(w) - sum(w^2) / sum(w)
    tau2 <- if (C > 0) max(0, (Q - df) / C) else 0
    wr <- 1 / (z$se^2 + tau2)
    b <- sum(wr * z$b) / sum(wr)
    se <- sqrt(1 / sum(wr))
    data.table(meta_b = b, meta_se = se, meta_lci = b - 1.96 * se, meta_uci = b + 1.96 * se,
               meta_p = 2 * pnorm(abs(b / se), lower.tail = FALSE), Q = Q, df = df,
               I2_percent = if (Q > 0) max(0, (Q - df) / Q) * 100 else 0, tau2 = tau2)
  }
  meta <- focus[, random_effects(.SD), by = .(outcome_group, biomarker_accession, endpoint_key)]
  meta <- merge(meta, focus[, .(target = target[1], endpoint = endpoint[1], min_source_p = min(pval)),
                           by = .(outcome_group, biomarker_accession, endpoint_key)],
                by = c("outcome_group", "biomarker_accession", "endpoint_key"), all.x = TRUE)
  meta[, `:=`(meta_effect_ratio = exp(meta_b), meta_lci_ratio = exp(meta_lci), meta_uci_ratio = exp(meta_uci))]
  setorder(meta, outcome_group, meta_p)
  fwrite(meta, file.path(out_dir, "priority_targets_random_effects_meta_analysis.csv"))

  make_meta_forest <- function(x, m, filename) {
    pdf(filename, width = 13, height = 10, onefile = TRUE)
    for (grp in unique(x$outcome_group)) {
      sx <- x[outcome_group == grp]
      sm <- m[outcome_group == grp]
      keys <- sm[, paste(biomarker_accession, endpoint_key)]
      sx[, key := paste(biomarker_accession, endpoint_key)]
      sx[, source_order := match(database, names(db_cols))]
      sx[, group_order := match(key, keys)]
      setorder(sx, group_order, source_order)
      sm[, key := paste(biomarker_accession, endpoint_key)]
      n_groups <- nrow(sm)
      # Four rows per result: a heading, three source estimates, then pooled RE.
      y_top <- seq(from = n_groups * 5, by = -5, length.out = n_groups)
      group_y <- setNames(y_top, sm$key)
      sx[, y := group_y[key] - source_order]
      sm[, y := group_y[key] - 4]
      xs <- c(sx$exp_lci95, sx$exp_uci95, sm$meta_lci_ratio, sm$meta_uci_ratio)
      xs <- xs[is.finite(xs) & xs > 0]
      xlim <- range(xs)
      xlim <- exp(log(xlim) + c(-0.12, 0.12) * diff(log(xlim)))
      par(mar = c(4.5, 18, 3.8, 1))
      plot(NA, xlim = xlim, ylim = c(0, n_groups * 5 + 1), log = "x", yaxt = "n",
           xlab = if (grp == "incidence") "Odds ratio (95% CI; log scale)" else "Hazard ratio (95% CI; log scale)",
           ylab = "", main = paste(group_titles[[grp]], "- random-effects meta-analysis"))
      abline(v = 1, lty = 2, col = "grey55")
      for (i in seq_len(n_groups)) {
        mm <- sm[i]
        text(xlim[1], group_y[mm$key], labels = paste0(mm$target, " - ", mm$endpoint,
             "  [I2=", formatC(mm$I2_percent, format = "f", digits = 0), "%]"), pos = 4, font = 2, cex = 0.75)
      }
      segments(sx$exp_lci95, sx$y, sx$exp_uci95, sx$y, col = db_cols[sx$database], lwd = 1)
      points(sx$exp_b, sx$y, pch = 16, col = db_cols[sx$database], cex = 0.8)
      segments(sm$meta_lci_ratio, sm$y, sm$meta_uci_ratio, sm$y, col = "black", lwd = 1.4)
      points(sm$meta_effect_ratio, sm$y, pch = 18, col = "black", cex = 1)
      axis_labels <- c(sx[, .(y, label = database)], sm[, .(y, label = "Random effects")])
      axis(2, at = axis_labels$y, labels = axis_labels$label, las = 1, cex.axis = 0.65)
      legend("topright", legend = c(names(db_cols), "Random effects"), col = c(db_cols, "black"),
             pch = c(16, 16, 16, 18), bty = "n", cex = 0.75)
      mtext("Exploratory pooling; assumes comparable per-SD protein scaling across platforms", side = 3, line = 0.35, cex = 0.7)
    }
    dev.off()
  }
  make_meta_forest(focus, meta, file.path(out_dir, "priority_targets_random_effects_meta_forest.pdf"))

  # A two-lane Sankey avoids an invalid single chain: FDR and three-database
  # replication are not nested events.  The upper lane is database-specific
  # testing; the lower lane is matched protein-endpoint replication.
  make_sankey <- function(all, shared, filename) {
    pdf(filename, width = 14, height = 8)
    par(mar = c(1, 1, 3, 1), xpd = NA)
    plot.new(); plot.window(xlim = c(0, 1.3), ylim = c(0, 1))
    rect_node <- function(x, y, h, text, col) {
      rect(x - 0.06, y - h / 2, x + 0.06, y + h / 2, col = col, border = "white")
      text(x, y, text, cex = 0.8)
    }
    flow <- function(x1, y1, h1, x2, y2, h2, col) {
      polygon(c(x1 + .06, x2 - .06, x2 - .06, x1 + .06),
              c(y1 - h1/2, y2 - h2/2, y2 + h2/2, y1 + h1/2), col = col, border = NA)
    }
    cols <- c("#B3CDE3", "#6497B1", "#005B96", "#41AB5D", "#238B45")
    for (i in seq_along(group_titles)) {
      grp <- names(group_titles)[i]
      yy <- 0.84 - (i - 1) * 0.28
      d <- all[outcome_group == grp]
      # upper lane: database-specific association tests
      n_test <- nrow(d); n_nom <- nrow(d[pval < .05]); n_fdr <- nrow(d[fdr_significant == TRUE])
      # lower lane: matched protein-endpoint pairs measured by all three sources
      pres <- d[!is.na(biomarker_accession), .(n_db = uniqueN(database)), by = .(biomarker_accession, endpoint_key)]
      n_common <- nrow(pres[n_db == 3])
      reps <- unique(shared[outcome_group == grp, .(biomarker_accession, endpoint_key)])
      n_rep <- nrow(reps)
      repdat <- merge(d, reps, by = c("biomarker_accession", "endpoint_key"))
      n_concordant <- repdat[, .(same_sign = all(b > 0) || all(b < 0)), by = .(biomarker_accession, endpoint_key)][same_sign == TRUE, .N]
      h1 <- 0.10; h2 <- max(0.015, h1 * n_nom / max(n_test, 1)); h3 <- max(0.008, h1 * n_fdr / max(n_test, 1))
      flow(.14, yy + .04, h1, .43, yy + .04, h2, adjustcolor(cols[2], .55))
      flow(.43, yy + .04, h2, .72, yy + .04, h3, adjustcolor(cols[3], .55))
      rect_node(.14, yy + .04, h1, paste0("Tests\n", n_test), cols[1])
      rect_node(.43, yy + .04, h2, paste0("Nominal\n", n_nom), cols[2])
      rect_node(.72, yy + .04, h3, paste0("FDR\n", n_fdr), cols[3])
      h4 <- .075; h5 <- max(.012, h4 * n_rep / max(n_common, 1)); h6 <- max(.012, h4 * n_concordant / max(n_common, 1))
      flow(.14, yy - .09, h4, .43, yy - .09, h5, adjustcolor(cols[4], .55))
      flow(.43, yy - .09, h5, .72, yy - .09, h6, adjustcolor(cols[5], .55))
      rect_node(.14, yy - .09, h4, paste0("Tested in all 3\n", n_common), "#C7E9C0")
      rect_node(.43, yy - .09, h5, paste0("Nominal in all 3\n", n_rep), cols[4])
      rect_node(.72, yy - .09, h6, paste0("Same direction\n", n_concordant), cols[5])
      text(.83, yy, group_titles[[grp]], adj = c(0, .5), font = 2, cex = .85)
    }
    legend("topleft", legend = c("Database-specific evidence", "Matched three-database replication"),
           fill = c(cols[2], cols[4]), border = NA, bty = "n", cex = 0.75)
    title("Inflammatory biomarker MR evidence flow")
    mtext("Top lane: individual MR tests. Bottom lane: UniProt accession + matched endpoint pairs. FDR and replication are separate evidence dimensions.", side = 1, line = -1.5, cex = .7)
    invisible(dev.off())
  }
  make_sankey(all_results, shared_all, file.path(out_dir, "mr_evidence_flow_sankey.pdf"))
}

cat("Wrote comparison tables and plots to: ", normalizePath(out_dir), "\n", sep = "")
