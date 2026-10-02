#!/usr/bin/env Rscript

# T2D positive-control MR for the BCAC 2021 survival summary statistics.
#
# Exposure: Mahajan et al. (2018b) DIAGRAM European T2D GWAS, unadjusted for BMI.
# Outcomes: raw BCAC 2021 survival summary statistics and Slope-Hunter-adjusted
#           BCAC survival summary statistics.
#
# Run from the project root, for example:
#   Rscript script/mr_t2d_survival_positive_control.R \
#     --plink /path/to/plink
#
# PLINK 1.9 is required for LD clumping. The supplied EUR reference panel is
# GRCh37, which matches the Mahajan 2018b and BCAC 2021 positions.

suppressPackageStartupMessages(library(data.table))

usage <- paste(
  "Usage: Rscript script/mr_t2d_survival_positive_control.R --plink /path/to/plink [options]",
  "",
  "Options:",
  "  --exposure FILE       Default: T2D/Mahajan.NatGenet2018b.T2D.European.txt",
  "  --raw-dir DIR         Default: outcome/survival_2021_gwas_results_files_with_rsid_no_chr23",
  "  --slopehunter-dir DIR Default: slopehunter/slopehunter_out_gzfile_formr",
  "  --ref-bfile PREFIX    Default: reference/EUR",
  "  --out-dir DIR         Default: T2D/Mahajan.NatGenet2018b.T2D.European_results",
  "  --clump-r2 VALUE      Default: 0.001",
  "  --clump-kb VALUE      Default: 10000",
  sep = "\n"
)

parse_args <- function(args) {
  opts <- list(
    exposure = "T2D/Mahajan.NatGenet2018b.T2D.European.txt",
    raw_dir = "outcome/survival_2021_gwas_results_files_with_rsid_no_chr23",
    slopehunter_dir = "slopehunter/slopehunter_out_gzfile_formr",
    ref_bfile = "reference/EUR",
    out_dir = "T2D/Mahajan.NatGenet2018b.T2D.European_results",
    plink = Sys.getenv("PLINK_BIN", unset = "plink"),
    clump_r2 = 0.001,
    clump_kb = 10000
  )
  if (!length(args)) return(opts)
  if ("--help" %in% args || "-h" %in% args) {
    cat(usage, "\n")
    quit(status = 0)
  }
  if (length(args) %% 2L != 0L) stop(usage, call. = FALSE)
  keys <- sub("^--", "", args[seq(1, length(args), 2)])
  values <- args[seq(2, length(args), 2)]
  valid <- names(opts)
  if (any(!keys %in% valid)) stop("Unknown option(s): ", paste(keys[!keys %in% valid], collapse = ", "), "\n", usage, call. = FALSE)
  for (i in seq_along(keys)) opts[[keys[i]]] <- values[i]
  opts$clump_r2 <- as.numeric(opts$clump_r2)
  opts$clump_kb <- as.numeric(opts$clump_kb)
  opts
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
must_exist <- c(opt$exposure, paste0(opt$ref_bfile, c(".bed", ".bim", ".fam")), opt$raw_dir, opt$slopehunter_dir)
if (any(!file.exists(must_exist))) stop("Missing input(s): ", paste(must_exist[!file.exists(must_exist)], collapse = ", "), call. = FALSE)
if (!nzchar(Sys.which(opt$plink)) && !file.exists(opt$plink)) {
  stop("PLINK 1.9 was not found: ", opt$plink, ". Supply --plink /absolute/path/to/plink or set PLINK_BIN.", call. = FALSE)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

complement <- function(x) chartr("ACGT", "TGCA", toupper(x))
pair_key <- function(a, b) {
  a <- toupper(a); b <- toupper(b)
  ifelse(a < b, paste(a, b, sep = "/"), paste(b, a, sep = "/"))
}
is_palindromic <- function(a, b) pair_key(a, b) %in% c("A/T", "C/G")

# Map Mahajan's CHR:POS IDs to rsIDs. Mapping only a genomic position is not
# enough: we additionally require the alleles to match the EUR reference panel
# on either strand.
message("Reading and mapping DIAGRAM T2D instruments to rsIDs...")
required_exposure <- c("SNP", "Chr", "Pos", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue", "Neff")
expo <- fread(opt$exposure, select = required_exposure)
if (!all(required_exposure %in% names(expo))) stop("Unexpected exposure columns. Required: ", paste(required_exposure, collapse = ", "), call. = FALSE)
expo <- expo[Pvalue < 5e-8 & is.finite(Beta) & is.finite(SE) & SE > 0]
expo[, `:=`(Chr = as.integer(Chr), Pos = as.integer(Pos), EA = toupper(EA), NEA = toupper(NEA))]
expo[, pair := pair_key(EA, NEA)]

bim <- fread(paste0(opt$ref_bfile, ".bim"), header = FALSE,
             col.names = c("Chr", "rsid", "cm", "Pos", "A1_ref", "A2_ref"))
bim[, `:=`(Chr = as.integer(Chr), Pos = as.integer(Pos), A1_ref = toupper(A1_ref), A2_ref = toupper(A2_ref))]
bim[, pair := pair_key(A1_ref, A2_ref)]
bim[, pair_complement := pair_key(complement(A1_ref), complement(A2_ref))]
mapped <- merge(expo, bim[, .(Chr, Pos, rsid, pair, pair_complement)], by = c("Chr", "Pos"), allow.cartesian = TRUE)
mapped <- mapped[pair.x == pair.y | pair.x == pair_complement]
setnames(mapped, c("pair.x", "pair.y"), c("pair", "reference_pair"))
mapped <- mapped[!duplicated(rsid)]
if (!nrow(mapped)) stop("No genome-wide-significant T2D variants mapped to EUR rsIDs.", call. = FALSE)

candidate_file <- file.path(opt$out_dir, "t2d_gws_candidates_mapped_to_rsid.tsv")
fwrite(mapped[, .(rsid, Chr, Pos, EA, NEA, EAF, Beta, SE, Pvalue, Neff)], candidate_file, sep = "\t")
message("Mapped ", nrow(mapped), " genome-wide-significant T2D variants. LD clumping with PLINK...")

clump_prefix <- file.path(opt$out_dir, "t2d_diagram_clump")
plink_args <- c("--bfile", opt$ref_bfile, "--clump", candidate_file,
                "--clump-snp-field", "rsid", "--clump-field", "Pvalue",
                "--clump-p1", "5e-8", "--clump-r2", as.character(opt$clump_r2),
                "--clump-kb", as.character(opt$clump_kb), "--out", clump_prefix)
plink_status <- system2(opt$plink, plink_args)
if (plink_status != 0L) stop("PLINK clumping failed; inspect the PLINK output in ", opt$out_dir, call. = FALSE)
clumped_file <- paste0(clump_prefix, ".clumped")
if (!file.exists(clumped_file) || file.info(clumped_file)$size == 0) stop("PLINK did not produce a clumped-instrument file.", call. = FALSE)
clumped <- fread(clumped_file, fill = TRUE)
if (!"SNP" %in% names(clumped) || !nrow(clumped)) stop("No instruments remain after LD clumping.", call. = FALSE)
inst <- mapped[rsid %in% clumped$SNP]
setorder(inst, Pvalue)
fwrite(inst[, .(rsid, Chr, Pos, EA, NEA, EAF, Beta, SE, Pvalue, Neff)],
       file.path(opt$out_dir, "t2d_diagram_independent_instruments.tsv"), sep = "\t")
message("Retained ", nrow(inst), " independent T2D instruments.")

harmonise <- function(exposure, outcome, has_outcome_eaf) {
  d <- merge(exposure[, .(SNP = rsid, ea_x = EA, oa_x = NEA, eaf_x = EAF, bx = Beta, bxse = SE, px = Pvalue)], outcome, by = "SNP")
  if (!nrow(d)) return(list(data = d, exclusions = data.table(reason = "no_matching_rsids", n = 0L)))
  d[, `:=`(ea_y = toupper(ea_y), oa_y = toupper(oa_y), ea_x = toupper(ea_x), oa_x = toupper(oa_x))]
  direct <- d$ea_x == d$ea_y & d$oa_x == d$oa_y
  reverse <- d$ea_x == d$oa_y & d$oa_x == d$ea_y
  strand_direct <- d$ea_x == complement(d$ea_y) & d$oa_x == complement(d$oa_y)
  strand_reverse <- d$ea_x == complement(d$oa_y) & d$oa_x == complement(d$ea_y)
  d[, alignment := fifelse(direct | strand_direct, 1L, fifelse(reverse | strand_reverse, -1L, NA_integer_))]
  d[, palindromic := is_palindromic(ea_x, oa_x)]
  d[, ambiguous_palindrome := FALSE]
  if (has_outcome_eaf) {
    d[palindromic == TRUE & (!is.finite(eaf_x) | !is.finite(eaf_y) | pmin(eaf_x, 1 - eaf_x, eaf_y, 1 - eaf_y) > 0.42), ambiguous_palindrome := TRUE]
  } else {
    d[palindromic == TRUE, ambiguous_palindrome := TRUE]
  }
  d[, by := alignment * by]
  excluded <- rbindlist(list(
    data.table(reason = "unmatched_alleles", n = sum(is.na(d$alignment))),
    data.table(reason = "ambiguous_palindromes", n = sum(d$ambiguous_palindrome)),
    data.table(reason = "invalid_outcome_estimate", n = sum(!is.finite(d$by) | !is.finite(d$byse) | d$byse <= 0))
  ))
  keep <- !is.na(d$alignment) & !d$ambiguous_palindrome & is.finite(d$by) & is.finite(d$byse) & d$byse > 0
  list(data = d[keep], exclusions = excluded)
}

ivw_result <- function(d, outcome_id, family) {
  n <- nrow(d)
  if (!n) return(data.table(outcome_id = outcome_id, outcome_family = family, method = "none", nsnp = 0L))
  if (n == 1L) {
    b <- d$by / d$bx; se <- d$byse / abs(d$bx)
    return(data.table(outcome_id = outcome_id, outcome_family = family, method = "Wald ratio", nsnp = 1L,
                      b = b, se = se, pval = 2 * pnorm(-abs(b / se)), Q = NA_real_, Q_pval = NA_real_))
  }
  w <- 1 / d$byse^2
  b <- sum(w * d$bx * d$by) / sum(w * d$bx^2)
  se_fe <- sqrt(1 / sum(w * d$bx^2))
  q <- sum(w * (d$by - b * d$bx)^2)
  q_df <- n - 1L
  # Multiplicative random-effects IVW, as used by default in contemporary MR workflows.
  se_mre <- se_fe * sqrt(max(1, q / q_df))
  data.table(outcome_id = outcome_id, outcome_family = family, method = "IVW multiplicative random effects",
             nsnp = n, b = b, se = se_mre, pval = 2 * pnorm(-abs(b / se_mre)), Q = q, Q_pval = pchisq(q, q_df, lower.tail = FALSE))
}

egger_result <- function(d, outcome_id, family) {
  if (nrow(d) < 3L) return(NULL)
  fit <- tryCatch(lm(by ~ bx, weights = 1 / byse^2, data = d), error = function(e) NULL)
  if (is.null(fit) || nrow(coef(summary(fit))) < 2L) return(NULL)
  cf <- coef(summary(fit))
  data.table(outcome_id = outcome_id, outcome_family = family, method = "MR-Egger",
             nsnp = nrow(d), b = cf["bx", "Estimate"], se = cf["bx", "Std. Error"], pval = cf["bx", "Pr(>|t|)"],
             egger_intercept = cf["(Intercept)", "Estimate"], egger_intercept_pval = cf["(Intercept)", "Pr(>|t|)"])
}

read_raw_outcome <- function(path) {
  required <- c("rsid", "a0", "a1", "Beta", "SE", "P")
  x <- fread(path, select = required)
  if (!all(required %in% names(x))) stop("Unexpected raw BCAC columns in ", basename(path))
  x[, .(SNP = rsid, ea_y = a1, oa_y = a0, eaf_y = NA_real_, by = Beta, byse = SE, py = P)]
}

read_slopehunter_outcome <- function(path) {
  required <- c("SNP", "EA.prognosis", "OA.prognosis", "EAF.prognosis", "ybeta.adj", "yse.adj", "yp.adj")
  x <- fread(path, select = required)
  if (!all(required %in% names(x))) stop("Unexpected Slope-Hunter columns in ", basename(path))
  # Use the corrected prognosis association estimates, rather than ybeta/yse.
  x[, .(SNP, ea_y = `EA.prognosis`, oa_y = `OA.prognosis`, eaf_y = as.numeric(`EAF.prognosis`),
        by = `ybeta.adj`, byse = `yse.adj`, py = `yp.adj`)]
}

run_family <- function(files, family, reader, has_outcome_eaf) {
  result_list <- list(); exclusion_list <- list(); harmonised_dir <- file.path(opt$out_dir, "harmonised", family)
  dir.create(harmonised_dir, recursive = TRUE, showWarnings = FALSE)
  for (path in files) {
    outcome_id <- sub("\\.slopehunter\\.tsv\\.gz$|\\.txt$", "", basename(path))
    message("[", family, "] ", outcome_id)
    out <- tryCatch(reader(path), error = function(e) e)
    if (inherits(out, "error")) {
      exclusion_list[[length(exclusion_list) + 1L]] <- data.table(outcome_id, outcome_family = family, reason = conditionMessage(out), n = NA_integer_)
      next
    }
    h <- harmonise(inst, out, has_outcome_eaf)
    ex <- copy(h$exclusions); ex[, `:=`(outcome_id = outcome_id, outcome_family = family)]
    exclusion_list[[length(exclusion_list) + 1L]] <- ex
    fwrite(h$data, file.path(harmonised_dir, paste0(outcome_id, "_harmonised.tsv")), sep = "\t")
    result_list[[length(result_list) + 1L]] <- ivw_result(h$data, outcome_id, family)
    eg <- egger_result(h$data, outcome_id, family)
    if (!is.null(eg)) result_list[[length(result_list) + 1L]] <- eg
  }
  list(results = rbindlist(result_list, fill = TRUE), exclusions = rbindlist(exclusion_list, fill = TRUE))
}

raw_files <- list.files(opt$raw_dir, pattern = "\\.txt$", full.names = TRUE)
slope_files <- list.files(opt$slopehunter_dir, pattern = "\\.slopehunter\\.tsv\\.gz$", full.names = TRUE)
if (!length(raw_files) || !length(slope_files)) stop("No raw and/or Slope-Hunter outcome files found.", call. = FALSE)
raw_run <- run_family(raw_files, "raw_bcac_2021", read_raw_outcome, FALSE)
slope_run <- run_family(slope_files, "slopehunter_adjusted", read_slopehunter_outcome, TRUE)
results <- rbindlist(list(raw_run$results, slope_run$results), fill = TRUE)
results[, `:=`(HR = exp(b), HR_lci95 = exp(b - 1.96 * se), HR_uci95 = exp(b + 1.96 * se))]
# Correct within outcome family and MR method, avoiding a combined correction
# across raw and Slope-Hunter versions of the same endpoints.
results[, qval_bh := p.adjust(pval, method = "BH"), by = .(outcome_family, method)]
fwrite(results, file.path(opt$out_dir, "t2d_survival_mr_results.tsv"), sep = "\t")
fwrite(rbindlist(list(raw_run$exclusions, slope_run$exclusions), fill = TRUE),
       file.path(opt$out_dir, "t2d_survival_harmonisation_exclusions.tsv"), sep = "\t")
writeLines(c(
  "Exposure: Mahajan et al. 2018b DIAGRAM European T2D GWAS, unadjusted for BMI.",
  paste("Independent instruments:", nrow(inst)),
  paste("Raw BCAC outcomes:", length(raw_files)),
  paste("Slope-Hunter outcomes:", length(slope_files)),
  "Slope-Hunter MR uses ybeta.adj, yse.adj, and yp.adj."
), file.path(opt$out_dir, "README.txt"))
message("Done. Results written to: ", opt$out_dir)
