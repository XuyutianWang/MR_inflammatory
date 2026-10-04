# Tumour-expression survival analysis

## Inputs

- Biomarker list: `exposure_all/download_UKBPPP/inflammatory_biomarkers_uniprot_ids.csv`
- Use `matched_gene` as the tumour-expression gene for the primary analysis (206 biomarkers currently have a value).
- Before analysis, inspect all `parent_biomarker`, `biomarker`, `uniprot_id`, `gene_names`, and `matched_gene` mappings. Save approved mappings as `data/candidate_gene_map.csv`; record excluded biomarkers and reasons.
- Do not use a tumour-gene proxy for biomarkers that are not meaningfully represented by tumour expression (for example composite blood-cell ratios). Flag these as “not tested in tumour data”.

## 1. Set up

```text
third-method/
  data/raw/{tcga,scanb,metabric}/
  data/derived/
  scripts/
  results/
  logs/
```

Save source URL, download date, release/version, manifest/query, and checksum for every download. Do not commit raw patient-level data.

Pre-specify before viewing results:

- Primary model: Cox model with expression as a continuous, within-cohort z-score.
- Discovery cohort: TCGA-BRCA.
- Validation cohorts: SCAN-B, then METABRIC if usable.
- Primary outcome: TCGA progression-free interval (PFI); SCAN-B/METABRIC recurrence/relapse-free endpoint as supplied.
- Secondary outcome: overall survival (OS).
- Multiplicity: Benjamini-Hochberg FDR across the final approved biomarker list.

## 2. Download TCGA-BRCA (discovery)

### Portal download

1. Open the [GDC Data Portal](https://portal.gdc.cancer.gov/).
2. Filter: project `TCGA-BRCA`; data category `Transcriptome Profiling`; data type `Gene Expression Quantification`; workflow `STAR - Counts`; sample type `Primary Tumor`; access `open`.
3. Add files to cart and download the manifest.
4. Download with the GDC Data Transfer Tool:

```bash
gdc-client download -m manifest.txt -d third-method/data/raw/tcga
```

5. From the same project, download open clinical and biospecimen annotations.
6. Download the TCGA Clinical Data Resource (TCGA-CDR) supplemental endpoint table. Use its PFI time/event variables; do not recreate PFI from raw clinical fields.

## 3. Download SCAN-B (primary validation)

1. Download processed RNA-seq expression from [GEO GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) into `data/raw/scanb/`.
2. Download the matching SCAN-B clinical/follow-up file from the release linked by the GEO record/publication.
3. Confirm that expression and clinical files use the same patient/sample IDs. Create and save an ID-mapping table.
4. Record endpoint definition, time origin, units, censoring rule, and covariate coding from the data dictionary.

If using the full SCAN-B release instead, download the expression and clinical/follow-up files from the same release only; do not mix SCAN-B releases.

## 4. Download METABRIC (optional validation)

1. Open [METABRIC in cBioPortal](https://www.cbioportal.org/study?id=brca_metabric).
2. Select **Download** and save clinical and mRNA-expression files in `data/raw/metabric/`.
3. Record the cBioPortal study name, platform, download date, endpoint definitions, and sample/patient IDs.

## 5. Build each analysis cohort

For every cohort:

1. Keep one primary-tumour specimen per invasive breast-cancer patient; exclude normal-adjacent, metastatic/recurrent, and non-tumour samples.
2. If multiple primary samples exist, use one deterministic selection rule and save selected IDs.
3. Join expression and clinical files by stable patient ID; exclude invalid/missing outcome data and check zero-day events.
4. Save a flow table: raw samples → unique patients → matched patients → outcome eligible → final N/events.
5. Extract approved `matched_gene` values; check missing/duplicated genes and distributions.

## 6. Prepare variables

### Expression

- TCGA: normalise STAR counts with edgeR TMM then use log2-CPM (with a prior count), or use another documented count-appropriate transformation.
- SCAN-B/METABRIC: use the supplied processed expression values; do not apply TCGA count normalisation.
- Within each cohort, standardise each analysed gene to mean 0 and SD 1.
- Analyse cohorts separately. Never pool raw RNA-seq and microarray expression values.

### Outcomes

- TCGA primary: TCGA-CDR PFI.
- SCAN-B/METABRIC primary: source-defined RFS/relapse-free survival. Write the exact source definition into `data/derived/endpoint_definitions.md`.
- Secondary: OS.

### Covariates

Create harmonised cohort-specific fields for age at diagnosis, stage, nodal status, grade, ER status, HER2 status, systemic therapy, radiotherapy, and surgery where available. Save a coding table.

Use:

- Model 0: gene expression only.
- Model 1 (primary): expression + age + stage + nodal status + grade + ER + HER2.
- Model 2: Model 1 + treatment variables available in the cohort.

For Model 1, use complete cases and report missingness, final N, and event count. If covariate missingness is material, add a sensitivity analysis using multiple imputation of covariates only.

## 7. Run survival models

Use R packages `survival`, `broom`, and `metafor`; save `sessionInfo()`.

For each gene and each endpoint, fit the model below with the cohort-specific field names:

```r
fit <- coxph(
  Surv(time, event) ~ expression_z + age + stage + nodes + grade + ER + HER2,
  data = analysis_data
)
```

For every gene, save:

- HR per 1-SD higher expression, 95% CI, P value, FDR-adjusted P value;
- N analysed, number of events, and median follow-up;
- Model 0, Model 1, and Model 2 estimates;
- output from `cox.zph(fit)`.

Checks:

1. Check proportional hazards with `cox.zph(fit)`; use a time-varying effect or documented alternative when violated.
2. Check non-linearity with restricted cubic splines when events allow.
3. Use Kaplan-Meier plots only for presentation; keep the continuous Cox model primary.
4. Run subtype analyses only if event counts allow and test interaction first.

## 8. Validate and synthesise

1. Run the same gene list, transformations, model sequence, and FDR procedure in SCAN-B.
2. Repeat in METABRIC if downloaded.
3. Keep cohort-specific models separate. Do not select genes based on TCGA P values before validation.
4. For the same endpoint definition in at least two cohorts, meta-analyse log(HR) and standard error with a random-effects model. Report `I²`, `tau²`, and a forest plot.
5. If endpoints differ (e.g. PFI versus RFS), show results side by side and do not pool them as one survival endpoint.

## 9. Report results

For each biomarker, report biomarker name, UniProt ID, approved gene mapping, cohort, outcome, model, HR (95% CI), P/FDR, N, events, and endpoint definition. In the triangulation table, label this evidence “tumour-expression prognosis”; do not call it circulating-protein evidence or make causal/treatment-effect claims.
