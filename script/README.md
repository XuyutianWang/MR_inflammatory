MR Analysis Workflow
This document outlines the step-by-step pipeline for preparing exposure and outcome data for Mendelian Randomization (MR) analysis using UKB PPP data.

1. Data Download
Retrieve the required biomarker and GWAS summary data.

Step	Description	Script Path
1.1	Retrieve UniProt IDs for prioritized biomarkers.	MR/exposure/download_UKBPPP/find_uniprot_ids.py
1.2	Download GWAS summary statistics from Synapse.	MR/exposure/download_UKBPPP/synapse_download.py
1.3 Unzip the tar files

2. Concatenate Chromosome Files
For each protein (e.g., ADIPOQ), merge the per-chromosome GWAS summary files (chr1–chr22) into a single complete dataset.

Script: MR/script/concatenate_ukbppp.py

3. Preprocess Exposure Data
Prepare the exposure GWAS data for instrumental variable selection.

Step	Description	Script Path
3.1	Filter instruments using genome-wide significance threshold (LOG10P > 7.30103).	MR/script/filterp.R
3.2	Convert LOG10P values to standard P-values.	MR/script/p_value_log10p.py
3.3	Match rsIDs using SNPlocs.Hsapiens.dbSNP155.GRCh38 reference.	MR/script/exposure_rsid.R

3.4 seperating cis and trans pqtls.
labelling cis and trans: script/label_pqtl_cis_trans.py
remove trans pqtls: script/filter_cis_lead_snps.R

3.5 ld clumping
script/run_clump.sh

3.6 remove weak instruments


4. Preprocess Survival Outcome Data
Prepare the outcome GWAS data for harmonization with exposure data.

Step	Description	Script Path
4.1	Remove chromosome 23 (sex chromosome) from the outcome data.	MR/script/remove_23.R
4.2	Match rsIDs using SNPlocs.Hsapiens.dbSNP155.GRCh37 reference.	MR/script/outcome_rsid.R


5. Preprocess incidence data

Step    Seperating different subtypes
5.1  In icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt, the subtype data are mixed 
SNP--var_name
beta--*_log_or_meta
se--*_se_meta
effect_allele--Effect.Meta
other_allele--Baseline.Meta
pval--Calculated using β and SE 
chr--chr.iCOGs
pos--Position.iCOGs

*The p-value for the subtypes are not provided, "FTOP_p_value_meta" and "MTOP_p_value_meta" are p-values for different statistcal methods. So we calculated the p-value using beta and SE.
  
5.2  Matching rsid for incidence data  MR/script/incidence_rsid.R

5.3  remove chr 23 
MR/script/remove_23_incidence.R

6. run mr

Step    harmonisation
6.1 


7. cis and trans 
[protein_gene_map.csv](/home/df24679/projects/MR/exposure_all/protein_gene_map.csv) contains all 120 exposure IDs with gene symbol, GRCh38 chromosome, and TSS.
[lead_snps_cis_trans](/home/df24679/projects/MR/exposure_all/lead_snps_cis_trans) contains 120 copied exposure files, each with a new qtl_type column (cis or trans).
[label_pqtl_cis_trans.py](/home/df24679/projects/MR/script/label_pqtl_cis_trans.py) reproduces the process.



Dependencies
R (with relevant packages)

Python 3

Synapse client for data download

Bioconductor packages: SNPlocs.Hsapiens.dbSNP155.GRCh38, SNPlocs.Hsapiens.dbSNP155.GRCh37