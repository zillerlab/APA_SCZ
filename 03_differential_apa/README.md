# 03 Differential APA

Differential 3'UTR alternative polyadenylation (dAPA) using the two-PAS reference of
module 02.

| Script | Content | Manuscript |
|---|---|---|
| `01_dapa_iN.R` | SCZ, BD and MDD vs HC in iNs | Fig. 2c; Ext. Data Fig. 6b; Suppl. Table 3 |
| `02_dapa_postmortem.R` | SCZ vs controls, postmortem DLPFC (CommonMind) | Fig. 2e-g; Suppl. Table 3 |
| `03_dapa_rbp_knockdown.R` | PTBP2, PCBP2 and CPSF6 knockdown vs scramble control | Fig. 5e-h,j; Suppl. Table 3 |
| `04_fig5ef_rbp_knockdown.R` | Overlap of knockdown dAPA sets; knockdown vs cohort effect; overlap with SCZ dAPA genes | Fig. 5e,f |
| `05_fig5h_rbp_variance_explained.R` | APA variance explained by RBP expression; run twice: `--genes all` (targets vs other genes, Fig. 5h left) and default (SCZ dAPA genes, unique contribution per RBP, Fig. 5h right) | Fig. 5h |
| `06_fig5g_syngo_rbp_knockdown.R` | SynGO enrichment of knockdown dAPA genes | Fig. 5g |

## Run

```bash
Rscript 03_differential_apa/01_dapa_iN.R --cores 8
Rscript 03_differential_apa/02_dapa_postmortem.R --cmc <dir> --cores 8   # controlled access
Rscript 03_differential_apa/03_dapa_rbp_knockdown.R --cores 8
Rscript 03_differential_apa/04_fig5ef_rbp_knockdown.R
Rscript 03_differential_apa/05_fig5h_rbp_variance_explained.R
Rscript 03_differential_apa/05_fig5h_rbp_variance_explained.R --genes all
Rscript 03_differential_apa/06_fig5g_syngo_rbp_knockdown.R
```

Scripts 04 and 05 use the outputs of 01, 03 and `05_differential_expression/01`.
Script 02 needs the CommonMind RNA-seq and covariates (Synapse syn2759792, controlled
access); the required files are listed in the script header.

## Inputs

| File | Content |
|---|---|
| `resources/dapa/utr_counts_iN_cohort.txt` | Read counts on proximal and distal UTR regions, iN cohort |
| `resources/rnaseq/rnaseq_iN_cohort_metadata.txt` | Library metadata (diagnosis, sex, site, size factor, technical PCs) |
| `resources/rbp_knockdown/utr_counts_rbp_knockdown.txt.gz` | UTR read counts and segment lengths, 20 knockdown and control libraries |
| `resources/rbp_knockdown/rbp_knockdown_samples.txt` | Knockdown sample table: donor, condition, batch, knockdown arm(s), size factor, knockdown efficiency |
| `resources/rbp_knockdown/syngo/syngo_1.3_dAPA_<RBP>_bh1e-4.zip` | SynGO portal output (release 1.3, default settings, brain-expressed background) for `results/dapa/dAPA_RBP_KD_genes_<RBP>.txt` |

## Outputs

| File | Content |
|---|---|
| `results/dapa/dAPA_iN_<SCZ,BD,MDD>_vs_Ctrl.txt` | dAPA per gene, iN cohort |
| `results/dapa/apa_ratio_iN_cohort.txt` | Distal usage ratio per gene and library |
| `results/dapa/dAPA_DLPFC_SCZ_vs_Control*.txt` | dAPA per gene, DLPFC |
| `results/dapa/dAPA_RBP_KD_<PTBP2,PCBP2,CPSF6>.txt` | dAPA per gene, knockdowns |
| `results/dapa/dAPA_RBP_KD_overlap.txt` | Knockdown dAPA overlaps |
| `results/dapa/dAPA_RBP_KD_genes_<RBP>.txt`, `dAPA_RBP_KD_genes_background.txt` | Knockdown dAPA gene lists and tested genes |
| `results/figures/Fig5{e,f,g,h}_*` | Fig. 5e-h and source data |
| `results/figures/Fig5_dAPA_KD_vs_cohort_*` | Overlap of knockdown dAPA and SCZ dAPA genes |

Columns of the cohort and DLPFC dAPA tables: `ID`, `method`, `n_ref`, `n_comp`, `beta`,
`se` (condition effect, logit scale), `p_wald`, `p_lrt`, `delta`, `mean_expr`, `qvalue`,
`padj_BH`. Columns of the knockdown dAPA tables: `ID`, `beta`, `se` (condition effect,
logit scale), `p_wald`, `padj_BH`, `dAPA`, `R_ctrl`, `R_kd` (fitted distal usage per
group), `deltaR`, `p_lrt`, `phi`, `n_ref`, `n_comp`, `meanExpr`.

## Method

Distal usage ratio R = distal / (distal + proximal) normalised counts. Beta regression
(logit link) per gene: iNs `R ~ condition + sex + TechnicalPC1` (dAPA: q <= 0.05);
DLPFC with institution, sex, age, PMI, RIN, RIN², library cluster, ancestry PCs 1-5 and
technical PCs 1-3 (q <= 0.1). q-values: Storey.

Knockdowns (3-5 knockdown and 5-7 control libraries per arm): beta regression
`R ~ condition + donor` on length-normalised read densities (`R/dapa_v2_functions.R`),
Wald test, Benjamini-Hochberg adjustment; dAPA: adjusted P <= 1e-4. With this number of
libraries the Wald test is anticonservative (the precision is estimated per gene from few
residual degrees of freedom), so the cutoff is stricter than for the cohort.

## Dependencies

R >= 4.0; betareg, ggplot2, readxl (CRAN); qvalue, limma, DESeq2 (Bioconductor).
