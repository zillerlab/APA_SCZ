# 05 Differential expression

DESeq2 analyses of gene and microRNA expression, and power of the cohort analysis.

| Script | Content | Manuscript |
|---|---|---|
| `01_deseq_iN_cohort.R` | SCZ, BD and MDD vs HC, iN cohort | Fig. 1d; Ext. Data Fig. 4f,g; Suppl. Table 3 |
| `02_deseq_rbp_knockdown.R` | PTBP2, PCBP2 and CPSF6 knockdown vs scramble control | Fig. 5d; Suppl. Table 3 |
| `03_deseq_gria1_3utr.R` | GRIA1 long 3'UTR deletion vs WT | Ext. Data Fig. 7d |
| `04_fig5d_rbp_knockdown.R` | Overlap of knockdown DEG sets | Fig. 5d |
| `05_deseq_microrna_iN_cohort.R` | microRNA expression, SCZ vs HC; microRNA PCA | Fig. 1c,e; Suppl. Table 3 |
| `06_power_iN_cohort.R` | Power: model-based and donor resampling | Ext. Data Fig. 4c,d |

## Run

```bash
Rscript 05_differential_expression/01_deseq_iN_cohort.R
Rscript 05_differential_expression/02_deseq_rbp_knockdown.R
Rscript 05_differential_expression/03_deseq_gria1_3utr.R
Rscript 05_differential_expression/04_fig5d_rbp_knockdown.R
Rscript 05_differential_expression/05_deseq_microrna_iN_cohort.R
Rscript 05_differential_expression/06_power_iN_cohort.R --cores 4   # uses the output of 01
```

## Inputs

| File | Content |
|---|---|
| `resources/rnaseq/gene_counts_iN_cohort.txt.gz` | Gene counts (GENCODE v27), iN cohort, one library per line |
| `resources/rnaseq/rnaseq_iN_cohort_metadata.txt` | Diagnosis, sex, differentiation site, technical covariates |
| `resources/rnaseq/gene_annotation_ensembl.txt` | Ensembl gene ID and symbol |
| `resources/rbp_knockdown/gene_counts_rbp_knockdown.txt.gz` | Gene counts, 20 knockdown and control libraries |
| `resources/rbp_knockdown/rbp_knockdown_samples.txt` | Donor, condition, batch, libraries used per knockdown |
| `resources/gria1_3utr/gene_counts_gria1_3utr.txt.gz`, `gria1_3utr_samples.txt` | Gene counts and samples, GRIA1 3'UTR KO and WT |
| `resources/microrna/mature_counts_iN_cohort.txt.gz` | Mature microRNA counts (miRBase v22) |
| `resources/microrna/microrna_iN_cohort_metadata.txt` | Line ID, diagnosis, sex, site, library used |

## Outputs

| File | Content |
|---|---|
| `results/deg/DEG_iN_<SCZ,BD,MDD>_vs_Ctrl.txt`, `DEG_iN_overlap.txt` | iN cohort DE |
| `results/deg/DEG_RBP_KD_<PTBP2,PCBP2,CPSF6>.txt`, `DEG_RBP_KD_overlap.txt` | Knockdown DE |
| `results/deg/DEG_GRIA1_3UTR_KO_vs_WT.txt` | GRIA1 3'UTR KO DE |
| `results/deg/DEmiR_iN_SCZ_vs_Ctrl.txt` | microRNA DE |
| `results/deg/Power_*.txt` | Power analyses |
| `results/figures/` | Fig. 1c,e, 5d; Ext. Data Fig. 4c,d and source data |

In the cohort and microRNA tables, `log2FoldChange` is log2(HC / case).

## Method

| Analysis | Design | Differential |
|---|---|---|
| iN cohort | `~ sex + site + DE_TechnicalPC1 + diagnosis` | padj <= 0.05, \|log2FC\| >= 0.4, baseMean >= 50 |
| RBP knockdown | `~ donor + condition` | padj <= 0.05, \|log2FC\| >= 0.4 |
| GRIA1 3'UTR KO | `~ genotype` | q <= 0.05, \|log2FC\| >= 1, baseMean >= 50 |
| microRNA | `~ sex + site + diagnosis` | padj <= 0.05 |

Power: `RnaSeqSampleSize` (dispersion from the SCZ vs HC DEGs); donor resampling with
limma on variance-stabilised data, 100 draws per group size.

## Dependencies

R >= 4.0; DESeq2, apeglm, qvalue, limma, RnaSeqSampleSize (Bioconductor).
