# 08 Chromatin accessibility (ATAC-seq)

Chromatin accessibility in iNs across the cohort.

| Script | Content | Manuscript |
|---|---|---|
| `01_diff_atac_iN_cohort.R` | PCA; differentially accessible peaks, SCZ vs HC | Fig. 1c; Ext. Data Fig. 4e |

## Run

```bash
Rscript 08_chromatin_accessibility/01_diff_atac_iN_cohort.R
```

## Inputs (`resources/atac/`)

| File | Content |
|---|---|
| `atac_peak_counts_iN_cohort.txt.gz` | Fragments per peak of the union peak set (hg19), one column per library |
| `atac_iN_cohort_metadata.txt` | Line ID, diagnosis, sex, site, FRiP, TSS enrichment |

## Outputs

| File | Content |
|---|---|
| `results/atac/diffATAC_iN_SCZ_vs_Ctrl.txt` | Per peak: logFC (SCZ / HC), logCPM, F, p, FDR, q |
| `results/atac/RUV_W1.txt` | Unwanted-variation factor per library |
| `results/figures/Fig1c_ATAC_PCA.*`, `ExtFig4e_diffATAC_volcano.*` | Fig. 1c; Ext. Data Fig. 4e |

## Method

edgeR (upper-quartile normalisation, quasi-likelihood F-test),
`~ W_1 + batch + FRiP + diagnosis`, with W_1 from RUVg (k = 1; negative controls: peaks
outside the 10,000 most significant of an initial fit). Differential: Storey q <= 0.1.

## Dependencies

R >= 4.0; edgeR, qvalue (Bioconductor), ggplot2 (CRAN).
