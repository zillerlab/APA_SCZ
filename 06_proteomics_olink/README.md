# 06 Proteomics (Olink)

Protein abundance (Olink Explore HT) after RBP knockdown in iNs.

| Script | Content | Manuscript |
|---|---|---|
| `01_dep_olink_rbp_knockdown.R` | PTBP2 and PCBP2 knockdown vs control; permutation test | Fig. 5i |
| `02_fig5j_dapa_vs_protein.R` | APA change vs protein change upon knockdown | Fig. 5j |

## Run

```bash
Rscript 06_proteomics_olink/01_dep_olink_rbp_knockdown.R --cores 8
Rscript 06_proteomics_olink/02_fig5j_dapa_vs_protein.R   # uses 03_differential_apa/03
```

## Inputs (`resources/olink/`)

| File | Content |
|---|---|
| `olink_npx_rbp_knockdown.txt.gz` | NPX values per protein and sample; `in_APA_set` marks the test set |
| `olink_samples_rbp_knockdown.txt` | Sample, line ID, condition, batch, QC |

## Outputs

| File | Content |
|---|---|
| `results/proteomics/DEP_OLINK_RBP_KD.txt` | Per protein and knockdown: logFC, t, p, q |
| `results/proteomics/DEP_OLINK_permutation.txt` | Permutation test |
| `results/figures/Fig5i_*`, `Fig5j_*` | Fig. 5i,j and source data |

## Method

limma, `NPX ~ condition + batch + mean sample NPX`; Benjamini-Hochberg over the proteins
encoded by genes of the 3'UTR reference (936). Number of proteins with p <= 0.05 compared
with all within-batch label permutations. Fig. 5j: Spearman correlation between knockdown
dAPA effect (condition coefficient of the beta regression; dAPA genes at BH <= 1e-4)
and protein logFC (MAD outliers removed).

## Dependencies

R >= 4.0; limma (Bioconductor), ggplot2 (CRAN).
