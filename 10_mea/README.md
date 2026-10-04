# 10 Multielectrode array (MEA)

Neuronal network activity of iNs (Axion well-averaged parameters).

| Script | Content | Manuscript |
|---|---|---|
| `01_mea_rbp_knockdown.R` | CPSF6, PCBP2 and PTBP2 knockdown vs scramble control | Fig. 6f,g |
| `02_mea_gria1_3utr.R` | GRIA1 long 3'UTR deletion vs WT | Fig. 7d |

## Run

```bash
Rscript 10_mea/01_mea_rbp_knockdown.R
Rscript 10_mea/02_mea_gria1_3utr.R
```

## Inputs (`resources/mea/`)

| File | Content |
|---|---|
| `mea_rbp_knockdown_wells.txt` | One row per well: plate, batch, line ID, condition, MEA parameters |
| `mea_gria1_3utr_wells.txt` | One row per well: plate, differentiation batch, line ID, genotype, clone, MEA parameters |

## Outputs

| File | Content |
|---|---|
| `results/mea/MEA_RBP_KD_LMM.txt` | Per parameter and knockdown: estimate, SE, t, p, q |
| `results/mea/MEA_GRIA1_3UTR_LM.txt` | Per parameter: estimate (KO - WT), SE, t, p, q |
| `results/figures/Fig6f_*`, `Fig6g_*`, `Fig7d_*` | Figures and source data |

## Method

- RBP knockdown: wells with >= 1 spike and >= 3 active electrodes; 1.5 x IQR outlier
  removal per parameter; `value ~ condition + (1 | line) + (1 | batch)` (lme4,
  Satterthwaite t-tests); Benjamini-Hochberg over parameters per knockdown.
- GRIA1 3'UTR KO: wells with >= 3 active electrodes and electrode resistance >= 25 kOhm;
  1.5 x IQR outlier removal on batch-scaled values; `value ~ batch + electrode
  resistance + genotype`; Benjamini-Hochberg over parameters.

## Dependencies

R >= 4.0; lme4, lmerTest, ggplot2 (CRAN); limma (Bioconductor).
