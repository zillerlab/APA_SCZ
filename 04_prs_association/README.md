# 04 PRS association

Association of polygenic risk scores (SCZ, BD, MDD) with cumulative differential APA.

| Script | Content | Manuscript |
|---|---|---|
| `01_prs_apa_iN.R` | PRS vs cumulative dAPA in iNs; permutation null | Fig. 4b; Ext. Data Fig. 6c-e,g |
| `02_plot_prs_iN.R` | Scatter plots | Fig. 4b; Ext. Data Fig. 6c-e |
| `03_prs_apa_dlpfc.R` | SCZ-PRS vs cumulative dAPA, postmortem DLPFC | Ext. Data Fig. 6i |
| `04_plot_prs_dlpfc.R` | Plots | Ext. Data Fig. 6i |
| `05_prs_extremes_iN.R` | dAPA between top and bottom 10% SCZ-PRS donors | Fig. 4d,e; Ext. Data Fig. 6f,h |

## Data access

PRS are individual-level data and are not included; they are available upon request to the authors and subject to proper ethics approval and DTA. Scripts read them with `--prs <file>` (columns `sample` = line ID, `PRS_SCZ`,
`PRS_BD`, `PRS_MDD`) and write per-donor tables

## Run

```bash
Rscript 04_prs_association/01_prs_apa_iN.R --prs <prs_file> --cores 8
Rscript 04_prs_association/02_plot_prs_iN.R
Rscript 04_prs_association/05_prs_extremes_iN.R --prs <prs_file> --cores 8

# DLPFC (CommonMind, controlled access)
Rscript 03_differential_apa/02_dapa_postmortem.R --cmc <dir> --ethnicity Caucasian \
        --prs <cmc_prs_file> --tag SCZ_vs_Control_EUR_PRS --cores 8
Rscript 03_differential_apa/02_dapa_postmortem.R --cmc <dir> --cases SCZ,BP,AFF \
        --tag allCases_vs_Control --residuals-out <private>/dlpfc_residuals.txt --cores 8
Rscript 04_prs_association/03_prs_apa_dlpfc.R --prs <cmc_prs_file> \
        --residuals <private>/dlpfc_residuals.txt --cmc <dir>
Rscript 04_prs_association/04_plot_prs_dlpfc.R
```

## Outputs (`results/prs/`)

| File | Content |
|---|---|
| `prs_apa_association_iN.txt` | Association statistics, Fig. 4b and Ext. Data Fig. 6c-e |
| `prs_apa_permutation_null_iN.txt` | Permutation null (Ext. Data Fig. 6g) |
| `dAPA_iN_PRShigh_vs_PRSlow.txt` | dAPA, high vs low SCZ-PRS donors |
| `prs_extremes_beta_iN.txt`, `prs_extremes_concordance_iN.txt` | Fig. 4e; Ext. Data Fig. 6f,h |
| `prs_apa_association_dlpfc.txt` | Ext. Data Fig. 6i |

## Method

Covariate-adjusted APA (residuals of `R ~ sex + TechnicalPC1`); per donor, mean |z| over
the dAPA genes; Tukey outlier removal; linear model cumulative APA ~ PRS; empirical p
from 1,000 permutations.

## Dependencies

R >= 4.0; betareg, ggplot2 (CRAN); qvalue (Bioconductor).
