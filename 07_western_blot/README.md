# 07 Western blot

Quantification of the RBP knockdowns and of GRIA1 in GRIA1 long 3'UTR KO iNs.

| Script | Content | Manuscript |
|---|---|---|
| `01_wb_rbp_knockdown.R` | Knockdown vs scramble control | Fig. 5c |
| `02_wb_gria1_3utr.R` | GRIA1 in GRIA1 long 3'UTR KO vs WT | Fig. 7e |

## Run

```bash
Rscript 07_western_blot/01_wb_rbp_knockdown.R
Rscript 07_western_blot/02_wb_gria1_3utr.R
```

## Input

`resources/western_blot/wb_rbp_knockdown.txt`: one row per lane; protein, line ID,
batch, condition (CTRL / KD), NRI (target / neurofilament L peak area).

`resources/western_blot/wb_gria1_3utr.txt`: one row per capillary; differentiation batch,
replicate (independent wells and extraction), blot date, line ID, sample (WT, KO clone),
condition (WT / KO), NRI (GRIA1 120 + 150 kDa / neurofilament L peak area).

## Outputs (`results/figures/`)

`Fig5c_WB_knockdown.pdf`, `Fig5c_WB_knockdown.txt` (source data),
`Fig5c_WB_knockdown_tests.txt`; `Fig7e_WB_GRIA1.pdf`, `Fig7e_WB_GRIA1.txt` (source data),
`Fig7e_WB_GRIA1_test.txt`.

## Method

Paired test: one-sample t-test of log2(NRI_KD / NRI_CTRL) per knockdown-control pair
(same line and batch), one-sided (knockdown lower than control). Figure values are
NRI divided by the geometric mean of the pair (display only).

GRIA1: both KO clones pooled; linear mixed model `log2(NRI) ~ condition + (1 | batch)`
(lme4), two-sided Satterthwaite t-test (lmerTest).

## Dependencies

R >= 4.0; ggplot2; lme4, lmerTest (Fig. 7e).
