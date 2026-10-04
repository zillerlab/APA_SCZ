# 11 GRIA1 3'UTR qPCR

GRIA1 mRNA levels in GRIA1 long 3'UTR KO iNs, and GFP reporter mRNA levels with the GRIA1
long vs short 3'UTR in HEK293T cells.

| Script | Content | Manuscript |
|---|---|---|
| `01_gria1_utr_reporter_qpcr.R` | qPCR of GFP-GRIA1 3'UTR reporters, long vs short | Fig. 7h |
| `02_gria1_3utr_ko_qpcr.R` | qPCR of GRIA1 in GRIA1 long 3'UTR KO vs WT iNs | Fig. 7c |

## Run

```bash
Rscript 11_gria1_reporter/01_gria1_utr_reporter_qpcr.R
Rscript 11_gria1_reporter/02_gria1_3utr_ko_qpcr.R
```

## Inputs (`resources/gria1_reporter/`)

| File | Content |
|---|---|
| `gria1_utr_reporter_qpcr.txt` | One row per transfection: batch, experiment (long and short 3'UTR transfected in parallel), 3'UTR, plasmid, GFP reporter expression relative to RTF2 adjusted for transfection efficiency, and its propagated SD |
| `gria1_3utr_ko_qpcr.txt` | One row per qPCR reaction: differentiation batch, line ID, genotype, sample, technical replicate, GRIA1 expression relative to EID2 |

## Outputs (`results/figures/`)

`Fig7h_GRIA1_reporter.pdf`, `Fig7h_GRIA1_reporter.txt` (source data),
`Fig7h_GRIA1_reporter_test.txt`; `Fig7c_qPCR_GRIA1.pdf`, `Fig7c_qPCR_GRIA1.txt` (source data),
`Fig7c_qPCR_GRIA1_test.txt`.

## Method

- Reporter: linear mixed model `rel_expr_adj ~ utr + (1 | experiment) + (1 | batch)` (lme4),
  two-sided Satterthwaite t-test of the 3'UTR coefficient (lmerTest). Figure values are
  normalised to the mean of the long 3'UTR transfections per batch (display only).
- GRIA1 3'UTR KO: technical replicates averaged per sample on the log2 scale; linear mixed
  model `log2(expression) ~ genotype + (1 | batch)` on the sample means, two-sided
  Satterthwaite t-test. Figure values are batch-normalised sample means (per batch, the
  mean of the WT and KO means on the log2 scale is removed; display only).

## Dependencies

R >= 4.0; lme4, lmerTest, ggplot2.
