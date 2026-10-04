# 01 EB scorecard

Germ-layer and pluripotency signature scores of embryoid bodies (EBs) from each iPSC
line, tested against undifferentiated iPSC reference samples.

**Manuscript:** Ext. Data Fig. 1g,h.

## Run

```bash
Rscript 01_eb_scorecard/EB_score_pipeline.R
```

## Inputs (`resources/eb_scorecard/`)

| File | Content |
|---|---|
| `EB_counts_clean.txt` | 3' RNA-seq counts, 20 reference iPSC and 153 EB samples |
| `EB_metadata_clean.txt` | Sample, line ID, condition (iPSC / EB), sequencing batch, reference flag |
| `EB_marker_panel.txt` | Marker genes (ectoderm, mesoderm, endoderm, pluripotency) with weights |
| `EB_marker_selection_audit.txt` | Candidate markers and selection criteria |

## Outputs (`results/eb_scorecard/`)

| File | Content |
|---|---|
| `EB_scores_per_sample.txt` | Score, t, p and FDR per signature and sample |
| `EB_scores_per_line.txt` | Pass / fail per line |
| `ExtDataFig1gh_EB_scores.txt` | Source data of Ext. Data Fig. 1g,h |

## Method

log2(CPM + 1); per-gene robust z-score against the reference iPSCs; weighted mean per
signature; one-sided Crawford-Howell single-case t-test against the reference scores
(reference samples scored leave-one-out); Benjamini-Hochberg per signature. An EB passes
if all germ-layer scores are increased and the pluripotency score is decreased
(FDR < 0.05); a line passes if at least one EB passes.

## Dependencies

R >= 4.0 (base R only).
