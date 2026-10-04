# 09 PTBP2 CLIP

PTBP2 binding sites in iNs (quick-irCLIP) and their overlap with 3'UTRs showing
differential APA.

| Script | Content | Manuscript |
|---|---|---|
| `01_clip_enrichment.R` | Peaks enriched in PTBP2 IP over input | Fig. 6a; Suppl. Table 1 |
| `02_fig6a_clip_dapa_overlap.R` | dAPA genes vs PTBP2 binding in three CLIP datasets | Fig. 6a |

## Run

```bash
Rscript 09_ptbp2_clip/01_clip_enrichment.R
Rscript 09_ptbp2_clip/02_fig6a_clip_dapa_overlap.R   # uses modules 02 and 03
```

## Inputs (`resources/clip/`)

| File | Content |
|---|---|
| `clip_PTBP2_iN_peak_counts.txt.gz` | Reads per peak (union of MACS3 peaks of the IP libraries), 8 libraries |
| `clip_PTBP2_iN_libraries.txt` | Library, mapped reads, donor sex, fraction (IP / input), replicate |
| `clip_PTBP2_iN_peaks_hg38.bed` | Peak coordinates |
| `eclip_PTBP2_postmortem_PFC_peaks_hg38.bed.gz` | Published PTBP2 eCLIP peaks, postmortem PFC (Dawicki-McKenna et al. 2023) |
| `eclip_PTBP2_iPSC_neurons_peaks_hg38.bed.gz` | Published PTBP2 eCLIP peaks, iPSC-derived neurons (same study) |

## Outputs

| File | Content |
|---|---|
| `results/clip/CLIP_PTBP2_iN_peak_stats.txt` | Tested peaks: logFC (IP / input), p, q, enriched |
| `results/clip/CLIP_PTBP2_iN_enriched_peaks_hg38.bed` | Enriched peaks |
| `results/clip/Fig6a_CLIP_dAPA_overlap.txt` | Fisher's exact tests per dAPA set and CLIP dataset |
| `results/clip/Fig6a_CLIP_binding_<iN,DLPFC>.txt` | Per gene: dAPA and binding flags |
| `results/figures/Fig6a_CLIP_dAPA_overlap.*` | Fig. 6a |

## Method

edgeR with total mapped reads as library size (no further normalisation), peaks with
>= 15 reads, `~ donor_sex + fraction`, quasi-likelihood F-test; Benjamini-Hochberg over
peaks with logFC > 0, enriched: q < 0.1. A 3'UTR (proximal start to distal end) is bound
if it overlaps a peak; two-sided Fisher's exact test among the genes tested for dAPA.

## Dependencies

R >= 4.0; edgeR, GenomicRanges (Bioconductor); ggplot2, scales (CRAN).
