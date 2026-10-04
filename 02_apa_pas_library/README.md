# 02 APA PAS library

Detection of polyadenylation sites (PAS) in iNs and construction of the 3'UTR reference
with one dominant proximal and one distal region per gene, used by all APA analyses.

**Manuscript:** Fig. 2a,b; Suppl. Table 3 (`APA_sites`).

## Workflow

| Step | Script | Output |
|---|---|---|
| 1 | `01_detect_3prime_ends.R` | 3' ends from merged 3' RNA-seq (transcriptR) |
| 2 | `qapa build` | 3'UTR isoform library (QAPA) |
| 3 | `03_segment_utrs.R` | Non-overlapping 3'UTR segments |
| 4 | `04_count_segments.R` | Mean normalised read count per segment |
| 5 | `05_select_dominant_pas.R` | Proximal and distal UTR regions |

Steps 1, 2 and 4 require the sequencing data (see Data availability); their outputs are
provided in `resources/apa_pas_library/`. Steps 3 and 5 run from these:

```bash
Rscript 02_apa_pas_library/03_segment_utrs.R
Rscript 02_apa_pas_library/05_select_dominant_pas.R
```

Steps 1, 2 and 4:

```bash
Rscript 02_apa_pas_library/01_detect_3prime_ends.R --bam <merged_3prime.bam>
qapa build --db ensembl_identifiers.txt -o transcriptR_3prime_ends_hg38.bed \
           gencode.basic.txt > qapa_3utrs_iN_hg38.bed
Rscript 02_apa_pas_library/04_count_segments.R --samples samples.tsv
```

## Inputs (`resources/apa_pas_library/`)

| File | Content |
|---|---|
| `gencode.v27.annotation_genes.bed` | GENCODE v27 gene intervals |
| `qapa_3utrs_iN_hg38.bed` | QAPA 3'UTR isoform library (step 2) |
| `segment_mean_normalized_counts_hg38.txt` | Mean size-factor-normalised count per segment (step 4) |

## Outputs (`results/apa_pas_library/`)

| File | Content |
|---|---|
| `qapa_3utrs_iN_segments_hg38.{bed,gtf}` | 3'UTR segments |
| `dominant_pas_segments_hg38.{txt,bed}` | Selected proximal and distal segments |
| `dominant_pas_utrs_hg38.{txt,bed,gtf}` | Proximal and distal UTR regions (Suppl. Table 3) |

## Method

For genes with more than two PAS, the dominant proximal PAS is placed at the largest
change in mean normalised read density between adjacent segments upstream of the most
distal segment. Read counting settings (featureCounts) are in
`R/apa_pas_functions.R::count_utr_features()`.

## Dependencies

R >= 4.0; transcriptR, GenomicRanges, Rsamtools (step 1) and Rsubread (step 4), all
Bioconductor; QAPA (step 2).
