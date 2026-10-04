#!/usr/bin/env Rscript
## ===========================================================================
## 05_select_dominant_pas.R
##
## Reduces each gene with more than one PAS to two 3'UTR regions: a proximal
## region (CDS end to the dominant proximal PAS) and a distal region (proximal
## PAS to the most distal PAS). This two-PAS model is the reference for all
## differential APA analyses (Supplementary Table 3, sheet APA_sites).
##
## Dominant proximal PAS: among the segments upstream of the most distal PAS,
## the PAS at the largest change in mean normalised read density (reads per bp)
## between adjacent segments, using the full-length RNA-seq of the iN cohort.
## See select_dominant_pas() in R/apa_pas_functions.R.
##
## USAGE
##   Rscript 02_apa_pas_library/05_select_dominant_pas.R [--counts <file>] [--out <dir>]
##
## INPUT
##   --counts  per-segment mean normalised counts (step 4)
##             (default resources/apa_pas_library/segment_mean_normalized_counts_hg38.txt)
##
## OUTPUTS
##   dominant_pas_segments_hg38.{txt,bed}  the selected proximal and distal segments
##   dominant_pas_utrs_hg38.{txt,bed,gtf}  proximal and distal UTR regions;
##                                         the .gtf is the annotation quantified in
##                                         all downstream APA analyses
##
## Dependencies: base R
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))

COUNTS  <- need_file(get_arg("--counts", file.path(REPO, "resources", "apa_pas_library",
                                                   "segment_mean_normalized_counts_hg38.txt")), "--counts")
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "apa_pas_library"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cnt  <- read.delim(COUNTS, stringsAsFactors = FALSE)
ann  <- cnt[, c("ID", "Chr", "Start", "End", "Strand", "Length")]
dens <- cnt$mean_norm_count / (cnt$End - cnt$Start)          # normalised reads per bp

seg <- select_dominant_pas(ann, dens)
utr <- extend_to_utr(seg)

wr <- function(x, stem) {
  write.table(x, file.path(OUT_DIR, paste0(stem, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  bed <- x[, c("Chr", "Start", "End", "ID", "Length", "Strand")]
  write.table(bed, file.path(OUT_DIR, paste0(stem, ".bed")), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  invisible(bed)
}
wr(seg, "dominant_pas_segments_hg38")
bed <- wr(utr, "dominant_pas_utrs_hg38")
bed_to_gtf(bed, file.path(OUT_DIR, "dominant_pas_utrs_hg38.gtf"))   # GTF start = Start + 1

cat(sprintf("%d genes, %d proximal + %d distal UTR regions -> %s\n",
            length(unique(utr$geneId)), sum(utr$class == "proximal"),
            sum(utr$class == "distal"), OUT_DIR))
