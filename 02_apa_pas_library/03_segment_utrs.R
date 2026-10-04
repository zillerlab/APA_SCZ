#!/usr/bin/env Rscript
## ===========================================================================
## 03_segment_utrs.R
##
## Splits the nested 3'UTR isoforms reported by QAPA into non-overlapping
## segments, one per interval between consecutive PAS, so that reads can be
## assigned to the region between two PAS unambiguously.
##
## USAGE
##   Rscript 02_apa_pas_library/03_segment_utrs.R [--qapa <bed>] [--out <dir>]
##
## INPUT
##   --qapa  QAPA 3'UTR library built with the de novo 3' ends of step 1
##           (default resources/apa_pas_library/qapa_3utrs_iN_hg38.bed)
##
## OUTPUTS
##   qapa_3utrs_iN_segments_hg38.bed   non-overlapping segments (BED)
##   qapa_3utrs_iN_segments_hg38.gtf   same, as GTF for featureCounts (step 4)
##
## Dependencies: base R
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))

QAPA    <- need_file(get_arg("--qapa", file.path(REPO, "resources", "apa_pas_library",
                                                 "qapa_3utrs_iN_hg38.bed")), "--qapa")
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "apa_pas_library"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

qapa <- read.table(QAPA, sep = "\t", stringsAsFactors = FALSE, quote = "")
seg  <- segment_utrs(qapa)

bed <- file.path(OUT_DIR, "qapa_3utrs_iN_segments_hg38.bed")
write.table(seg, bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
bed_to_gtf(seg, sub("\\.bed$", ".gtf", bed))

g <- table(qapa[[7]])
cat(sprintf("%d 3'UTR isoforms from %d genes (%d with >1 PAS) -> %d segments\n",
            nrow(qapa), length(g), sum(g > 1), nrow(seg)))
cat("wrote", bed, "and .gtf\n")
