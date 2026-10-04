#!/usr/bin/env Rscript
## ===========================================================================
## 04_count_segments.R
##
## Counts full-length RNA-seq reads on every 3'UTR segment across the iN cohort
## and summarises them as the mean size-factor-normalised count per segment.
## Only this per-segment mean enters the choice of the dominant proximal PAS
## (step 5), so the summary is distributed in resources/ and step 5 can be run
## without access to individual-level data.
##
## USAGE
##   Rscript 02_apa_pas_library/04_count_segments.R --samples <sheet>
##          [--gtf <segments gtf>] [--out <dir>] [--threads <n>]
##
## INPUTS
##   --samples  tab-separated sheet with columns
##                sample       sample id
##                bam          GRCh38 alignment (paired-end, duplicate-marked)
##                size_factor  DESeq2 size factor from the gene-level analysis
##              (one library per donor, the samples of the differential
##               expression analysis)
##   --gtf      default results/apa_pas_library/qapa_3utrs_iN_segments_hg38.gtf
##
## OUTPUTS
##   segment_counts_per_sample.rds               featureCounts output (not distributed)
##   segment_mean_normalized_counts_hg38.txt     ID, Chr, Start, End, Strand, Length,
##                                               n_samples, mean_norm_count
##                                               -> copy to resources/apa_pas_library/
##
## Dependencies: Rsubread (Bioconductor)
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))

SHEET   <- need_file(get_arg("--samples"), "--samples")
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "apa_pas_library"))
GTF     <- need_file(get_arg("--gtf", file.path(OUT_DIR, "qapa_3utrs_iN_segments_hg38.gtf")), "--gtf")
THREADS <- as.integer(get_arg("--threads", "4"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ss <- read.delim(SHEET, stringsAsFactors = FALSE)
stopifnot(all(c("sample", "bam", "size_factor") %in% names(ss)), !anyDuplicated(ss$sample))

fc <- count_utr_features(ss$bam, GTF, paired = TRUE, threads = THREADS)
colnames(fc$counts) <- ss$sample
saveRDS(fc, file.path(OUT_DIR, "segment_counts_per_sample.rds"))

norm <- sweep(fc$counts, 2, ss$size_factor, "/")
ann  <- fc$annotation
names(ann)[names(ann) == "GeneID"] <- "ID"
out  <- data.frame(ann[, c("ID", "Chr", "Start", "End", "Strand", "Length")],
                   n_samples = ncol(norm), mean_norm_count = rowMeans(norm))
f <- file.path(OUT_DIR, "segment_mean_normalized_counts_hg38.txt")
write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d segments x %d samples -> %s\n", nrow(out), ncol(norm), f))
