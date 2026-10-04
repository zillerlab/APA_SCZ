#!/usr/bin/env Rscript
## ===========================================================================
## 01_clip_enrichment.R
##
## PTBP2 CLIP-seq (quick-irCLIP) in iNs: peaks enriched in the PTBP2 immunoprecipitate
## (IP) over the input.
##
## Manuscript: Fig. 6a (iN CLIP binding), Supplementary Table 1 (CLIP-Seq).
##
## USAGE
##   Rscript 09_ptbp2_clip/01_clip_enrichment.R [--out <dir>]
##
## INPUTS (resources/clip/)
##   clip_PTBP2_iN_peak_counts.txt.gz   reads per peak of the union peak set (MACS3 peaks of
##                                      the IP libraries), 30,892 peaks x 8 libraries
##   clip_PTBP2_iN_libraries.txt        library, mapped_reads (total mapped reads), donor_sex
##                                      (F / M: one female and one male donor), fraction
##                                      (IP / input), replicate
##   clip_PTBP2_iN_peaks_hg38.bed       coordinates of the peaks (chrom, start, end, peak_id)
##
## OUTPUTS (results/clip/)
##   CLIP_PTBP2_iN_peak_stats.txt    one row per tested peak: peak_id, logFC (log2 IP / input),
##                                   logCPM, F, PValue, q_directional (BH over peaks with
##                                   logFC > 0), enriched (q_directional < 0.10)
##   CLIP_PTBP2_iN_enriched_peaks_hg38.bed   enriched peaks (chrom, start, end, peak_id)
##
## METHOD
##   edgeR. Library size = total mapped reads, normalisation factors 1 (the peaks were
##   called on the IP libraries, so the IP/input ratio within peaks is the signal and
##   must not be normalised away). Peaks with >= 15 reads summed over the 8 libraries are
##   tested; library sizes stay the total mapped reads after filtering. Design
##   ~ donor_sex + fraction; robust
##   dispersion and quasi-likelihood fit, QL F-test IP vs input. One-sided test for
##   enrichment: Benjamini-Hochberg over the two-sided p-values of peaks with logFC > 0;
##   enriched: q < 0.10.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
suppressPackageStartupMessages(library(edgeR))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "clip"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
MIN_TOTAL <- 15; Q <- 0.10

## ------------------------------------------------------------ inputs
rs   <- file.path(REPO, "resources", "clip")
cnt  <- read.delim(need_file(file.path(rs, "clip_PTBP2_iN_peak_counts.txt.gz"), "resources"), row.names = 1, check.names = FALSE)
libs <- read.delim(need_file(file.path(rs, "clip_PTBP2_iN_libraries.txt"), "resources"), stringsAsFactors = FALSE)
bed  <- read.delim(need_file(file.path(rs, "clip_PTBP2_iN_peaks_hg38.bed"), "resources"), header = FALSE,
                   col.names = c("chrom", "start", "end", "peak_id"), stringsAsFactors = FALSE)
libs <- libs[match(colnames(cnt), libs$library), ]
grp  <- factor(libs$fraction, levels = c("input", "IP"))
sex  <- factor(libs$donor_sex)
design <- model.matrix(~ sex + grp)

## ------------------------------------------------------------ edgeR
y <- DGEList(as.matrix(cnt), group = grp, lib.size = libs$mapped_reads)
y$samples$norm.factors <- 1
keep <- rowSums(y$counts) >= MIN_TOTAL
y <- y[keep, , keep.lib.sizes = TRUE]           # keep total mapped reads as library size
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
tt  <- topTags(glmQLFTest(fit, coef = "grpIP"), n = Inf, sort.by = "none")$table
tt  <- data.frame(peak_id = rownames(tt), tt[, c("logFC", "logCPM", "F", "PValue")], row.names = NULL)
up  <- tt$logFC > 0
tt$q_directional <- NA_real_; tt$q_directional[up] <- p.adjust(tt$PValue[up], "BH")
tt$enriched <- !is.na(tt$q_directional) & tt$q_directional < Q
tt <- tt[order(tt$PValue), ]
write.table(tt, file.path(OUT_DIR, "CLIP_PTBP2_iN_peak_stats.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

enr <- bed[bed$peak_id %in% tt$peak_id[tt$enriched], ]
stopifnot(nrow(enr) == sum(tt$enriched))           # coordinates available for every enriched peak
enr <- enr[order(enr$chrom, enr$start), ]
write.table(enr, file.path(OUT_DIR, "CLIP_PTBP2_iN_enriched_peaks_hg38.bed"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cat(sprintf("%d peaks, %d tested (>= %d reads in total); common dispersion %.4f (BCV %.3f)\n",
            nrow(cnt), sum(keep), MIN_TOTAL, y$common.dispersion, sqrt(y$common.dispersion)))
cat(sprintf("logFC > 0: %d peaks; enriched (q < %.2f): %d (q < 0.05: %d)\n",
            sum(up), Q, sum(tt$enriched), sum(tt$q_directional < 0.05, na.rm = TRUE)))
