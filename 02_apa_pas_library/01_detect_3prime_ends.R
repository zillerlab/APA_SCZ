#!/usr/bin/env Rscript
## ===========================================================================
## 01_detect_3prime_ends.R
##
## De novo detection of polyadenylated 3' transcript ends in iNs from 3'
## RNA-seq (QuantSeq), using transcriptR. The resulting 3' ends are the custom
## PAS set supplied to QAPA in step 2.
##
## Input BAM: coordinate-sorted merge (Picard MergeSamFiles) of the
## duplicate-filtered, GRCh38-aligned 3' RNA-seq libraries of HC and SCZ iNs.
## The BAM is not distributed (EGA); this step documents how the PAS set in
## resources/apa_pas_library/ was produced.
##
## USAGE
##   Rscript 02_apa_pas_library/01_detect_3prime_ends.R --bam <merged.bam>
##          [--genes <gencode genes BED>] [--out <dir>]
##
## INPUTS
##   --bam     merged 3' RNA-seq BAM (+ .bai)
##   --genes   GENCODE v27 gene intervals, BED6
##             (default resources/apa_pas_library/gencode.v27.annotation_genes.bed)
##
## OUTPUT
##   transcriptR_3prime_ends_hg38.bed   BED6: chr, start, end, id, fragments,
##                                      strand; input for qapa build (-o)
##
## PARAMETERS (Supplementary Methods)
##   fragment size 250, strand swapped on import (see note below),
##   background FDR 0.01, gap distance 75 bp, transcriptR minimum length
##   50 bp and FPKM 0.5, then kept if within 100 bp of a GENCODE v27 gene,
##   >= 10 fragments, FPKM >= 2 and width <= 2 kb.
##
## NOTE on strand: reads are imported with swap.strand = TRUE and the strand of
## the detected 3' ends is swapped back when writing the output, so the output
## strand is the strand of the reads.
##
## Dependencies: transcriptR, GenomicRanges, Rsamtools (Bioconductor)
## ===========================================================================

suppressPackageStartupMessages({
  library(transcriptR)
  library(GenomicRanges)
  library(Rsamtools)
})
.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))

BAM     <- need_file(get_arg("--bam"), "--bam")
GENES   <- need_file(get_arg("--genes", file.path(REPO, "resources", "apa_pas_library",
                                                  "gencode.v27.annotation_genes.bed")), "--genes")
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "apa_pas_library"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FRAGMENT_SIZE <- 250
BG_FDR        <- 0.01
GAP_DIST      <- 75
MIN_LENGTH    <- 50
MIN_FPKM_TR   <- 0.5
GENE_MAXGAP   <- 100
MIN_FRAGMENTS <- 10
MIN_FPKM      <- 2
MAX_WIDTH     <- 2000
CHROMS        <- c(paste0("chr", 1:22), "chrX", "chrY")

## ---------------------------------------------------------- regions: chr1-22, X, Y
len <- seqlengths(BamFile(BAM))[CHROMS]
if (anyNA(len)) stop("BAM header lacks: ", paste(CHROMS[is.na(len)], collapse = ", "))
region <- GRanges(CHROMS, IRanges(1, len))

## ---------------------------------------------------------- transcriptR
tds <- constructTDS(file = BAM, region = region, fragment.size = FRAGMENT_SIZE,
                    unique = FALSE, paired.end = FALSE, swap.strand = TRUE)
levels(seqnames(tds@fragments)) <- levels(droplevels(seqnames(tds@fragments)))
estimateBackground(tds, fdr.cutoff = BG_FDR)
detectTranscripts(tds, estimate.params = TRUE, gap.dist = GAP_DIST)
trx <- getTranscripts(tds, min.length = MIN_LENGTH, min.fpkm = MIN_FPKM_TR)

## ---------------------------------------------------------- keep genic 3' ends
genes <- read.table(GENES, sep = "\t", quote = "")[, 1:6]
names(genes) <- BED6
genes <- makeGRangesFromDataFrame(genes)
hit   <- findOverlaps(trx, genes, maxgap = GENE_MAXGAP, select = "first", ignore.strand = TRUE)
tx    <- as.data.frame(trx[!is.na(hit)])
# columns: seqnames start end width strand id length bases.covered coverage fragments fpkm

keep  <- tx$fragments >= MIN_FRAGMENTS & tx$fpkm >= MIN_FPKM & tx$width <= MAX_WIDTH
tx    <- tx[keep, ]
out   <- data.frame(chr = tx$seqnames, start = tx$start, end = tx$end,
                    name = paste(tx$seqnames, tx$start, tx$end, sep = "_"),
                    score = tx$fragments,
                    strand = ifelse(tx$strand == "-", "+", "-"))   # undo swap.strand
f <- file.path(OUT_DIR, "transcriptR_3prime_ends_hg38.bed")
write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("%d transcripts detected, %d genic, %d after filtering -> %s\n",
            length(trx), sum(!is.na(hit)), nrow(out), f))
