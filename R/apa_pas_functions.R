## ===========================================================================
## apa_pas_functions.R
##
## Shared functions for building the 3'UTR / polyadenylation-site (PAS)
## reference and quantifying reads on it. Sourced by the scripts in
## 02_apa_pas_library/ and by downstream APA modules.
##
##   segment_utrs()          QAPA 3'UTR isoforms -> non-overlapping segments
##   bed_to_gtf()            BED (0-based start) -> GTF (1-based start)
##   select_dominant_pas()   dominant proximal + distal segment per gene
##   extend_to_utr()         segments -> proximal UTR and distal extension
##   count_utr_features()    Rsubread::featureCounts with the APA settings
##   parse_utr_id()          transcript / gene / UTR coordinates from a QAPA id
## ===========================================================================

BED6 <- c("chr", "start", "end", "name", "score", "strand")

## ------------------------------------------------------------ QAPA identifiers
## QAPA 3'UTR ids have the form
##   <ENST>_<gene>_<species>_<chr>_<start>_<end>_<strand>_utr_<start>_<end>
## where the first coordinate pair spans the 3'UTR of that isoform
## (CDS end to PAS, BED coordinates). When QAPA collapses several transcripts
## with the same 3'UTR, the leading <ENST>_<gene> part is a comma-separated list;
## transcript and gene are then taken from its first entry.
parse_utr_id <- function(id) {
  tg <- strsplit(sub(",.*", "", id), "_", fixed = TRUE)   # first <ENST>_<gene>
  p  <- strsplit(id, "_", fixed = TRUE)                    # full id for coordinates
  k  <- vapply(p, function(r) grep("chr", r, fixed = TRUE)[1], 1L)
  data.frame(transcriptId = vapply(tg, `[`, "", 1),
             geneId       = vapply(tg, `[`, "", 2),
             utr_start    = as.numeric(mapply(`[`, p, k + 1)),
             utr_end      = as.numeric(mapply(`[`, p, k + 2)),
             stringsAsFactors = FALSE)
}

## ------------------------------------------------------------ segmentation
## All 3'UTR isoforms of a gene share the CDS end and differ in their PAS.
## Nested isoforms are cut into non-overlapping segments, each spanning the
## region between two consecutive PAS. Rows of the QAPA output are ordered by
## position within each gene. Genes with a single isoform are kept unchanged.
##   qapa   data frame read from the QAPA BED (column 4 = id, column 6 =
##          strand, column 7 = gene)
segment_utrs <- function(qapa) {
  gene <- qapa[[7]]
  multi <- gene %in% names(which(table(gene) > 1))
  segs <- lapply(split(qapa[multi, ], factor(gene[multi], levels = unique(gene[multi]))),
                 function(tmp) {
                   n <- nrow(tmp)
                   if (tmp[1, 6] == "+") {
                     tmp[2:n, 2] <- tmp[1:(n - 1), 3] + 1   # start after previous PAS
                   } else {
                     tmp[1:(n - 1), 3] <- tmp[2:n, 2] - 1   # end before next PAS
                   }
                   tmp
                 })
  segs <- do.call(rbind, segs)
  out  <- rbind(segs, qapa[!qapa[[4]] %in% segs[[4]], ])
  rownames(out) <- NULL
  out[order(out[[1]], out[[2]]), ]
}

## ------------------------------------------------------------ BED -> GTF
## One "transcript" feature per interval; featureCounts uses transcript_id as
## the feature id (GTF.attrType = "transcript_id").
bed_to_gtf <- function(bed, file, source = "qapa", shift_start = 1L) {
  id  <- bed[[4]]
  gtf <- data.frame(bed[[1]], source, "transcript",
                    format(bed[[2]] + shift_start, scientific = FALSE, trim = TRUE),
                    format(bed[[3]], scientific = FALSE, trim = TRUE),
                    ".", bed[[6]], ".",
                    sprintf('gene_id "%s"; transcript_id "%s";', id, id))
  write.table(gtf, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(file)
}

## ------------------------------------------------------------ dominant PAS
## For every gene with more than one segment, keep the most distal segment
## (distal PAS) and one dominant proximal segment. With two segments these are
## the two segments. With more, the proximal PAS is placed at the largest change
## in mean read density (normalised reads per bp) between adjacent segments,
## ignoring the step into the distal segment:
##   largest drop  (density falls 5'->3')  -> proximal = 5' segment of the pair
##   largest rise  (density rises 5'->3')  -> proximal = 3' segment of the pair
##
##   ann   featureCounts annotation of the segments: ID, Chr, Start, End, Strand
##   dens  matrix of normalised read density per segment (rows = ann), one
##         column per sample or a single column of the mean across samples
select_dominant_pas <- function(ann, dens) {
  dens <- as.matrix(dens)
  ids  <- parse_utr_id(ann$ID)
  ann$transcriptId <- ids$transcriptId
  ann$geneId       <- ids$geneId
  o    <- order(ann$Chr, ann$Start)
  ann  <- ann[o, ]; dens <- dens[o, , drop = FALSE]

  g   <- ann$geneId
  sel <- unique(g)[table(g)[unique(g)] > 1]
  message(sprintf("identified %d genes with more than one PAS", length(sel)))

  pick <- t(vapply(sel, function(X) {
    i <- which(g == X); n <- length(i); d <- dens[i, , drop = FALSE]
    m <- ann[i, ]
    if (n == 2) {
      if (m$Strand[1] == "+") c(m$ID[1], m$ID[n]) else c(m$ID[n], m$ID[1])
    } else if (m$Strand[1] == "+") {
      delta <- rowMeans(d[1:(n - 2), , drop = FALSE] - d[2:(n - 1), , drop = FALSE])
      id <- which.max(abs(delta)); if (sign(delta[id]) < 0) id <- id + 1
      c(m$ID[id], m$ID[n])
    } else {
      delta <- rowMeans(d[2:(n - 1), , drop = FALSE] - d[3:n, , drop = FALSE])
      id <- which.max(abs(delta)); if (sign(delta[id]) < 0) id <- id + 1
      c(m$ID[id + 1], m$ID[1])
    }
  }, c(proximal = "", distal = "")))

  keep <- ann[ann$geneId %in% sel, ]
  out  <- rbind(cbind(keep[keep$ID %in% pick[, "proximal"], ], class = "proximal"),
                cbind(keep[keep$ID %in% pick[, "distal"],   ], class = "distal"))
  out  <- out[order(out$Chr, out$Start), ]
  rownames(out) <- NULL
  out
}

## ------------------------------------------------------------ full UTRs
## Proximal region: CDS end to proximal PAS (the shared part of both isoforms).
## Distal region:   proximal PAS to distal PAS (present only in the long isoform).
extend_to_utr <- function(seg) {
  out <- do.call(rbind, lapply(split(seg, factor(seg$geneId, levels = unique(seg$geneId))),
    function(tmp) {
      u <- parse_utr_id(tmp$ID[1]); p <- tmp$class == "proximal"
      if (tmp$Strand[1] == "+") {
        tmp$Start[p]  <- u$utr_start
        tmp$Start[!p] <- tmp$End[p] + 1
      } else {
        tmp$End[p]    <- u$utr_end
        tmp$End[!p]   <- tmp$Start[p] - 1
      }
      tmp
    }))
  out <- out[order(out$Chr, out$Start), ]
  rownames(out) <- NULL
  out
}

## ------------------------------------------------------------ quantification
## Reads per UTR feature, used for all APA quantifications in the study
## (reference building, cohort, RBP knockdown). Paired-end libraries require
## both mates to be mapped; single-end libraries use paired = FALSE.
count_utr_features <- function(bams, gtf, paired = TRUE, multimap = FALSE, threads = 4) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) stop("package Rsubread is required")
  Rsubread::featureCounts(files = bams, annot.ext = gtf, isGTFAnnotationFile = TRUE,
                          GTF.featureType = "transcript", GTF.attrType = "transcript_id",
                          useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                          countMultiMappingReads = multimap, primaryOnly = TRUE,
                          ignoreDup = FALSE, minOverlap = 10,
                          isPairedEnd = paired, requireBothEndsMapped = paired,
                          nthreads = threads)
}

## ------------------------------------------------------------ paths
## Command-line helpers: --flag <value> overrides a default location.
get_arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (!is.na(i) && i < length(a)) a[i + 1] else default
}
need_file <- function(p, flag) {
  if (is.null(p) || !file.exists(p)) stop("input not found: ", p,
                                          "\n  set it with ", flag, " <file>")
  p
}
