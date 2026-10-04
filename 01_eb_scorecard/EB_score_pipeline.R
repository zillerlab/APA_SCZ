#!/usr/bin/env Rscript
## ===========================================================================
## EB_score_pipeline.R
##
## Germ-layer and pluripotency signature scoring of iPSC-derived embryoid
## bodies (EBs), used as tri-lineage differentiation QC of all iPSC lines.
## Produces the data shown in Extended Data Fig. 1g,h.
##
## Raabe, Atella, Hausruckinger, Gagliardi, Almeida et al.
## "Polygenic risk for schizophrenia is associated with alternative
##  polyadenylation as molecular mechanism contributing to synaptic impairment"
##
## Dependencies: base R (>= 4.0) and stats. No additional packages.
##
## USAGE (from the repository root or any other directory)
##   Rscript 01_eb_scorecard/EB_score_pipeline.R
##   Rscript 01_eb_scorecard/EB_score_pipeline.R --resources <dir> --out <dir>
##
##   Defaults:  --resources  <repo>/resources/eb_scorecard
##              --out        <repo>/results/eb_scorecard
##
## INPUTS (in the resources directory)
##   EB_counts_clean.txt    raw counts, 16,638 genes x 173 samples
##                          (20 reference iPSC + 153 EB); first column geneName
##   EB_metadata_clean.txt  one row per count column: sample, cell_line,
##                          condition (iPSC|EB), batch, is_reference,
##                          library_size, n_genes_detected
##   EB_marker_panel.txt    layer, gene, weight, ref_median_log2CPM, ref_MAD,
##                          ref_batch_span, direction (176 genes)
##
## OUTPUTS (in the output directory)
##   EB_scores_per_sample.txt      score / t / p / FDR per signature, pass calls
##   EB_scores_per_line.txt        per-line verdict
##   ExtDataFig1gh_EB_scores.txt   reference iPSCs plus one representative EB
##                                 per R1-passing line (Extended Data Fig. 1g,h)
##   sessionInfo.txt
##
## METHOD
##   1  normalise   counts -> CPM -> log2(CPM + 1). Depth normalisation only;
##                  scores are contrasts against the reference iPSC samples.
##   2  robust z    per gene, against the reference iPSC samples:
##                    z = (x - median_ref) / max(MAD_ref * 1.4826, 0.25),
##                  winsorised at +/-10 so that no single gene dominates.
##   3  aggregate   weighted mean of z over the genes of each signature,
##                    S = sum(w_g * z_g) / sum(|w_g|),
##                  with w_g = 1 / (1 + MAD_ref) scaled to mean 1 per layer
##                  (from the panel file).
##   4  test        Crawford & Howell (1998) single case vs control group:
##                    t = (S - mean(S_ref)) / (sd(S_ref) * sqrt(1 + 1/n)),
##                  df = n - 1, one-sided in the direction given by the panel
##                  (germ layers up, pluripotency down). Reference samples are
##                  scored leave-one-out (df = n - 2). Benjamini-Hochberg FDR
##                  across EB samples, per signature.
##   5  call        R2 = all three germ-layer signatures FDR < ALPHA
##                  R1 = R2 and pluripotency signature FDR < ALPHA
##                  A cell line passes if at least one of its EBs passes.
##   6  figure      For each R1-passing line, the EB with the smallest maximum
##                  FDR across the four signatures is shown in Extended Data
##                  Fig. 1g,h together with all reference iPSC samples.
##
## EXPECTED RESULT (printed at the end)
##   EB samples  R2 130/153   R1 119/153
##   cell lines  R2  93/103   R1  85/103
##   reference iPSC passing: 0
## ===========================================================================

ALPHA   <- 0.05   # FDR threshold for signature calls
Z_FLOOR <- 0.25   # floor on the reference MAD, in log2 units
WINSOR  <- 10     # cap on per-gene |z|
MIN_REF <- 8      # minimum number of reference samples for the test

## ------------------------------------------------------------------- paths
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)   # source()
  if (!is.null(of)) return(dirname(normalizePath(of)))
  getwd()
}
get_arg <- function(flag, default) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (!is.na(i) && i < length(a)) a[i + 1] else default
}
REPO    <- normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
RES_DIR <- get_arg("--resources", file.path(REPO, "resources", "eb_scorecard"))
OUT_DIR <- get_arg("--out",       file.path(REPO, "results",   "eb_scorecard"))

in_file <- function(f) {
  p <- file.path(RES_DIR, f)
  if (!file.exists(p)) stop("input not found: ", p,
                            "\n  set the location with --resources <dir>")
  p
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("resources:", normalizePath(RES_DIR), "\noutput:   ", normalizePath(OUT_DIR), "\n")

## ------------------------------------------------------------- read inputs
cnt   <- as.matrix(read.delim(in_file("EB_counts_clean.txt"), row.names = 1,
                              check.names = FALSE))
meta  <- read.delim(in_file("EB_metadata_clean.txt"), stringsAsFactors = FALSE)
panel <- read.delim(in_file("EB_marker_panel.txt"),   stringsAsFactors = FALSE)

stopifnot(!anyDuplicated(meta$sample), setequal(meta$sample, colnames(cnt)))
meta <- meta[match(colnames(cnt), meta$sample), ]
rownames(meta) <- meta$sample
meta$is_reference <- as.logical(meta$is_reference)

missing <- setdiff(panel$gene, rownames(cnt))
if (length(missing)) {
  warning(sprintf("%d panel genes absent from the count matrix and skipped: %s",
                  length(missing), paste(missing, collapse = ", ")))
  panel <- panel[panel$gene %in% rownames(cnt), ]
}
refi <- meta$sample[meta$is_reference]
ebi  <- meta$sample[meta$condition == "EB"]
n    <- length(refi)
if (n < MIN_REF) stop("fewer than ", MIN_REF, " reference samples; ",
                      "the single-case test is not usable")
LAYERS <- unique(panel$layer)
SIDE   <- vapply(split(panel$direction, panel$layer), function(x) x[1], "")[LAYERS]
cat(sprintf("%d samples: %d reference iPSC, %d EB | panel %s\n",
            ncol(cnt), n, length(ebi),
            paste(sprintf("%s=%d", LAYERS, as.vector(table(panel$layer)[LAYERS])),
                  collapse = " ")))

## -------------------------------------------------------------- 1 normalise
cpm <- sweep(cnt, 2, colSums(cnt), "/") * 1e6
lg  <- log2(cpm + 1)

## -------------------------------------------------------------- 2 robust z
ref_med <- apply(lg[, refi, drop = FALSE], 1, median)
ref_mad <- apply(abs(lg[, refi, drop = FALSE] - ref_med), 1, median) * 1.4826
Z <- pmin(pmax((lg - ref_med) / pmax(ref_mad, Z_FLOOR), -WINSOR), WINSOR)

## consistency check against the reference statistics stored in the panel
if (all(c("ref_median_log2CPM", "ref_MAD") %in% names(panel))) {
  d <- max(abs(ref_med[panel$gene] - panel$ref_median_log2CPM),
           abs(ref_mad[panel$gene] - panel$ref_MAD))
  cat(sprintf("reference statistics reproduce the panel file to %.3f log2 units\n", d))
  if (d > 0.001) warning("reference statistics differ from the panel file; ",
                         "has the reference set changed since the panel was built?")
}

## -------------------------------------------------------------- 3 aggregate
S <- sapply(LAYERS, function(L) {
  g <- panel$gene[panel$layer == L]
  w <- panel$weight[panel$layer == L]
  as.numeric(crossprod(Z[g, , drop = FALSE], w) / sum(abs(w)))
})
rownames(S) <- colnames(lg)

## -------------------------------------------------------------- 4 test
## Crawford & Howell single case vs control group; leave-one-out for controls
crawford <- function(s, refi, side) {
  r  <- s[refi]; nn <- length(refi)
  tt <- (s - mean(r)) / (sd(r) * sqrt(1 + 1 / nn))
  for (k in refi) {
    o <- r[setdiff(refi, k)]
    tt[k] <- (s[k] - mean(o)) / (sd(o) * sqrt(1 + 1 / length(o)))
  }
  df <- ifelse(names(s) %in% refi, nn - 2, nn - 1)
  p  <- if (side == "up") pt(tt, df, lower.tail = FALSE) else pt(tt, df, lower.tail = TRUE)
  list(t = tt, p = setNames(p, names(s)))
}

res <- data.frame(sample = colnames(lg), cell_line = meta$cell_line,
                  condition = meta$condition, batch = meta$batch,
                  is_reference = meta$is_reference,
                  row.names = colnames(lg), stringsAsFactors = FALSE)
for (L in LAYERS) {
  ch <- crawford(setNames(S[, L], rownames(S)), refi, SIDE[[L]])
  res[[paste0(L, "_score")]] <- round(S[, L], 4)
  res[[paste0(L, "_t")]]     <- round(ch$t, 3)
  res[[paste0(L, "_p")]]     <- ch$p
  q <- rep(NA_real_, nrow(res))
  q[res$condition == "EB"] <- p.adjust(ch$p[ebi], "BH")
  res[[paste0(L, "_fdr")]]   <- q
}

## -------------------------------------------------------------- 5 call
GERM  <- LAYERS[SIDE == "up"]
PLURI <- LAYERS[SIDE == "down"]
fdrG  <- as.matrix(res[, paste0(GERM, "_fdr"), drop = FALSE])
isEB  <- res$condition == "EB"
res$n_layers_significant <- ifelse(isEB, rowSums(fdrG < ALPHA), NA)
res$pass_R2_3layers <- isEB & rowSums(fdrG < ALPHA) == length(GERM)
res$pass_R1_3layers_and_pluri <- res$pass_R2_3layers &
  apply(as.matrix(res[, paste0(PLURI, "_fdr"), drop = FALSE]) < ALPHA, 1, all)

e  <- res[isEB, ]
ln <- data.frame(cell_line    = names(table(e$cell_line)),
                 n_EB         = as.vector(table(e$cell_line)),
                 n_EB_pass_R1 = as.vector(tapply(e$pass_R1_3layers_and_pluri, e$cell_line, sum)),
                 n_EB_pass_R2 = as.vector(tapply(e$pass_R2_3layers,           e$cell_line, sum)),
                 stringsAsFactors = FALSE)
for (L in LAYERS) ln[[paste0("best_", L, "_fdr")]] <-
  round(as.vector(tapply(e[[paste0(L, "_fdr")]], e$cell_line, min)), 4)
ln$LINE_PASS_R1 <- ln$n_EB_pass_R1 > 0
ln$LINE_PASS_R2 <- ln$n_EB_pass_R2 > 0

## -------------------------------------------------------------- 6 figure table
## one representative EB per R1-passing line: smallest maximum FDR across the
## four signatures (first in input order on ties), plus all reference iPSCs
fdrAll  <- as.matrix(e[, paste0(LAYERS, "_fdr")])
e$maxFDR <- apply(fdrAll, 1, max)
ePass   <- e[e$pass_R1_3layers_and_pluri, ]
repEB   <- vapply(split(seq_len(nrow(ePass)), ePass$cell_line),
                  function(i) ePass$sample[i[which.min(ePass$maxFDR[i])]], "")
figCols <- c("sample", "cell_line", "condition", "is_reference",
             paste0(rep(c("Ectoderm", "Mesoderm", "Endoderm", "Pluri"), each = 3),
                    c("_score", "_t", "_p")))
fig <- res[res$sample %in% c(refi, repEB), figCols]

## -------------------------------------------------------------- write
wt <- function(x, f) write.table(x, file.path(OUT_DIR, f), sep = "\t",
                                 quote = FALSE, row.names = FALSE)
wt(res, "EB_scores_per_sample.txt")
wt(ln,  "EB_scores_per_line.txt")
wt(fig, "ExtDataFig1gh_EB_scores.txt")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat(sprintf("\nEB samples  R2 %d/%d   R1 %d/%d\n", sum(e$pass_R2_3layers), nrow(e),
            sum(e$pass_R1_3layers_and_pluri), nrow(e)))
cat(sprintf("cell lines  R2 %d/%d   R1 %d/%d\n", sum(ln$LINE_PASS_R2), nrow(ln),
            sum(ln$LINE_PASS_R1), nrow(ln)))
cat(sprintf("reference iPSC passing (must be 0): %d\n",
            sum(res$pass_R2_3layers[!isEB], na.rm = TRUE)))
fail <- ln$cell_line[!ln$LINE_PASS_R2]
if (length(fail)) cat("lines with no tri-lineage EB:", paste(fail, collapse = ", "), "\n")
cat(sprintf("Extended Data Fig. 1g,h table: %d reference iPSC + %d EB\n",
            sum(fig$is_reference), sum(!fig$is_reference)))
