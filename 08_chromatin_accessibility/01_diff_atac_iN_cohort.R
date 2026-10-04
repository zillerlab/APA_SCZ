#!/usr/bin/env Rscript
## ===========================================================================
## 01_diff_atac_iN_cohort.R
##
## Chromatin accessibility (ATAC-seq) in iNs (day 49): PCA of all libraries and
## differentially accessible peaks between SCZ and HC.
##
## Manuscript: Fig. 1c (open chromatin PCA), Extended Data Fig. 4e, main text.
##
## USAGE
##   Rscript 08_chromatin_accessibility/01_diff_atac_iN_cohort.R [--out <dir>]
##
## INPUTS (resources/atac/)
##   atac_peak_counts_iN_cohort.txt.gz  fragments per peak of the union peak set (hg19;
##                                      peaks called in >= 6 libraries), 109,114 peaks x
##                                      106 libraries (columns <PluriCore ID>_<replicate>)
##   atac_iN_cohort_metadata.txt        library, line (PluriCore ID), diagnosis, sex, site
##                                      (differentiation site), FRiP, TSSE
##
## OUTPUTS
##   results/atac/diffATAC_iN_SCZ_vs_Ctrl.txt   one row per peak: peak_id, chr, start, end,
##       length, logFC (log2 SCZ / HC), logCPM, F, PValue, FDR (BH), qvalue (Storey),
##       diff (qvalue <= 0.1)
##   results/atac/RUV_W1.txt                    unwanted-variation factor per library
##   results/figures/Fig1c_ATAC_PCA.pdf / .txt  PCA (Source Data Fig. 1c)
##   results/figures/ExtFig4e_diffATAC_volcano.pdf
##
## METHOD
##   PCA: counts scaled to the mean library size, log2(x + 1); peaks with a value above
##   the mean (over libraries) of the 99th percentile in any library are removed;
##   prcomp, centred and scaled.
##   Differential accessibility (edgeR): upper-quartile normalisation; design
##   ~0 + batch + FRiP + diagnosis (all four diagnoses, HC = reference); common and
##   tagwise NB dispersion; quasi-likelihood F-test of the SCZ coefficient. Peaks
##   outside the 10,000 most significant are used as negative controls for RUVg
##   (k = 1; log(count + 1), centred, first left singular vector). The model is refitted
##   with W_1 as covariate (~0 + W_1 + batch + FRiP + diagnosis). Storey q-values;
##   differential: q <= 0.1.
##   All 106 libraries (96 lines; 10 lines with two libraries) enter as independent
##   samples.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(edgeR); library(qvalue) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "atac"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
N_TOP_EXCLUDED <- 10000; Q <- 0.1
LEGACY <- if (packageVersion("edgeR") >= "4.0.0") list(legacy = TRUE) else list()   # edgeR 4: keep the classic QL fit

## ------------------------------------------------------------ inputs
rs   <- file.path(REPO, "resources", "atac")
meta <- read.delim(need_file(file.path(rs, "atac_iN_cohort_metadata.txt"), "resources"), stringsAsFactors = FALSE)
x    <- read.delim(need_file(file.path(rs, "atac_peak_counts_iN_cohort.txt.gz"), "resources"),
                   check.names = FALSE, colClasses = c(chr = "character"))
peaks <- x[, c("peak_id", "chr", "start", "end", "length")]
cnt  <- as.matrix(x[, meta$library]); rownames(cnt) <- peaks$peak_id
an   <- data.frame(batch = factor(meta$site), FRiP = meta$FRiP,
                   diagnosis = relevel(factor(meta$diagnosis), ref = "Ctrl"), row.names = meta$library)
cat(sprintf("%d peaks x %d libraries (%d lines): %s\n", nrow(cnt), ncol(cnt), length(unique(meta$line)),
            paste(names(table(an$diagnosis)), table(an$diagnosis), collapse = ", ")))

## ------------------------------------------------------------ PCA (Fig. 1c)
lnorm <- log2(t(t(cnt) / (colSums(cnt) / mean(colSums(cnt)))) + 1)
cut   <- mean(apply(lnorm, 2, quantile, 0.99))
pca   <- prcomp(t(lnorm[rowSums(lnorm > cut) == 0, ]), center = TRUE, scale. = TRUE)
ve    <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pcs   <- data.frame(library = meta$library, line = meta$line, diagnosis = meta$diagnosis,
                    PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write.table(pcs, file.path(FIG_DIR, "Fig1c_ATAC_PCA.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(pcs, aes(PC1, PC2, colour = diagnosis)) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(Ctrl = "grey55", SCZ = "#c8321f", BD = "#1f78c8", MDD = "#e0a020")) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL,
       title = sprintf("Open chromatin (ATAC-seq, n = %d)", ncol(cnt))) +
  theme_classic(base_size = 10) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(FIG_DIR, "Fig1c_ATAC_PCA.pdf"), p, width = 3.6, height = 3)

## ------------------------------------------------------------ edgeR QL fit
ql_fit <- function(design) {
  y <- calcNormFactors(DGEList(counts = cnt), method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design); y <- estimateGLMTagwiseDisp(y, design)
  f <- do.call(glmQLFit, c(list(y, design), LEGACY))
  topTags(glmQLFTest(f, coef = "diagnosisSCZ"), n = Inf, sort.by = "PValue")$table
}
ruvg_w1 <- function(counts, ctrl) {                 # RUVSeq::RUVg(k = 1), matrix method
  Y <- t(log(counts + 1)); Y <- sweep(Y, 2, colMeans(Y))
  svd(Y[, ctrl])$u[, 1]
}
r0 <- ql_fit(model.matrix(~ 0 + batch + FRiP + diagnosis, an))
an$W_1 <- ruvg_w1(cnt, !rownames(cnt) %in% rownames(r0)[seq_len(N_TOP_EXCLUDED)])
res <- ql_fit(model.matrix(~ 0 + W_1 + batch + FRiP + diagnosis, an))
res$qvalue <- qvalue(res$PValue)$qvalue
res <- data.frame(peaks[match(rownames(res), peaks$peak_id), ], res[, c("logFC", "logCPM", "F", "PValue", "FDR", "qvalue")],
                  diff = res$qvalue <= Q, row.names = NULL)
write.table(res, file.path(OUT_DIR, "diffATAC_iN_SCZ_vs_Ctrl.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(library = meta$library, W_1 = an$W_1), file.path(OUT_DIR, "RUV_W1.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ volcano (Ext. Data Fig. 4e)
res$class <- ifelse(!res$diff, "NO", ifelse(res$logFC > 0, "Up", "Down"))
vd <- res[order(res$class != "NO"), ]
p <- ggplot(vd, aes(logFC, -log10(PValue), colour = class)) +
  geom_point(size = 0.5, stroke = 0) +
  scale_colour_manual(values = c(NO = "grey80", Up = "#c8321f", Down = "#1f78c8"), breaks = c("Down", "Up"),
                      labels = c(Down = sprintf("less accessible in SCZ (%d)", sum(res$class == "Down")),
                                 Up = sprintf("more accessible in SCZ (%d)", sum(res$class == "Up")))) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "log2 fold change (SCZ / HC)", y = "-log10 p", colour = NULL,
       caption = sprintf("%d of %d peaks q <= %.1f", sum(res$diff), nrow(res), Q)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = FIG_SEC), legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.key.height = unit(0.35, "cm"), legend.text = element_text(size = 8),
        plot.caption = element_text(size = 7, colour = FIG_SEC, hjust = 0))
ggsave(file.path(FIG_DIR, "ExtFig4e_diffATAC_volcano.pdf"), p, width = 4.6, height = 3.4)
ggsave(file.path(FIG_DIR, "ExtFig4e_diffATAC_volcano.png"), p, width = 4.6, height = 3.4, dpi = 300)

cat(sprintf("PCA on %d peaks; differential peaks (q <= %.1f): %d (%d up, %d down in SCZ); BH FDR <= 0.1: %d\n",
            nrow(pca$rotation), Q, sum(res$diff), sum(res$diff & res$logFC > 0), sum(res$diff & res$logFC < 0), sum(res$FDR <= 0.1)))
