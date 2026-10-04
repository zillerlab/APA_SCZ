#!/usr/bin/env Rscript
## ===========================================================================
## 05_deseq_microrna_iN_cohort.R
##
## microRNA expression (small RNA-seq) in iNs (day 49): PCA and
## differential expression between SCZ and HC.
##
## Manuscript: Fig. 1c (microRNA PCA), Fig. 1e, Supplementary Table 3 (DiffMicroRNAs).
##
## USAGE
##   Rscript 05_differential_expression/05_deseq_microrna_iN_cohort.R [--out <dir>]
##
## INPUTS (resources/microrna/)
##   mature_counts_iN_cohort.txt.gz    raw counts of mature microRNAs (miRBase v22),
##                                     2,657 microRNAs x 101 QC-passing libraries
##                                     (columns <PluriCore ID>_<replicate>)
##   microrna_iN_cohort_metadata.txt   library, line (PluriCore ID), diagnosis, sex, site
##                                     (differentiation site), in_analysis (one library per
##                                     line, 80 libraries)
##
## OUTPUTS
##   results/deg/DEmiR_iN_SCZ_vs_Ctrl.txt   one row per tested microRNA: microRNA, baseMean,
##       log2FoldChange, lfcSE, stat, pvalue, padj, class
##       log2FoldChange = log2(HC / SCZ), i.e. positive = lower in SCZ (orientation of
##       Fig. 1d,e); class = Up / Down in SCZ for padj <= 0.01, |log2FC| >= 0.5,
##       baseMean >= 10 (highlighted in Fig. 1e), otherwise NO
##   results/figures/Fig1c_microRNA_PCA.pdf / .txt
##   results/figures/Fig1e_microRNA_volcano.pdf / .png
##
## METHOD
##   0  One library per line (in_analysis).
##   1  Size factors (DESeq2 median-of-ratios) on all microRNAs.
##   2  microRNAs with a total count > 100 over all libraries are tested.
##   3  PCA: variance-stabilising transformation (blind), prcomp centred and scaled.
##   4  DESeq2 ~ sex + site + diagnosis (SCZ, HC, BD, MDD), Wald test HC vs SCZ,
##      Benjamini-Hochberg FDR. Differential: padj <= 0.05.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages(library(DESeq2))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "deg"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
MIN_TOTAL <- 100; FDR <- 0.05
FIG_FDR <- 0.01; FIG_LFC <- 0.5; FIG_MEAN <- 10

## ------------------------------------------------------------ inputs
rs   <- file.path(REPO, "resources", "microrna")
meta <- read.delim(need_file(file.path(rs, "microrna_iN_cohort_metadata.txt"), "resources"), stringsAsFactors = FALSE)
meta <- meta[as.logical(meta$in_analysis), ]
stopifnot(!anyDuplicated(meta$line))
cnt  <- read.delim(need_file(file.path(rs, "mature_counts_iN_cohort.txt.gz"), "resources"), row.names = 1, check.names = FALSE)
cnt  <- as.matrix(cnt[, meta$library])
cd   <- data.frame(sex = factor(meta$sex), site = factor(meta$site), diagnosis = factor(meta$diagnosis),
                   row.names = meta$library)
cat(sprintf("%d libraries (%d lines): %s\n", ncol(cnt), length(unique(meta$line)),
            paste(names(table(cd$diagnosis)), table(cd$diagnosis), collapse = ", ")))

## ------------------------------------------------------------ 1-2 size factors, filter
dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ sex + site + diagnosis)
dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(cnt) > MIN_TOTAL, ]
cat(sprintf("%d of %d microRNAs tested\n", nrow(dds), nrow(cnt)))

## ------------------------------------------------------------ 3 PCA (Fig. 1c)
vsd <- assay(vst(dds, nsub = 300, blind = TRUE))
pca <- prcomp(t(vsd), center = TRUE, scale. = TRUE); ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pcs <- data.frame(library = meta$library, line = meta$line, diagnosis = meta$diagnosis, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write.table(pcs, file.path(FIG_DIR, "Fig1c_microRNA_PCA.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(pcs, aes(PC1, PC2, colour = diagnosis)) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(Ctrl = "grey55", SCZ = "#c8321f", BD = "#1f78c8", MDD = "#e0a020")) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL,
       title = sprintf("microRNA (small RNA-seq, n = %d)", ncol(cnt))) +
  theme_classic(base_size = 10) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(FIG_DIR, "Fig1c_microRNA_PCA.pdf"), p, width = 3.6, height = 3)

## ------------------------------------------------------------ 4 DESeq2
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("diagnosis", "Ctrl", "SCZ")))
res <- data.frame(microRNA = rownames(res), res, row.names = NULL)
hit <- !is.na(res$padj) & res$padj <= FIG_FDR & abs(res$log2FoldChange) >= FIG_LFC & res$baseMean >= FIG_MEAN
res$class <- ifelse(!hit, "NO", ifelse(res$log2FoldChange > 0, "Down", "Up"))   # direction in SCZ
res <- res[order(res$pvalue), ]
write.table(res, file.path(OUT_DIR, "DEmiR_iN_SCZ_vs_Ctrl.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ------------------------------------------------------------ volcano (Fig. 1e)
vd  <- res[!is.na(res$pvalue), ]
vd$grp <- ifelse(vd$class != "NO", vd$class, ifelse(!is.na(vd$padj) & vd$padj <= FDR, "FDR", "NO"))
vd  <- vd[order(match(vd$grp, c("NO", "FDR", "Down", "Up"))), ]
lab <- vd[vd$class != "NO", ]; lab$name <- sub("^hsa-", "", lab$microRNA)
xl  <- 3   # x axis; microRNAs beyond (low counts, not significant) are drawn at the border
vd$log2FoldChange_plot <- pmax(pmin(vd$log2FoldChange, xl), -xl); lab$log2FoldChange_plot <- lab$log2FoldChange
p <- ggplot(vd, aes(log2FoldChange_plot, -log10(pvalue), colour = grp)) +
  geom_vline(xintercept = c(-FIG_LFC, FIG_LFC), colour = FIG_GRID, linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 1, stroke = 0) +
  ggrepel::geom_text_repel(data = lab, aes(label = name), size = 2.2, show.legend = FALSE, max.overlaps = Inf,
                           segment.size = 0.2, segment.colour = FIG_MUT, min.segment.length = 0, box.padding = 0.25, seed = 1) +
  scale_colour_manual(values = c(NO = "grey80", FDR = "grey50", Up = "#c8321f", Down = "#1f78c8"),
                      breaks = c("Down", "Up", "FDR"),
                      labels = c(Down = "lower in SCZ", Up = "higher in SCZ",
                                 FDR = sprintf("padj <= %.2f", FDR))) +
  scale_x_continuous(limits = c(-xl, xl)) +
  labs(x = "log2 fold change (HC / SCZ)", y = "-log10 p", colour = NULL,
       caption = sprintf("%d of %d microRNAs padj <= %.2f; coloured: padj <= %.2f, |log2FC| >= %.1f",
                         sum(res$padj <= FDR, na.rm = TRUE), nrow(res), FDR, FIG_FDR, FIG_LFC)) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = FIG_SEC), legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.key.height = unit(0.35, "cm"), legend.text = element_text(size = 8),
        plot.caption = element_text(size = 7, colour = FIG_SEC, hjust = 0))
ggsave(file.path(FIG_DIR, "Fig1e_microRNA_volcano.pdf"), p, width = 4.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Fig1e_microRNA_volcano.png"), p, width = 4.6, height = 3.6, dpi = 300)

cat(sprintf("differential (padj <= %.2f): %d (%d lower, %d higher in SCZ); highlighted in Fig. 1e: %d\n", FDR,
            sum(res$padj <= FDR, na.rm = TRUE), sum(res$padj <= FDR & res$log2FoldChange > 0, na.rm = TRUE),
            sum(res$padj <= FDR & res$log2FoldChange < 0, na.rm = TRUE), sum(hit)))
