#!/usr/bin/env Rscript
## ===========================================================================
## 02_deseq_rbp_knockdown.R
##
## Differential gene expression after shRNA knockdown of PTBP2, PCBP2 and CPSF6
## vs scramble control in iNs, and the overlap of the resulting DEG sets.
##
## Manuscript: Fig. 5d; Supplementary Table 3.
##
## USAGE
##   Rscript 05_differential_expression/02_deseq_rbp_knockdown.R [--out <dir>]
##
## INPUTS (resources/rbp_knockdown/)
##   gene_counts_rbp_knockdown.txt.gz   gene-level read counts (GENCODE v27), 20 libraries
##   rbp_knockdown_samples.txt          donor (iPSC line ID), condition, batch,
##                                      and the knockdown arm(s) each library belongs to
##
## OUTPUTS (results/deg/)
##   DEG_RBP_KD_<PTBP2|PCBP2|CPSF6>.txt  one row per tested gene:
##       gene_id, gene_name, baseMean, log2FoldChange (knockdown / control), lfcShrink
##       (apeglm), lfcSE, stat, pvalue, padj (BH), qvalue (Storey), DEG (Up / Down / NO)
##   DEG_RBP_KD_overlap.txt              DEG counts per knockdown and pairwise overlaps
##                                       (Fig. 5d)
##
## METHOD
##   Each knockdown is compared with the scramble controls of the same donors and
##   experiment (column in_<RBP>; controls can serve more than one arm). The same libraries
##   are used for differential APA (module 03). block = donor.
##   1  Genes with >= 10 reads in at least as many libraries as the smaller group.
##   2  DESeq2 ~ block + condition, Wald test, knockdown vs control; BH-adjusted
##      p-values (alpha = 0.05 for independent filtering); apeglm-shrunken log2FC and
##      Storey q-values are reported alongside.
##   3  DEG: padj <= 0.05 and |log2FC| >= 0.4.
##   4  Pairwise overlap of DEG sets, Fisher's exact test over genes tested in at
##      least one knockdown (as for dAPA, module 03). Venn diagram: 04_fig5d_rbp_knockdown.R.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
suppressPackageStartupMessages({ library(DESeq2); library(qvalue) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "deg"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RBPS <- c(PTBP2 = "ENSG00000117569", PCBP2 = "ENSG00000197111", CPSF6 = "ENSG00000111605")
MIN_READS <- 10; FDR <- 0.05; LFC <- 0.4

## ------------------------------------------------------------ inputs
rs  <- file.path(REPO, "resources", "rbp_knockdown")
si0 <- read.delim(need_file(file.path(rs, "rbp_knockdown_samples.txt"), "resources"), stringsAsFactors = FALSE)
cm  <- read.delim(need_file(file.path(rs, "gene_counts_rbp_knockdown.txt.gz"), "resources"),
                  check.names = FALSE, stringsAsFactors = FALSE)
cts <- as.matrix(cm[, si0$sample]); rownames(cts) <- cm$gene_id
gene_name <- setNames(cm$gene_name, cm$gene_id)

## ------------------------------------------------------------ one knockdown
run_kd <- function(A) {
  si <- si0[si0[[paste0("in_", A)]], ]
  rownames(si) <- si$sample
  si$condition <- relevel(factor(si$condition), ref = "CTRL")
  si$block <- factor(si$donor)
  cat(sprintf("\n%s knockdown: %d libraries\n", A, nrow(si)))
  print(table(condition = si$condition, block = si$block))
  dds <- suppressMessages(DESeqDataSetFromMatrix(cts[, si$sample], si[, c("condition", "block")],
                                                 design = ~ block + condition))
  dds <- dds[rowSums(counts(dds) >= MIN_READS) >= min(table(si$condition)), ]
  dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  r  <- as.data.frame(results(dds, contrast = c("condition", A, "CTRL"), alpha = FDR))
  sh <- lfcShrink(dds, coef = paste0("condition_", A, "_vs_CTRL"), type = "apeglm", quiet = TRUE)
  r$lfcShrink <- sh$log2FoldChange[match(rownames(r), rownames(sh))]
  ok <- !is.na(r$pvalue); r$qvalue <- NA_real_; r$qvalue[ok] <- qvalue(r$pvalue[ok])$qvalues
  sig <- which(r$padj <= FDR & abs(r$log2FoldChange) >= LFC)
  deg <- rep("NO", nrow(r)); deg[sig] <- ifelse(r$log2FoldChange[sig] > 0, "Up", "Down")
  out <- data.frame(gene_id = rownames(r), gene_name = gene_name[rownames(r)],
                    r[, c("baseMean", "log2FoldChange", "lfcShrink", "lfcSE", "stat", "pvalue", "padj", "qvalue")],
                    DEG = deg)
  out <- out[order(out$pvalue), ]
  write.table(out, file.path(OUT_DIR, sprintf("DEG_RBP_KD_%s.txt", A)), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  tg <- out[out$gene_id == RBPS[[A]], ]
  cat(sprintf("%d genes tested; %d DEGs (%d up, %d down); %s itself: log2FC %.2f, padj %.2g\n",
              nrow(out), length(sig), sum(deg == "Up"), sum(deg == "Down"), A, tg$log2FoldChange, tg$padj))
  out
}
res <- lapply(setNames(names(RBPS), names(RBPS)), run_kd)

## ------------------------------------------------------------ Fig. 5d overlaps
univ <- Reduce(union, lapply(res, `[[`, "gene_id"))   # genes tested in at least one knockdown
degs <- lapply(res, function(r) intersect(r$gene_id[r$DEG != "NO"], univ))
pairs <- combn(names(RBPS), 2, simplify = FALSE)
ov <- do.call(rbind, lapply(pairs, function(p) {
  a <- univ %in% degs[[p[1]]]; b <- univ %in% degs[[p[2]]]
  ft <- fisher.test(table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))))
  data.frame(set_1 = p[1], set_2 = p[2], n_1 = sum(a), n_2 = sum(b), n_overlap = sum(a & b),
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}))
ov <- rbind(ov, data.frame(set_1 = "all three", set_2 = "", n_1 = NA, n_2 = NA,
                           n_overlap = length(Reduce(intersect, degs)), odds_ratio = NA, p_value = NA))
ov$n_universe <- length(univ)
write.table(ov, file.path(OUT_DIR, "DEG_RBP_KD_overlap.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("\nDEGs among %d genes tested in at least one knockdown: %s\n", length(univ),
            paste(names(degs), lengths(degs), collapse = ", ")))
print(ov[, 1:7], row.names = FALSE, digits = 3)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_deseq_rbp_knockdown.txt"))
