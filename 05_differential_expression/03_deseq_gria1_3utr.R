#!/usr/bin/env Rscript
## ===========================================================================
## 03_deseq_gria1_3utr.R
##
## Differential gene expression in iNs carrying a deletion of the GRIA1 long 3'UTR
## (region between proximal and distal PAS; two clones) vs unedited iNs.
##
## Manuscript: Extended Data Fig. 7d; main text (97 DEGs).
##
## USAGE
##   Rscript 05_differential_expression/03_deseq_gria1_3utr.R [--out <dir>]
##
## INPUTS (resources/gria1_3utr/)
##   gene_counts_gria1_3utr.txt.gz   gene-level read counts (GENCODE v27), 5 libraries
##   gria1_3utr_samples.txt          line (PluriCore ID) and genotype (KO / WT)
##
## OUTPUTS (results/deg/)
##   DEG_GRIA1_3UTR_KO_vs_WT.txt     one row per tested gene: gene_id, gene_name, baseMean,
##       log2FoldChange (KO / WT), lfcSE, stat, pvalue, padj (BH), qvalue (Storey), DEG
##
## METHOD
##   KO: two edited clones of PSYLMUi002-A. WT: the unedited parental line and two
##   scramble-shRNA control libraries (PSYLMUi002-A, MPIPi013-A) of the RBP knockdown
##   experiment.
##   1  Size factors on all genes.
##   2  Genes tested if their total count is > 1 and < 5,000; MAP2 is kept.
##   3  DESeq2 ~ genotype, Wald test KO vs WT; Storey q-values over all tested genes.
##   4  DEG: q <= 0.05, |log2FC| >= 1 and baseMean >= 50.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
suppressPackageStartupMessages({ library(DESeq2); library(qvalue) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "deg"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
MIN_TOTAL <- 1; MAX_TOTAL <- 5000; KEEP_GENES <- "ENSG00000078018"   # MAP2
QV <- 0.05; LFC <- 1; MIN_MEAN <- 50
GRIA1 <- "ENSG00000155511"

## ------------------------------------------------------------ inputs
rs  <- file.path(REPO, "resources", "gria1_3utr")
si  <- read.delim(need_file(file.path(rs, "gria1_3utr_samples.txt"), "resources"), stringsAsFactors = FALSE)
cm  <- read.delim(need_file(file.path(rs, "gene_counts_gria1_3utr.txt.gz"), "resources"),
                  check.names = FALSE, stringsAsFactors = FALSE)
cts <- as.matrix(cm[, si$sample]); rownames(cts) <- cm$gene_id
gene_name <- setNames(cm$gene_name, cm$gene_id)
cd  <- data.frame(genotype = factor(si$genotype), row.names = si$sample)   # alphabetical factor levels (the reference
                                                                        # level affects low-count genes)
print(si[, c("sample", "line", "genotype")], row.names = FALSE)

## ------------------------------------------------------------ 1-2 size factors, gene filter
dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cts, cd, design = ~ genotype))
tot  <- rowSums(cts[, si$sample != "PSYLMUi002-A_GRIA1dUTR_cl13"])
keep <- union(rownames(cts)[tot > MIN_TOTAL & tot < MAX_TOTAL], KEEP_GENES)
dds  <- dds[keep, ]
cat(sprintf("%d genes tested\n", nrow(dds)))

## ------------------------------------------------------------ 3 DESeq2
dds <- DESeq(dds, quiet = TRUE)
r <- as.data.frame(results(dds, contrast = c("genotype", "KO", "WT"), test = "Wald"))
r$qvalue <- NA_real_; ok <- !is.na(r$pvalue); r$qvalue[ok] <- qvalue(r$pvalue[ok])$qvalues

## ------------------------------------------------------------ 4 DEGs
sig <- which(r$qvalue <= QV & abs(r$log2FoldChange) >= LFC & r$baseMean >= MIN_MEAN)
deg <- rep("NO", nrow(r)); deg[sig] <- ifelse(r$log2FoldChange[sig] > 0, "Up", "Down")
out <- data.frame(gene_id = rownames(r), gene_name = gene_name[rownames(r)],
                  r[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "qvalue")], DEG = deg)
out <- out[order(out$pvalue), ]
write.table(out, file.path(OUT_DIR, "DEG_GRIA1_3UTR_KO_vs_WT.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
g <- out[out$gene_id == GRIA1, ]
cat(sprintf("%d DEGs (%d up, %d down in KO); GRIA1: log2FC %.2f, q %.2g\n",
            length(sig), sum(deg == "Up"), sum(deg == "Down"), g$log2FoldChange, g$qvalue))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_deseq_gria1_3utr.txt"))
