#!/usr/bin/env Rscript
## ===========================================================================
## 01_deseq_iN_cohort.R
##
## Differential gene expression between iNs (day 49) from individuals with SCZ,
## BD or MDD and healthy controls (HC).
##
## Manuscript: Fig. 1d (SCZ vs HC), Extended Data Fig. 4f,g (SCZ, BD, MDD vs HC),
##             Supplementary Table 3.
##
## USAGE
##   Rscript 05_differential_expression/01_deseq_iN_cohort.R [--out <dir>]
##
## INPUTS
##   resources/rnaseq/gene_counts_iN_cohort.txt.gz   gene-level read counts (GENCODE v27,
##                                                   featureCounts), 95 iN libraries
##   resources/rnaseq/rnaseq_iN_cohort_metadata.txt  diagnosis, sex, site, DE_TechnicalPC1
##   resources/rnaseq/gene_annotation_ensembl.txt    Ensembl gene ID -> gene symbol
##
## OUTPUTS (results/deg/)
##   DEG_iN_<SCZ|BD|MDD>_vs_Ctrl.txt  one row per tested gene:
##       gene_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, DEG
##       log2FoldChange = log2(HC / case), i.e. positive = lower expression in cases
##       (orientation of Fig. 1d and Source Data); DEG = Up / Down in cases, or NO
##   DEG_iN_overlap.txt               genes differential in at least one contrast
##                                    (Ext. Data Fig. 4f,g)
##
## METHOD
##   1  Size factors (DESeq2 median-of-ratios) on all quantified genes.
##   2  Genes kept if total count over all libraries > 3 (20th percentile) and
##      < 707,707.7 (99th percentile), annotated in Ensembl and not mitochondrial.
##      MAP2 (above the upper cut-off) is kept as neuronal marker.
##   3  DESeq2 negative binomial GLM ~ sex + site + DE_TechnicalPC1 + diagnosis,
##      Wald test for each diagnosis vs HC, Benjamini-Hochberg FDR.
##   4  DEG: padj <= 0.05, |log2FC| >= 0.4, baseMean >= 50.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
suppressPackageStartupMessages(library(DESeq2))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "deg"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_TOTAL  <- 3            # total count across libraries must exceed this
MAX_TOTAL  <- 707707.7     # ... and stay below this
KEEP_GENES <- "ENSG00000078018"   # MAP2
FDR <- 0.05; LFC <- 0.4; MIN_MEAN <- 50
DIAGNOSES <- c("SCZ", "BD", "MDD")

## ------------------------------------------------------------ inputs
rs   <- file.path(REPO, "resources", "rnaseq")
meta <- read.delim(need_file(file.path(rs, "rnaseq_iN_cohort_metadata.txt"), "resources"), stringsAsFactors = FALSE)
cnt  <- read.delim(need_file(file.path(rs, "gene_counts_iN_cohort.txt.gz"), "resources"),
                   row.names = 1, check.names = FALSE)
ann  <- read.delim(need_file(file.path(rs, "gene_annotation_ensembl.txt"), "resources"),
                   stringsAsFactors = FALSE, na.strings = "")
cnt  <- as.matrix(cnt[, meta$sample])
cd   <- data.frame(sex = factor(meta$sex), site = factor(meta$site),
                   DE_TechnicalPC1 = meta$DE_TechnicalPC1,
                   diagnosis = factor(meta$diagnosis),   # alphabetical factor levels; the reference level
                   row.names = meta$sample)               # affects the fit of low-count genes
cat(sprintf("%d libraries: %s\n", ncol(cnt), paste(names(table(cd$diagnosis)), table(cd$diagnosis), collapse = ", ")))

## ------------------------------------------------------------ 1-2 size factors, gene filter
dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ sex + site + DE_TechnicalPC1 + diagnosis)
dds <- estimateSizeFactors(dds)
tot  <- rowSums(cnt)
mito <- ann$gene_id[grepl("^MT-", ann$gene_name)]
keep <- rownames(cnt)[tot > MIN_TOTAL & tot < MAX_TOTAL & rownames(cnt) %in% ann$gene_id & !rownames(cnt) %in% mito]
keep <- union(keep, KEEP_GENES)
dds  <- dds[keep, ]
cat(sprintf("%d of %d genes tested\n", nrow(dds), nrow(cnt)))

## ------------------------------------------------------------ 3 DESeq2
dds <- DESeq(dds, quiet = TRUE)

## ------------------------------------------------------------ 4 results
res <- lapply(setNames(DIAGNOSES, DIAGNOSES), function(D) {
  r <- as.data.frame(results(dds, contrast = c("diagnosis", "Ctrl", D), test = "Wald"))
  sig <- which(r$padj <= FDR & abs(r$log2FoldChange) >= LFC & r$baseMean >= MIN_MEAN)
  deg <- rep("NO", nrow(r)); deg[sig] <- ifelse(r$log2FoldChange[sig] > 0, "Down", "Up")
  out <- data.frame(gene_id = rownames(r), gene_name = ann$gene_name[match(rownames(r), ann$gene_id)],
                    r[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")], DEG = deg)
  out <- out[order(out$gene_id), ]
  write.table(out, file.path(OUT_DIR, sprintf("DEG_iN_%s_vs_Ctrl.txt", D)), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  cat(sprintf("%-3s vs HC: %d DEGs (%d up, %d down in %s)\n", D, length(sig), sum(deg == "Up"), sum(deg == "Down"), D))
  out
})

## Ext. Data Fig. 4f,g: overlap and log2FC of genes differential in any contrast
sig <- sapply(res, function(r) r$DEG != "NO")
lfc <- sapply(res, function(r) r$log2FoldChange)
any <- rowSums(sig) > 0
ov  <- data.frame(gene_id = res$SCZ$gene_id, gene_name = res$SCZ$gene_name,
                  setNames(as.data.frame(sig), paste0("DEG_", DIAGNOSES)),
                  setNames(as.data.frame(lfc), paste0("log2FC_", DIAGNOSES)))[any, ]
write.table(ov, file.path(OUT_DIR, "DEG_iN_overlap.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n%d genes differential in at least one contrast; Pearson r of log2FC across them:\n", sum(any)))
print(round(cor(lfc[any, ]), 3))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_deseq_iN_cohort.txt"))
