#!/usr/bin/env Rscript
## ===========================================================================
## 06_power_iN_cohort.R
##
## Sensitivity of the SCZ vs HC differential expression analysis in the iN cohort:
##   (1) model-based power curves (negative binomial, RnaSeqSampleSize) with the
##       dispersion estimated from the cohort;
##   (2) donor-level resampling: differential expression in random subsets of donors,
##       compared with the full cohort.
##
## Manuscript: Extended Data Fig. 4c (model-based), 4d (resampling); main text.
##
## USAGE
##   Rscript 05_differential_expression/06_power_iN_cohort.R [--cores <n>] [--out <dir>]
##
## INPUTS
##   resources/rnaseq/gene_counts_iN_cohort.txt.gz, rnaseq_iN_cohort_metadata.txt,
##   gene_annotation_ensembl.txt                       (as for 01_deseq_iN_cohort.R)
##   results/deg/DEG_iN_SCZ_vs_Ctrl.txt                01_deseq_iN_cohort.R
##
## OUTPUTS
##   results/deg/Power_model_based.txt        power per FDR, log2 fold change and group size
##   results/deg/Power_resampling.txt         one row per resampling iteration
##   results/deg/Power_resampling_summary.txt medians per group size
##   results/figures/ExtFig4c_power_model.pdf / .png
##   results/figures/ExtFig4d_power_resampling.pdf / .png
##
## METHOD
##   Libraries: the SCZ and HC libraries of the differential expression analysis (one per
##   donor; 36 SCZ, 33 HC).
##   (1) Dispersion: edgeR common dispersion (RnaSeqSampleSize::est_count_dispersion, all
##   libraries, two groups) of the raw counts of the SCZ vs HC DEGs (01_deseq_iN_cohort.R).
##   Power curves: RnaSeqSampleSize::est_power_curve, n = 1-80 per group, mean count of
##   200 reads (lambda0), log2 fold change 0.4, 0.5 and 1, FDR 0.05, 0.1 and 0.01.
##   (2) Variance-stabilised expression (DESeq2 vst, blind) of all cohort libraries,
##   genes with > 20 and < 707,707.7 reads in total, annotated and not mitochondrial (as
##   for the DE analysis). Full cohort: limma ~ sex + site + DE_TechnicalPC1 + diagnosis;
##   DEGs: BH q <= 0.05, |log2FC| >= 0.3, mean expression > 6. Resampling: for n = 5-30
##   donors per group, 100 random draws of n SCZ and n HC donors (seed fixed per draw);
##   limma ~ site + diagnosis; DEGs as above. Per draw: DEGs found, fraction of the
##   full-cohort DEGs recovered, Jaccard index, Spearman correlation of log2FC with the
##   full cohort (all genes; full-cohort DEGs), and the fraction of full-cohort DEGs with
##   the same direction.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(DESeq2); library(limma); library(RnaSeqSampleSize) })

CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "deg"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
LOG2FC <- c(0.4, 0.5, 1); FDRS <- c(0.05, 0.1, 0.01); LAMBDA0 <- 200; N_MAX <- 80
SIZES <- c(5, 10, 15, 20, 25, 30); N_ITER <- 100
Q_DEG <- 0.05; LFC_DEG <- 0.3; MIN_EXPR <- 6

## ------------------------------------------------------------ inputs
rs   <- file.path(REPO, "resources", "rnaseq")
meta <- read.delim(need_file(file.path(rs, "rnaseq_iN_cohort_metadata.txt"), "resources"), stringsAsFactors = FALSE)
cnt  <- as.matrix(read.delim(need_file(file.path(rs, "gene_counts_iN_cohort.txt.gz"), "resources"),
                             row.names = 1, check.names = FALSE)[, meta$sample])
ann  <- read.delim(need_file(file.path(rs, "gene_annotation_ensembl.txt"), "resources"), stringsAsFactors = FALSE, na.strings = "")
deg  <- read.delim(need_file(file.path(OUT_DIR, "DEG_iN_SCZ_vs_Ctrl.txt"), "01_deseq_iN_cohort.R"), stringsAsFactors = FALSE)
cs   <- meta$diagnosis %in% c("SCZ", "Ctrl")
cat(sprintf("SCZ %d, HC %d libraries (one per donor)\n", sum(meta$diagnosis == "SCZ"), sum(meta$diagnosis == "Ctrl")))

## ------------------------------------------------------------ (1) model-based power
dg   <- deg$gene_id[deg$DEG != "NO"]
disp <- est_count_dispersion(cnt[rownames(cnt) %in% dg, cs], group = meta$diagnosis[cs], subSampleNum = sum(cs))
phi0 <- disp$common.dispersion
cat(sprintf("dispersion of the %d DEGs: %.4f (BCV %.3f)\n", length(dg), phi0, sqrt(phi0)))
pw <- do.call(rbind, lapply(LOG2FC, function(l) do.call(rbind, lapply(FDRS, function(f) {
  x <- est_power_curve(n = N_MAX, f = f, rho = 2^l, lambda0 = LAMBDA0, phi0 = phi0)$process
  data.frame(FDR = f, log2FC = l, n_per_group = x[, 1], power = x[, 2])
}))))
pw$dispersion <- phi0
write.table(pw, file.path(OUT_DIR, "Power_model_based.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
pd <- pw; pd$FDR <- factor(pd$FDR, levels = FDRS); pd$log2FC <- factor(sprintf("log2FC %.1f", pd$log2FC), levels = sprintf("log2FC %.1f", LOG2FC))
p <- ggplot(pd, aes(n_per_group, power, colour = FDR)) + geom_line(linewidth = 0.6) + geom_point(size = 0.9) +
  facet_wrap(~ log2FC, nrow = 1) + scale_colour_manual(values = c("0.05" = "#1f78c8", "0.1" = "#e0a020", "0.01" = "#c8321f")) +
  scale_y_continuous(limits = c(0, 1)) + labs(x = "Donors per group", y = "Power", colour = "FDR") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(FIG_DIR, "ExtFig4c_power_model.pdf"), p, width = 7, height = 2.6)
ggsave(file.path(FIG_DIR, "ExtFig4c_power_model.png"), p, width = 7, height = 2.6, dpi = 300)

## ------------------------------------------------------------ (2) resampling
tot  <- rowSums(cnt); mito <- ann$gene_id[grepl("^MT-", ann$gene_name)]
keep <- union(rownames(cnt)[tot > 20 & tot < 707707.7 & rownames(cnt) %in% ann$gene_id & !rownames(cnt) %in% mito],
              "ENSG00000078018")   # MAP2
dds  <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt, data.frame(row.names = meta$sample, all = factor(rep(1, ncol(cnt)))), ~ 1))
expr <- assay(vst(dds[keep, ], blind = TRUE))[, cs]
md   <- meta[cs, ]; md$diagnosis <- relevel(factor(md$diagnosis), ref = "Ctrl"); md$site <- factor(md$site)
coefn <- "diagnosisSCZ"
call_degs <- function(tt) rownames(tt)[tt$adj.P.Val <= Q_DEG & abs(tt$logFC) >= LFC_DEG & tt$AveExpr > MIN_EXPR]
full <- topTable(eBayes(lmFit(expr, model.matrix(~ sex + site + DE_TechnicalPC1 + diagnosis, md))), coef = coefn, number = Inf, sort.by = "none")
full_deg <- call_degs(full); full_lfc <- setNames(full$logFC, rownames(full))
cat(sprintf("%d genes; full-cohort limma DEGs: %d\n", nrow(expr), length(full_deg)))

scz <- unique(md$sample[md$diagnosis == "SCZ"]); hc <- unique(md$sample[md$diagnosis == "Ctrl"])
one_draw <- function(n, it) {
  set.seed(100000 + n * 1000 + it)
  k <- md$sample %in% c(sample(scz, n), sample(hc, n)); m <- droplevels(md[k, ])
  X <- model.matrix(~ site + diagnosis, m); if (qr(X)$rank < ncol(X)) return(NULL)
  tt <- topTable(eBayes(lmFit(expr[, k], X)), coef = coefn, number = Inf, sort.by = "none")
  sig <- call_degs(tt); lfc <- setNames(tt$logFC, rownames(tt))
  data.frame(n_per_group = n, iter = it, n_deg = length(sig), n_recovered = length(intersect(sig, full_deg)),
             frac_recovered = length(intersect(sig, full_deg)) / length(full_deg),
             jaccard = length(intersect(sig, full_deg)) / length(union(sig, full_deg)),
             lfc_cor_all = cor(lfc, full_lfc[names(lfc)], method = "spearman"),
             lfc_cor_full_degs = cor(lfc[full_deg], full_lfc[full_deg], method = "spearman"),
             direction_concordance_full_degs = mean(sign(lfc[full_deg]) == sign(full_lfc[full_deg])))
}
grid <- expand.grid(n = SIZES, it = seq_len(N_ITER))
rs_res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) one_draw(grid$n[i], grid$it[i]), mc.cores = CORES))
rs_res <- rs_res[order(rs_res$n_per_group, rs_res$iter), ]
write.table(rs_res, file.path(OUT_DIR, "Power_resampling.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
sm <- aggregate(rs_res[, -(1:2)], list(n_per_group = rs_res$n_per_group), median)
sm$n_draws <- as.vector(table(rs_res$n_per_group))
write.table(sm, file.path(OUT_DIR, "Power_resampling_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

long <- rbind(data.frame(metric = "Full-cohort DEGs recovered", n = rs_res$n_per_group, value = rs_res$frac_recovered),
              data.frame(metric = "log2FC correlation (full-cohort DEGs)", n = rs_res$n_per_group, value = rs_res$lfc_cor_full_degs),
              data.frame(metric = "Direction concordance (full-cohort DEGs)", n = rs_res$n_per_group, value = rs_res$direction_concordance_full_degs))
long$metric <- factor(long$metric, levels = unique(long$metric))
p <- ggplot(long, aes(factor(n), value)) + geom_boxplot(outlier.size = 0.4, linewidth = 0.3, fill = "grey90") +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") + labs(x = "Donors per group", y = NULL) +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(FIG_DIR, "ExtFig4d_power_resampling.pdf"), p, width = 7, height = 2.6)
ggsave(file.path(FIG_DIR, "ExtFig4d_power_resampling.png"), p, width = 7, height = 2.6, dpi = 300)

print(sm, row.names = FALSE, digits = 3)
print(pw[pw$n_per_group %in% c(20, 30, 35, 40), c("FDR", "log2FC", "n_per_group", "power")], row.names = FALSE)
