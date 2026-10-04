#!/usr/bin/env Rscript
## ===========================================================================
## 01_prs_apa_iN.R
##
## Association of polygenic risk (SCZ, BD, MDD PRS; PRS-CS) with cumulative
## differential APA per iN donor.
##
## Manuscript: Fig. 4b; Extended Data Fig. 6c-e, 6g.
##
## The PRS are individual-level genetic data and are not distributed with this
## repository (controlled access, see Data availability). Supply them with
## --prs: a tab-separated file with columns
##     sample   PluriCore ID (as in resources/rnaseq/rnaseq_iN_cohort_metadata.txt)
##     PRS_SCZ, PRS_BD, PRS_MDD
##
## USAGE
##   Rscript 04_prs_association/01_prs_apa_iN.R --prs <file> [--cores <n>]
##          [--out <dir>] [--private-out <dir>]
##
## OUTPUTS
##   <out>/prs_apa_association_iN.txt         statistics of Fig. 4b and ED Fig. 6c-e (public)
##   <out>/prs_apa_permutation_null_iN.txt    1,000 permuted Pearson r (public; ED Fig. 6g)
##   <private-out>/prs_apa_per_donor_iN.txt   cumulative APA and PRS per library
##                                            (contains PRS; do not publish)
##
## METHOD
##   1  Covariate-adjusted APA: residuals of the covariate-only beta regression
##      R ~ sex + TechnicalPC1, fitted per gene across all 93 iN libraries
##      (R/dapa_functions.R; controls vs all cases for the testability rule).
##   2  Cumulative APA per library: mean |z| across the dAPA genes of a contrast,
##      z = residual standardised per gene across the libraries with a PRS.
##   3  Outlier libraries removed by Tukey's fences on the mean |residual|.
##   4  Linear model cumulative APA ~ PRS (t-test of the slope), Pearson r.
##   5  Empirical p: 1,000 random re-assignments of cumulative APA to donors
##      (seed 42).
##   Gene sets: SCZ dAPA genes (q <= 0.05; Fig. 4b, ED Fig. 6c,d) and BD genes
##   with p <= 0.05 in the BD vs HC comparison (ED Fig. 6e).
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "dapa_functions.R"))
source(file.path(REPO, "R", "prs_functions.R"))

PRS     <- need_file(get_arg("--prs"), "--prs")
CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "prs"))
PRIV    <- get_arg("--private-out", file.path(REPO, "results", "private"))
for (d in c(OUT_DIR, PRIV)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

COVARIATES <- c("sex", "TechnicalPC1")
SCORES     <- c("PRS_SCZ", "PRS_BD", "PRS_MDD")

## ------------------------------------------------------------ inputs
meta <- read.delim(file.path(REPO, "resources", "rnaseq", "rnaseq_iN_cohort_metadata.txt"),
                   stringsAsFactors = FALSE)
meta <- meta[meta$dapa, ]                      # the 93 libraries of the APA analysis
rownames(meta) <- meta$sample
cnt  <- as.matrix(read.delim(file.path(REPO, "resources", "dapa", "utr_counts_iN_cohort.txt"),
                             row.names = 1, check.names = FALSE))[, meta$sample]
ann  <- read.delim(file.path(REPO, "results", "apa_pas_library", "dominant_pas_utrs_hg38.txt"),
                   stringsAsFactors = FALSE)
dSCZ <- read.delim(file.path(REPO, "results", "dapa", "dAPA_iN_SCZ_vs_Ctrl.txt"), stringsAsFactors = FALSE)
dBD  <- read.delim(file.path(REPO, "results", "dapa", "dAPA_iN_BD_vs_Ctrl.txt"),  stringsAsFactors = FALSE)
prs  <- read.delim(PRS, stringsAsFactors = FALSE)
stopifnot(all(c("sample", SCORES) %in% names(prs)), !anyDuplicated(prs$sample))

GENE_SETS <- list(SCZ_dAPA_q05 = dSCZ$ID[dSCZ$qvalue <= 0.05],
                  BD_dAPA_p05  = dBD$ID[dBD$p_wald <= 0.05])

## ------------------------------------------------------------ 1 residuals
ap <- apa_ratios(cnt, ann, setNames(meta$size_factor, meta$sample))
si <- data.frame(condition = ifelse(meta$diagnosis == "Ctrl", "Ctrl", "Case"),
                 meta[, COVARIATES], row.names = meta$sample)
d  <- diff_apa(ap$R, si, COVARIATES, c("Ctrl", "Case"), min_per_group = 10, cores = CORES)
resid <- attr(d, "residuals")

## ------------------------------------------------------------ 2-4 association
libs <- meta$sample[meta$sample %in% prs$sample]          # library order of the metadata
P    <- prs[match(libs, prs$sample), SCORES]
cat(sprintf("%d of %d libraries with PRS
", length(libs), nrow(meta)))

CS <- lapply(GENE_SETS, function(g) cumulative_score(resid, g, libs))
PANELS <- data.frame(panel    = c("Fig4b", "EDFig6c", "EDFig6d", "EDFig6e"),
                     gene_set = c("SCZ_dAPA_q05", "SCZ_dAPA_q05", "SCZ_dAPA_q05", "BD_dAPA_p05"),
                     score    = c("PRS_SCZ", "PRS_BD", "PRS_MDD", "PRS_BD"))
stats <- do.call(rbind, lapply(seq_len(nrow(PANELS)), function(j) {
  cs <- CS[[PANELS$gene_set[j]]]; k <- cs$tukey_keep
  data.frame(PANELS[j, ], n_genes = length(GENE_SETS[[PANELS$gene_set[j]]]),
             t(prs_association(cs$mean_abs_z[k], P[[PANELS$score[j]]][k])))
}))

## 5 permutation null for Fig. 4b (Ext. Data Fig. 6g)
cs   <- CS$SCZ_dAPA_q05; k <- cs$tukey_keep
null <- permutation_null(cs$mean_abs_z[k], P$PRS_SCZ[k], n_perm = 1000, seed = 42)
stats$empirical_p <- NA
stats$empirical_p[stats$panel == "Fig4b"] <- mean(null > stats$pearson_r[stats$panel == "Fig4b"])

perdonor <- data.frame(sample = libs, diagnosis = meta[libs, "diagnosis"], P)
for (gs in names(CS)) {
  perdonor[[paste0(gs, "_mean_abs_z")]] <- CS[[gs]]$mean_abs_z
  perdonor[[paste0(gs, "_tukey_keep")]] <- CS[[gs]]$tukey_keep
}
write.table(stats, file.path(OUT_DIR, "prs_apa_association_iN.txt"), sep = "	", quote = FALSE, row.names = FALSE)
write.table(data.frame(Pearson = null), file.path(OUT_DIR, "prs_apa_permutation_null_iN.txt"),
            sep = "	", quote = FALSE, row.names = FALSE)
write.table(perdonor, file.path(PRIV, "prs_apa_per_donor_iN.txt"), sep = "	", quote = FALSE, row.names = FALSE)
print(stats[, c("panel", "gene_set", "score", "n", "pearson_r", "p_value", "empirical_p")], row.names = FALSE, digits = 3)
