#!/usr/bin/env Rscript
## ===========================================================================
## 05_prs_extremes_iN.R
##
## Differential APA between iN donors with the highest and lowest SCZ polygenic
## risk, independent of diagnosis, and its concordance with SCZ vs HC dAPA.
##
## Manuscript: Fig. 4d,e; Extended Data Fig. 6f,h.
##
## USAGE
##   Rscript 04_prs_association/05_prs_extremes_iN.R --prs <file> [--cores <n>]
##          [--out <dir>] [--private-out <dir>]
##   --prs   as for 01_prs_apa_iN.R (sample = PluriCore ID, PRS_SCZ, ...); private
##
## OUTPUTS
##   <out>/dAPA_iN_PRShigh_vs_PRSlow.txt            gene-level dAPA, high vs low SCZ-PRS
##   <out>/prs_extremes_concordance_iN.txt          statistics of Fig. 4e, ED Fig. 6f,h
##   <out>/prs_extremes_beta_iN.txt                 per-gene effects plotted in Fig. 4e
##   <private-out>/prs_extremes_groups_iN.txt       PRS and group per library (Fig. 4d)
##
## METHOD
##   1  Extreme groups: libraries with SCZ-PRS below the 10th ("Low") or above the
##      90th percentile ("High") of the 80 iN donors with PRS, all diagnoses.
##   2  dAPA High vs Low with the model of 03_differential_apa (R ~ group + sex +
##      TechnicalPC1), >= 3 libraries per group, expression filter mean max
##      normalised count >= 50 over the tested libraries, Storey q-values.
##   3  Fig. 4e: correlation of the condition coefficients (High vs Low; SCZ vs HC)
##      across genes tested in both.
##   4  ED Fig. 6h: overlap of genes with nominal p <= 0.05 (High vs Low) and genes
##      with SCZ vs HC dAPA (q <= 0.05), Fisher's exact test across genes tested in
##      the extreme-group analysis.
##   5  ED Fig. 6f: cumulative APA (mean |residual| across High-vs-Low genes with
##      q <= 0.1) vs SCZ-PRS across all 80 donors.
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

## ------------------------------------------------------------ inputs
meta <- read.delim(file.path(REPO, "resources", "rnaseq", "rnaseq_iN_cohort_metadata.txt"), stringsAsFactors = FALSE)
meta <- meta[meta$dapa, ]                      # the 93 libraries of the APA analysis
rownames(meta) <- meta$sample
cnt  <- as.matrix(read.delim(file.path(REPO, "resources", "dapa", "utr_counts_iN_cohort.txt"),
                             row.names = 1, check.names = FALSE))[, meta$sample]
ann  <- read.delim(file.path(REPO, "results", "apa_pas_library", "dominant_pas_utrs_hg38.txt"), stringsAsFactors = FALSE)
dSCZ <- read.delim(file.path(REPO, "results", "dapa", "dAPA_iN_SCZ_vs_Ctrl.txt"), stringsAsFactors = FALSE)
prs  <- read.delim(PRS, stringsAsFactors = FALSE)
ap   <- apa_ratios(cnt, ann, setNames(meta$size_factor, meta$sample))

## ------------------------------------------------------------ 1 extreme groups
libs <- meta$sample[meta$sample %in% prs$sample]
x    <- prs$PRS_SCZ[match(libs, prs$sample)]
q    <- stats::quantile(x, probs = c(0.1, 0.9))
grp  <- ifelse(x < q[1], "Low", ifelse(x > q[2], "High", "None"))
write.table(data.frame(sample = libs, diagnosis = meta[libs, "diagnosis"], PRS_SCZ = x, group = grp),
            file.path(PRIV, "prs_extremes_groups_iN.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
s  <- libs[grp != "None"]
si <- data.frame(condition = grp[grp != "None"], meta[s, COVARIATES], row.names = s)
cat(sprintf("extreme groups: %d Low, %d High of %d donors with PRS\n", sum(grp == "Low"), sum(grp == "High"), length(libs)))

## ------------------------------------------------------------ 2 dAPA High vs Low
d <- diff_apa(ap$R[, s], si, COVARIATES, c("Low", "High"), min_per_group = 3, cores = CORES)
r <- finalize_dapa(d, ap$expr, si, min_expr = 50)
r <- r[order(r$p_wald), ]
write.table(r, file.path(OUT_DIR, "dAPA_iN_PRShigh_vs_PRSlow.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ 3 Fig. 4e
g  <- intersect(r$ID, dSCZ$ID)
bp <- r$beta[match(g, r$ID)]; bs <- dSCZ$beta[match(g, dSCZ$ID)]
ok <- is.finite(bp) & is.finite(bs)
write.table(data.frame(ID = g[ok], betaPRS = bp[ok], betaSCZ = bs[ok]),
            file.path(OUT_DIR, "prs_extremes_beta_iN.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
ct <- cor.test(bp[ok], bs[ok])

## ------------------------------------------------------------ 4 ED Fig. 6h
inPRS <- r$ID %in% r$ID[r$p_wald <= 0.05]
inSCZ <- r$ID %in% dSCZ$ID[dSCZ$qvalue <= 0.05]
ft <- fisher.test(table(factor(inPRS, c(TRUE, FALSE)), factor(inSCZ, c(TRUE, FALSE))))

## ------------------------------------------------------------ 5 ED Fig. 6f
siA <- data.frame(condition = ifelse(meta$diagnosis == "Ctrl", "Ctrl", "Case"), meta[, COVARIATES],
                  row.names = meta$sample)
dA  <- diff_apa(ap$R, siA, COVARIATES, c("Ctrl", "Case"), min_per_group = 10, cores = CORES)
res <- attr(dA, "residuals")
cumA <- colMeans(abs(res[intersect(r$ID[r$qvalue <= 0.1], rownames(res)), libs, drop = FALSE]), na.rm = TRUE)
a6f <- prs_association(cumA, x)

stats <- data.frame(
  panel = c("Fig4e", "EDFig6h", "EDFig6f"),
  test  = c("Pearson, beta High vs Low vs beta SCZ vs HC", "Fisher, p<=0.05 High vs Low vs q<=0.05 SCZ vs HC",
            "linear model, cumulative APA (mean |residual|, High-vs-Low genes q<=0.1) ~ SCZ-PRS"),
  n = c(sum(ok), length(inPRS), a6f[["n"]]),
  estimate = c(unname(ct$estimate), unname(ft$estimate), a6f[["pearson_r"]]),
  p_value  = c(ct$p.value, ft$p.value, a6f[["p_value"]]),
  n_overlap = c(NA, sum(inPRS & inSCZ), NA))
write.table(stats, file.path(OUT_DIR, "prs_extremes_concordance_iN.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("High vs Low: %d genes tested, %d at p <= 0.05, %d at q <= 0.1\n", nrow(r), sum(r$p_wald <= 0.05), sum(r$qvalue <= 0.1)))
print(stats[, -2], row.names = FALSE, digits = 4)
