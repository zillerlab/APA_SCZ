#!/usr/bin/env Rscript
## ===========================================================================
## 01_dapa_iN.R
##
## Differential 3'UTR alternative polyadenylation (dAPA) between iNs from
## individuals with SCZ, BD or MDD and healthy controls (HC).
##
## Manuscript: Fig. 2c (SCZ vs HC), Extended Data Fig. 6b (BD, MDD vs HC),
##             Supplementary Table 3 (DiffAPA_iPSC_PM, iN columns).
##
## USAGE
##   Rscript 03_differential_apa/01_dapa_iN.R [--cores <n>] [--out <dir>]
##
## INPUTS
##   resources/dapa/utr_counts_iN_cohort.txt          read counts on the proximal and
##                                                    distal UTR regions, 93 iN libraries
##   resources/rnaseq/rnaseq_iN_cohort_metadata.txt   sample metadata, libraries with dapa = TRUE (diagnosis, sex,
##                                                    TechnicalPC1, size factors, ...)
##   results/apa_pas_library/dominant_pas_utrs_hg38.txt  UTR region annotation (module 02)
##
## OUTPUTS (results/dapa/)
##   dAPA_iN_<SCZ|BD|MDD>_vs_Ctrl.txt   one row per tested gene:
##       ID, method, n_ref, n_comp, beta, se, p_wald, p_lrt, delta, mean_expr,
##       qvalue, padj_BH
##       p_wald / qvalue are the values reported in the manuscript (q <= 0.05);
##       delta > 0 = higher distal 3'UTR usage in controls
##   apa_ratio_iN_cohort.txt            distal usage ratio R per gene and sample
##
## MODEL (R/dapa_functions.R)
##   R = distal / (distal + proximal) normalised reads per gene and sample;
##   beta regression R ~ condition + sex + TechnicalPC1 (Wald test), genes with
##   >= 10 usable samples per group; quasibinomial GLM if only one group has
##   >= 10; genes kept at mean max(proximal, distal) normalised count >= 50;
##   Storey q-values.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "dapa_functions.R"))

CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "dapa"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COVARIATES    <- c("sex", "TechnicalPC1")
MIN_PER_GROUP <- 10
MIN_EXPR      <- 50
CONTRASTS     <- list(c("Ctrl", "SCZ"), c("Ctrl", "BD"), c("Ctrl", "MDD"))

## ------------------------------------------------------------ inputs
cnt  <- read.delim(need_file(file.path(REPO, "resources", "dapa", "utr_counts_iN_cohort.txt"), "resources"),
                   row.names = 1, check.names = FALSE)
meta <- read.delim(need_file(file.path(REPO, "resources", "rnaseq", "rnaseq_iN_cohort_metadata.txt"), "resources"),
                   stringsAsFactors = FALSE)
meta <- meta[meta$dapa, ]                      # the 93 libraries of the APA analysis
ann  <- read.delim(need_file(file.path(REPO, "results", "apa_pas_library", "dominant_pas_utrs_hg38.txt"),
                             "module 02_apa_pas_library"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
stopifnot(setequal(colnames(cnt), meta$sample))
cnt <- as.matrix(cnt[, meta$sample])

ap <- apa_ratios(cnt, ann, setNames(meta$size_factor, meta$sample))
write.table(data.frame(ID = rownames(ap$R), round(ap$R, 6), check.names = FALSE),
            file.path(OUT_DIR, "apa_ratio_iN_cohort.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ contrasts
for (g in CONTRASTS) {
  s  <- meta$sample[meta$diagnosis %in% g]
  si <- data.frame(condition = meta[s, "diagnosis"], meta[s, COVARIATES], row.names = s)
  d  <- diff_apa(ap$R[, s], si, COVARIATES, g, min_per_group = MIN_PER_GROUP, cores = CORES)
  r  <- finalize_dapa(d, ap$expr, si, min_expr = MIN_EXPR)
  r  <- r[order(r$p_wald), ]
  f  <- file.path(OUT_DIR, sprintf("dAPA_iN_%s_vs_%s.txt", g[2], g[1]))
  write.table(r, f, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s vs %s: %d genes tested, %d at q <= 0.05 -> %s\n", g[2], g[1], nrow(r),
              sum(r$qvalue <= 0.05), f))
}
