#!/usr/bin/env Rscript
## ===========================================================================
## 03_prs_apa_dlpfc.R
##
## Association of SCZ polygenic risk with cumulative differential APA per donor
## in postmortem DLPFC (CommonMind Consortium).
##
## Manuscript: Extended Data Fig. 6i.
##
## Individual-level inputs are controlled access and are not distributed:
##   --prs        CMC PRS, tab-separated: sample (RNA-seq sample ID), PRS_SCZ,
##                (PRS_BD, PRS_MDD may be present; European donors with imputed genotypes)
##   --residuals  per-library residuals of the covariate-only model, written by
##                03_differential_apa/02_dapa_postmortem.R
##                  --cases SCZ,BP,AFF --tag allCases_vs_Control --residuals-out <file>
##   --cmc        CMC data directory (clinical table)
## Gene set (public): results/dapa/dAPA_DLPFC_SCZ_vs_Control_EUR_PRS.txt, from
##   02_dapa_postmortem.R --ethnicity Caucasian --prs <file> --tag SCZ_vs_Control_EUR_PRS
##   (dAPA at q <= 0.1 between 150 control and 173 SCZ donors of European
##   ancestry with PRS)
##
## USAGE
##   Rscript 04_prs_association/03_prs_apa_dlpfc.R --prs <file> --residuals <file> --cmc <dir>
##          [--out <dir>] [--private-out <dir>]
##
## OUTPUTS
##   <out>/prs_apa_association_dlpfc.txt        statistics of ED Fig. 6i (public)
##   <private-out>/prs_apa_per_donor_dlpfc.txt  per-donor values (private)
##
## METHOD
##   As for iNs (01_prs_apa_iN.R): cumulative APA = mean |z| across the dAPA genes,
##   z standardised per gene across all donors with PRS (all diagnoses), Tukey
##   filter on the mean |residual|, linear model cumulative APA ~ SCZ-PRS (left
##   panel); SCZ-PRS in SCZ vs controls, Wilcoxon test (right panel).
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "prs_functions.R"))

PRS     <- need_file(get_arg("--prs"), "--prs")
RESID   <- need_file(get_arg("--residuals"), "--residuals")
CMC     <- get_arg("--cmc"); if (is.null(CMC)) stop("set the CMC data directory with --cmc <dir>")
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "prs"))
PRIV    <- get_arg("--private-out", file.path(REPO, "results", "private"))
for (d in c(OUT_DIR, PRIV)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

Q_MAX  <- 0.1
SCORES <- "PRS_SCZ"

## ------------------------------------------------------------ inputs
res  <- as.matrix(read.delim(RESID, row.names = 1, check.names = FALSE))
prs  <- read.delim(PRS, stringsAsFactors = FALSE)
clin <- read.delim(need_file(file.path(CMC, "CMC_MSSM-Penn-Pitt_Clinical.txt"), "--cmc"), stringsAsFactors = FALSE)
dg   <- read.delim(file.path(REPO, "results", "dapa", "dAPA_DLPFC_SCZ_vs_Control_EUR_PRS.txt"), stringsAsFactors = FALSE)
genes <- dg$ID[dg$qvalue <= Q_MAX]

libs <- colnames(res)[colnames(res) %in% prs$sample]     # order of the residual matrix
P    <- prs[match(libs, prs$sample), SCORES, drop = FALSE]
cl   <- clin[match(libs, clin$DLPFC_RNA_Sequencing_Sample_ID), ]
Dx   <- factor(cl$Dx)
cat(sprintf("%d libraries with PRS | %d dAPA genes (q <= %.2f)\n", length(libs), length(genes), Q_MAX))

## ------------------------------------------------------------ association
cs <- cumulative_score(res, genes, libs); k <- cs$tukey_keep
a  <- prs_association(cs$mean_abs_z[k], P$PRS_SCZ[k])
w  <- wilcox.test(P$PRS_SCZ[Dx == "SCZ"], P$PRS_SCZ[Dx == "Control"])$p.value
stats <- data.frame(panel = c("EDFig6i_left", "EDFig6i_right"),
                    test  = c("cumulative APA ~ SCZ-PRS (linear model)", "SCZ-PRS, SCZ vs Control (Wilcoxon)"),
                    n_genes = c(length(genes), NA), n = c(a[["n"]], sum(Dx %in% c("SCZ", "Control"))),
                    pearson_r = c(a[["pearson_r"]], NA), slope = c(a[["slope"]], NA), p_value = c(a[["p_value"]], w))
write.table(stats, file.path(OUT_DIR, "prs_apa_association_dlpfc.txt"), sep = "	", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = libs, diagnosis = cl$Dx, PRS_SCZ = P$PRS_SCZ, mean_abs_z = cs$mean_abs_z,
                       tukey_keep = k),
            file.path(PRIV, "prs_apa_per_donor_dlpfc.txt"), sep = "	", quote = FALSE, row.names = FALSE)
print(stats[, -2], row.names = FALSE, digits = 3)
