#!/usr/bin/env Rscript
## ===========================================================================
## 03_dapa_rbp_knockdown.R
##
## Differential APA after shRNA knockdown of PTBP2, PCBP2 and CPSF6 in iNs vs scramble
## control.
##
## Manuscript: Fig. 5e-h,j (input), Supplementary Table 3.
##
## USAGE
##   Rscript 03_differential_apa/03_dapa_rbp_knockdown.R [--cores <n>] [--out <dir>]
##
## INPUTS (resources/rbp_knockdown/)
##   utr_counts_rbp_knockdown.txt.gz  read counts on the proximal and distal UTR segment
##                                    (featureID, geneId, class, Length, 20 libraries)
##   rbp_knockdown_samples.txt        donor, condition, batch, knockdown arm(s) (in_<RBP>),
##                                    size factor
##
## OUTPUTS (results/dapa/)
##   dAPA_RBP_KD_<PTBP2|PCBP2|CPSF6>.txt   one row per tested gene:
##       ID, beta (knockdown vs control, logit scale), se, p_wald, padj_BH, dAPA
##       (lengthening / shortening / NO), R_ctrl, R_kd (fitted distal usage per group),
##       deltaR, p_lrt, phi (precision), n_ref, n_comp, meanExpr
##   dAPA_RBP_KD_overlap.txt               dAPA counts and pairwise overlaps (Fisher's exact
##                                         test over the genes tested in at least one knockdown)
##   dAPA_RBP_KD_genes_<RBP>.txt           dAPA gene lists (input for SynGO, Fig. 5g)
##   dAPA_RBP_KD_genes_background.txt      genes tested in at least one knockdown
##
## METHOD (R/dapa_v2_functions.R)
##   R = distal / (distal + proximal) length-normalised read density per gene and library.
##   Per knockdown arm (libraries with in_<RBP> = TRUE: knockdowns and the scramble controls
##   of the same donors), beta regression R ~ condition + block, block = donor, >= 3
##   libraries per group; genes with mean normalised count >= 10; Wald test;
##   Benjamini-Hochberg adjustment. dAPA: adjusted P <= 1e-4. With 8-12 libraries per arm
##   the Wald test is anticonservative, and the cutoff is therefore stricter than for the
##   iN cohort.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "dapa_v2_functions.R"))

CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "dapa"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RBPS  <- c("PTBP2", "PCBP2", "CPSF6")
BH_KD <- 1e-4                                           # knockdown dAPA: BH-adjusted Wald P <= 1e-4

## ------------------------------------------------------------ inputs
rs  <- file.path(REPO, "resources", "rbp_knockdown")
si0 <- read.delim(need_file(file.path(rs, "rbp_knockdown_samples.txt"), "resources"), stringsAsFactors = FALSE)
uc  <- read.delim(need_file(file.path(rs, "utr_counts_rbp_knockdown.txt.gz"), "resources"),
                  check.names = FALSE, stringsAsFactors = FALSE)
cts <- as.matrix(uc[, si0$sample]); rownames(cts) <- uc$featureID
ap  <- build_apa_ratios_v2(cts, uc[, c("featureID", "geneId", "class", "Length")],
                           setNames(si0$size_factor, si0$sample))

## ------------------------------------------------------------ per knockdown
res <- lapply(setNames(RBPS, RBPS), function(A) {
  s  <- si0[si0[[paste0("in_", A)]], ]
  si <- data.frame(condition = s$condition, block = s$donor, row.names = s$sample)
  d  <- diff_apa_v2(ap$R[, rownames(si)], si, "block", c("CTRL", A), min_per_group = 3, cores = CORES)
  r  <- finalize_dapa_v2(d, ap, si, min_expr = 10)
  r$dAPA <- ifelse(r$padj_wald <= BH_KD, ifelse(r$beta > 0, "lengthening", "shortening"), "NO")
  r <- data.frame(ID = r$ID, beta = r$beta, se = r$se, p_wald = r$p_wald, padj_BH = r$padj_wald, dAPA = r$dAPA,
                  R_ctrl = r$fitted_ref, R_kd = r$fitted_comp, deltaR = r$deltaRatio_adj, p_lrt = r$p_lrt,
                  phi = r$phi, n_ref = r$n_ref, n_comp = r$n_comp, meanExpr = r$meanExpr)
  write.table(r, file.path(OUT_DIR, sprintf("dAPA_RBP_KD_%s.txt", A)), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  writeLines(r$ID[r$dAPA != "NO"], file.path(OUT_DIR, sprintf("dAPA_RBP_KD_genes_%s.txt", A)))
  r
})

## ------------------------------------------------------------ Fig. 5e overlaps
univ <- sort(Reduce(union, lapply(res, `[[`, "ID")))    # genes tested in at least one knockdown
writeLines(univ, file.path(OUT_DIR, "dAPA_RBP_KD_genes_background.txt"))
sets <- lapply(res, function(r) r$ID[r$dAPA != "NO"])
ov <- do.call(rbind, lapply(combn(RBPS, 2, simplify = FALSE), function(p) {
  a <- univ %in% sets[[p[1]]]; b <- univ %in% sets[[p[2]]]
  ft <- fisher.test(table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))))
  data.frame(set_1 = p[1], set_2 = p[2], n_1 = sum(a), n_2 = sum(b), n_overlap = sum(a & b),
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}))
ov <- rbind(ov, data.frame(set_1 = "all three", set_2 = "", n_1 = NA, n_2 = NA,
                           n_overlap = length(Reduce(intersect, sets)), odds_ratio = NA, p_value = NA))
ov$n_universe <- length(univ)
write.table(ov, file.path(OUT_DIR, "dAPA_RBP_KD_overlap.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("\ndAPA (BH <= %g): %s\n", BH_KD, paste(RBPS, lengths(sets), collapse = ", ")))
cat(sprintf("shortening: %s\n", paste(RBPS, sapply(res, function(r) sum(r$dAPA == "shortening")), collapse = ", ")))
cat(sprintf("overlaps among %d genes tested in at least one knockdown:\n", length(univ)))
print(ov[, 1:7], row.names = FALSE, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_dapa_rbp_knockdown.txt"))
