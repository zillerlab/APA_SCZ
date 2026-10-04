#!/usr/bin/env Rscript
## ===========================================================================
## 05_fig5h_rbp_variance_explained.R
##
## Variance in SCZ-associated APA across the iN cohort explained by the expression
## of PTBP2, PCBP2 and CPSF6.
##
## Manuscript: Fig. 5h. The panel combines the two runs: left (targets vs other genes),
## all tested genes (--genes all); right (unique contribution of each RBP), SCZ dAPA
## genes (--genes scz).
##
## USAGE
##   Rscript 03_differential_apa/05_fig5h_rbp_variance_explained.R [--genes scz|all] [--cores <n>] [--out <dir>]
##     --genes scz (default): SCZ vs HC dAPA genes; output prefix Fig5h (Fig. 5h, right)
##     --genes all: all genes tested in the cohort dAPA analysis; output prefix Fig5h_allgenes
##                  (Fig. 5h, left)
##
## INPUTS
##   results/dapa/apa_ratio_iN_cohort.txt             01_dapa_iN.R (distal usage ratio R, 93 iNs)
##   results/dapa/dAPA_iN_SCZ_vs_Ctrl.txt             01_dapa_iN.R (SCZ dAPA genes)
##   results/dapa/dAPA_RBP_KD_<RBP>.txt               03_dapa_rbp_knockdown.R (RBP targets)
##   resources/rnaseq/gene_counts_iN_cohort.txt.gz    gene counts, 95 iNs
##   resources/rnaseq/rnaseq_iN_cohort_metadata.txt   diagnosis, sex, site, technical PCs
##
## OUTPUTS (results/figures/)
##   Fig5h_variance_explained.pdf / .png
##   Fig5h_varExp_all.txt   per SCZ dAPA gene: R2 of diagnosis, of diagnosis + the three
##                          RBPs, partial R2 of the RBPs beyond diagnosis, R2 loss when
##                          each RBP is dropped (and its F-test p), target class
##                          (Source Data Fig. 5h)
##   Fig5h_summary.txt      tests
##
## METHOD
##   Genes: SCZ vs HC dAPA genes in iNs (q <= 0.05), or all genes tested in the cohort
##   dAPA analysis (--genes all).
##   APA: per gene, residuals (response scale) of the covariate-only beta regression
##   R ~ sex + TechnicalPC1 over the iNs with 0 < R < 1 (the model of the cohort dAPA
##   analysis without diagnosis).
##   RBP expression: DESeq2 variance-stabilising transformation (blind) of the gene
##   counts of all 95 iNs, with site, sex and DE_TechnicalPC1 removed
##   (limma::removeBatchEffect, diagnosis protected; the covariates of the DE model).
##   Per gene: base model lm(APA ~ diagnosis), full model
##   lm(APA ~ diagnosis + PTBP2 + PCBP2 + CPSF6).
##   Left: partial R2 of the three RBPs beyond diagnosis,
##   (R2_full - R2_base) / (1 - R2_base), for RBP targets vs other genes; RBP target =
##   dAPA upon knockdown of at least one of the three RBPs (BH <= 1e-4; the gene sets of
##   Fig. 5e). One-sided Wilcoxon rank-sum test (targets greater); knockdowns with fewer
##   than 5 targets among the SCZ dAPA genes are not tested.
##   Right: unique contribution of each RBP, R2_full - R2 of the full model without that
##   RBP; F-test of the nested models per gene; pairwise two-sided paired Wilcoxon tests
##   between RBPs.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(limma); library(DESeq2); library(betareg) })

CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
GENESET <- match.arg(get_arg("--genes", "scz"), c("scz", "all"))
PREFIX  <- if (GENESET == "scz") "Fig5h" else "Fig5h_allgenes"
GLAB    <- if (GENESET == "scz") "SCZ dAPA genes" else "tested genes"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RBPS <- c(PTBP2 = "ENSG00000117569", PCBP2 = "ENSG00000197111", CPSF6 = "ENSG00000111605")
P <- names(RBPS)
Q_COHORT <- 0.05; BH_KD <- 1e-4; MIN_LIB <- 20; MIN_TARGETS <- 5
DAPA <- file.path(REPO, "results", "dapa"); RS <- file.path(REPO, "resources", "rnaseq")

## ------------------------------------------------------------ inputs
R <- as.matrix(read.delim(need_file(file.path(DAPA, "apa_ratio_iN_cohort.txt"), "01_dapa_iN.R"),
                          check.names = FALSE, row.names = 1))
meta <- read.delim(need_file(file.path(RS, "rnaseq_iN_cohort_metadata.txt"), "resources"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
coh <- read.delim(need_file(file.path(DAPA, "dAPA_iN_SCZ_vs_Ctrl.txt"), "01_dapa_iN.R"), stringsAsFactors = FALSE)
kd  <- lapply(setNames(P, P), function(A)
  read.delim(need_file(file.path(DAPA, sprintf("dAPA_RBP_KD_%s.txt", A)), "03_dapa_rbp_knockdown.R"),
             stringsAsFactors = FALSE))
cm  <- read.delim(need_file(file.path(RS, "gene_counts_iN_cohort.txt.gz"), "resources"), check.names = FALSE,
                  stringsAsFactors = FALSE)
cm  <- cm[!duplicated(cm$gene_id), ]
cts <- as.matrix(cm[, meta$sample]); rownames(cts) <- cm$gene_id

genes <- intersect(coh$ID[!is.na(coh$qvalue) & (GENESET == "all" | coh$qvalue <= Q_COHORT)], rownames(R))
m93 <- meta[colnames(R), ]
dx  <- factor(m93$diagnosis)

## ------------------------------------------------------------ APA: covariate-only beta-regression residuals
resid_beta <- function(y) {
  k <- !is.na(y) & y > 0 & y < 1; out <- rep(NA_real_, length(y))
  if (sum(k) < MIN_LIB) return(out)
  d <- data.frame(ratio = y[k], sex = m93$sex[k], TechnicalPC1 = m93$TechnicalPC1[k])
  f <- tryCatch(betareg(ratio ~ sex + TechnicalPC1, data = d, link = "logit", type = "ML", link.phi = "identity"),
                error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(f)) out[k] <- d$ratio - fitted(f)
  out
}
A <- do.call(rbind, parallel::mclapply(genes, function(g) resid_beta(R[g, ]), mc.cores = CORES))
dimnames(A) <- list(genes, colnames(R))

## ------------------------------------------------------------ RBP expression
dds <- estimateSizeFactors(suppressWarnings(suppressMessages(
  DESeqDataSetFromMatrix(cts, meta[, c("diagnosis", "sex")], ~ diagnosis))))
V <- assay(vst(dds, blind = TRUE))
E <- removeBatchEffect(V, batch = meta$site, batch2 = meta$sex, covariates = meta$DE_TechnicalPC1,
                       design = model.matrix(~ diagnosis, meta))[RBPS, colnames(R)]
rownames(E) <- P

## ------------------------------------------------------------ per-gene models
fit_gene <- function(y) {
  ok <- !is.na(y); if (sum(ok) < MIN_LIB) return(NULL)
  d <- data.frame(y = y[ok], dx = droplevels(dx[ok]), t(E[, ok, drop = FALSE]))
  base <- lm(y ~ dx, d); full <- lm(y ~ dx + PTBP2 + PCBP2 + CPSF6, d)
  r0 <- summary(base)$r.squared; r1 <- summary(full)$r.squared
  out <- c(r2_diagnosis = r0, r2_full = r1, partial_r2_rbp = (r1 - r0) / (1 - r0),
           p_rbp = anova(base, full)$`Pr(>F)`[2])
  for (g in P) {
    red <- lm(reformulate(c("dx", setdiff(P, g)), "y"), d)
    out[paste0("r2_loss_", g)] <- r1 - summary(red)$r.squared
    out[paste0("p_", g)] <- anova(red, full)$`Pr(>F)`[2]
  }
  out
}
res <- lapply(genes, function(g) fit_gene(A[g, ]))
keep <- !vapply(res, is.null, logical(1))
res <- data.frame(gene = genes[keep], do.call(rbind, res[keep]), check.names = FALSE)
targets <- Reduce(union, lapply(kd, function(r) r$ID[r$padj_BH <= BH_KD]))
res$class <- factor(ifelse(res$gene %in% targets, "RBPtarget", "NoTarget"), levels = c("NoTarget", "RBPtarget"))
for (A_ in P) res[[paste0("KD_dAPA_", A_)]] <- res$gene %in% kd[[A_]]$ID[kd[[A_]]$padj_BH <= BH_KD]
write.table(res, file.path(OUT_DIR, paste0(PREFIX, "_varExp_all.txt")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ------------------------------------------------------------ tests
tg <- res$class == "RBPtarget"
wt <- wilcox.test(res$partial_r2_rbp[tg], res$partial_r2_rbp[!tg], alternative = "greater")
summ <- data.frame(test = "partial R2 beyond diagnosis, RBP targets > non-targets (one-sided Wilcoxon rank-sum)",
                   n_1 = sum(tg), n_2 = sum(!tg), median_1 = median(res$partial_r2_rbp[tg]),
                   median_2 = median(res$partial_r2_rbp[!tg]), p_value = wt$p.value)
for (A_ in P) {
  x <- res[[paste0("KD_dAPA_", A_)]]
  summ <- rbind(summ, data.frame(test = sprintf("partial R2 beyond diagnosis, %s knockdown targets > other genes (one-sided Wilcoxon rank-sum)", A_),
                                 n_1 = sum(x), n_2 = sum(!x), median_1 = if (any(x)) median(res$partial_r2_rbp[x]) else NA,
                                 median_2 = median(res$partial_r2_rbp[!x]),
                                 p_value = if (sum(x) >= MIN_TARGETS) wilcox.test(res$partial_r2_rbp[x], res$partial_r2_rbp[!x], alternative = "greater")$p.value else NA))
}
for (p in combn(P, 2, simplify = FALSE)) {
  x <- res[[paste0("r2_loss_", p[1])]]; y <- res[[paste0("r2_loss_", p[2])]]
  summ <- rbind(summ, data.frame(test = sprintf("R2 loss when dropped, %s vs %s (two-sided paired Wilcoxon)", p[1], p[2]),
                                 n_1 = length(x), n_2 = length(y), median_1 = median(x), median_2 = median(y),
                                 p_value = wilcox.test(x, y, paired = TRUE)$p.value))
}
for (A_ in P) summ <- rbind(summ, data.frame(test = sprintf("genes with a significant unique contribution of %s (nested F-test p < 0.05)", A_),
                                             n_1 = sum(res[[paste0("p_", A_)]] < 0.05), n_2 = nrow(res),
                                             median_1 = NA, median_2 = NA, p_value = NA))
write.table(summ, file.path(OUT_DIR, paste0(PREFIX, "_summary.txt")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ------------------------------------------------------------ figure
BLUE <- "#1f78c8"; GREY <- "#b9b8b2"
fmt_p <- function(p) if (p < 0.001) formatC(p, format = "e", digits = 1) else as.character(signif(p, 2))
vtheme <- theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text = element_text(colour = FIG_SEC, size = 9),
        axis.title = element_text(colour = FIG_SEC), plot.subtitle = element_text(hjust = 0.5, colour = FIG_SEC, size = 8.5))
nl <- table(res$class); ymax <- max(res$partial_r2_rbp)
p1 <- ggplot(res, aes(class, partial_r2_rbp, fill = class)) +
  geom_violin(colour = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.14, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = c(NoTarget = GREY, RBPtarget = BLUE)) +
  scale_x_discrete(labels = sprintf("%s\nn = %s", names(nl), format(nl, big.mark = ","))) +
  annotate("segment", x = 1, xend = 2, y = ymax * 1.04, yend = ymax * 1.04, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ymax * 1.09, size = 3, label = sprintf("p = %s", fmt_p(wt$p.value))) +
  labs(x = NULL, y = expression("Partial "*R^2*"("*Delta*"APA) by RBP expression"),
       subtitle = sprintf("%s %s, beyond diagnosis", format(nrow(res), big.mark = ","), GLAB)) + vtheme
dl <- data.frame(RBP = factor(rep(P, each = nrow(res)), levels = P),
                 loss = unlist(res[, paste0("r2_loss_", P)], use.names = FALSE))
p2 <- ggplot(dl, aes(RBP, loss)) +
  geom_violin(fill = BLUE, colour = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.14, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  labs(x = NULL, y = expression("Unique "*R^2*"("*Delta*"APA)"),
       subtitle = sprintf("%s %s, loss when dropped", format(nrow(res), big.mark = ","), GLAB)) + vtheme
lay <- grid::grid.layout(2, 2, heights = grid::unit(c(0.1, 0.9), "npc"), widths = grid::unit(c(0.5, 0.5), "npc"))
draw <- function() {
  grid::pushViewport(grid::viewport(layout = lay))
  grid::grid.text(expression(Delta*"APA variance explained by RBP expression"),
                  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = grid::gpar(fontsize = 12))
  print(p1, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid::popViewport()
}
pdf(file.path(OUT_DIR, paste0(PREFIX, "_variance_explained.pdf")), width = 7, height = 3.6); grid::grid.newpage(); draw(); invisible(dev.off())
png(file.path(OUT_DIR, paste0(PREFIX, "_variance_explained.png")), width = 7, height = 3.6, units = "in", res = 300)
grid::grid.newpage(); draw(); invisible(dev.off())

## ------------------------------------------------------------ summary
cat(sprintf("%d %s (%d RBP targets)\n", nrow(res), GLAB, sum(tg)))
print(summ, row.names = FALSE, digits = 3)
