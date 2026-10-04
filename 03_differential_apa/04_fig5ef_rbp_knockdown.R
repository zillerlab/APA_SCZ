#!/usr/bin/env Rscript
## ===========================================================================
## 04_fig5ef_rbp_knockdown.R
##
## Fig. 5e: overlap of the knockdown dAPA gene sets (PTBP2, PCBP2, CPSF6).
## Fig. 5f: correlation of the knockdown APA effect with the SCZ vs HC APA effect
##          in the iN cohort.
## Overlap of each knockdown dAPA gene set with the SCZ vs HC dAPA genes (Venn diagrams).
##
## USAGE
##   Rscript 03_differential_apa/04_fig5ef_rbp_knockdown.R [--dapa <dir>] [--deg <dir>] [--out <dir>]
##
## INPUTS
##   results/dapa/dAPA_RBP_KD_<PTBP2|PCBP2|CPSF6>.txt   03_dapa_rbp_knockdown.R
##   results/dapa/dAPA_iN_SCZ_vs_Ctrl.txt               01_dapa_iN.R
##   results/deg/DEG_iN_SCZ_vs_Ctrl.txt                 05_differential_expression/01_deseq_iN_cohort.R
##                                                      (direction of each RBP in SCZ)
##
## OUTPUTS (results/figures/)
##   Fig5e_dAPA_overlap.pdf / .png          Venn diagram, knockdown dAPA genes (BH <= 1e-4)
##   Fig5e_dAPA_membership.txt              gene x knockdown membership (Source Data Fig. 5e)
##   Fig5f_KD_vs_cohort.pdf / .png          forest plot
##   Fig5f_KD_vs_cohort.txt                 correlation statistics (Source Data Fig. 5f)
##   Fig5f_KD_vs_cohort_scatter.pdf / .png / .txt  the same relation as scatter plots, one per knockdown
##   Fig5_dAPA_KD_vs_cohort_overlap.pdf / .png   Venn diagrams, knockdown dAPA vs SCZ dAPA genes
##   Fig5_dAPA_KD_vs_cohort_overlap.txt          counts and Fisher's exact tests
##   Fig5_dAPA_KD_vs_cohort_membership.txt       gene x set membership (source data)
##
## METHOD
##   Fig. 5e: universe = genes tested in at least one knockdown; knockdown dAPA =
##   Benjamini-Hochberg-adjusted Wald P <= 1e-4.
##   Fig. 5f: both effects are positive for higher distal 3'UTR usage.
##     knockdown:  KD - Ctrl        = deltaR of dAPA_RBP_KD_<RBP>.txt (difference in
##                                    fitted R, knockdown - control)
##     cohort:     SCZ - HC         = -delta of dAPA_iN_SCZ_vs_Ctrl.txt
##                                    (delta is the difference in R, HC - SCZ)
##   rho > 0: the knockdown shifts 3'UTR usage in the same direction as SCZ (mimics);
##   rho < 0: in the opposite direction (opposes).
##   Expected sign: a factor higher in SCZ (log2FC HC/SCZ < 0 in the iN DE analysis)
##   should oppose, a factor lower in SCZ should mimic.
##   Spearman rho over the cohort dAPA genes (q <= 0.05) tested in the knockdown.
##   Two-sided p (asymptotic t approximation); 95% CI from the Fisher z transform with
##   the Bonett-Wright SE, sqrt((1 + rho^2 / 2) / (n - 3)); percentile bootstrap CI
##   (2,000 resamples) in addition.
##   Overlap: per knockdown, genes tested in both the knockdown and the cohort; knockdown
##   dAPA (BH <= 1e-4) vs SCZ dAPA (q <= 0.05); two-sided Fisher's exact test.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))    # plot_venn3()

DAPA_DIR <- get_arg("--dapa", file.path(REPO, "results", "dapa"))
DEG_DIR  <- get_arg("--deg",  file.path(REPO, "results", "deg"))
OUT_DIR  <- get_arg("--out",  file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RBPS <- c("PTBP2", "PCBP2", "CPSF6")
Q <- 0.05; BH_KD <- 1e-4; N_BOOT <- 2000; SEED <- 1; MIN_N <- 10   # cohort dAPA q <= 0.05, knockdown dAPA BH <= 1e-4

INK <- "#0b0b0b"; SEC <- "#52514e"; MUT <- "#8a8a86"; GRID <- "#e4e3df"
NEG <- "#2a78d6"; POS <- "#d1344e"; NEU <- "#b9b8b2"

## ------------------------------------------------------------ inputs
kd <- lapply(setNames(RBPS, RBPS), function(A)
  read.delim(need_file(file.path(DAPA_DIR, sprintf("dAPA_RBP_KD_%s.txt", A)), "03_dapa_rbp_knockdown.R"),
             stringsAsFactors = FALSE))
coh <- read.delim(need_file(file.path(DAPA_DIR, "dAPA_iN_SCZ_vs_Ctrl.txt"), "01_dapa_iN.R"), stringsAsFactors = FALSE)
deg <- read.delim(need_file(file.path(DEG_DIR, "DEG_iN_SCZ_vs_Ctrl.txt"), "05_differential_expression"),
                  stringsAsFactors = FALSE)
coh$d_scz <- -coh$delta                                  # SCZ - HC

## ------------------------------------------------------------ Fig. 5e
univ <- sort(Reduce(union, lapply(kd, `[[`, "ID")))
mem  <- data.frame(gene = univ, sapply(kd, function(r) as.integer(univ %in% r$ID[r$padj_BH <= BH_KD])))
mem$n_RBPs <- rowSums(mem[, RBPS])
write.table(mem, file.path(OUT_DIR, "Fig5e_dAPA_membership.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

reg <- venn3_counts(mem[, RBPS])
pe  <- plot_venn3(mem[, RBPS], expression(paste(Delta, "APA (RBP-KD vs. Ctrl)")))
ggsave(file.path(OUT_DIR, "Fig5e_dAPA_overlap.pdf"), pe, width = 3.4, height = 3.4)
ggsave(file.path(OUT_DIR, "Fig5e_dAPA_overlap.png"), pe, width = 3.4, height = 3.4, dpi = 300)

## ------------------------------------------------------------ Fig. 5f statistics
dir_scz <- setNames(deg$log2FoldChange[match(RBPS, deg$gene_name)], RBPS)   # log2(HC / SCZ)
stopifnot(!anyNA(dir_scz))
cor_row <- function(A, subset, ids) {
  r <- kd[[A]]; ids <- intersect(ids, intersect(r$ID, coh$ID))
  x <- coh$d_scz[match(ids, coh$ID)]; y <- r$deltaR[match(ids, r$ID)]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < MIN_N) return(NULL)
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate); se <- sqrt((1 + rho^2 / 2) / (n - 3))
  set.seed(SEED)
  bs <- replicate(N_BOOT, { i <- sample.int(n, n, replace = TRUE); cor(x[i], y[i], method = "spearman") })
  expect <- if (dir_scz[[A]] < 0) "opposes" else "mimics"
  data.frame(arm = A, RBP_in_SCZ = if (dir_scz[[A]] < 0) "up" else "down",
             log2FC_HC_vs_SCZ = round(dir_scz[[A]], 3), expectation = expect, subset = subset, n = n,
             rho = rho, ci_low = tanh(atanh(rho) - 1.96 * se), ci_high = tanh(atanh(rho) + 1.96 * se),
             boot_low = unname(quantile(bs, 0.025)), boot_high = unname(quantile(bs, 0.975)),
             p_value = ct$p.value,
             pct_expected_sign = 100 * mean(if (expect == "opposes") x * y < 0 else x * y > 0))
}
cohort_sig <- coh$ID[!is.na(coh$qvalue) & coh$qvalue <= Q]
st <- do.call(rbind, lapply(RBPS, function(A) cor_row(A, "cohort q<=0.05", cohort_sig)))
st$row_order <- rev(seq_len(nrow(st)))                  # 1 = bottom row
write.table(st, file.path(OUT_DIR, "Fig5f_KD_vs_cohort.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ Fig. 5f forest plot
d <- st[order(st$row_order), ]
d$y   <- d$row_order
d$sig <- d$p_value < 0.05
d$col <- ifelse(!d$sig, NEU, ifelse(d$rho < 0, NEG, POS))
fmt_p <- function(p) ifelse(p >= 0.05, "n.s.", sub("e-0", "e-", sprintf("p=%.1e", p)))
d$lab_left  <- sprintf("%s | n=%s", d$subset, formatC(d$n, width = 5, big.mark = ","))
d$lab_right <- sprintf("%+.3f  %s", d$rho, fmt_p(d$p_value))
armY <- aggregate(y ~ arm, d, mean); armY$exp <- d$expectation[match(armY$arm, d$arm)]
armY$dir <- d$RBP_in_SCZ[match(armY$arm, d$arm)]
seps <- head(cumsum(rle(d$arm)$lengths), -1) + 0.5
XL <- -0.5; XR <- 0.3; YT <- nrow(d) + 1.1

pf <- ggplot(d, aes(rho, y)) +
  geom_hline(yintercept = seps, colour = GRID, linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.4, linetype = "dashed") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, colour = I(col)), orientation = "y", width = 0.3, linewidth = 0.8) +
  geom_point(aes(fill = I(col)), shape = 21, colour = "white", size = 2.6, stroke = 0.6) +
  geom_text(aes(x = XL - 0.01, label = lab_left), hjust = 1, size = 2.6, colour = SEC) +
  geom_text(aes(x = XR + 0.02, label = lab_right, colour = I(ifelse(sig, INK, MUT))), hjust = 0, size = 2.6) +
  geom_text(data = armY, aes(x = XL - 0.37, y = y + 0.14, label = arm), hjust = 1, size = 3.6, colour = INK) +
  geom_text(data = armY, aes(x = XL - 0.37, y = y - 0.2, label = sprintf("%s in SCZ", dir)),
            hjust = 1, size = 2.3, colour = MUT) +
  annotate("segment", x = -0.03, xend = -0.25, y = YT, yend = YT, colour = NEG, linewidth = 0.8,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  annotate("segment", x = 0.03, xend = 0.25, y = YT, yend = YT, colour = POS, linewidth = 0.8,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  annotate("text", x = -0.14, y = YT + 0.75, label = "KD OPPOSES\nSCZ shift", colour = NEG, size = 3, lineheight = 0.9) +
  annotate("text", x = 0.14, y = YT + 0.75, label = "KD MIMICS\nSCZ shift", colour = POS, size = 3, lineheight = 0.9) +
  scale_x_continuous(breaks = seq(-0.5, 0.3, 0.1), labels = function(x) sprintf("%.1f", x)) +
  coord_cartesian(xlim = c(XL, XR), ylim = c(0.5, nrow(d) + 2), clip = "off") +
  labs(x = expression(paste("Spearman ", rho, "  ", Delta, "R (SCZ - HC) vs. ", Delta, "R (KD - Ctrl)  [95% CI]")),
       y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.line.x = element_line(colour = MUT, linewidth = 0.4), axis.ticks.x = element_line(colour = MUT),
        axis.text.x = element_text(colour = SEC), axis.title.x = element_text(colour = SEC, size = 8.5),
        plot.margin = margin(10, 95, 8, 175))
hf <- 1.6 + 0.4 * nrow(d)
ggsave(file.path(OUT_DIR, "Fig5f_KD_vs_cohort.pdf"), pf, width = 7, height = hf)
ggsave(file.path(OUT_DIR, "Fig5f_KD_vs_cohort.png"), pf, width = 7, height = hf, dpi = 300, bg = "white")

## ------------------------------------------------------------ Fig. 5f scatter plots (one per knockdown)
sc <- do.call(rbind, lapply(RBPS, function(A) {
  r <- kd[[A]]; ids <- intersect(cohort_sig, intersect(r$ID, coh$ID))
  data.frame(arm = A, gene = ids, dR_SCZ_minus_HC = coh$d_scz[match(ids, coh$ID)], dR_KD_minus_Ctrl = r$deltaR[match(ids, r$ID)])
}))
sc <- sc[is.finite(sc$dR_SCZ_minus_HC) & is.finite(sc$dR_KD_minus_Ctrl), ]
write.table(sc, file.path(OUT_DIR, "Fig5f_KD_vs_cohort_scatter.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
sc$arm <- factor(sc$arm, RBPS)
lab <- setNames(sprintf("%s\nrho = %+.2f, p = %s, n = %d", st$arm, st$rho,
                        ifelse(st$p_value < 0.001, sub("e-0", "e-", sprintf("%.1e", st$p_value)), sprintf("%.2f", st$p_value)), st$n), st$arm)
ps <- ggplot(sc, aes(dR_SCZ_minus_HC, dR_KD_minus_Ctrl)) +
  geom_hline(yintercept = 0, colour = GRID) + geom_vline(xintercept = 0, colour = GRID) +
  geom_point(size = 0.9, alpha = 0.6, colour = "grey35") +
  geom_smooth(method = "lm", formula = y ~ x, colour = INK, fill = "grey80", linewidth = 0.5) +
  facet_wrap(~ arm, nrow = 1, scales = "free_y", labeller = as_labeller(lab)) +
  labs(x = expression(paste(Delta, "R (SCZ - HC), SCZ dAPA genes")), y = expression(paste(Delta, "R (KD - Ctrl)"))) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 8.5, lineheight = 1.1),
        axis.text = element_text(colour = SEC), axis.title = element_text(colour = SEC))
ggsave(file.path(OUT_DIR, "Fig5f_KD_vs_cohort_scatter.pdf"), ps, width = 7.5, height = 2.8)
ggsave(file.path(OUT_DIR, "Fig5f_KD_vs_cohort_scatter.png"), ps, width = 7.5, height = 2.8, dpi = 300, bg = "white")

## ------------------------------------------------------------ knockdown dAPA vs cohort dAPA
ovl <- list(); memc <- list(); pv <- list()
for (A in RBPS) {
  r <- kd[[A]]; u <- intersect(r$ID, coh$ID[!is.na(coh$qvalue)])
  a <- u %in% r$ID[r$padj_BH <= BH_KD]; b <- u %in% cohort_sig
  ft <- fisher.test(table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))))
  ovl[[A]] <- data.frame(knockdown = A, n_tested_in_both = length(u), n_KD_dAPA = sum(a), n_SCZ_dAPA = sum(b),
                         n_overlap = sum(a & b), n_expected = sum(a) * sum(b) / length(u),
                         odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
                         p_value = ft$p.value)
  memc[[A]] <- data.frame(knockdown = A, gene = u, KD_dAPA = as.integer(a), SCZ_dAPA = as.integer(b))
  pv[[A]] <- plot_venn2(a, b, c(sprintf("%s KD\ndAPA", A), "SCZ vs HC\ndAPA"), title = A,
                        subtitle = sprintf("OR = %.2f, p = %s\n%s genes tested in both", unname(ft$estimate),
                                           signif(ft$p.value, 2), format(length(u), big.mark = ",")))
}
ovl <- do.call(rbind, ovl); memc <- do.call(rbind, memc)
write.table(ovl, file.path(OUT_DIR, "Fig5_dAPA_KD_vs_cohort_overlap.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(memc, file.path(OUT_DIR, "Fig5_dAPA_KD_vs_cohort_membership.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
draw_ov <- function() {
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(RBPS))))
  for (i in seq_along(RBPS)) print(pv[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
  grid::popViewport()
}
pdf(file.path(OUT_DIR, "Fig5_dAPA_KD_vs_cohort_overlap.pdf"), width = 8.4, height = 2.9); grid::grid.newpage(); draw_ov(); invisible(dev.off())
png(file.path(OUT_DIR, "Fig5_dAPA_KD_vs_cohort_overlap.png"), width = 8.4, height = 2.9, units = "in", res = 300)
grid::grid.newpage(); draw_ov(); invisible(dev.off())

## ------------------------------------------------------------ summary
cat(sprintf("Fig. 5e: %d genes tested in >= 1 knockdown; dAPA %s\n", length(univ),
            paste(RBPS, colSums(mem[, RBPS]), collapse = ", ")))
cat("Venn regions:", paste(names(reg), reg, collapse = ", "), "\n")
cat("Knockdown dAPA vs SCZ dAPA:\n"); print(ovl, row.names = FALSE, digits = 3)
cat("Fig. 5f:\n"); print(st[, c("arm", "expectation", "subset", "n", "rho", "ci_low", "ci_high", "p_value")],
                        row.names = FALSE, digits = 3)
