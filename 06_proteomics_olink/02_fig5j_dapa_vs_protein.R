#!/usr/bin/env Rscript
## ===========================================================================
## 02_fig5j_dapa_vs_protein.R
##
## APA change vs protein change upon PTBP2 and PCBP2 knockdown.
##
## Manuscript: Fig. 5j (PTBP2); main text (PTBP2, PCBP2).
##
## USAGE
##   Rscript 06_proteomics_olink/02_fig5j_dapa_vs_protein.R [--out <dir>]
##
## INPUTS
##   results/dapa/dAPA_RBP_KD_<PTBP2|PCBP2|CPSF6>.txt   03_differential_apa/03_dapa_rbp_knockdown.R
##   results/proteomics/DEP_OLINK_RBP_KD.txt           01_dep_olink_rbp_knockdown.R
##
## OUTPUTS (results/figures/)
##   Fig5j_dAPA_vs_protein.pdf / .png     scatter, PTBP2 knockdown
##   Fig5j_dAPA_vs_protein.txt            genes plotted: dAPA beta, protein logFC (Source Data Fig. 5j)
##   Fig5j_dAPA_vs_protein_stats.txt      correlation per knockdown
##
## METHOD
##   Genes: on the Olink panel, tested in all three knockdown dAPA analyses, and dAPA upon
##   knockdown of the RBP (BH <= 1e-4). x = knockdown effect on distal usage (condition
##   coefficient of the beta regression, logit scale); y = protein logFC of the same knockdown.
##   Genes that are outliers on either axis (|MAD z| > 3 for the dAPA beta or for the
##   protein logFC) are excluded. Spearman correlation, two-sided p; knockdowns with fewer
##   than 10 genes are not tested. The line is a least-squares fit with its 95% confidence
##   band.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ARMS <- c("PTBP2", "PCBP2", "CPSF6"); KDS <- c("PTBP2", "PCBP2"); FIG_KD <- "PTBP2"
BH_KD <- 1e-4; MAX_Z <- 3; MIN_N <- 10

kd <- lapply(setNames(ARMS, ARMS), function(A) {
  r <- read.delim(need_file(file.path(REPO, "results", "dapa", sprintf("dAPA_RBP_KD_%s.txt", A)), "03_dapa_rbp_knockdown.R"),
                  stringsAsFactors = FALSE)
  r[is.finite(r$beta), ]
})
prot <- read.delim(need_file(file.path(REPO, "results", "proteomics", "DEP_OLINK_RBP_KD.txt"), "01_dep_olink_rbp_knockdown.R"),
                   stringsAsFactors = FALSE)
univ <- intersect(Reduce(intersect, lapply(kd, `[[`, "ID")), prot$protein)
madz <- function(x) { m <- median(x); d <- median(abs(x - m)); if (d == 0) 0 * x else 0.6745 * (x - m) / d }

pts <- list(); st <- list()
for (A in KDS) {
  r <- kd[[A]]; g <- intersect(r$ID[r$padj_BH <= BH_KD], univ)
  d <- data.frame(knockdown = A, gene = g, dAPA_beta = r$beta[match(g, r$ID)],
                  protein_logFC = prot[[paste0(A, "_logFC")]][match(g, prot$protein)])
  d <- d[is.finite(d$protein_logFC), ]
  n_sel <- nrow(d); d <- d[abs(madz(d$protein_logFC)) <= MAX_Z & abs(madz(d$dAPA_beta)) <= MAX_Z, ]
  c2 <- if (nrow(d) >= MIN_N) suppressWarnings(cor.test(d$dAPA_beta, d$protein_logFC, method = "spearman", exact = FALSE)) else
          list(estimate = NA_real_, p.value = NA_real_)
  st[[A]] <- data.frame(knockdown = A, n_dAPA_on_panel = n_sel, n = nrow(d), rho = unname(c2$estimate),
                        p_value = c2$p.value)
  pts[[A]] <- d
}
st <- do.call(rbind, st); pts <- do.call(rbind, pts)
write.table(st, file.path(OUT_DIR, "Fig5j_dAPA_vs_protein_stats.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pts, file.path(OUT_DIR, "Fig5j_dAPA_vs_protein.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

d <- pts[pts$knockdown == FIG_KD, ]; s <- st[st$knockdown == FIG_KD, ]
p <- ggplot(d, aes(dAPA_beta, protein_logFC)) +
  geom_hline(yintercept = 0, colour = FIG_GRID) + geom_vline(xintercept = 0, colour = FIG_GRID) +
  geom_smooth(method = "lm", formula = y ~ x, colour = FIG_INK, fill = "grey80", linewidth = 0.6) +
  geom_point(size = 1.4, alpha = 0.8, colour = "grey20") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 3, colour = FIG_SEC,
           label = sprintf("rho = %+.3f\np = %.3f\nn = %d", s$rho, s$p_value, s$n)) +
  labs(x = sprintf("dAPA beta (%s KD)", FIG_KD), y = sprintf("log2FC protein (%s KD)", FIG_KD),
       title = sprintf("%s KD: dAPA vs protein", FIG_KD)) +
  theme_classic(base_size = 10) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(OUT_DIR, "Fig5j_dAPA_vs_protein.pdf"), p, width = 3.4, height = 3.2)
ggsave(file.path(OUT_DIR, "Fig5j_dAPA_vs_protein.png"), p, width = 3.4, height = 3.2, dpi = 300)

cat(sprintf("genes on the panel tested in all three knockdowns: %d\n", length(univ)))
print(st, row.names = FALSE, digits = 3)
