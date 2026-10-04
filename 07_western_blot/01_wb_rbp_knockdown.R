#!/usr/bin/env Rscript
## ===========================================================================
## 01_wb_rbp_knockdown.R
##
## Western blot quantification of the RBP knockdowns in iNs.
##
## Manuscript: Fig. 5c.
##
## USAGE
##   Rscript 07_western_blot/01_wb_rbp_knockdown.R [--out <dir>]
##
## INPUT (resources/western_blot/)
##   wb_rbp_knockdown.txt   protein (antibody), donor (PluriCore ID), batch, condition
##                          (CTRL = scramble shRNA, KD = shRNA against the protein), NRI
##                          (target peak area / neurofilament L peak area of the same sample)
##
## OUTPUTS (results/figures/)
##   Fig5c_WB_knockdown.pdf / .png
##   Fig5c_WB_knockdown.txt         per band: NRI and batch-normalised NRI (Source Data Fig. 5c)
##   Fig5c_WB_knockdown_tests.txt   paired test per protein
##
## METHOD
##   Each knockdown lane is paired with the scramble lane of the same donor and batch.
##   Test: ratio-paired t-test, i.e. one-sample t-test of log2(NRI_KD / NRI_CTRL) against 0,
##   one-sided (knockdown lower than control); the two-sided p-value is reported alongside.
##   Figure: batch-normalised NRI, i.e. NRI divided by the geometric mean of its pair (display
##   only); lines connect the lanes of a pair.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ORDER <- c("PTBP2", "CPSF6", "PCBP2")

d <- read.delim(need_file(file.path(REPO, "resources", "western_blot", "wb_rbp_knockdown.txt"), "resources"),
                stringsAsFactors = FALSE)
d$pair <- paste(d$protein, d$donor, d$batch, sep = "_")
stopifnot(all(table(d$pair) == 2), all(tapply(d$condition, d$pair, function(x) setequal(x, c("CTRL", "KD")))))
d$NRI_batchNorm <- 2^(log2(d$NRI) - ave(log2(d$NRI), d$pair))
prots <- intersect(ORDER, unique(d$protein))

tests <- do.call(rbind, lapply(prots, function(A) {
  s <- d[d$protein == A, ]
  kd <- s$NRI[s$condition == "KD"][order(s$pair[s$condition == "KD"])]
  ct <- s$NRI[s$condition == "CTRL"][order(s$pair[s$condition == "CTRL"])]
  lr <- log2(kd / ct); tt <- t.test(lr); t1 <- t.test(lr, alternative = "less")
  data.frame(protein = A, n_pairs = length(lr), n_donors = length(unique(s$donor)),
             mean_ratio_KD_CTRL = 2^mean(lr), ratio_lo95 = 2^tt$conf.int[1], ratio_hi95 = 2^tt$conf.int[2],
             t = unname(tt$statistic), p_value = t1$p.value, p_two_sided = tt$p.value)
}))
write.table(d[, c("protein", "donor", "batch", "condition", "NRI", "NRI_batchNorm")],
            file.path(OUT_DIR, "Fig5c_WB_knockdown.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tests, file.path(OUT_DIR, "Fig5c_WB_knockdown_tests.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

d$protein <- factor(d$protein, levels = prots)
d$condition <- factor(d$condition, levels = c("CTRL", "KD"))
lab <- data.frame(protein = factor(prots, levels = prots), x = 1.5, y = max(d$NRI_batchNorm) * 1.12,
                  label = sprintf("p = %s", signif(tests$p_value[match(prots, tests$protein)], 2)))
p <- ggplot(d, aes(condition, NRI_batchNorm)) +
  geom_boxplot(aes(fill = condition), width = 0.55, outlier.shape = NA, linewidth = 0.4, colour = FIG_INK) +
  geom_line(aes(group = pair), colour = FIG_MUT, linewidth = 0.3) +
  geom_point(size = 1.6, colour = FIG_INK) +
  geom_text(data = lab, aes(x, y, label = label), size = 3, colour = FIG_SEC) +
  annotate("segment", x = 1, xend = 2, y = max(d$NRI_batchNorm) * 1.06, yend = max(d$NRI_batchNorm) * 1.06, linewidth = 0.3) +
  facet_wrap(~ protein, nrow = 1) +
  scale_fill_manual(values = c(CTRL = "grey80", KD = "#1f78c8"), guide = "none") +
  scale_x_discrete(labels = function(x) ifelse(x == "KD", "shRNA", "Ctrl")) +
  scale_y_continuous(limits = c(0, max(d$NRI_batchNorm) * 1.18)) +
  labs(x = NULL, y = "Protein / NF-L [WB], batch-normalised") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(OUT_DIR, "Fig5c_WB_knockdown.pdf"), p, width = 1.9 * length(prots) + 0.6, height = 3)
ggsave(file.path(OUT_DIR, "Fig5c_WB_knockdown.png"), p, width = 1.9 * length(prots) + 0.6, height = 3, dpi = 300)
print(tests, row.names = FALSE, digits = 3)
