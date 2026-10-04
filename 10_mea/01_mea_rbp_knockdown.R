#!/usr/bin/env Rscript
## ===========================================================================
## 01_mea_rbp_knockdown.R
##
## Neuronal network activity (multielectrode array, MEA) of iNs after shRNA knockdown
## of CPSF6, PCBP2 or PTBP2, compared with the scramble shRNA control.
##
## Manuscript: Fig. 6f,g.
##
## USAGE
##   Rscript 10_mea/01_mea_rbp_knockdown.R [--out <dir>]
##
## INPUT (resources/mea/)
##   mea_rbp_knockdown_wells.txt   one row per well: well, plate, batch (experimental
##                                 batch, B1-B5), line (PluriCore ID), condition (CTRL = scramble
##                                 shRNA, CPSF6, PCBP2, PTBP2), then the well-averaged MEA
##                                 parameters (Axion Navigator)
##
## OUTPUTS
##   results/mea/MEA_RBP_KD_LMM.txt       per parameter and knockdown: estimate (knockdown -
##                                        control, analysis scale), SE, df, t, p, q (BH over the
##                                        parameters, per knockdown), wells used
##   results/figures/Fig6f_MEA_heatmap.pdf / .png / .txt
##   results/figures/Fig6g_MEA_features.pdf / .png / .txt
##
## METHOD
##   Wells with at least one spike and >= 3 active electrodes. Per parameter: spike, burst
##   and network-burst counts log2-transformed; values outside 1.5 x IQR (over all wells)
##   removed; linear mixed model value ~ condition + (1 | line) + (1 | batch) (lme4, REML),
##   reference = CTRL, Satterthwaite t-tests (lmerTest). Benjamini-Hochberg over the 27
##   parameters, separately per knockdown.
##   Fig. 6f: -log10 p per parameter and knockdown; stars: q <= 0.05 (*), 0.01 (**),
##   0.001 (***). Fig. 6g: selected parameters, values divided by the mean of their
##   batch relative to the first batch (display only; the model uses the values
##   before this scaling).
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "mea"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
KDS <- c("CPSF6", "PCBP2", "PTBP2")
LOG2 <- c("Number.of.Spikes", "Number.of.Bursts", "Number.of.Network.Bursts")
EXCLUDE <- c("Number.of.Active.Electrodes", "Number.of.Elecs.Participating.in.Burst...Avg")   # not tested
FIG_G <- c("Burst.Percentage...Avg", "Network.Burst.Frequency")

## ------------------------------------------------------------ input
w <- read.delim(need_file(file.path(REPO, "resources", "mea", "mea_rbp_knockdown_wells.txt"), "resources"),
                stringsAsFactors = FALSE)
w <- w[!is.na(w$Number.of.Spikes) & w$Number.of.Spikes > 0 & w$Number.of.Active.Electrodes >= 3, ]
w$condition <- relevel(factor(w$condition), ref = "CTRL")
w$line <- factor(w$line); w$batch <- factor(w$batch)
pars <- names(w)[(which(names(w) == "condition") + 1):ncol(w)]
pars <- setdiff(pars[!grepl("\\.Std", pars)], EXCLUDE)
cat(sprintf("%d wells, %d lines, %d batches; %d parameters\n", nrow(w), nlevels(w$line), nlevels(w$batch), length(pars)))
print(table(w$condition))

prep <- function(p) {
  v <- w[[p]]; if (p %in% LOG2) v <- log2(v)
  lo <- quantile(v, 0.25, na.rm = TRUE) - 1.5 * IQR(v, na.rm = TRUE)
  hi <- quantile(v, 0.75, na.rm = TRUE) + 1.5 * IQR(v, na.rm = TRUE)
  v[!(is.finite(v) & v >= lo & v <= hi)] <- NA
  d <- data.frame(w[, c("line", "batch", "condition")], value = v)
  d[!is.na(d$value), ]
}

## ------------------------------------------------------------ mixed models
res <- do.call(rbind, lapply(pars, function(p) {
  d  <- prep(p)
  co <- summary(suppressMessages(lmer(value ~ condition + (1 | line) + (1 | batch), d)))$coefficients
  data.frame(parameter = p, knockdown = KDS, scale = ifelse(p %in% LOG2, "log2", "linear"),
             estimate = co[paste0("condition", KDS), "Estimate"], se = co[paste0("condition", KDS), "Std. Error"],
             df = co[paste0("condition", KDS), "df"], t = co[paste0("condition", KDS), "t value"],
             p_value = co[paste0("condition", KDS), "Pr(>|t|)"], n_wells = nrow(d))
}))
res$q_value <- ave(res$p_value, res$knockdown, FUN = function(x) p.adjust(x, "BH"))
write.table(res, file.path(OUT_DIR, "MEA_RBP_KD_LMM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ Fig. 6f
lab <- function(x) trimws(gsub("\\s+", " ", gsub("\\.", " ", sub("\\.\\.\\.Avg.*$|\\.\\.Hz\\.$|\\.\\.sec\\.$", "", x))))
res$stars <- ifelse(res$q_value <= 0.001, "***", ifelse(res$q_value <= 0.01, "**", ifelse(res$q_value <= 0.05, "*", "")))
fd <- res; fd$par <- factor(lab(fd$parameter), levels = rev(unique(lab(pars))))
write.table(res[, c("parameter", "knockdown", "estimate", "p_value", "q_value", "stars")],
            file.path(FIG_DIR, "Fig6f_MEA_heatmap.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(fd, aes(knockdown, par, fill = -log10(p_value))) +
  geom_tile(colour = "white") + geom_text(aes(label = stars), size = 3, colour = FIG_INK) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "-log10 p") +
  labs(x = NULL, y = NULL, caption = "q (BH over parameters): * <= 0.05, ** <= 0.01, *** <= 0.001") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = FIG_INK),
        plot.caption = element_text(size = 6.5, colour = FIG_SEC, hjust = 0), plot.caption.position = "plot")
ggsave(file.path(FIG_DIR, "Fig6f_MEA_heatmap.pdf"), p, width = 4.8, height = 6.2)
ggsave(file.path(FIG_DIR, "Fig6f_MEA_heatmap.png"), p, width = 4.8, height = 6.2, dpi = 300)

## ------------------------------------------------------------ Fig. 6g
g <- do.call(rbind, lapply(FIG_G, function(p) {
  d <- prep(p); bf <- tapply(d$value, d$batch, mean); bf <- bf / bf[as.character(d$batch[1])]
  data.frame(parameter = lab(p), d[, c("line", "batch", "condition")], value = d$value,
             value_batch_scaled = d$value / bf[as.character(d$batch)])
}))
write.table(g, file.path(FIG_DIR, "Fig6g_MEA_features.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
ql <- res[res$parameter %in% FIG_G, ]; ql$parameter <- lab(ql$parameter)
ql$label <- sprintf("q = %.2g", ql$q_value); ql$condition <- factor(ql$knockdown, levels = levels(w$condition))
ymax <- tapply(g$value_batch_scaled, g$parameter, max); ql$y <- ymax[ql$parameter] * 1.08
p <- ggplot(g, aes(condition, value_batch_scaled)) +
  geom_violin(fill = "grey90", colour = FIG_MUT, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3) +
  geom_text(data = ql, aes(condition, y, label = label), size = 2.6, colour = FIG_SEC) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Value (scaled to batch mean)") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(FIG_DIR, "Fig6g_MEA_features.pdf"), p, width = 5.6, height = 2.8)
ggsave(file.path(FIG_DIR, "Fig6g_MEA_features.png"), p, width = 5.6, height = 2.8, dpi = 300)

cat("parameters with q <= 0.05 (decreased / increased vs control):\n")
for (k in KDS) { r <- res[res$knockdown == k & res$q_value <= 0.05, ]
  cat(sprintf("  %s: %d (%d / %d)\n", k, nrow(r), sum(r$estimate < 0), sum(r$estimate > 0))) }
print(res[res$parameter %in% FIG_G, c("parameter", "knockdown", "estimate", "p_value", "q_value")], row.names = FALSE, digits = 3)
