#!/usr/bin/env Rscript
## ===========================================================================
## 02_gria1_3utr_ko_qpcr.R
##
## GRIA1 mRNA levels (qPCR) in iNs with a deletion of the GRIA1 long 3'UTR vs the unedited
## parental line, day 49.
##
## Manuscript: Fig. 7c.
##
## USAGE
##   Rscript 11_gria1_reporter/02_gria1_3utr_ko_qpcr.R [--out <dir>]
##
## INPUT (resources/gria1_reporter/)
##   gria1_3utr_ko_qpcr.txt   one row per qPCR reaction: batch (differentiation batch), line
##                            (PluriCore ID), genotype (WT / KO), sample, tech_rep (technical
##                            replicate), rel_expr (GRIA1 relative to EID2, 2^-dCt)
##
## OUTPUTS (results/figures/)
##   Fig7c_qPCR_GRIA1.pdf / .png
##   Fig7c_qPCR_GRIA1.txt         per reaction: values, sample mean and batch-normalised
##                                sample mean (Source Data Fig. 7c)
##   Fig7c_qPCR_GRIA1_test.txt    linear mixed model result
##
## METHOD
##   Technical replicates are averaged per sample on the log2 scale. Linear mixed model
##   log2(expression) ~ genotype + (1 | batch) on the sample means (lme4), two-sided t-test of
##   the genotype coefficient with Satterthwaite degrees of freedom (lmerTest).
##   Figure: batch-normalised sample means (display only). Per batch, a normalisation factor
##   is the mean of the WT and KO means on the log2 scale; it is subtracted and the mean factor
##   over batches added back, so that all batches are on a common scale.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(need_file(file.path(REPO, "resources", "gria1_reporter", "gria1_3utr_ko_qpcr.txt"), "resources"),
                stringsAsFactors = FALSE)
d$genotype <- factor(d$genotype, levels = c("WT", "KO"))

## ------------------------------------------------------------ sample means (technical replicates)
s <- aggregate(list(log2_expr = log2(d$rel_expr)), d[, c("batch", "genotype", "sample")], mean)
gm <- tapply(s$log2_expr, list(s$batch, s$genotype), mean)       # batch x genotype means
stopifnot(!anyNA(gm))
bf <- rowMeans(gm)                                                # batch normalisation factor (log2)
s$batch_norm <- 2^(s$log2_expr - bf[s$batch] + mean(bf))

## ------------------------------------------------------------ test
m  <- suppressMessages(lmer(log2_expr ~ genotype + (1 | batch), data = s))
co <- summary(m)$coefficients["genotypeKO", ]
test <- data.frame(model = "log2(expression) ~ genotype + (1 | batch), sample means",
                   n_WT = sum(s$genotype == "WT"), n_KO = sum(s$genotype == "KO"), n_batches = length(unique(s$batch)),
                   log2FC_KO_vs_WT = co[["Estimate"]], se = co[["Std. Error"]],
                   df = co[["df"]], t = co[["t value"]], p_value = co[["Pr(>|t|)"]])

i <- match(d$sample, s$sample)
out <- data.frame(d[, c("batch", "line", "genotype", "sample", "tech_rep", "rel_expr")],
                  sample_mean = 2^s$log2_expr[i], sample_mean_batchNorm = s$batch_norm[i])
write.table(out, file.path(OUT_DIR, "Fig7c_qPCR_GRIA1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(test, file.path(OUT_DIR, "Fig7c_qPCR_GRIA1_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ figure
ymax <- max(s$batch_norm)
p <- ggplot(s, aes(genotype, batch_norm)) +
  geom_boxplot(aes(fill = genotype), width = 0.55, outlier.shape = NA, linewidth = 0.4, colour = FIG_INK) +
  geom_point(size = 1.6, colour = FIG_INK) +
  annotate("segment", x = 1, xend = 2, y = ymax * 1.06, yend = ymax * 1.06, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ymax * 1.12, label = sprintf("p = %.3f", test$p_value), size = 3, colour = FIG_SEC) +
  scale_fill_manual(values = c(WT = "grey80", KO = "#1f78c8"), guide = "none") +
  scale_y_continuous(limits = c(0, ymax * 1.18)) +
  labs(x = NULL, y = "GRIA1 mRNA vs EID2\n(batch-normalised)") +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(OUT_DIR, "Fig7c_qPCR_GRIA1.pdf"), p, width = 2, height = 3)
ggsave(file.path(OUT_DIR, "Fig7c_qPCR_GRIA1.png"), p, width = 2, height = 3, dpi = 300)
print(test, row.names = FALSE, digits = 3)
