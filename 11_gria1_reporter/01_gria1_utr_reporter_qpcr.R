#!/usr/bin/env Rscript
## ===========================================================================
## 01_gria1_utr_reporter_qpcr.R
##
## GFP reporter mRNA levels with the GRIA1 long vs short 3'UTR in HEK293T cells (qPCR of
## sorted GFP-positive cells).
##
## Manuscript: Fig. 7h.
##
## USAGE
##   Rscript 11_gria1_reporter/01_gria1_utr_reporter_qpcr.R [--out <dir>]
##
## INPUT (resources/gria1_reporter/)
##   gria1_utr_reporter_qpcr.txt   one row per transfection: batch, experiment (a long and a
##                                 short 3'UTR transfection performed in parallel), label,
##                                 utr (Long / Short), plasmid, gfp_sort, rel_expr_adj (GFP
##                                 relative to RTF2, adjusted for transfection efficiency),
##                                 rel_expr_sd (propagated SD of technical replicates)
##
## OUTPUTS (results/figures/)
##   Fig7h_GRIA1_reporter.pdf / .png
##   Fig7h_GRIA1_reporter.txt         per transfection: values and batch-normalised values
##                                    (Source Data Fig. 7h)
##   Fig7h_GRIA1_reporter_test.txt    linear mixed model result
##
## METHOD
##   Linear mixed model rel_expr_adj ~ utr + (1 | experiment) + (1 | batch) (lme4), two-sided
##   t-test of the utr coefficient with Satterthwaite degrees of freedom (lmerTest).
##   Figure: values divided by the mean of the long 3'UTR transfections of the same batch
##   (display only); lines connect the transfections of an experiment.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(need_file(file.path(REPO, "resources", "gria1_reporter", "gria1_utr_reporter_qpcr.txt"), "resources"),
                stringsAsFactors = FALSE)
stopifnot(all(table(d$experiment) == 2),
          all(tapply(d$utr, d$experiment, function(x) setequal(x, c("Long", "Short")))),
          all(d$gfp_sort == "pos"))
d$utr <- factor(d$utr, levels = c("Long", "Short"))

## ------------------------------------------------------------ test
m  <- lmer(rel_expr_adj ~ utr + (1 | experiment) + (1 | batch), data = d)
co <- summary(m)$coefficients["utrShort", ]
test <- data.frame(model = "rel_expr_adj ~ utr + (1 | experiment) + (1 | batch)",
                   n_experiments = length(unique(d$experiment)), n_batches = length(unique(d$batch)),
                   estimate_Short_minus_Long = co[["Estimate"]], se = co[["Std. Error"]],
                   df = co[["df"]], t = co[["t value"]], p_value = co[["Pr(>|t|)"]])

## ------------------------------------------------------------ batch normalisation (display)
long_mean <- tapply(d$rel_expr_adj[d$utr == "Long"], d$batch[d$utr == "Long"], mean)
d$rel_expr_batchNorm <- d$rel_expr_adj / long_mean[d$batch]

write.table(d[, c("batch", "experiment", "label", "utr", "plasmid", "rel_expr_adj", "rel_expr_sd", "rel_expr_batchNorm")],
            file.path(OUT_DIR, "Fig7h_GRIA1_reporter.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(test, file.path(OUT_DIR, "Fig7h_GRIA1_reporter_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ figure
ymax <- max(d$rel_expr_batchNorm)
p <- ggplot(d, aes(utr, rel_expr_batchNorm)) +
  geom_boxplot(aes(fill = utr), width = 0.55, outlier.shape = NA, linewidth = 0.4, colour = FIG_INK) +
  geom_line(aes(group = experiment), colour = FIG_MUT, linewidth = 0.3) +
  geom_point(size = 1.6, colour = FIG_INK) +
  annotate("segment", x = 1, xend = 2, y = ymax * 1.06, yend = ymax * 1.06, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ymax * 1.12, label = sprintf("p = %.3f", test$p_value), size = 3, colour = FIG_SEC) +
  scale_fill_manual(values = c(Long = "grey80", Short = "#1f78c8"), guide = "none") +
  scale_y_continuous(limits = c(0, ymax * 1.18)) +
  labs(x = "GRIA1 3'UTR", y = "GFP-GRIA1 3'UTR mRNA\n(relative to long 3'UTR)") +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(OUT_DIR, "Fig7h_GRIA1_reporter.pdf"), p, width = 2.2, height = 3)
ggsave(file.path(OUT_DIR, "Fig7h_GRIA1_reporter.png"), p, width = 2.2, height = 3, dpi = 300)
print(test, row.names = FALSE, digits = 3)
