#!/usr/bin/env Rscript
## ===========================================================================
## 02_wb_gria1_3utr.R
##
## Western blot quantification of GRIA1 in iNs with a deletion of the GRIA1 long 3'UTR
## (CRISPR, two clones) vs the unedited parental line, day 49.
##
## Manuscript: Fig. 7e.
##
## USAGE
##   Rscript 07_western_blot/02_wb_gria1_3utr.R [--out <dir>]
##
## INPUT (resources/western_blot/)
##   wb_gria1_3utr.txt   one row per capillary: batch (differentiation batch), replicate
##                       (independent wells and protein extraction of that batch), wb_date, blots,
##                       capillary, line (PluriCore ID), sample (WT, KO_Cl13, KO_Cl19), condition
##                       (WT / KO), NRI (GRIA1 120 + 150 kDa peak area / neurofilament L peak area
##                       of the same sample)
##
## OUTPUTS (results/figures/)
##   Fig7e_WB_GRIA1.pdf / .png
##   Fig7e_WB_GRIA1.txt         per capillary: areas and NRI (Source Data Fig. 7e)
##   Fig7e_WB_GRIA1_test.txt    linear mixed model result
##
## METHOD
##   Both KO clones are pooled as KO. Linear mixed model log2(NRI) ~ condition + (1 | batch)
##   (lme4), two-sided t-test of the condition coefficient with Satterthwaite degrees of
##   freedom (lmerTest).
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(need_file(file.path(REPO, "resources", "western_blot", "wb_gria1_3utr.txt"), "resources"),
                stringsAsFactors = FALSE)
stopifnot(all(tapply(d$condition, paste(d$batch, d$replicate), function(x) "WT" %in% x && "KO" %in% x)))
d$condition <- factor(d$condition, levels = c("WT", "KO"))

## ------------------------------------------------------------ test
m  <- suppressMessages(lmer(log2(NRI) ~ condition + (1 | batch), data = d))
co <- summary(m)$coefficients["conditionKO", ]
test <- data.frame(model = "log2(NRI) ~ condition + (1 | batch)",
                   n_WT = sum(d$condition == "WT"), n_KO = sum(d$condition == "KO"),
                   n_replicates = length(unique(paste(d$batch, d$replicate))), n_batches = length(unique(d$batch)),
                   log2FC_KO_vs_WT = co[["Estimate"]], se = co[["Std. Error"]],
                   df = co[["df"]], t = co[["t value"]], p_value = co[["Pr(>|t|)"]])

write.table(d[, c("batch", "replicate", "wb_date", "capillary", "line", "sample", "condition", "NRI")],
            file.path(OUT_DIR, "Fig7e_WB_GRIA1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(test, file.path(OUT_DIR, "Fig7e_WB_GRIA1_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ figure
ymax <- max(d$NRI)
p <- ggplot(d, aes(condition, NRI)) +
  geom_boxplot(aes(fill = condition), width = 0.55, outlier.shape = NA, linewidth = 0.4, colour = FIG_INK) +
  geom_point(size = 1.6, colour = FIG_INK) +
  annotate("segment", x = 1, xend = 2, y = ymax * 1.06, yend = ymax * 1.06, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ymax * 1.12, label = sprintf("p = %.3f", test$p_value), size = 3, colour = FIG_SEC) +
  scale_fill_manual(values = c(WT = "grey80", KO = "#1f78c8"), guide = "none") +
  scale_y_continuous(limits = c(0, ymax * 1.18)) +
  labs(x = NULL, y = "GRIA1 / NF-L [WB]") +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(OUT_DIR, "Fig7e_WB_GRIA1.pdf"), p, width = 2, height = 3)
ggsave(file.path(OUT_DIR, "Fig7e_WB_GRIA1.png"), p, width = 2, height = 3, dpi = 300)
print(test, row.names = FALSE, digits = 3)
