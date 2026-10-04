#!/usr/bin/env Rscript
## ===========================================================================
## 04_plot_prs_dlpfc.R
##
## Extended Data Fig. 6i, from the per-donor table of 03_prs_apa_dlpfc.R
## (contains PRS; private):
##   left   SCZ-PRS vs cumulative APA across DLPFC dAPA transcripts
##   right  SCZ-PRS by diagnosis (Wilcoxon test, SCZ vs Control)
## Output: <private-out>/figures/EDFig6i_{left,right}.{pdf,png}
##
## USAGE   Rscript 04_prs_association/04_plot_prs_dlpfc.R [--private-out <dir>]
## Dependencies: ggplot2
## ===========================================================================

suppressPackageStartupMessages(library(ggplot2))
.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))

PRIV <- get_arg("--private-out", file.path(REPO, "results", "private"))
d    <- read.delim(need_file(file.path(PRIV, "prs_apa_per_donor_dlpfc.txt"), "--private-out"),
                   stringsAsFactors = FALSE)
FIG  <- file.path(PRIV, "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
save <- function(g, stem, w = 3.2, h = 3.2) {
  ggsave(file.path(FIG, paste0(stem, ".pdf")), g, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(stem, ".png")), g, width = w, height = h, dpi = 300, type = "cairo")
}

## left: SCZ-PRS vs cumulative APA (Tukey-filtered donors, linear fit with SE)
df  <- d[d$tukey_keep, ]
ft  <- summary(lm(mean_abs_z ~ PRS_SCZ, df))$coefficients
lab <- sprintf("n = %d\nr = %.3f\np = %.3f", nrow(df), cor(df$mean_abs_z, df$PRS_SCZ), ft[2, 4])
g <- ggplot(df, aes(PRS_SCZ, mean_abs_z)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2b6cb0", fill = "grey75", linewidth = 0.7) +
  geom_point(size = 1.2, colour = "black", alpha = 0.7) +
  annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 3.2) +
  labs(x = "SCZ-PRS", y = "Cumulative APA (mean |z|)",
       title = parse(text = "DLPFC~dAPA~transcripts~(italic(q) <= 0.1)")[[1]]) +
  theme_classic(base_size = 10) + theme(plot.title = element_text(size = 9, face = "plain"))
save(g, "EDFig6i_left")

## right: SCZ-PRS by diagnosis
d$diagnosis <- factor(d$diagnosis, levels = intersect(c("Control", "SCZ", "BP", "AFF"), unique(d$diagnosis)))
p <- wilcox.test(d$PRS_SCZ[d$diagnosis == "SCZ"], d$PRS_SCZ[d$diagnosis == "Control"])$p.value
n <- table(d$diagnosis)
g <- ggplot(d, aes(diagnosis, PRS_SCZ)) +
  geom_violin(fill = "grey90", colour = "grey40", linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.4) +
  scale_x_discrete(labels = sprintf("%s\n(n = %d)", names(n), as.vector(n))) +
  labs(x = NULL, y = "SCZ-PRS",
       title = sprintf("SCZ vs Control: Wilcoxon p = %s", format(signif(p, 2), scientific = TRUE))) +
  theme_classic(base_size = 10) + theme(plot.title = element_text(size = 9, face = "plain"))
save(g, "EDFig6i_right")
cat(gsub("\n", ", ", lab), "| Wilcoxon SCZ vs Control p =", signif(p, 3), "\n")
