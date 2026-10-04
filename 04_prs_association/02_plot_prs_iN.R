#!/usr/bin/env Rscript
## ===========================================================================
## 02_plot_prs_iN.R
##
## Scatter plots of PRS vs cumulative APA per iN donor, all with the same
## procedure (mean |z| across the gene set, Tukey-filtered donors, linear fit
## with standard error):
##   Fig. 4b        SCZ-PRS vs cumulative SCZ dAPA
##   ED Fig. 6c     BD-PRS  vs cumulative SCZ dAPA
##   ED Fig. 6d     MDD-PRS vs cumulative SCZ dAPA
##   ED Fig. 6e     BD-PRS  vs cumulative BD APA (nominal p <= 0.05 genes)
##
## Input: the per-donor table written by 01_prs_apa_iN.R (contains PRS; private).
## Output: PDF panels in <private-out>/figures/.
##
## USAGE
##   Rscript 04_prs_association/02_plot_prs_iN.R [--private-out <dir>]
## Dependencies: ggplot2
## ===========================================================================

suppressPackageStartupMessages(library(ggplot2))
.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()

PRIV <- get_arg("--private-out", file.path(REPO, "results", "private"))
d    <- read.delim(need_file(file.path(PRIV, "prs_apa_per_donor_iN.txt"), "--private-out"),
                   stringsAsFactors = FALSE)
FIG  <- file.path(PRIV, "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

PANELS <- list(
  Fig4b    = list(set = "SCZ_dAPA_q05", score = "PRS_SCZ", x = "SCZ-PRS", t = 'SCZ~dAPA~transcripts~(italic(q) <= 0.05)'),
  EDFig6c  = list(set = "SCZ_dAPA_q05", score = "PRS_BD",  x = "BD-PRS",  t = 'SCZ~dAPA~transcripts~(italic(q) <= 0.05)'),
  EDFig6d  = list(set = "SCZ_dAPA_q05", score = "PRS_MDD", x = "MDD-PRS", t = 'SCZ~dAPA~transcripts~(italic(q) <= 0.05)'),
  EDFig6e  = list(set = "BD_dAPA_p05",  score = "PRS_BD",  x = "BD-PRS",  t = 'BD~APA~transcripts~(nominal~italic(p) <= 0.05)'))

for (nm in names(PANELS)) {
  p  <- PANELS[[nm]]
  df <- data.frame(prs = d[[p$score]], apa = d[[paste0(p$set, "_mean_abs_z")]],
                   keep = d[[paste0(p$set, "_tukey_keep")]])
  df <- df[df$keep, ]
  ft <- summary(lm(apa ~ prs, df))$coefficients
  lab <- sprintf("n = %d\nr = %.3f\np = %.3f", nrow(df), cor(df$apa, df$prs), ft[2, 4])
  g <- ggplot(df, aes(prs, apa)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#2b6cb0", fill = "grey75", linewidth = 0.7) +
    geom_point(size = 1.8, colour = "black") +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.1, vjust = 1.2, size = 3.2) +
    labs(x = p$x, y = "Cumulative APA (mean |z|)", title = parse(text = p$t)[[1]]) +
    theme_classic(base_size = 10) + theme(plot.title = element_text(size = 9, face = "plain"))
  ggsave(file.path(FIG, paste0(nm, "_PRS_cumulativeAPA.pdf")), g, width = 3.2, height = 3.2, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(nm, "_PRS_cumulativeAPA.png")), g, width = 3.2, height = 3.2, dpi = 300, type = "cairo")
  cat(sprintf("%-8s %s\n", nm, gsub("\n", ", ", lab)))
}
