#!/usr/bin/env Rscript
## ===========================================================================
## 02_mea_gria1_3utr.R
##
## Neuronal network activity (MEA) of iNs with a deletion of the GRIA1 long 3'UTR
## (CRISPR, two clones) vs the unedited parental line, day 49.
##
## Manuscript: Fig. 7d.
##
## USAGE
##   Rscript 10_mea/02_mea_gria1_3utr.R [--out <dir>]
##
## INPUT (resources/mea/)
##   mea_gria1_3utr_wells.txt   one row per well: well, plate, batch (differentiation batch),
##                              line (PluriCore ID), genotype (WT / GRIA1_3UTR_KO), clone
##                              (WT, KO_Cl13, KO_Cl19), then the well-averaged MEA parameters
##                              (Axion Navigator)
##
## OUTPUTS
##   results/mea/MEA_GRIA1_3UTR_LM.txt        per parameter: estimate (KO - WT, analysis scale),
##                                            SE, t, p, q (BH over parameters), wells used
##   results/figures/Fig7d_MEA_GRIA1.pdf / .png / .txt
##
## METHOD
##   Wells with >= 3 active electrodes and weighted mean electrode resistance >= 25 kOhm.
##   Number of spikes, bursts and network bursts and the burst peak rate are log2-transformed.
##   Parameters: all well averages except standard deviations, electrode counts, resistances
##   and parameters missing in > 5 wells (35 parameters). Per parameter, outliers are removed
##   at 1.5 x IQR of the values scaled to their batch mean; linear model
##   value ~ batch + electrode resistance + genotype (both KO clones vs WT), t-test of the
##   genotype coefficient; Benjamini-Hochberg over the parameters.
##   Figure: parameters with q <= 0.1; values adjusted for batch and electrode resistance
##   (limma::removeBatchEffect, design ~ genotype), display only.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages(library(limma))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "mea"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
LOG2 <- c("Number_of_Spikes", "Number_of_Bursts", "Number_of_Network_Bursts", "Burst_Peak_Max_Spikes_per_sec_")
EXCLUDE <- c("Number_of_Active_Electrodes", "Weighted_Mean_Resistance_kOhms_", "Resistance_Avg_kOhms_",
             "Number_of_Bursting_Electrodes", "Number_of_Covered_Electrodes", "Start_Electrode")
RES <- "Weighted_Mean_Resistance_kOhms_"; Q_FIG <- 0.1

## ------------------------------------------------------------ input
w <- read.delim(need_file(file.path(REPO, "resources", "mea", "mea_gria1_3utr_wells.txt"), "resources"),
                stringsAsFactors = FALSE, check.names = FALSE)
meas <- names(w)[(which(names(w) == "clone") + 1):ncol(w)]
for (p in intersect(LOG2, meas)) w[[p]] <- log2(w[[p]])
pars <- setdiff(meas[!grepl("Std", meas)], EXCLUDE)
pars <- pars[colSums(is.na(w[, pars])) <= 5]
w <- w[which(w$Number_of_Active_Electrodes >= 3 & w[[RES]] >= 25), ]
w$genotype <- relevel(factor(w$genotype), ref = "WT"); w$batch <- factor(w$batch)
cat(sprintf("%d wells (%s), %d batches, %d parameters\n", nrow(w),
            paste(names(table(w$clone)), table(w$clone), collapse = ", "), nlevels(w$batch), length(pars)))

## ------------------------------------------------------------ linear models
keep_wells <- function(p) {
  v <- w[[p]]; bm <- tapply(v, w$batch, mean, na.rm = TRUE)
  vn <- v / bm[as.character(w$batch)]
  lo <- quantile(vn, 0.25, na.rm = TRUE) - 1.5 * IQR(vn, na.rm = TRUE)
  hi <- quantile(vn, 0.75, na.rm = TRUE) + 1.5 * IQR(vn, na.rm = TRUE)
  !is.na(vn) & is.finite(vn) & vn >= lo & vn <= hi
}
res <- do.call(rbind, lapply(pars, function(p) {
  k <- keep_wells(p); d <- data.frame(value = w[[p]], batch = w$batch, res = w[[RES]], genotype = w$genotype)[k, ]
  co <- summary(lm(value ~ batch + res + genotype, d))$coefficients["genotypeGRIA1_3UTR_KO", ]
  data.frame(parameter = p, scale = ifelse(p %in% LOG2, "log2", "linear"), estimate = co[[1]], se = co[[2]],
             t = co[[3]], p_value = co[[4]], n_wells = nrow(d))
}))
res$q_value <- p.adjust(res$p_value, "BH")
res <- res[order(res$p_value), ]
write.table(res, file.path(OUT_DIR, "MEA_GRIA1_3UTR_LM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ Fig. 7d
sel <- res$parameter[res$q_value <= Q_FIG]
lab <- function(x) paste0(trimws(gsub("_+", " ", sub("_Avg(_sec_)?$", "", x))), ifelse(x %in% LOG2, " (log2)", ""))
fd <- do.call(rbind, lapply(sel, function(p) {
  k <- keep_wells(p) & !is.na(w[[p]])
  adj <- removeBatchEffect(matrix(w[[p]][k], nrow = 1), batch = w$batch[k], covariates = w[[RES]][k],
                           design = model.matrix(~ genotype, data.frame(genotype = w$genotype[k])))
  data.frame(parameter = lab(p), well = w$well[k], plate = w$plate[k], batch = w$batch[k], clone = w$clone[k],
             genotype = w$genotype[k], value = w[[p]][k], value_adjusted = as.numeric(adj))
}))
write.table(fd, file.path(FIG_DIR, "Fig7d_MEA_GRIA1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
st <- res[res$parameter %in% sel, ]; st$parameter <- lab(st$parameter)
st$label <- sprintf("p = %.2g\nq = %.2g", st$p_value, st$q_value)
ymax <- tapply(fd$value_adjusted, fd$parameter, max); ymin <- tapply(fd$value_adjusted, fd$parameter, min)
st$y <- ymax[st$parameter] + 0.12 * (ymax[st$parameter] - ymin[st$parameter])
fd$parameter <- factor(fd$parameter, levels = lab(sel)); st$parameter <- factor(st$parameter, levels = lab(sel))
p <- ggplot(fd, aes(genotype, value_adjusted)) +
  geom_violin(fill = "grey90", colour = FIG_MUT, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3) +
  geom_text(data = st, aes(x = 1.5, y = y, label = label), size = 2.5, colour = FIG_SEC, lineheight = 0.9) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3, labeller = label_wrap_gen(28)) +
  scale_x_discrete(labels = c(WT = "WT", GRIA1_3UTR_KO = "GRIA1\n3'UTR KO")) +
  labs(x = NULL, y = "Value (adjusted for batch and electrode resistance)") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), axis.text = element_text(colour = FIG_SEC))
ggsave(file.path(FIG_DIR, "Fig7d_MEA_GRIA1.pdf"), p, width = 6.2, height = 2.9 * ceiling(length(sel) / 3))
ggsave(file.path(FIG_DIR, "Fig7d_MEA_GRIA1.png"), p, width = 6.2, height = 2.9 * ceiling(length(sel) / 3), dpi = 300)

print(head(res[, c("parameter", "estimate", "p_value", "q_value", "n_wells")], 8), row.names = FALSE, digits = 3)
