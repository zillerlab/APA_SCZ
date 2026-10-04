#!/usr/bin/env Rscript
## ===========================================================================
## 02_fig6a_clip_dapa_overlap.R
##
## PTBP2 binding of 3'UTRs with SCZ-associated differential APA: overlap of dAPA genes
## (iN cohort; postmortem DLPFC) with PTBP2 CLIP peaks from three datasets.
##
## Manuscript: Fig. 6a.
##
## USAGE
##   Rscript 09_ptbp2_clip/02_fig6a_clip_dapa_overlap.R [--out <dir>]
##
## INPUTS
##   results/apa_pas_library/dominant_pas_utrs_hg38.txt   02_apa_pas_library (3'UTR reference)
##   results/dapa/dAPA_iN_SCZ_vs_Ctrl.txt                 03_differential_apa/01
##   results/dapa/dAPA_DLPFC_SCZ_vs_Control.txt           03_differential_apa/02
##   results/clip/CLIP_PTBP2_iN_enriched_peaks_hg38.bed   01_clip_enrichment.R
##   resources/clip/eclip_PTBP2_postmortem_PFC_peaks_hg38.bed.gz   published PTBP2 eCLIP peaks
##   resources/clip/eclip_PTBP2_iPSC_neurons_peaks_hg38.bed.gz     (manuscript ref. 44)
##
## OUTPUTS
##   results/clip/Fig6a_CLIP_dAPA_overlap.txt        per dAPA set x CLIP dataset: genes tested,
##       dAPA, bound, dAPA and bound, odds ratio, 95% CI, Fisher p (two-sided)
##   results/clip/Fig6a_CLIP_binding_<iN|DLPFC>.txt  per tested gene: dAPA flag and binding flag
##       per CLIP dataset (Source Data Fig. 6a)
##   results/figures/Fig6a_CLIP_dAPA_overlap.pdf / .png
##
## METHOD
##   3'UTR of a gene: from the start of its proximal to the end of its distal UTR segment
##   (two-PAS reference, 4,506 genes). A UTR is bound if it overlaps a CLIP peak (same
##   strand for the stranded eCLIP peaks; iN peaks are unstranded). Coordinates are used
##   as given in the files. Universe: genes of the reference tested in the respective dAPA
##   analysis. dAPA: q <= 0.05 (iN), q <= 0.1 (DLPFC), as in module 03. Fisher's exact
##   test (two-sided) of dAPA vs bound within the universe.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages(library(GenomicRanges))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "clip"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
DAPA <- list(iN    = list(file = "dAPA_iN_SCZ_vs_Ctrl.txt",       q = 0.05, label = "dAPA iN"),
             DLPFC = list(file = "dAPA_DLPFC_SCZ_vs_Control.txt", q = 0.10, label = "dAPA DLPFC"))
CLIP <- c(iN = "PTBP2 iN (this study)", iPSC_neurons = "PTBP2 iPSC neurons", postmortem_PFC = "PTBP2 postmortem PFC")

## ------------------------------------------------------------ 3'UTRs per gene
ref <- read.delim(need_file(file.path(REPO, "results", "apa_pas_library", "dominant_pas_utrs_hg38.txt"),
                            "02_apa_pas_library/05_select_dominant_pas.R"), stringsAsFactors = FALSE)
utr <- do.call(rbind, lapply(split(ref, factor(ref$geneId, levels = unique(ref$geneId))), function(t)
  data.frame(chrom = t$Chr[1], start = min(t$Start), end = max(t$End), strand = t$Strand[1], gene = t$geneId[1])))
gu  <- makeGRangesFromDataFrame(utr); names(gu) <- utr$gene

## ------------------------------------------------------------ CLIP peaks -> bound genes
read_bed <- function(p, stranded) {
  x <- read.delim(p, header = FALSE, stringsAsFactors = FALSE)
  makeGRangesFromDataFrame(data.frame(chrom = x[[1]], start = x[[2]], end = x[[3]],
                                      strand = if (stranded) x[[6]] else "*"))
}
rs <- file.path(REPO, "resources", "clip")
peaks <- list(
  iN             = read_bed(need_file(file.path(REPO, "results", "clip", "CLIP_PTBP2_iN_enriched_peaks_hg38.bed"), "01_clip_enrichment.R"), FALSE),
  iPSC_neurons   = read_bed(need_file(file.path(rs, "eclip_PTBP2_iPSC_neurons_peaks_hg38.bed.gz"), "resources"), TRUE),
  postmortem_PFC = read_bed(need_file(file.path(rs, "eclip_PTBP2_postmortem_PFC_peaks_hg38.bed.gz"), "resources"), TRUE))
bound <- sapply(peaks, function(g) suppressWarnings(overlapsAny(gu, g)))
rownames(bound) <- utr$gene

## ------------------------------------------------------------ Fisher tests
stats <- list()
for (d in names(DAPA)) {
  r <- read.delim(need_file(file.path(REPO, "results", "dapa", DAPA[[d]]$file), "03_differential_apa"), stringsAsFactors = FALSE)
  r <- r[!is.na(r$qvalue) & r$ID %in% utr$gene, ]
  tab <- data.frame(gene = r$ID, dAPA = r$qvalue <= DAPA[[d]]$q, bound[r$ID, , drop = FALSE], row.names = NULL)
  write.table(tab, file.path(OUT_DIR, sprintf("Fig6a_CLIP_binding_%s.txt", d)), sep = "\t", quote = FALSE, row.names = FALSE)
  for (cl in names(CLIP)) {
    f <- fisher.test(table(factor(tab$dAPA, c(FALSE, TRUE)), factor(tab[[cl]], c(FALSE, TRUE))))
    stats[[paste(d, cl)]] <- data.frame(dAPA_set = d, CLIP = cl, n_tested = nrow(tab), n_dAPA = sum(tab$dAPA),
      n_bound = sum(tab[[cl]]), n_dAPA_bound = sum(tab$dAPA & tab[[cl]]), odds_ratio = unname(f$estimate),
      or_lo95 = f$conf.int[1], or_hi95 = f$conf.int[2], p_value = f$p.value)
  }
}
stats <- do.call(rbind, stats)
write.table(stats, file.path(OUT_DIR, "Fig6a_CLIP_dAPA_overlap.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ Fig. 6a
fd <- stats
fd$x <- factor(sapply(DAPA, `[[`, "label")[fd$dAPA_set], levels = sapply(DAPA, `[[`, "label"))
fd$y <- factor(CLIP[fd$CLIP], levels = rev(CLIP))
fd$stars <- cut(fd$p_value, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
fd$lab <- sprintf("%d/%d\nOR %.2f %s", fd$n_dAPA_bound, fd$n_dAPA, fd$odds_ratio, fd$stars)
p <- ggplot(fd, aes(x, y, fill = log2(odds_ratio))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = lab), size = 2.6, colour = FIG_INK, lineheight = 0.9) +
  scale_fill_gradient2(low = "#1f78c8", mid = "white", high = "#c8321f", midpoint = 0, limits = c(-1, 1),
                       oob = scales::squish, name = "log2 OR") +
  labs(x = NULL, y = NULL, caption = "dAPA and bound / dAPA genes\nFisher's exact test: * p < 0.05, ** p < 0.01, *** p < 0.001") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = FIG_INK),
        plot.caption = element_text(size = 6.5, colour = FIG_SEC, hjust = 0), plot.caption.position = "plot")
ggsave(file.path(FIG_DIR, "Fig6a_CLIP_dAPA_overlap.pdf"), p, width = 4.6, height = 2.8)
ggsave(file.path(FIG_DIR, "Fig6a_CLIP_dAPA_overlap.png"), p, width = 4.6, height = 2.8, dpi = 300)

cat(sprintf("%d genes in the UTR reference; bound: %s\n", nrow(utr),
            paste(names(CLIP), colSums(bound), sep = " ", collapse = ", ")))
print(stats[, c("dAPA_set", "CLIP", "n_tested", "n_dAPA", "n_bound", "n_dAPA_bound", "odds_ratio", "p_value")],
      row.names = FALSE, digits = 3)
