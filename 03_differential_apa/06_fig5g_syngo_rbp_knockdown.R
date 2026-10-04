#!/usr/bin/env Rscript
## ===========================================================================
## 06_fig5g_syngo_rbp_knockdown.R
##
## SynGO enrichment of the knockdown dAPA genes (PTBP2, PCBP2, CPSF6): dot plot.
##
## Manuscript: Fig. 5g.
##
## USAGE
##   Rscript 03_differential_apa/06_fig5g_syngo_rbp_knockdown.R [--out <dir>]
##
## INPUTS (resources/rbp_knockdown/syngo/)
##   syngo_1.3_dAPA_<RBP>_bh1e-4.zip SynGO portal (https://www.syngoportal.org, release 1.3)
##                                  output for the knockdown dAPA genes of each arm
##                                  (BH <= 1e-4; results/dapa/dAPA_RBP_KD_genes_<RBP>.txt
##                                  of 03_dapa_rbp_knockdown.R), default settings,
##                                  "brain expressed genes" background. An arm without
##                                  portal output is shown as not tested. Used: readme.txt
##                                  (input and background size) and
##                                  syngo_ontologies_with_annotations_matching_user_input.xlsx
##                                  (per-term counts, GSEA p-value and FDR).
##
## OUTPUTS (results/figures/)
##   Fig5g_syngo_dotplot.pdf / .png
##   Fig5g_syngo_effectsizes.txt    one row per RBP x SynGO term (Source Data Fig. 5g);
##                                  arms without portal output carry no statistics
##   Fig5g_syngo_terms_shown.txt    terms in the figure
##
## METHOD
##   p-values and FDR are the portal's. Effect sizes from the portal's counts, input
##   list (n genes) vs background (N genes): k = input genes in the term, K = background
##   genes in the term; fold = (k / n) / (K / N); odds ratio and 95% CI from Fisher's
##   exact test on the 2 x 2 table (k, n - k, K - k, N - n - K + k).
##   Figure: union of the 10 most significant terms (lowest FDR; FDR <= 0.05 on >= 3
##   genes) of each knockdown; all terms are in Fig5g_syngo_effectsizes.txt.
##   Colour = fold enrichment (comparable across knockdowns); size = -log10 FDR, which
##   depends on the size of the input list and is not comparable across knockdowns.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages(library(readxl))

OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RBPS <- c("PTBP2", "PCBP2", "CPSF6")
SYNGO  <- file.path(REPO, "resources", "rbp_knockdown", "syngo")
PORTAL <- RBPS[file.exists(file.path(SYNGO, sprintf("syngo_1.3_dAPA_%s_bh1e-4.zip", RBPS)))]   # arms with portal output
if (!length(PORTAL)) stop("no SynGO portal output in ", SYNGO, " (syngo_1.3_dAPA_<RBP>_bh1e-4.zip)")
cfg <- list(fdr_cut = 0.05, min_k = 3, top_n_per_kd = 10, domains = c("CC", "BP"),
            size_range = c(1.6, 6.0), size_breaks = c(1.3, 3, 6, 8), size_cap = 8,
            fold_breaks = c(0.25, 1, 4, 16), pdf_width = 7.4, row_height = 0.175, base_height = 1.6)

## ------------------------------------------------------------ portal output -> effect sizes
read_portal <- function(A) {
  z <- need_file(file.path(SYNGO, sprintf("syngo_1.3_dAPA_%s_bh1e-4.zip", A)), "resources")
  tmp <- file.path(tempdir(), A); unzip(z, exdir = tmp)
  rd <- readLines(file.path(tmp, "readme.txt"))
  n <- as.integer(sub(".* / ([0-9]+) genes from your gene list.*", "\\1", grep("genes from your gene list", rd, value = TRUE)))
  N <- as.integer(sub(".*It contains ([0-9]+) unique genes in total.*", "\\1", grep("unique genes in total", rd, value = TRUE)))
  x <- as.data.frame(read_excel(file.path(tmp, "syngo_ontologies_with_annotations_matching_user_input.xlsx")))
  d <- data.frame(RBP = A, domain = x[["GO domain"]], GO_ID = x[["GO term ID"]], term = x[["GO term name"]],
                  k = as.integer(x[["GSEA count foreground/input"]]), term_size = as.integer(x[["GSEA count background"]]),
                  p_portal = suppressWarnings(as.numeric(x[["GSEA p-value"]])),
                  FDR_portal = suppressWarnings(as.numeric(x[["GSEA FDR corrected p-value"]])),
                  n_input = n, n_background = N, genes = x[["genes - hgnc_symbol"]], stringsAsFactors = FALSE)
  d$fold <- ifelse(d$k > 0, (d$k / n) / (d$term_size / N), NA_real_)
  or <- t(vapply(seq_len(nrow(d)), function(i) {
    if (is.na(d$p_portal[i])) return(rep(NA_real_, 3))
    k <- d$k[i]; K <- d$term_size[i]
    ft <- fisher.test(matrix(c(k, n - k, K - k, N - n - K + k), 2))
    c(unname(ft$estimate), ft$conf.int)
  }, numeric(3)))
  d$OR <- or[, 1]; d$OR_lo95 <- or[, 2]; d$OR_hi95 <- or[, 3]
  d
}
d <- do.call(rbind, lapply(PORTAL, read_portal))
## arms without portal output: every term not tested
for (A in setdiff(RBPS, PORTAL)) {
  g <- readLines(need_file(file.path(REPO, "results", "dapa", sprintf("dAPA_RBP_KD_genes_%s.txt", A)), "03_dapa_rbp_knockdown.R"))
  e <- d[d$RBP == PORTAL[1], ]
  e$RBP <- A; e$n_input <- length(g); e$genes <- NA_character_
  e[, c("k", "p_portal", "FDR_portal", "n_background", "fold", "OR", "OR_lo95", "OR_hi95")] <- NA
  d <- rbind(d, e)
}
d <- d[order(match(d$RBP, RBPS)), ]
d$tested <- !is.na(d$FDR_portal)
d$significant <- d$tested & !is.na(d$k) & d$FDR_portal <= cfg$fdr_cut & d$k >= cfg$min_k
d$mlogFDR <- -log10(d$FDR_portal)
write.table(d[, c("RBP", "domain", "GO_ID", "term", "k", "term_size", "n_input", "n_background", "fold", "OR", "OR_lo95",
                  "OR_hi95", "p_portal", "FDR_portal", "mlogFDR", "significant", "genes")],
            file.path(OUT_DIR, "Fig5g_syngo_effectsizes.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ------------------------------------------------------------ terms shown
d <- d[d$domain %in% cfg$domains, ]
d$thin <- d$significant & !is.na(d$k) & d$k < 5
key <- paste(d$domain, d$GO_ID, sep = "|")
keep <- do.call(rbind, lapply(split(d, key), function(s) data.frame(
  domain = s$domain[1], GO_ID = s$GO_ID[1], term = s$term[1], n_sig = sum(s$significant),
  max_fold = if (any(s$significant)) max(s$fold[s$significant]) else NA_real_,
  min_fdr = if (any(s$significant)) min(s$FDR_portal[s$significant]) else NA_real_)))
top <- unlist(lapply(split(d[d$significant, ], d$RBP[d$significant]), function(s)
  head(paste(s$domain, s$GO_ID, sep = "|")[order(s$FDR_portal, -s$fold)], cfg$top_n_per_kd)))
keep <- keep[paste(keep$domain, keep$GO_ID, sep = "|") %in% top, ]
keep <- keep[order(-keep$max_fold), ]
if (!nrow(keep)) stop("no SynGO term passes the filter")
write.table(keep, file.path(OUT_DIR, "Fig5g_syngo_terms_shown.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
## rows: strongest fold at the top of each domain block
lev <- keep$term[order(factor(keep$domain, c("CC", "BP")), keep$max_fold, decreasing = c(TRUE, FALSE), method = "radix")]
p <- d[paste(d$domain, d$GO_ID) %in% paste(keep$domain, keep$GO_ID), ]
p$term <- factor(p$term, levels = unique(lev))
p$RBP <- factor(p$RBP, levels = RBPS)
p$domain <- factor(p$domain, levels = c("CC", "BP"), labels = c("Cellular component", "Biological process"))
p$log2fold <- log2(p$fold)

## ------------------------------------------------------------ plot
pal_div <- c("#9e2726", "#c93736", "#de5a59", "#e88080", "#f2a9a8", "#fbd9d8",
             "#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b")
vlim <- max(abs(p$log2fold[p$tested]), na.rm = TRUE)
nin <- tapply(d$n_input, d$RBP, `[`, 1)[RBPS]
gg <- ggplot(p, aes(x = RBP, y = term)) +
  geom_blank() +                                      # trains the y scale on the full term order
  geom_point(data = p[!p$tested, ], shape = 4, size = 1.3, colour = "#d8d7d1", stroke = 0.5) +
  geom_point(data = p[p$tested, ], aes(size = pmin(mlogFDR, cfg$size_cap), fill = log2fold,
                                       colour = significant, stroke = significant), shape = 21) +
  geom_point(data = p[p$thin, ], shape = 21, size = 0.45, fill = "white", colour = "white", stroke = 0) +
  scale_colour_manual(values = c(`TRUE` = "#0b0b0b", `FALSE` = "#b9b8b0"), guide = "none") +
  scale_discrete_manual("stroke", values = c(`TRUE` = 0.55, `FALSE` = 0.3), guide = "none") +
  scale_fill_gradientn(colours = pal_div, limits = c(-vlim, vlim), breaks = log2(cfg$fold_breaks),
                       labels = cfg$fold_breaks, name = "fold\nenrichment",
                       guide = guide_colourbar(barwidth = 0.45, barheight = 6, order = 1)) +
  scale_size_continuous(range = cfg$size_range, limits = c(0, cfg$size_cap), breaks = cfg$size_breaks,
                        labels = cfg$size_breaks, name = expression(-log[10]~FDR),
                        guide = guide_legend(override.aes = list(fill = "#9ec5f4", colour = "#0b0b0b", stroke = 0.55), order = 2)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.6), drop = FALSE) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "SynGO enrichment of RBP knockdown dAPA genes",
       subtitle = sprintf("%d terms: union of the %d most significant terms per knockdown (FDR <= %.2g on >= %d genes)",
                          nrow(keep), cfg$top_n_per_kd, cfg$fdr_cut, cfg$min_k),
       caption = paste0("Colour is fold enrichment (observed / expected overlap) and is comparable across columns. Size is -log10 FDR and is not:\n",
                        sprintf("the input sets differ in size (%s knockdown dAPA genes, BH <= 1e-4), so FDR tracks power as much as effect.\n",
                                paste(sprintf("%s n = %s", RBPS, format(nin, big.mark = ",", trim = TRUE)), collapse = ", ")),
                        "Black outline = significant; grey outline = not; white centre = fewer than 5 overlapping genes; x = not tested (fewer than 3 matching genes).\n",
                        "SynGO v1.3, default settings, brain-expressed background."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), panel.grid.major.y = element_line(colour = "#f2f1ed", linewidth = 0.25),
        panel.spacing.y = unit(6, "pt"), strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 0, size = 7, colour = "#52514e"),
        axis.text.x.top = element_text(face = "bold", size = 8.5, colour = "#0b0b0b"),
        axis.text.y = element_text(size = 7, colour = "#0b0b0b"), axis.ticks = element_blank(),
        plot.title = element_text(face = "bold", size = 10, colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 7.2, colour = "#52514e", margin = margin(b = 8)),
        plot.caption = element_text(size = 6, colour = "#52514e", hjust = 0, lineheight = 1.4, margin = margin(t = 10)),
        plot.caption.position = "plot", plot.title.position = "plot",
        legend.title = element_text(size = 7.2, colour = "#52514e"), legend.text = element_text(size = 6.6, colour = "#52514e"),
        legend.key.height = unit(11, "pt"), plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 10, 8, 10))
h <- cfg$base_height + cfg$row_height * nrow(keep)
ggsave(file.path(OUT_DIR, "Fig5g_syngo_dotplot.pdf"), gg, width = cfg$pdf_width, height = h, device = grDevices::pdf, bg = "white")
ggsave(file.path(OUT_DIR, "Fig5g_syngo_dotplot.png"), gg, width = cfg$pdf_width, height = h, dpi = 300, bg = "white")

## ------------------------------------------------------------ summary
cat(sprintf("input genes: %s\n", paste(RBPS, nin, collapse = ", ")))
cat(sprintf("terms significant (FDR <= %.2g, k >= %d): %s; shown: %d\n", cfg$fdr_cut, cfg$min_k,
            paste(RBPS, tapply(d$significant, factor(d$RBP, RBPS), sum), collapse = ", "), nrow(keep)))
