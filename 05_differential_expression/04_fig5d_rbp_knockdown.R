#!/usr/bin/env Rscript
## ===========================================================================
## 04_fig5d_rbp_knockdown.R
##
## Fig. 5d: overlap of the DEG sets after PTBP2, PCBP2 and CPSF6 knockdown.
##
## USAGE
##   Rscript 05_differential_expression/04_fig5d_rbp_knockdown.R [--deg <dir>] [--out <dir>]
##
## INPUTS
##   results/deg/DEG_RBP_KD_<PTBP2|PCBP2|CPSF6>.txt   02_deseq_rbp_knockdown.R
##
## OUTPUTS (results/figures/)
##   Fig5d_DEG_overlap.pdf / .png     Venn diagram of the DEGs (padj <= 0.05, |log2FC| >= 0.4)
##   Fig5d_DEG_membership.txt         gene x knockdown membership, log2FC, p and padj per
##                                    knockdown (Source Data Fig. 5d); genes tested in at
##                                    least one knockdown
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))    # plot_venn3()

DEG_DIR <- get_arg("--deg", file.path(REPO, "results", "deg"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "figures"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RBPS <- c("PTBP2", "PCBP2", "CPSF6")

de <- lapply(setNames(RBPS, RBPS), function(A)
  read.delim(need_file(file.path(DEG_DIR, sprintf("DEG_RBP_KD_%s.txt", A)), "02_deseq_rbp_knockdown.R"),
             stringsAsFactors = FALSE))
univ <- sort(Reduce(union, lapply(de, `[[`, "gene_id")))
mem <- data.frame(gene_id = univ,
                  gene_name = Reduce(function(x, r) ifelse(is.na(x), r$gene_name[match(univ, r$gene_id)], x),
                                     de, rep(NA_character_, length(univ))),
                  sapply(de, function(r) as.integer(univ %in% r$gene_id[r$DEG != "NO"])))
mem$n_RBPs <- rowSums(mem[, RBPS])
for (A in RBPS) { r <- de[[A]]; i <- match(univ, r$gene_id)
  mem[[paste0("log2FoldChange_", A)]] <- r$log2FoldChange[i]
  mem[[paste0("pvalue_", A)]] <- r$pvalue[i]; mem[[paste0("padj_", A)]] <- r$padj[i] }
write.table(mem, file.path(OUT_DIR, "Fig5d_DEG_membership.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")

p <- plot_venn3(mem[, RBPS], expression(paste(Delta, "mRNA (RBP-KD vs. Ctrl)")))
ggsave(file.path(OUT_DIR, "Fig5d_DEG_overlap.pdf"), p, width = 3.4, height = 3.4)
ggsave(file.path(OUT_DIR, "Fig5d_DEG_overlap.png"), p, width = 3.4, height = 3.4, dpi = 300)

reg <- venn3_counts(mem[, RBPS])
cat(sprintf("Fig. 5d: %d genes tested in >= 1 knockdown; DEGs %s\n", length(univ),
            paste(RBPS, colSums(mem[, RBPS]), collapse = ", ")))
cat("Venn regions:", paste(names(reg), reg, collapse = ", "), "\n")
