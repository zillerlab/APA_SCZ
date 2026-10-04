#!/usr/bin/env Rscript
## ===========================================================================
## 01_dep_olink_rbp_knockdown.R
##
## Differential protein abundance (Olink Explore HT) after shRNA knockdown of PTBP2 and
## PCBP2 vs lentivirus-infected (scramble) control in iNs, and the number of
## differentially abundant proteins compared with a within-batch label permutation null.
##
## Manuscript: Fig. 5i; main text (proteomics).
##
## USAGE
##   Rscript 06_proteomics_olink/01_dep_olink_rbp_knockdown.R [--cores <n>] [--out <dir>]
##
## INPUTS (resources/olink/)
##   olink_npx_rbp_knockdown.txt.gz    NPX, 5,415 proteins x 9 samples; in_APA_set = protein
##                                     is one of the 936 proteins of the pre-specified test set
##                                     (genes with two PAS in the 3'UTR reference)
##   olink_samples_rbp_knockdown.txt   donor (PluriCore ID), condition, batch, Olink QC flag
##
## OUTPUTS
##   results/proteomics/DEP_OLINK_RBP_KD.txt          one row per protein: mean NPX; per
##       knockdown logFC (knockdown - control, NPX), moderated t, p, BH q over the 936
##       proteins of the test set (q_APAset) and over all proteins (q_allProteins)
##   results/proteomics/DEP_OLINK_permutation.txt     proteins with p <= 0.05 among the 936,
##       observed vs permutation null
##   results/figures/Fig5i_DEP_counts.pdf / .png / .txt
##
## METHOD
##   limma: NPX ~ condition + batch + mean sample NPX (centred; absorbs global loading
##   differences), reference = control; moderated t (empirical Bayes). Counts: proteins of
##   the 936-protein test set with p <= 0.05 (expected by chance: 936 x 0.05 = 47).
##   Permutation null: the condition labels are permuted within batch over all 540
##   distinct arrangements, the model is refitted and the counts recomputed;
##   p = fraction of arrangements with a count >= the observed one (the observed
##   arrangement included).
##   Samples: 9 (control 3, PTBP2 3, PCBP2 3; two donors, two batches). The CPSF6
##   knockdown samples and their batch are not used.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "figure_functions.R"))
suppressPackageStartupMessages(library(limma))

CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "proteomics"))
FIG_DIR <- file.path(REPO, "results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
KDS <- c("PTBP2", "PCBP2"); FIG_KDS <- c("PCBP2", "PTBP2"); P_CUT <- 0.05

## ------------------------------------------------------------ inputs
rs  <- file.path(REPO, "resources", "olink")
x   <- read.delim(need_file(file.path(rs, "olink_npx_rbp_knockdown.txt.gz"), "resources"), check.names = FALSE,
                  stringsAsFactors = FALSE)
si  <- read.delim(need_file(file.path(rs, "olink_samples_rbp_knockdown.txt"), "resources"), stringsAsFactors = FALSE)
npx <- as.matrix(x[, si$sample]); rownames(npx) <- x$protein
inset <- as.logical(x$in_APA_set)
si$condition <- relevel(factor(si$condition), ref = "CTRL")
si$batch <- factor(si$batch)
si$mean_npx <- colMeans(npx); si$mean_npx <- si$mean_npx - mean(si$mean_npx)
print(table(condition = si$condition, batch = si$batch))

## ------------------------------------------------------------ model
fit_counts <- function(cond, full = FALSE) {
  d <- data.frame(condition = relevel(factor(cond), ref = "CTRL"), batch = si$batch, mean_npx = si$mean_npx)
  f <- eBayes(lmFit(npx, model.matrix(~ condition + batch + mean_npx, d)))
  tt <- lapply(setNames(KDS, KDS), function(k) topTable(f, coef = paste0("condition", k), number = Inf, sort.by = "none"))
  if (full) return(list(fit = f, tt = tt))
  vapply(tt, function(t) sum(t$P.Value[inset] <= P_CUT), numeric(1))
}
obs <- fit_counts(as.character(si$condition), full = TRUE)
res <- data.frame(protein = rownames(npx), in_APA_set = inset, mean_NPX = obs$fit$Amean)
for (k in KDS) {
  t <- obs$tt[[k]]
  q <- rep(NA_real_, nrow(t)); q[inset] <- p.adjust(t$P.Value[inset], "BH")
  res[[paste0(k, "_logFC")]] <- t$logFC; res[[paste0(k, "_t")]] <- t$t; res[[paste0(k, "_p")]] <- t$P.Value
  res[[paste0(k, "_q_APAset")]] <- q; res[[paste0(k, "_q_allProteins")]] <- p.adjust(t$P.Value, "BH")
}
write.table(res, file.path(OUT_DIR, "DEP_OLINK_RBP_KD.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
n_obs <- vapply(KDS, function(k) sum(res[[paste0(k, "_p")]][inset] <= P_CUT), numeric(1))

## ------------------------------------------------------------ within-batch permutation null
perms <- function(v) if (length(v) <= 1) list(v) else
  unique(do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(r) c(v[i], r)))))
idx <- split(seq_len(nrow(si)), si$batch)
pb  <- lapply(idx, function(ix) perms(as.character(si$condition[ix])))
grid <- expand.grid(lapply(pb, seq_along))
null <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(r) {
  cond <- as.character(si$condition)
  for (b in names(idx)) cond[idx[[b]]] <- pb[[b]][[grid[r, b]]]
  fit_counts(cond)
}, mc.cores = CORES))
perm <- data.frame(knockdown = KDS, n_test_set = sum(inset), n_p_le_0.05 = n_obs,
                   expected_by_chance = sum(inset) * P_CUT, null_mean = colMeans(null),
                   null_95th = apply(null, 2, quantile, 0.95),
                   p_permutation = colMeans(sweep(null, 2, n_obs, ">=")), n_permutations = nrow(null))
write.table(perm, file.path(OUT_DIR, "DEP_OLINK_permutation.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------ Fig. 5i
fd <- perm[match(FIG_KDS, perm$knockdown), ]
fd$knockdown <- factor(fd$knockdown, levels = FIG_KDS)
write.table(fd, file.path(FIG_DIR, "Fig5i_DEP_counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(fd, aes(knockdown, n_p_le_0.05)) +
  geom_col(fill = "#1f78c8", width = 0.8) +
  geom_hline(yintercept = sum(inset) * P_CUT, linetype = "dashed", colour = FIG_MUT, linewidth = 0.4) +
  geom_text(aes(label = sprintf("p = %.3f", p_permutation)), vjust = -0.5, size = 3, colour = FIG_SEC) +
  scale_y_continuous(limits = c(0, 200), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Differentially abundant proteins (KD - Ctrl)",
       caption = sprintf("p <= 0.05, %d proteins\ndashed: expected by chance\np: within-batch permutation", sum(inset))) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = FIG_SEC), plot.caption = element_text(size = 7, colour = FIG_SEC, hjust = 0))
ggsave(file.path(FIG_DIR, "Fig5i_DEP_counts.pdf"), p, width = 2.6, height = 3.4)
ggsave(file.path(FIG_DIR, "Fig5i_DEP_counts.png"), p, width = 2.6, height = 3.4, dpi = 300)

## ------------------------------------------------------------ summary
cat(sprintf("%d samples, %d proteins (%d in the test set), residual df %d, prior df %.2f\n", ncol(npx), nrow(npx),
            sum(inset), ncol(npx) - ncol(model.matrix(~ condition + batch + mean_npx, si)), obs$fit$df.prior))
print(perm, row.names = FALSE, digits = 3)
