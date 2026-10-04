## ===========================================================================
## prs_functions.R
##
## Association of polygenic risk scores (PRS) with the cumulative
## covariate-adjusted APA (or expression) of a gene set, per donor.
##
##   cumulative_score()  per-library cumulative score over a gene set
##   prs_association()   linear model / Pearson correlation with a PRS
##   permutation_null()  PRS-label permutation null of the correlation
## ===========================================================================

#' Cumulative score of a gene set per library.
#' @param M        genes x libraries matrix of covariate-adjusted values
#'                 (residuals of the covariate-only model); NA allowed
#' @param genes    gene set (rownames of M)
#' @param libs     libraries to include (the z-scores are computed across these)
#' @return data.frame(sample, mean_abs_resid, mean_abs_z, tukey_keep)
#'   mean_abs_resid mean |value| across the genes (NA ignored)
#'   mean_abs_z     mean |z| across the genes, z = per-gene standardisation
#'                  across `libs`; genes with any NA are left out
#'   tukey_keep     mean_abs_resid within Tukey's fences (1.5 x IQR)
cumulative_score <- function(M, genes, libs) {
  X <- M[rownames(M) %in% genes, libs, drop = FALSE]
  a <- colMeans(abs(X), na.rm = TRUE)
  Z <- t(apply(X, 1, function(x) (x - mean(x)) / stats::sd(x)))
  z <- colMeans(abs(Z), na.rm = TRUE)
  q <- stats::quantile(a); iqr <- q[4] - q[2]
  data.frame(sample = libs, mean_abs_resid = a, mean_abs_z = z,
             tukey_keep = a >= q[2] - 1.5 * iqr & a <= q[4] + 1.5 * iqr,
             n_genes_resid = nrow(X), n_genes_z = sum(stats::complete.cases(Z)),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Linear model cumulative score ~ PRS (t-test of the slope) and Pearson r.
prs_association <- function(score, prs) {
  ok <- is.finite(score) & is.finite(prs)
  f <- summary(stats::lm(score[ok] ~ prs[ok]))$coefficients
  c(n = sum(ok), pearson_r = stats::cor(score[ok], prs[ok]),
    slope = f[2, 1], p_value = f[2, 4])
}

#' Null distribution of the Pearson correlation under random assignment of the
#' cumulative scores to donors (seeded).
permutation_null <- function(score, prs, n_perm = 1000, seed = 42) {
  set.seed(seed)
  n <- length(score)
  vapply(seq_len(n_perm), function(i) stats::cor(score[sample.int(n, n)], prs), numeric(1))
}
