## ===========================================================================
## dapa_functions.R
##
## Differential alternative polyadenylation (dAPA) from read counts on the
## proximal and distal 3'UTR regions of each gene (reference built in
## 02_apa_pas_library).
##
##   apa_ratios()   distal usage ratio R and expression per gene and sample
##   diff_apa()     per-gene regression of R on condition + covariates
##   finalize_dapa() expression filter, q-values, residual-based delta
##
## MODEL
##   R  = d / (d + p), with d and p the size-factor-normalised read counts on
##        the distal and proximal region (size factors cancel in the ratio).
##   Per gene, samples with R in (0, 1) are modelled by beta regression (logit
##   link, ML, constant precision):  R ~ condition + covariates.
##   The condition effect is tested by its Wald z-test (p_wald, used in the
##   manuscript); a likelihood-ratio test against R ~ covariates is reported
##   alongside (p_lrt).
##   If only one of the two groups has >= min_per_group usable samples, the
##   gene is modelled by a quasibinomial GLM on all samples instead
##   (method = "quasibinomial"; t-test on the condition coefficient).
##   Fits that return an error or a warning are treated as failed.
##
## EFFECT SIZE
##   delta = mean(residual | reference) - mean(residual | comparison), with the
##   residuals of the covariate-only model (response scale for beta
##   regression, working scale for the GLM fallback). Positive delta = higher
##   distal 3'UTR usage in the reference group (e.g. controls).
##
## Dependencies: betareg, qvalue (CRAN / Bioconductor); parallel (base R)
## ===========================================================================

suppressPackageStartupMessages({
  library(betareg)
  library(qvalue)
})

## ------------------------------------------------------------ ratios
#' @param cts   counts, UTR regions x samples (rownames = region ID)
#' @param ann   region annotation with columns ID, geneId, class (proximal/distal)
#' @param size_factors named numeric vector (one per column of cts)
#' @return list(R, expr): genes x samples. R is 0 where a gene has no reads
#'         in a sample (such samples are excluded from the beta regression).
#'         expr = max(normalised proximal, normalised distal count).
apa_ratios <- function(cts, ann, size_factors) {
  cts <- as.matrix(cts)
  ann <- ann[match(rownames(cts), ann$ID), ]
  if (anyNA(ann$ID)) stop("count rows without annotation")
  if (!all(colnames(cts) %in% names(size_factors))) stop("size factors missing for some samples")
  norm <- sweep(cts, 2, size_factors[colnames(cts)], "/")

  genes <- unique(ann$geneId)
  iP <- match(paste(genes, "proximal"), paste(ann$geneId, ann$class))
  iD <- match(paste(genes, "distal"),   paste(ann$geneId, ann$class))
  ok <- !is.na(iP) & !is.na(iD)
  genes <- genes[ok]; P <- norm[iP[ok], , drop = FALSE]; D <- norm[iD[ok], , drop = FALSE]
  dimnames(P) <- dimnames(D) <- list(genes, colnames(cts))

  R <- D / (D + P); R[is.na(R)] <- 0
  list(R = R, expr = pmax(P, D))
}

## ------------------------------------------------------------ per-gene fit
.quiet_fit <- function(expr) tryCatch(expr, error = function(e) NULL, warning = function(w) NULL)

.varying <- function(d, covariates) covariates[vapply(covariates, function(v) {
  x <- d[[v]]
  if (is.factor(x) || is.character(x)) length(unique(x)) > 1 else stats::sd(x) > 0
}, logical(1))]

.fit_gene <- function(y, si, covariates, groups, min_per_group, ones_usable, fallback) {
  cond <- si$condition
  usable <- y > 0 & (if (ones_usable) y <= 1 else y < 1)
  n1 <- sum(usable & cond == groups[1]); n2 <- sum(usable & cond == groups[2])
  out <- list(method = NA_character_, beta = NA_real_, se = NA_real_, p_wald = NA_real_,
              p_lrt = NA_real_, n_ref = n1, n_comp = n2, resid = NULL)

  if (n1 >= min_per_group && n2 >= min_per_group) {
    keep <- y > 0 & y < 1
    d <- si[keep, , drop = FALSE]; d$ratio <- y[keep]
    cv <- .varying(d, covariates)
    fFull <- stats::reformulate(c("condition", cv), response = "ratio")
    fRed  <- stats::reformulate(if (length(cv)) cv else "1", response = "ratio")
    fit  <- .quiet_fit(betareg::betareg(fFull, data = d, link = "logit", type = "ML", link.phi = "identity"))
    fitR <- .quiet_fit(betareg::betareg(fRed,  data = d, link = "logit", type = "ML", link.phi = "identity"))
    if (is.null(fit) || is.null(fitR)) return(out)
    cf <- summary(fit)$coefficients$mean
    term <- paste0("condition", groups[2])
    out$method <- "betareg"
    out$beta <- cf[term, "Estimate"]; out$se <- cf[term, "Std. Error"]; out$p_wald <- cf[term, "Pr(>|z|)"]
    lr <- 2 * (as.numeric(stats::logLik(fit)) - as.numeric(stats::logLik(fitR)))
    out$p_lrt <- stats::pchisq(lr, df = 1, lower.tail = FALSE)
    out$resid <- setNames(d$ratio - stats::fitted(fitR), rownames(d))

  } else if (fallback == "quasibinomial" && (n1 >= min_per_group || n2 >= min_per_group)) {
    d <- si; d$ratio <- y
    cv <- .varying(d, covariates)
    fFull <- stats::reformulate(c("condition", cv), response = "ratio")
    fRed  <- stats::reformulate(if (length(cv)) cv else "1", response = "ratio")
    fit  <- .quiet_fit(stats::glm(fFull, data = d, family = stats::quasibinomial, maxit = 1000))
    fitR <- .quiet_fit(stats::glm(fRed,  data = d, family = stats::quasibinomial, maxit = 1000))
    if (is.null(fit) || is.null(fitR)) return(out)
    cf <- summary(fit)$coefficients
    term <- paste0("condition", groups[2])
    out$method <- "quasibinomial"
    out$beta <- cf[term, "Estimate"]; out$se <- cf[term, "Std. Error"]; out$p_wald <- cf[term, "Pr(>|t|)"]
    out$resid <- setNames(fitR$residuals, rownames(d))           # working residuals
  }
  out
}

## ------------------------------------------------------------ all genes
#' @param R          genes x samples ratio matrix (from apa_ratios)
#' @param sample_info data.frame, rownames = colnames(R), with `condition` and covariates
#' @param groups     c(reference, comparison); beta = comparison vs reference
#' @param ones_usable count libraries with R = 1 towards min_per_group (they are
#'                   excluded from the beta regression either way)
#' @param fallback   "quasibinomial": fit a quasibinomial GLM when only one group
#'                   reaches min_per_group; "none": leave such genes untested
diff_apa <- function(R, sample_info, covariates, groups, min_per_group = 10,
                     ones_usable = TRUE, fallback = c("quasibinomial", "none"), cores = 1) {
  fallback <- match.arg(fallback)
  R <- as.matrix(R)
  if (!identical(colnames(R), rownames(sample_info)))
    stop("colnames(R) must equal rownames(sample_info) (same samples, same order)")
  if (!all(sample_info$condition %in% groups)) stop("sample_info$condition outside `groups`")
  miss <- setdiff(covariates, names(sample_info))
  if (length(miss)) stop("covariates not in sample_info: ", paste(miss, collapse = ", "))
  si <- sample_info
  si$condition <- factor(si$condition, levels = groups)
  for (v in covariates) if (is.character(si[[v]])) si[[v]] <- factor(si[[v]])

  message(sprintf("diff_apa: %s vs %s | %d genes | %d / %d samples", groups[2], groups[1],
                  nrow(R), sum(si$condition == groups[1]), sum(si$condition == groups[2])))
  fits <- parallel::mclapply(seq_len(nrow(R)), function(i)
    .fit_gene(R[i, ], si, covariates, groups, min_per_group, ones_usable, fallback),
    mc.cores = cores)

  get <- function(f) vapply(fits, function(x) as.numeric(x[[f]]), numeric(1))
  res <- data.frame(ID = rownames(R),
                    method = vapply(fits, `[[`, "", "method"),
                    n_ref = get("n_ref"), n_comp = get("n_comp"),
                    beta = get("beta"), se = get("se"),
                    p_wald = get("p_wald"), p_lrt = get("p_lrt"),
                    stringsAsFactors = FALSE)
  resid <- matrix(NA_real_, nrow(R), ncol(R), dimnames = dimnames(R))
  for (i in seq_along(fits)) if (!is.null(fits[[i]]$resid))
    resid[i, names(fits[[i]]$resid)] <- fits[[i]]$resid
  attr(res, "residuals") <- resid
  attr(res, "groups") <- groups
  message(sprintf("  %d betareg, %d quasibinomial, %d not testable",
                  sum(res$method %in% "betareg"), sum(res$method %in% "quasibinomial"),
                  sum(is.na(res$method))))
  res
}

## ------------------------------------------------------------ filter + FDR
#' Keeps genes with a p-value and mean expression >= min_expr across the tested
#' samples, then computes q-values (Storey) and BH-adjusted p-values on the kept
#' genes, and the residual-based delta.
finalize_dapa <- function(res, expr, sample_info, min_expr = 50, p_col = "p_wald") {
  groups <- attr(res, "groups"); resid <- attr(res, "residuals")
  s  <- rownames(sample_info)
  g1 <- s[sample_info$condition == groups[1]]; g2 <- s[sample_info$condition == groups[2]]
  res$delta <- rowMeans(resid[res$ID, g1, drop = FALSE], na.rm = TRUE) -
               rowMeans(resid[res$ID, g2, drop = FALSE], na.rm = TRUE)
  res$mean_expr <- rowMeans(expr[res$ID, s, drop = FALSE])
  res <- res[!is.na(res[[p_col]]) & res$mean_expr >= min_expr, ]
  res$qvalue <- qvalue::qvalue(res[[p_col]])$qvalues
  res$padj_BH <- stats::p.adjust(res[[p_col]], "BH")
  message(sprintf("  %d genes pass mean expression >= %g | q <= 0.05: %d", nrow(res), min_expr,
                  sum(res$qvalue <= 0.05)))
  rownames(res) <- NULL
  res
}
