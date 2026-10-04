## ===========================================================================
## dapa_v2_functions.R
##
## Differential APA, version used for the RBP knockdown experiments.
##
## Differences to R/dapa_functions.R (iN cohort, postmortem DLPFC):
##   * R is computed from length-normalised read densities of the proximal and
##     distal UTR segment, R = dens_distal / (dens_distal + dens_proximal);
##     libraries without reads on either segment are missing (NA)
##   * R = 0 or 1 is squeezed into (0,1) (Smithson & Verkuilen 2006) instead of
##     being dropped
##   * covariates that do not vary among a gene's usable libraries are dropped
##   * no GLM fallback; fits that do not converge are not reported
##   * expression filter: mean size-factor-normalised count (proximal + distal)
##   * effect: condition coefficient (comparison vs reference, logit scale);
##     exp(beta) is the fold change in distal/proximal density ratio;
##     deltaRatio_adj = covariate-adjusted fitted R, comparison - reference
##
## Functions
##   build_apa_ratios_v2()  R, PDU and expression from segment counts
##   diff_apa_v2()          per-gene beta regression, Wald and LR test
##   finalize_dapa_v2()     expression filter, BH and Storey q-values
## ===========================================================================

suppressPackageStartupMessages({ library(betareg); library(qvalue) })

#' @param cts    segment counts, features x libraries
#' @param ann    feature annotation with featureID, geneId, class (proximal/distal), Length
#' @param size_factors named numeric vector (libraries); used for the expression filter
build_apa_ratios_v2 <- function(cts, ann, size_factors, min_len = 100) {
  ann <- ann[match(rownames(cts), ann$featureID), ]
  stopifnot(!anyNA(ann$featureID), all(names(size_factors) %in% colnames(cts)))
  sf <- size_factors[colnames(cts)]
  iP <- which(ann$class == "proximal"); iD <- which(ann$class == "distal")
  genes <- intersect(ann$geneId[iP], ann$geneId[iD])
  iP <- iP[match(genes, ann$geneId[iP])]; iD <- iD[match(genes, ann$geneId[iD])]
  ok <- ann$Length[iP] >= min_len & ann$Length[iD] >= min_len
  genes <- genes[ok]; iP <- iP[ok]; iD <- iD[ok]
  cP <- cts[iP, , drop = FALSE]; cD <- cts[iD, , drop = FALSE]
  dimnames(cP) <- dimnames(cD) <- list(genes, colnames(cts))
  nP <- t(t(cP) / sf); nD <- t(t(cD) / sf)
  dP <- nP / ann$Length[iP]; dD <- nD / ann$Length[iD]
  R <- dD / (dD + dP); PDU <- dD / dP
  R[(cP + cD) == 0 | !is.finite(R)] <- NA
  PDU[(cP + cD) == 0 | !is.finite(PDU) | PDU >= 1] <- NA
  list(R = R, PDU = PDU, expr = nP + nD)
}

.squeeze01 <- function(y) { n <- sum(!is.na(y)); if (n < 2) y else (y * (n - 1) + 0.5) / n }

.fit_gene_v2 <- function(y, si, covariates, groups, min_per_group, lr_test = TRUE) {
  keep <- !is.na(y)
  if (sum(keep) < 2 * min_per_group) return(NULL)
  d <- si[keep, , drop = FALSE]; d$ratio <- as.numeric(y[keep])
  d$condition <- factor(as.character(d$condition), levels = groups)
  n1 <- sum(d$condition == groups[1]); n2 <- sum(d$condition == groups[2])
  if (n1 < min_per_group || n2 < min_per_group) return(NULL)
  d <- droplevels(d)
  cov <- covariates[vapply(covariates, function(v) {
    x <- d[[v]]; if (is.factor(x) || is.character(x)) nlevels(factor(x)) > 1 else isTRUE(stats::sd(x) > 0)
  }, logical(1))]
  d$ratio <- .squeeze01(d$ratio)
  fF <- stats::reformulate(c("condition", cov), response = "ratio")
  fR <- stats::reformulate(if (length(cov)) cov else "1", response = "ratio")
  q <- function(e) tryCatch(suppressWarnings(e), error = function(e) NULL)
  fit  <- q(betareg::betareg(fF, data = d, link = "logit", type = "ML", link.phi = "identity"))
  fitR <- if (lr_test) q(betareg::betareg(fR, data = d, link = "logit", type = "ML", link.phi = "identity")) else NULL
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  if (lr_test && (is.null(fitR) || !isTRUE(fitR$converged))) return(NULL)
  cf <- summary(fit)$coefficients$mean; term <- paste0("condition", groups[2])
  if (!term %in% rownames(cf) || !is.finite(cf[term, 4])) return(NULL)
  lr <- if (lr_test) 2 * (as.numeric(stats::logLik(fit)) - as.numeric(stats::logLik(fitR))) else NA_real_
  ## covariate-adjusted fitted R per group (covariates at their most frequent level / mean)
  nd <- data.frame(condition = factor(groups, levels = groups))
  for (v in cov) { x <- d[[v]]; nd[[v]] <- if (is.factor(x)) factor(names(which.max(table(x))), levels = levels(x)) else mean(x) }
  pr <- tryCatch(suppressWarnings(as.numeric(stats::predict(fit, newdata = nd, type = "response"))),
                 error = function(e) rep(NA_real_, 2))
  list(beta = cf[term, 1], se = cf[term, 2], p_wald = cf[term, 4],
       fitted_ref = pr[1], fitted_comp = pr[2], deltaRatio_adj = pr[2] - pr[1],
       p_lrt = if (is.finite(lr) && lr >= 0) stats::pchisq(lr, 1, lower.tail = FALSE) else NA,
       phi = as.numeric(fit$coefficients$precision)[1], n1 = n1, n2 = n2, n = nrow(d),
       cov = paste(cov, collapse = ","),
       resid = setNames(if (lr_test) d$ratio - as.numeric(stats::fitted(fitR)) else rep(NA_real_, nrow(d)), rownames(d)))
}

#' @param R           genes x libraries, columns in the order of rownames(sample_info)
#' @param sample_info data.frame with `condition` and the covariates
#' @param groups      c(reference, comparison)
#' @param lr_test     also fit the model without condition (likelihood-ratio test, residual-based delta);
#'                    FALSE: Wald test only
diff_apa_v2 <- function(R, sample_info, covariates, groups, min_per_group = 3, cores = 1, lr_test = TRUE) {
  stopifnot(identical(colnames(R), rownames(sample_info)), all(sample_info$condition %in% groups))
  for (v in covariates) if (is.character(sample_info[[v]])) sample_info[[v]] <- factor(sample_info[[v]])
  message(sprintf("diff_apa_v2: %s vs %s | %d genes | %d / %d libraries", groups[2], groups[1], nrow(R),
                  sum(sample_info$condition == groups[1]), sum(sample_info$condition == groups[2])))
  fits <- parallel::mclapply(seq_len(nrow(R)), function(i)
    .fit_gene_v2(R[i, ], sample_info, covariates, groups, min_per_group, lr_test), mc.cores = cores)
  ok <- !vapply(fits, is.null, logical(1))
  get <- function(f) { v <- rep(NA_real_, length(fits)); v[ok] <- vapply(fits[ok], function(x) as.numeric(x[[f]]), numeric(1)); v }
  res <- data.frame(ID = rownames(R), beta = get("beta"), se = get("se"), p_wald = get("p_wald"),
                    p_lrt = get("p_lrt"), phi = get("phi"), fitted_ref = get("fitted_ref"),
                    fitted_comp = get("fitted_comp"), deltaRatio_adj = get("deltaRatio_adj"),
                    n_ref = get("n1"), n_comp = get("n2"),
                    n_used = get("n"), stringsAsFactors = FALSE)
  res$covUsed <- NA_character_; res$covUsed[ok] <- vapply(fits[ok], `[[`, character(1), "cov")
  resid <- matrix(NA_real_, nrow(R), ncol(R), dimnames = dimnames(R))
  for (i in which(ok)) resid[i, names(fits[[i]]$resid)] <- fits[[i]]$resid
  attr(res, "residuals") <- resid; attr(res, "groups") <- groups
  message(sprintf("  %d genes testable", sum(ok)))
  res
}

#' Expression filter (mean normalised proximal + distal count >= min_expr), BH and
#' Storey q-values for the Wald and LR tests, residual-based delta.
finalize_dapa_v2 <- function(res, ap, sample_info, min_expr = 10) {
  groups <- attr(res, "groups"); resid <- attr(res, "residuals")
  res$meanExpr <- rowMeans(ap$expr[res$ID, rownames(sample_info), drop = FALSE])
  keep <- res$meanExpr >= min_expr & !is.na(res$p_wald)
  out <- res[keep, ]
  g1 <- rownames(sample_info)[sample_info$condition == groups[1]]
  g2 <- rownames(sample_info)[sample_info$condition == groups[2]]
  out$deltaResid <- rowMeans(resid[out$ID, g2, drop = FALSE], na.rm = TRUE) -
                    rowMeans(resid[out$ID, g1, drop = FALSE], na.rm = TRUE)
  for (t in c("wald", "lrt")) {
    p <- out[[paste0("p_", t)]]; ok <- !is.na(p)
    if (!any(ok)) { out[[paste0("padj_", t)]] <- NA_real_; out[[paste0("qvalue_", t)]] <- NA_real_; next }
    out[[paste0("padj_", t)]] <- NA_real_; out[[paste0("qvalue_", t)]] <- NA_real_
    out[[paste0("padj_", t)]][ok] <- stats::p.adjust(p[ok], "BH")
    out[[paste0("qvalue_", t)]][ok] <- qvalue::qvalue(p[ok])$qvalues
  }
  out <- out[order(out$p_wald), ]
  message(sprintf("  %d genes pass mean expression >= %g | q <= 0.05 (Wald): %d", nrow(out), min_expr,
                  sum(out$qvalue_wald <= 0.05)))
  out
}
