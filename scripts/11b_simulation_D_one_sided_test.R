################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 11b, scenario D: power of the one-sided likelihood-ratio rule
##
##   The specification test is one-sided: with s the rotated residual variance
##   and theta the genetic variance,
##     D+ = 2 [ sup_{s >= 1, theta} l_R(s, theta) - sup_theta l_R(1, theta) ],
##   rejecting when D+ > 2.71, the 0.90 quantile of the chi-squared
##   distribution on one degree of freedom, which is the 0.95 quantile of the
##   50:50 mixture of chi-squared distributions on zero and one degrees of
##   freedom (Self and Liang 1987). D+ is zero when the free estimate of s is
##   at or below one. Script 11 measured detection with a threshold on the
##   residual estimate at the 95th percentile of the control; this script
##   applies D+ to the additive model in every scenario D dataset, each
##   regenerated from its own seed as in script 11, and reports the rejection
##   rate with its Monte Carlo standard error beside the threshold rule.
##
##   The additive model is fitted through the spectral form of the residual
##   likelihood, as in 10g. Check: the free residual reproduces the stored value
##   of script 11 for the same dataset.
##
##   Output: 11b_simD_one_sided.RData (read by F1 for Table 4)
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

suppressMessages({library(lme4); library(parallel)})
source("scripts/00_config.R")   # NREP_D, SEED_D, NCORE_SIM
OUTD  <- OUT
P7 <- new.env(); load(file.path(OUTD, "07_sim_primitives.RData"), envir = P7)
D1 <- new.env(); load(file.path(OUTD, "01_data.RData"),           envir = D1)
DD <- new.env(); load(file.path(OUTD, "11_simD.RData"),           envir = DD)

A <- D1$A; Q <- D1$Q; n <- nrow(A)
simulate_y <- P7$simulate_y; stage1_fast <- P7$stage1_fast; MARKERS <- P7$MARKERS

## make_transforms as in script 00
make_transforms <- function(Omega) {
  sv <- svd(Omega)
  if (min(sv$d) <= 0 || max(sv$d) / min(sv$d) > 1e10) stop("Omega not safely positive definite")
  list(rotated = sv$u %*% diag(1 / sqrt(sv$d)) %*% t(sv$u))
}

P <- list(pheno = D1$pheno, id_idx = P7$ID_IDX, yr_idx = P7$YR_IDX, yr_rep = P7$YR_REP,
          id_yr = P7$ID_YR, pars = P7$PARS, rootA = P7$ROOT_A, rootE = P7$ROOT_E,
          rootD = P7$ROOT_D, Gsub = P7$G_SUB, Q = Q, A = A, n = n, L = P7$L_mat)

NREP <- NREP_D                                          # from script 00
TRUTH  <- DD$TRUTH                                      # as script 11 stored it
CRIT   <- 2.71                                          # qchisq(0.90, 1) = 2.7055, rounded

## The residual likelihood of y = X b + u + e, Var = vg K + ve I, in the
## spectral form of 10g, for the additive model (rotated intercept and
## principal components, no marker). The free fit starts from two points and
## keeps the better; at ve = 1 only vg is optimised.
fit_one_sided <- function(s1) {
  Tm <- make_transforms(s1$Omega)$rotated
  y  <- as.numeric(Tm %*% s1$mu)
  X  <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q)
  K  <- Tm %*% A %*% t(Tm); K <- (K + t(K)) / 2
  S  <- diag(n) - X %*% solve(crossprod(X), t(X))
  ee <- eigen(S %*% K %*% S, symmetric = TRUE)
  keep <- seq_len(n - ncol(X))
  d   <- pmax(ee$values[keep], 0); eta <- as.numeric(crossprod(ee$vectors[, keep], y))
  dev <- function(vg, ve) sum(log(vg * d + ve)) + sum(eta^2 / (vg * d + ve))
  fr  <- lapply(list(c(log(60), 0), c(log(20), log(2))), function(st)
    optim(st, function(p) dev(exp(p[1]), exp(p[2])), method = "BFGS",
          control = list(reltol = 1e-12, maxit = 1000)))
  free <- fr[[which.min(vapply(fr, `[[`, numeric(1), "value"))]]
  one  <- optimize(function(lv) dev(exp(lv), 1), interval = log(c(0.01, 5000)), tol = 1e-10)
  ve   <- exp(free$par[2])
  cost <- one$objective - free$value
  c(vg_free = exp(free$par[1]), ve_free = ve, vg_one = exp(one$minimum),
    cost_of_one = cost, D_plus = if (ve > 1) max(cost, 0) else 0)
}

rep_one <- function(r) {
  out <- try({
    do.call(rbind, lapply(seq_len(nrow(TRUTH)), function(i) {
      src <- TRUTH$source[i]; s2 <- TRUTH$s2[i]
      set.seed(SEED_D + r * 10 + TRUTH$k[i])           # the seed of script 11 for this dataset
      b  <- setNames(rep(0, length(MARKERS)), MARKERS)
      s1 <- stage1_fast(simulate_y(b,
                                   s2_epi = if (src == "epi") s2 else 0,
                                   s2_dom = if (src == "het") s2 else 0,
                                   P = P), P = P)
      data.frame(rep = r, source = src, s2_true = s2, t(fit_one_sided(s1)))
    }))
  }, silent = TRUE)
  if (inherits(out, "try-error")) as.character(out) else out
}

ncore <- max(1, min(detectCores() - 1, NCORE_SIM))
cat("scenario D, one-sided rule:", NREP, "replicates x", nrow(TRUTH), "conditions,", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, suppressMessages(library(lme4))))
clusterExport(cl, c("P", "TRUTH", "MARKERS", "n", "Q", "A", "SEED_D", "simulate_y",
                    "stage1_fast", "make_transforms", "fit_one_sided"))
res <- parLapply(cl, seq_len(NREP), rep_one)
stopCluster(cl)
bad <- vapply(res, is.character, logical(1))
cat("replicates failed:", sum(bad), "\n"); if (any(bad)) cat(res[[which(bad)[1]]], "\n")
R1 <- do.call(rbind, res[!bad])

## ---------- check: the same datasets as script 11 ----------
m <- merge(R1, DD$D[, c("rep", "source", "s2_true", "Ve_additive")], by = c("rep", "source", "s2_true"))
dif <- abs(m$ve_free - m$Ve_additive)
away <- m$Ve_additive > 0.01
cat(sprintf("\nfree residual against script 11's stored value, %d datasets: max |difference| %.2e (residual above 0.01, n = %d), %.2e over all\n",
            nrow(m), max(dif[away]), sum(away), max(dif)))

## ---------- rejection rates ----------
q95 <- as.numeric(quantile(DD$D$Ve_additive[DD$D$source == "none"], 0.95))
pct <- function(x) 100 * mean(x)
mcse <- function(x) 100 * sqrt(mean(x) * (1 - mean(x)) / length(x))
lab <- c(none = "none", epi = "epistatic (0.56 with A)", het = "heterozygosity (0.20 with A)")
one_sided_tab <- do.call(rbind, lapply(split(m, list(m$source, m$s2_true), drop = TRUE), function(x) {
  rej <- x$D_plus > CRIT; thr <- x$Ve_additive > q95
  data.frame(structure = lab[x$source[1]], withheld = x$s2_true[1], n = nrow(x),
             one_sided_pct = pct(rej), one_sided_mcse = mcse(rej),
             threshold_pct = pct(thr), threshold_mcse = mcse(thr),
             two_sided_3.84_pct = pct(x$cost_of_one > 3.84),      # qchisq(0.95, 1) = 3.8415, rounded
             residual_above_one_pct = pct(x$ve_free > 1))
}))
one_sided_tab <- one_sided_tab[order(match(names(lab), names(lab))[match(one_sided_tab$structure, lab)],
                                     one_sided_tab$withheld), ]

options(width = 200)
cat(sprintf("\n=== 11b: one-sided rule D+ > %.4f against the threshold rule (residual above %.3f, the control's 95th percentile) ===\n\n", CRIT, q95))
print(format(one_sided_tab, digits = 4), row.names = FALSE)
cat(sprintf("\none_sided: D+ > %g, D+ = 0 when the free residual is at or below one; two_sided: cost_of_one > 3.84\n", CRIT))
cat("threshold: script 11's rule, reproduced from its stored residuals\n")
save(R1, one_sided_tab, q95, CRIT, file = file.path(OUTD, "11b_simD_one_sided.RData"))
