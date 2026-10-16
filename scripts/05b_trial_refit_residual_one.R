################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 05b: the refit with the rotated residual variance fixed at one,
##              on the observed trial
##
##   Markers close to the genome-wide threshold are fitted again with the marker
##   among the fixed effects, the genetic variance re-estimated and the rotated
##   residual variance fixed at one. The refit is applied to every marker of the
##   comparison set of script 06, for which the single-stage analysis exists,
##   and set against that analysis and against the rotated P3D scan of script
##   05; the markers within WINDOW units of -log10 p of the threshold are
##   reported as a subset.
##
##   It also gives the genetic variance explained by each marker, the drop in
##   the genetic variance when the marker is fitted, as a percentage of the
##   genetic variance without a marker; both fits have the residual at one.
##
##   Test: Wald t on n - p degrees of freedom with the standard error at the
##   fitted components; the single-stage reference of 06 uses the Kenward-Roger
##   F test.
##
##   Output: 05b_trial_refit_residual_one.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
load(file.path(OUT, "01_data.RData"))
e02 <- new.env(); load(file.path(OUT, "02_rotation.RData"), envir = e02)
e05 <- new.env(); load(file.path(OUT, "05_scan.RData"),     envir = e05)
e06 <- new.env(); load(file.path(OUT, "06_comparison.RData"), envir = e06)
t0 <- Sys.time()

WINDOW <- 1.5                                          # units of -log10 p on either side of the threshold
scan   <- e05$scans$rotated$res
thr    <- -log10(ALPHA / sum(!is.na(scan$neglog10p)))
near   <- scan$SNP[!is.na(scan$neglog10p) & abs(scan$neglog10p - thr) <= WINDOW]
sel    <- e06$ss$SNP[!is.na(e06$ss$neglog10p)]         # the comparison set, where the single-stage analysis exists
cat(sprintf("threshold %.3f; %d markers in the comparison set; %d within %g unit(s) of the threshold on the rotated scan (%d above, %d below)\n",
            thr, length(sel), length(near), WINDOW, sum(scan$neglog10p[match(near, scan$SNP)] > thr), sum(scan$neglog10p[match(near, scan$SNP)] <= thr)))
stopifnot(all(near %in% sel))                           # every marker near the threshold is in the comparison set

## ---------- the rotated inputs, as in scripts 02 and 05 ----------
W  <- e02$Tlist$rotated
y  <- as.numeric(W %*% e02$blues$BLUE)
X0 <- cbind(as.numeric(W %*% rep(1, n)), W %*% Q)
K  <- W %*% A %*% t(W); K <- (K + t(K)) / 2
DF <- n - ncol(X0) - 1

## ---------- one marker: marker as fixed effect, residual at one ----------
## The spectral form: S = I - X (X'X)^-1 X' and S K S = U diag(d) U'; minus
## twice the residual log-likelihood is sum log(vg d + ve) + sum eta^2 /
## (vg d + ve) over the n - p non-null directions. With ve = 1 only vg is
## optimized.
refit_one <- function(z) {
  X  <- cbind(X0, as.numeric(W %*% z))
  S  <- diag(n) - X %*% solve(crossprod(X), t(X))
  ee <- eigen(S %*% K %*% S, symmetric = TRUE)
  keep <- seq_len(n - ncol(X))
  d <- pmax(ee$values[keep], 0); eta <- as.numeric(crossprod(ee$vectors[, keep], y))
  dev <- function(lv) sum(log(exp(lv) * d + 1)) + sum(eta^2 / (exp(lv) * d + 1))
  o  <- optimize(dev, interval = log(c(1e-4, 1e6)), tol = 1e-10)
  stopifnot(o$minimum > log(1e-4) + 1e-3, o$minimum < log(1e6) - 1e-3)
  vg <- exp(o$minimum)
  Vi <- chol2inv(chol(vg * K + diag(n))); XtViX <- crossprod(X, Vi %*% X); C <- solve(XtViX)
  b  <- as.numeric(C %*% crossprod(X, Vi %*% y))[ncol(X)]; s <- sqrt(C[ncol(X), ncol(X)])
  c(vg = vg, beta = b, se = s, t = b / s, neglog10p = -log10(2 * pt(-abs(b / s), DF)))
}

refit <- do.call(rbind, lapply(sel, function(s) data.frame(SNP = s, t(refit_one(G[, s])), row.names = NULL)))

## ---------- the genetic variance explained by a marker ----------
## The same fit without a marker gives vg0, the genetic variance with no marker
## fitted. The variance explained by a marker is the drop from vg0 to the
## genetic variance of its own fit, as a percentage of vg0. Since A is the same
## in both fits, this equals the percentage drop in the average semivariance of
## the genetic part.
vg_at_one <- function(X) {
  S  <- diag(n) - X %*% solve(crossprod(X), t(X))
  ee <- eigen(S %*% K %*% S, symmetric = TRUE)
  keep <- seq_len(n - ncol(X))
  d <- pmax(ee$values[keep], 0); eta <- as.numeric(crossprod(ee$vectors[, keep], y))
  dev <- function(lv) sum(log(exp(lv) * d + 1)) + sum(eta^2 / (exp(lv) * d + 1))
  o  <- optimize(dev, interval = log(c(1e-4, 1e6)), tol = 1e-10)
  stopifnot(o$minimum > log(1e-4) + 1e-3, o$minimum < log(1e6) - 1e-3)
  exp(o$minimum)
}
vg0 <- vg_at_one(X0)

## ---------- against the scan and the single-stage analysis ----------
tab <- data.frame(SNP = sel,
                  scan_beta = scan$Beta[match(sel, scan$SNP)], scan_se = scan$SE[match(sel, scan$SNP)], scan_lp = scan$neglog10p[match(sel, scan$SNP)],
                  one_beta = refit$beta, one_se = refit$se, one_lp = refit$neglog10p, one_vg = refit$vg,
                  ss_beta = e06$ss$beta[match(sel, e06$ss$SNP)], ss_se = e06$ss$se[match(sel, e06$ss$SNP)],
                  ss_df = e06$ss$df[match(sel, e06$ss$SNP)], ss_lp = e06$ss$neglog10p[match(sel, e06$ss$SNP)])
tab$sig_scan <- tab$scan_lp > thr; tab$sig_one <- tab$one_lp > thr; tab$sig_ss <- tab$ss_lp > thr
tab$near <- tab$SNP %in% near
tab$one_explained <- 100 * (1 - tab$one_vg / vg0)   # percentage of vg0 the marker explains
tab <- tab[order(-tab$scan_lp), ]

summarise <- function(tab, label) data.frame(
  markers = label, n = nrow(tab), threshold = thr,
  significant_scan = sum(tab$sig_scan), significant_one = sum(tab$sig_one), significant_ss = sum(tab$sig_ss),
  decisions_one_vs_ss_agree = sum(tab$sig_one == tab$sig_ss), decisions_scan_vs_ss_agree = sum(tab$sig_scan == tab$sig_ss),
  max_abs_dbeta_one_vs_ss = max(abs(tab$one_beta - tab$ss_beta)), max_abs_dbeta_scan_vs_ss = max(abs(tab$scan_beta - tab$ss_beta)),
  se_ratio_one_vs_ss_min = min(tab$one_se / tab$ss_se), se_ratio_one_vs_ss_max = max(tab$one_se / tab$ss_se),
  se_ratio_scan_vs_ss_min = min(tab$scan_se / tab$ss_se), se_ratio_scan_vs_ss_max = max(tab$scan_se / tab$ss_se),
  max_abs_dlp_one_vs_ss = max(abs(tab$one_lp - tab$ss_lp)), max_abs_dlp_scan_vs_ss = max(abs(tab$scan_lp - tab$ss_lp)),
  vg_one_min = min(tab$one_vg), vg_one_max = max(tab$one_vg))
summary_tab <- rbind(summarise(tab, "comparison set"), summarise(tab[tab$near, ], sprintf("within %g units of the threshold", WINDOW)))

options(width = 220)
cat(sprintf("\n=== the markers within %g units of the threshold: rotated P3D scan, refit with the residual at one, single-stage (06) ===\n", WINDOW))
print(format(tab[tab$near, names(tab) != "near"], digits = 4), row.names = FALSE)
cat("\n=== summary, all 227 markers of the comparison set and the subset near the threshold ===\n")
print(format(summary_tab, digits = 4), row.names = FALSE)
w  <- function(x) 100 * max(abs(c(x$se_ratio_one_vs_ss_min,  x$se_ratio_one_vs_ss_max)  - 1))
w2 <- function(x) 100 * max(abs(c(x$se_ratio_scan_vs_ss_min, x$se_ratio_scan_vs_ss_max) - 1))
a <- summary_tab[1, ]; b <- summary_tab[2, ]
cat(sprintf("\nIn words: over the %d markers of the comparison set the refit agrees with the single-stage decision at %d, with estimates within %.3f,\n",
            a$n, a$decisions_one_vs_ss_agree, a$max_abs_dbeta_one_vs_ss))
cat(sprintf("standard errors within %.2f%% and -log10 p within %.3f of the single-stage values (the scan: within %.2f, %.1f%% and %.2f);\n",
            w(a), a$max_abs_dlp_one_vs_ss, a$max_abs_dbeta_scan_vs_ss, w2(a), a$max_abs_dlp_scan_vs_ss))
cat(sprintf("the %d markers within %g units of the threshold: %d significant on the scan, %d after the refit, %d single-stage; the refit and the single-stage analysis agree on %d of them.\n",
            b$n, WINDOW, b$significant_scan, b$significant_one, b$significant_ss, b$decisions_one_vs_ss_agree))
sig <- tab[tab$sig_scan, ]
cat(sprintf("\nGenetic variance with the residual at one and no marker: %.2f. The %d markers significant on the rotated scan explain %.1f%% to %.1f%% of it; the strongest, %s, explains %.1f%%.\n",
            vg0, nrow(sig), min(sig$one_explained), max(sig$one_explained), sig$SNP[1], sig$one_explained[1]))
save(tab, summary_tab, sel, near, thr, WINDOW, DF, vg0, file = file.path(OUT, "05b_trial_refit_residual_one.RData"))
cat(sprintf("\nsaved 05b_trial_refit_residual_one | %.2f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
