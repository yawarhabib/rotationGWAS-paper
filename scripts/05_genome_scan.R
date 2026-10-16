################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##                    Script 05: genome scan under each scheme
##
##   Single-marker scan of every marker under the three weighting schemes with
##   the additive kernel, which no extension fitted in 04 improves on clearly.
##
##   The variance components are estimated once on the model without a marker
##   and held fixed across the scan: the ratio of genetic to residual variance
##   is fixed, carried in Hinv, and the residual scale is re-estimated at each
##   marker. Section 6 measures the cost of this approximation (P3D, population
##   parameters previously determined).
##
##   The scan also records the standard error of every marker under every
##   scheme; section 3 compares the schemes' standard errors across the genome.
##
##   Output: 05_scan.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(rrBLUP)
library(parallel)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "02_rotation.RData"))
load(file.path(OUT, "04_variance.RData"))
t0 <- Sys.time()

sym <- function(M) (M + t(M)) / 2

## ---------- 1. scan under one transformation ----------
## mixed.solve tolerates a singular kernel, so the kernel is passed as it is,
## and its residual is the quantity that 04 tests against one.
## The scan is vectorised: for a partitioned model the residual sum of squares
## after fitting the marker equals y'Py - b z'Py, with P removing the fixed
## effects under Hinv, so every marker's estimate, standard error and p-value
## come from one matrix product.
scan_arm <- function(Tm) {
  y  <- as.numeric(Tm %*% blues$BLUE)
  X0 <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q)
  Zt <- Tm %*% G
  K  <- sym(Tm %*% A %*% t(Tm))

  nf <- mixed.solve(y = y, X = X0, Z = diag(n), K = K, return.Hinv = TRUE)
  Hi <- nf$Hinv

  XtHiX <- crossprod(X0, Hi %*% X0)
  P     <- Hi - Hi %*% X0 %*% solve(XtHiX) %*% crossprod(X0, Hi)

  PZ  <- P %*% Zt
  zPz <- colSums(Zt * PZ)
  zPy <- as.numeric(crossprod(PZ, y))
  yPy <- as.numeric(crossprod(y, P %*% y))

  p  <- ncol(X0) + 1
  v2 <- n - p

  ## a marker with no variance left after the transformation carries no test
  ok <- zPz > 1e-8 * max(zPz)

  beta <- se <- lp <- rep(NA_real_, ncol(Zt))
  b  <- zPy[ok] / zPz[ok]
  s2 <- (yPy - b * zPy[ok]) / v2
  sb <- sqrt(s2 / zPz[ok])
  beta[ok] <- b
  se[ok]   <- sb
  lp[ok]   <- -log10(pbeta(v2 / (v2 + (b / sb)^2), v2 / 2, 1 / 2))

  list(res = data.frame(SNP = map$SNP, Chromosome = map$Chromosome,
                        Position = map$Position, Beta = beta, SE = se,
                        neglog10p = lp),
       Vg = nf$Vu, Ve = nf$Ve)
}

cat("scanning", ncol(G), "markers under three weighting schemes ...\n")
scans <- lapply(Tlist, scan_arm)

## ---------- 2. inflation factor and significance ----------
lambda_of <- function(lp) {
  pv <- 10^(-lp); pv <- pv[!is.na(pv) & pv > 0 & pv <= 1]
  qchisq(median(pv), 1, lower.tail = FALSE) / qchisq(0.5, 1)
}

## One threshold for all three schemes; that they test the same markers is
## one of the checks.
dropped  <- sapply(scans, function(s) sum(is.na(s$res$neglog10p)))
n_tested <- nrow(map) - dropped[["rotated"]]
bonf     <- -log10(ALPHA / n_tested)

scan_tab <- do.call(rbind, lapply(names(scans), function(nm) {
  s <- scans[[nm]]
  data.frame(Scheme = nm,
             Vg = round(s$Vg, 3), Ve = round(s$Ve, 4),
             dropped = dropped[[nm]],
             lambda_P3D = round(lambda_of(s$res$neglog10p), 4),
             n_significant = sum(s$res$neglog10p > bonf, na.rm = TRUE))
}))

## The rotated scheme fits the model of 04 with a different solver; the two
## residual estimates are compared in the checks.
ve_04  <- fits$additive$vc[2]
ve_gap <- abs(scans$rotated$Ve - ve_04)

## ---------- 3. standard errors across the genome ----------
## The rotated scheme is the reference. A ratio below one at every marker means
## a scheme understates its standard errors throughout; a ratio scattered about
## one means the error is marker-specific. The share below one separates the
## two.
se_r <- scans$rotated$res$SE
ok3  <- !is.na(se_r) & !is.na(scans$unweighted$res$SE) & !is.na(scans$diagonal$res$SE)

pfreq <- colMeans(G) / 2
maf   <- pmin(pfreq, 1 - pfreq)

se_ratio <- list(unweighted = scans$unweighted$res$SE[ok3] / se_r[ok3],
                 diagonal   = scans$diagonal$res$SE[ok3]   / se_r[ok3])

se_tab <- do.call(rbind, lapply(names(se_ratio), function(nm) {
  r <- se_ratio[[nm]]
  data.frame(Scheme        = nm,
             Mean          = round(mean(r), 5),
             Median        = round(median(r), 5),
             Q2.5          = round(unname(quantile(r, 0.025)), 5),
             Q97.5         = round(unname(quantile(r, 0.975)), 5),
             Pct_below_one = round(100 * mean(r < 1), 2),
             Cor_with_MAF  = round(cor(r, maf[ok3]), 4))
}))

## ---------- 4. do the schemes select the same markers? ----------
sig <- lapply(scans, function(s) s$res$SNP[which(s$res$neglog10p > bonf)])
overlap_tab <- data.frame(
  Comparison = c("unweighted vs rotated", "diagonal vs rotated", "all three"),
  Shared = c(length(intersect(sig$unweighted, sig$rotated)),
             length(intersect(sig$diagonal, sig$rotated)),
             length(Reduce(intersect, sig))),
  Union  = c(length(union(sig$unweighted, sig$rotated)),
             length(union(sig$diagonal, sig$rotated)),
             length(Reduce(union, sig))))

## ---------- 5. the test marker ----------
test_tab <- do.call(rbind, lapply(names(scans), function(nm) {
  r <- scans[[nm]]$res
  r <- r[r$SNP == TEST_SNP, ]
  stopifnot(nrow(r) == 1)
  data.frame(Scheme = nm, Beta = round(r$Beta, 4), SE = round(r$SE, 4),
             neglog10p = round(r$neglog10p, 4))
}))

## ---------- 6. cost of the P3D approximation ----------
## The variance components are re-estimated per marker on a random sample of
## markers, and the inflation factor with and without the re-estimation is
## compared. The REML objective is that of 04, with the kernel W(A + cJ)W' as
## built there, so H is positive definite for every residual; the determinants
## come from Cholesky factors.
set.seed(SEED_MISC + 2)
VE_BOUNDARY <- 1e-3

valid <- which(!is.na(scans$rotated$res$neglog10p))
samp  <- sort(sample(valid, min(NREP_P3D, length(valid))))

y_r  <- y_rot
Xn_r <- cbind(ic_rot, Q_rot)
Zs   <- (W %*% G)[, samp, drop = FALSE]
Ka   <- K_add                    # W(A + cJ)W', as built in 04

neg2_1k <- function(pars, X) {
  ve <- exp(pars[2])
  H <- Ka * exp(pars[1]) + diag(n) * ve
  val <- try({
    R     <- chol(H)
    Hi    <- chol2inv(R)
    XtHiX <- crossprod(X, Hi %*% X)
    Rx    <- chol(XtHiX)
    P     <- Hi - Hi %*% X %*% chol2inv(Rx) %*% crossprod(X, Hi)
    2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) +
      as.numeric(crossprod(y_r, P %*% y_r))
  }, silent = TRUE)
  if (inherits(val, "try-error") || !is.finite(val)) return(1e10)
  val
}

wald_given <- function(X, vc) {
  H  <- Ka * vc[1] + diag(n) * vc[2]
  Hi <- chol2inv(chol(H))
  p  <- ncol(X); v2 <- n - p
  Wi <- solve(crossprod(X, Hi %*% X))
  b  <- Wi %*% crossprod(X, Hi %*% matrix(y_r, n, 1))
  r  <- matrix(y_r, n, 1) - X %*% b
  s2 <- as.double(crossprod(r, Hi %*% r)) / v2
  cb <- s2 * Wi
  -log10(pbeta(v2 / (v2 + b[p]^2 / cb[p, p]), v2 / 2, 1 / 2))
}

vc_null <- c(scans$rotated$Vg, scans$rotated$Ve)

cat("measuring the P3D correction on", length(samp), "markers ...\n")
ncore <- max(1, min(detectCores() - 1, 6))
cl <- makeCluster(ncore)
clusterExport(cl, c("Xn_r", "Zs", "y_r", "n", "Ka", "vc_null",
                    "neg2_1k", "wald_given"), envir = environment())
## the number of cores changes the timing, not the result
p3d <- parLapply(cl, seq_len(ncol(Zs)), function(j) {
  X <- cbind(Xn_r, Zs[, j])
  a <- wald_given(X, vc_null)
  par <- log(vc_null); dev <- Inf; conv <- NA_integer_
  for (i in 1:3) {
    o <- optim(par, neg2_1k, X = X, method = "Nelder-Mead",
               control = list(maxit = 2000, reltol = 1e-10))
    improved <- dev - o$value; par <- o$par; dev <- o$value; conv <- o$convergence
    if (improved < 1e-9) break
  }
  vc <- exp(par)
  c(p3d = a, exact = wald_given(X, vc), ve = vc[2], conv = conv, dev = dev)
})
stopCluster(cl)

p3d_raw  <- do.call(rbind, p3d)
lp_p3d   <- p3d_raw[, "p3d"]
lp_exact <- p3d_raw[, "exact"]

## as in 04, a residual near zero is a boundary solution, not a failure
p3d_boundary <- p3d_raw[, "ve"] < VE_BOUNDARY
p3d_failed   <- p3d_raw[, "conv"] != 0 & !p3d_boundary

lam_p3d    <- lambda_of(lp_p3d)
lam_exact  <- lambda_of(lp_exact)
lam_gw_rot <- lambda_of(scans$rotated$res$neglog10p)

p3d_tab <- data.frame(
  Quantity = c("lambda, sample, P3D", "lambda, sample, re-estimated",
               "paired difference", "lambda, genome-wide, P3D",
               "lambda, genome-wide, corrected",
               "markers at the boundary", "fits failed"),
  Value = c(round(c(lam_p3d, lam_exact, lam_exact - lam_p3d,
                    lam_gw_rot, lam_gw_rot + (lam_exact - lam_p3d)), 4),
            sum(p3d_boundary), sum(p3d_failed)))

## ---------- checks ----------
checks <- data.frame(
  Check = c("all schemes tested the same markers",
            "test marker present under every scheme",
            "rotated residual agrees with 04",
            "no P3D fit diverged",
            "every non-boundary P3D fit converged",
            "inflation factors all finite"),
  Result = c(length(unique(dropped)) == 1,
             nrow(test_tab) == 3 && !anyNA(test_tab$neglog10p),
             ve_gap < 0.05,
             all(is.finite(p3d_raw[, "dev"])) && min(p3d_raw[, "dev"]) > -1e6,
             !any(p3d_failed),
             all(is.finite(scan_tab$lambda_P3D))))

## ---------- report ----------
cat("\n=== 05 genome scan ===\n")
cat("\nby weighting scheme (Bonferroni threshold", round(bonf, 3),
    "on", n_tested, "markers)\n")
print(scan_tab, row.names = FALSE)
cat("   rotated residual here", round(scans$rotated$Ve, 4),
    "against", round(ve_04, 4), "from 04, difference",
    signif(ve_gap, 3), "\n")
cat("\nstandard error relative to the rotated scheme, across the genome\n")
print(se_tab, row.names = FALSE)
## ---------- significant markers, with where they are ----------
## A count of markers is not a count of loci: in a panel with linkage
## disequilibrium neighbouring markers carry the same signal. The positions show
## which.
sig_tab <- do.call(rbind, lapply(names(scans), function(nm) {
  r <- scans[[nm]]$res
  r <- r[!is.na(r$neglog10p) & r$neglog10p > bonf, ]
  if (!nrow(r)) return(NULL)
  r$Chromosome <- map$Chromosome[match(r$SNP, map$SNP)]
  r$Position   <- map$Position[match(r$SNP, map$SNP)]
  r <- r[order(-r$neglog10p), ]
  data.frame(Scheme = nm, SNP = r$SNP, Chromosome = r$Chromosome,
             Position = r$Position, neglog10p = round(r$neglog10p, 4),
             row.names = NULL)
}))

locus_tab <- do.call(rbind, lapply(split(sig_tab, sig_tab$Scheme), function(x)
  data.frame(Scheme = x$Scheme[1], n_markers = nrow(x),
             chromosomes = paste(sort(unique(x$Chromosome)), collapse = ", "),
             first = min(x$Position), last = max(x$Position),
             span_bp = diff(range(x$Position)), row.names = NULL)))

cat("\nmarkers passing the Bonferroni threshold\n")
print(sig_tab, row.names = FALSE)
cat("\nthe extent they cover\n")
print(locus_tab, row.names = FALSE)
cat("   a count of markers is not a count of loci; the span shows which\n")
## Markers significant under one scheme and not another, with their distance
## from the threshold under every scheme.
only_in <- setdiff(sig_tab$SNP[sig_tab$Scheme == "unweighted"],
                   sig_tab$SNP[sig_tab$Scheme == "rotated"])
margin_tab <- if (!length(only_in)) NULL else
  do.call(rbind, lapply(only_in, function(s) {
    data.frame(SNP = s,
               Position = map$Position[match(s, map$SNP)],
               Scheme = names(scans),
               neglog10p = round(vapply(names(scans), function(nm)
                 scans[[nm]]$res$neglog10p[scans[[nm]]$res$SNP == s], numeric(1)), 4),
               margin = round(vapply(names(scans), function(nm)
                 scans[[nm]]$res$neglog10p[scans[[nm]]$res$SNP == s], numeric(1)) - bonf, 4),
               row.names = NULL)
  }))

if (!is.null(margin_tab)) {
  cat("\nmarkers significant under one scheme and not another\n")
  print(margin_tab, row.names = FALSE)
  cat("   margin is the distance from the threshold; positive is above it\n")
}
cat("   the strongest marker overall is", sig_tab$SNP[1],
    "at", format(sig_tab$Position[1], big.mark = ","), "on chromosome",
    sig_tab$Chromosome[1], "\n")

cat("\nmarker sets selected\n");    print(overlap_tab, row.names = FALSE)
cat("\ntest marker\n");             print(test_tab, row.names = FALSE)
cat("\nP3D correction\n");          print(p3d_tab, row.names = FALSE)
cat("\nchecks\n");                  print(checks, row.names = FALSE)

if (!all(checks$Result)) stop("genome-scan checks failed")

save(scans, scan_tab, se_tab, se_ratio, overlap_tab, test_tab, p3d_tab, checks,
     sig_tab, locus_tab, margin_tab,
     bonf, n_tested, dropped, sig, samp, lp_p3d, lp_exact, p3d_raw,
     p3d_boundary, p3d_failed, lam_gw_rot, ve_gap, maf, ok3,
     file = file.path(OUT, "05_scan.RData"))

cat("\nsaved 05_scan.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
