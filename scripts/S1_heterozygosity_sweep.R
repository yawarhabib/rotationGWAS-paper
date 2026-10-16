################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##          Script S1: sensitivity to the heterozygosity filter
##
##   Quality control removed the markers whose observed heterozygosity, as
##   reported by PLINK --hardy on the 160 genotyped individuals, exceeded 0.10.
##   This script asks whether that threshold matters, by rebuilding the
##   analysis from the unfiltered marker set at a range of thresholds and
##   comparing the relationship matrix, the variance components and the
##   genomic inflation factor at each.
##
##   The first stage does not use the genotypes, so Omega, W and the rotated
##   response are the same under every threshold; only the relationship matrix
##   and what is derived from it change. The threshold 0.10 must reproduce the
##   analysed marker set exactly, which is checked.
##
##   Output: S1_het_sweep.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(genio)
library(data.table)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "02_rotation.RData"))
t0 <- Sys.time()

S05 <- file.path(ROOT, "data/genotypes_160_lines_65905_markers_unfiltered")
HWE <- file.path(ROOT, "data/marker_heterozygosity_65905_markers.txt")

## ---------- 1. the unfiltered marker set ----------
pl  <- read_plink(S05, verbose = FALSE)
G05 <- t(pl$X)
rownames(G05) <- as.character(pl$fam$id)
colnames(G05) <- as.character(pl$bim$id)

## the same accessions in the same order as script 01
G05 <- G05[rownames(A), , drop = FALSE]

## the observed heterozygosity from the PLINK output, the statistic quality
## control used
hwe <- fread(HWE)
ohet <- setNames(hwe$`O(HET)`, hwe$SNP)[colnames(G05)]
stopifnot(!anyNA(ohet))

## ---------- 2. one analysis at one threshold ----------
## Centring makes every relationship matrix singular, so each receives the
## matrix of ones before rotation, K + cJ, through with_ones in script 00, as
## in script 04.
Xf <- cbind(ic_rot, Q_rot)
p_f <- ncol(Xf)

## The REML objective, with the determinants from determinant() and a direct
## inverse of H.
neg2 <- function(pars, K, y) {
  H  <- K * exp(pars[1]) + diag(n) * exp(pars[2])
  Hi <- try(solve(H), silent = TRUE)
  if (inherits(Hi, "try-error")) return(1e10)
  XtHiX <- crossprod(Xf, Hi %*% Xf)
  P <- Hi - Hi %*% Xf %*% solve(XtHiX) %*% crossprod(Xf, Hi)
  as.numeric(determinant(H, logarithm = TRUE)$modulus +
             determinant(XtHiX, logarithm = TRUE)$modulus + crossprod(y, P %*% y))
}

## Nelder-Mead can stop short of the optimum, so it is restarted from its own
## answer until the deviance stops improving.
fit_reml <- function(K, y, start = c(70, 1)) {
  par <- log(start); dev <- Inf
  for (i in 1:20) {
    o <- optim(par, neg2, K = K, y = y, method = "Nelder-Mead",
               control = list(maxit = 5000, reltol = 1e-11))
    if (dev - o$value < 1e-9) { par <- o$par; dev <- o$value; break }
    par <- o$par; dev <- o$value
  }
  list(Vg = exp(par[1]), Ve = exp(par[2]), dev = dev, conv = o$convergence)
}

analyse <- function(thr) {
  keep <- ohet <= thr
  Gt   <- G05[, keep, drop = FALSE]

  ## imputation at the marker mean, using only the analysed accessions
  mm <- colMeans(Gt, na.rm = TRUE)
  ok <- !is.na(mm)                       # a marker missing in every accession
  Gt <- Gt[, ok, drop = FALSE]; mm <- mm[ok]
  ix <- which(is.na(Gt), arr.ind = TRUE)
  Gt[ix] <- mm[ix[, "col"]]

  ## the mean minor allele frequency, since a low observed heterozygosity also
  ## selects for a low allele frequency
  pf  <- colMeans(Gt) / 2
  maf <- pmin(pf, 1 - pf)

  ## sommer documents the coding -1, 0, 1, as in script 01
  At <- sommer::A.mat(Gt - 1)
  ev <- eigen(At, symmetric = TRUE, only.values = TRUE)$values
  Kt <- with_ones(At, W)

  fit <- fit_reml(Kt, y_rot)

  ## genome scan under the fitted variance components, all markers at once
  H  <- Kt * fit$Vg + diag(n) * fit$Ve
  Hi <- solve(H)
  XtHiX <- crossprod(Xf, Hi %*% Xf)
  P  <- Hi - Hi %*% Xf %*% solve(XtHiX) %*% crossprod(Xf, Hi)
  Zr <- W %*% Gt
  PZ <- P %*% Zr
  zPz <- colSums(Zr * PZ)
  zPy <- as.numeric(crossprod(PZ, y_rot))
  yPy <- as.numeric(crossprod(y_rot, P %*% y_rot))

  ## a marker with no variance left after rotation carries no test
  good <- zPz > 1e-8 * max(zPz)
  b  <- zPy[good] / zPz[good]
  v  <- n - p_f - 1
  s2 <- (yPy - b * zPy[good]) / v
  t2 <- (b / sqrt(s2 / zPz[good]))^2
  ## The inflation factor as script 05 defines it: each t statistic on v degrees
  ## of freedom is turned into its p-value and that into a chi-squared on one
  ## degree of freedom.
  pv  <- pbeta(v / (v + t2), v / 2, 1 / 2)
  chi <- qchisq(pv, 1, lower.tail = FALSE)

  list(row = data.frame(
         Threshold   = thr,
         Markers     = ncol(Gt),
         Dropped     = sum(!keep),
         Untestable  = sum(!good),
         Mean_MAF    = round(mean(maf), 4),
         Mean_diag   = round(mean(diag(At)), 4),
         F_implied   = round(mean(diag(At)) - 1, 4),
         Eff_rank    = sum(ev > 1e-6 * max(ev)),
         Vg          = round(fit$Vg, 3),
         Ve          = round(fit$Ve, 4),
         Lambda      = round(median(chi) / qchisq(0.5, 1), 4),
         Converged   = fit$conv == 0),
       A = At, Vg = fit$Vg, Ve = fit$Ve)
}

## ---------- 3. the sweep ----------
## 0.10 is the threshold actually applied; the rest bracket it, and 1.00 is no
## filter at all.
thresholds <- c(1.00, 0.50, 0.30, 0.20, 0.10, 0.05, 0.02)
res <- lapply(thresholds, analyse)

sweep_tab <- do.call(rbind, lapply(res, `[[`, "row"))

## ---------- 4. the threshold 0.10 must reproduce the analysis ----------
i10 <- which(thresholds == 0.10)
A_rebuilt <- res[[i10]]$A
dev_A <- max(abs(A_rebuilt - A))

## At 0.10 the sweep fits the model of script 04 with another objective, so the
## variance components must agree. Read into an environment of its own, since
## script 04 saves objects whose names are used here.
e04 <- new.env()
load(file.path(OUT, "04_variance.RData"), envir = e04)
vc_04 <- as.numeric(e04$fits$additive$vc[1:2])
rm(e04)
vc_S1  <- c(res[[i10]]$Vg, res[[i10]]$Ve)
dev_vc <- max(abs(vc_S1 - vc_04) / vc_04)

checks <- data.frame(
  Check = c("threshold 0.10 gives the analysed marker count",
            "threshold 0.10 reproduces the relationship matrix",
            "threshold 0.10 reproduces the variance components of script 04",
            "no filter leaves more markers than any filter",
            "all variance-component fits converged"),
  Result = c(sweep_tab$Markers[i10] == ncol(G),
             dev_A < 1e-10,
             dev_vc < 1e-3,
             sweep_tab$Markers[1] == max(sweep_tab$Markers),
             all(sweep_tab$Converged)))

cat("\n=== S1 heterozygosity-filter sweep ===\n")
cat("\nthreshold 0.10 is the filter actually applied\n\n")
print(sweep_tab, row.names = FALSE)
cat("\nmax |rebuilt A - pipeline A| at threshold 0.10:", signif(dev_A, 3), "\n")
cat("variance components at 0.10:", signif(vc_S1, 6), "| script 04:", signif(vc_04, 6),
    "| largest relative difference", signif(dev_vc, 3), "\n")
cat("\nchecks\n")
print(checks, row.names = FALSE)

## saved before the checks are enforced, so a failed check does not lose the run
save(sweep_tab, checks, thresholds, dev_A, dev_vc, vc_S1, vc_04,
     file = file.path(OUT, "S1_het_sweep.RData"))

cat("\nsaved S1_het_sweep.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

if (!all(checks$Result)) stop("sweep checks failed")
