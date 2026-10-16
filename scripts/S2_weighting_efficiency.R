################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script S2: why the weighting schemes agree on this trial
##
##   Scenario B increases the correlation among the adjusted means and the
##   three weighting schemes still return the same precision. This script
##   computes the reason exactly, without simulation.
##
##   Two quantities account for it. First, the exact variance of a marker
##   coefficient under each weighting, evaluated under the true Omega, with
##   generalized least squares as the reference: ignoring the correlation,
##   which is the one thing the rotation does that inverse-variance weighting
##   cannot, costs little on this design. The correlation at the observed block
##   variance does not come from shared replicate effects: it is the same with
##   the replicate variance at zero (section 3b) and comes from estimating the
##   fixed year effects in an unbalanced design.
##
##   Second, the share of an adjusted mean's variance that weighting can act
##   on at all. An adjusted mean carries genetic variance and first-stage
##   estimation error, and only the second is what Omega describes; on this
##   trial the genetic part dominates, so the inefficiency of a weighting
##   scheme is diluted before it reaches the marker test.
##
##   Two versions of the cost are computed. The Omega-only sandwich (section 2)
##   says how much each weighting loses on the estimation error alone. The
##   version with the kernel (section 5) takes the true covariance s2_g A +
##   Omega against each scheme's working covariance and is what the marker
##   test loses; it reproduces the ratios scenarios B and C measure. The
##   measured values are read from the output of script 09.
##
##   Output: S2_efficiency.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(lme4)
library(MASS)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
t0 <- Sys.time()

## ---------- 1. Omega at any block variance, analytically ----------
## V = s2_rep Zr Zr' + s2_ge Zge Zge' + s2_e I at the calibrated components,
## then Omega = L Var(beta) L' with L the averaging contrast from 07. No data
## are simulated.
d  <- pheno; d$ysim <- rnorm(nrow(d))     # only to obtain the design matrices
m0 <- lmer(ysim ~ Id + Year + (1|Year:Replication) + (1|Id:Year),
           data = d, REML = TRUE)
X    <- as.matrix(getME(m0, "X"))
Zr   <- model.matrix(~ 0 + YR_REP)
Zge  <- model.matrix(~ 0 + ID_YR)
N    <- nrow(pheno)
ZrZr <- tcrossprod(Zr); ZgeZge <- tcrossprod(Zge)

omega_at <- function(s2r) {
  V  <- s2r * ZrZr + PARS$s2_ge * ZgeZge + PARS$s2_e * diag(N)
  Vi <- solve(V)
  Vb <- ginv(crossprod(X, Vi %*% X))
  O  <- L_mat %*% Vb %*% t(L_mat)
  (O + t(O)) / 2
}

## ---------- 2. exact variance of a marker coefficient ----------
## For y = [1, Q, z] b + e with Var(e) = Omega, a scheme using weight matrix Wt
## has sandwich variance (X'WtX)^-1 X'Wt Omega Wt X (X'WtX)^-1. With Wt equal to
## Omega inverse this collapses to the generalized least squares variance, which
## is the minimum; every other choice is evaluated under the true Omega.
var_beta <- function(z, O, Wt) {
  Xm    <- cbind(1, Q, z)
  bread <- solve(crossprod(Xm, Wt %*% Xm))
  Vb    <- bread %*% crossprod(Xm, Wt %*% O %*% Wt %*% Xm) %*% bread
  Vb[ncol(Xm), ncol(Xm)]
}

S2R_GRID <- c(PARS$s2_rep, 2, 5, 10, 20, 40)     # the grid scenario B sweeps

eff_tab <- do.call(rbind, lapply(S2R_GRID, function(s2r) {
  O  <- omega_at(s2r)
  Oi <- solve(O); Di <- diag(1 / diag(O)); Ii <- diag(n)
  v  <- vapply(MARKERS, function(s) {
    z <- G[, s]
    c(gls = var_beta(z, O, Oi), dia = var_beta(z, O, Di), unw = var_beta(z, O, Ii))
  }, numeric(3))
  data.frame(s2_rep          = s2r,
             mean_corr       = mean(abs(cov2cor(O)[upper.tri(O)])),
             var_ratio       = max(diag(O)) / min(diag(O)),
             diag_vs_gls     = mean(v["dia", ] / v["gls", ]),
             diag_vs_gls_max = max(v["dia", ] / v["gls", ]),
             unw_vs_gls      = mean(v["unw", ] / v["gls", ]),
             unw_vs_gls_max  = max(v["unw", ] / v["gls", ]),
             row.names = NULL)
}))

## ---------- 3. are the markers orthogonal to the block structure? ----------
## The share of each marker's variance lying in the space of the replicate
## incidence, aggregated to accession level and centred; this space contains
## the year space. The share describes the design.
Bacc <- matrix(0, n, ncol(Zr), dimnames = list(levels(pheno$Id), colnames(Zr)))
for (j in seq_len(ncol(Zr))) Bacc[, j] <- tapply(Zr[, j], pheno$Id, sum)
Bacc <- Bacc[, colSums(Bacc) > 0, drop = FALSE]
Bc   <- scale(Bacc, center = TRUE, scale = FALSE)
Pb   <- Bc %*% ginv(crossprod(Bc)) %*% t(Bc)

blk_share <- vapply(MARKERS, function(s) {
  z <- scale(G[, s], center = TRUE, scale = FALSE)
  as.numeric(crossprod(z, Pb %*% z) / crossprod(z))
}, numeric(1))

## ---------- 3b. where the correlation comes from ----------
## Omega with the replicate variance set to zero, against the calibrated one.
O_rep0 <- omega_at(0); O_cal <- omega_at(PARS$s2_rep)
v_rep0 <- vapply(MARKERS, function(s) {
  z <- G[, s]
  c(gls = var_beta(z, O_rep0, solve(O_rep0)), dia = var_beta(z, O_rep0, diag(1 / diag(O_rep0))))
}, numeric(2))
corr_source <- data.frame(
  s2_rep        = c(0, PARS$s2_rep),
  mean_abs_corr = c(mean(abs(cov2cor(O_rep0)[upper.tri(O_rep0)])), mean(abs(cov2cor(O_cal)[upper.tri(O_cal)]))),
  max_abs_corr  = c(max(abs(cov2cor(O_rep0)[upper.tri(O_rep0)])),  max(abs(cov2cor(O_cal)[upper.tri(O_cal)]))),
  diag_vs_gls   = c(mean(v_rep0["dia", ] / v_rep0["gls", ]), eff_tab$diag_vs_gls[1]))

## the share expected by chance: the median share of a vector unrelated to a
## space of rank d in n - 1 dimensions is the median of Beta(d/2, (n-1-d)/2),
## and the mean is d / (n - 1)
chance_median <- function(d) qbeta(0.5, d / 2, (n - 1 - d) / 2)

## The share expected for a marker unrelated to the blocks is the rank of the
## projection over n - 1: a centred vector lies in n - 1 dimensions, and the
## aggregated replicate columns are collinear where every accession of a year
## is in every block of that year.
n_block_cols <- ncol(Bc)
blk_rank     <- qr(Bc)$rank
blk_expected <- blk_rank / (n - 1)

blk_tab <- data.frame(
  Quantity = c("block-space columns", "rank of the block space", "accessions",
               "share of marker variance in the block space, median",
               "the same, mean", "the same, maximum",
               "mean share expected by chance, rank / (n - 1)",
               "median share expected by chance",
               "markers' median over the chance median",
               "markers' mean over the chance mean",
               "the expectation as first coded, columns / n"),
  Value = c(format(n_block_cols), format(blk_rank), format(n),
            format(round(median(blk_share), 4), nsmall = 4),
            format(round(mean(blk_share), 4), nsmall = 4),
            format(round(max(blk_share), 4), nsmall = 4),
            format(round(blk_expected, 4), nsmall = 4),
            format(round(chance_median(blk_rank), 4), nsmall = 4),
            format(round(median(blk_share) / chance_median(blk_rank), 3), nsmall = 3),
            format(round(mean(blk_share) / blk_expected, 3), nsmall = 3),
            format(round(n_block_cols / n, 4), nsmall = 4)))

## ---------- 4. what share of the variance can weighting act on? ----------
## An adjusted mean carries genetic variance and first-stage estimation error,
## and Omega describes only the second. The share is computed on the mean
## diagonal and on a difference basis (Schmidt et al. 2019, Genetics 212, eq.
## 58): centred traces, so that anything common to all means, which the
## intercept absorbs, drops out; the centred trace divided by n - 1 is the
## average semivariance (ASV) of Piepho (2019, Biometrical Journal 61, Eq. 11),
## the basis the paper uses. The ratio of genetic to total on this basis is the
## heritability of the adjusted means on a genotype-difference basis.
Pc_ <- diag(n) - matrix(1 / n, n, n)
share_tab <- do.call(rbind, lapply(S2R_GRID, function(s2r) {
  O <- omega_at(s2r)
  gen <- PARS$s2_g * mean(diag(A))
  est <- mean(diag(O))
  gen_d <- PARS$s2_g * sum(diag(Pc_ %*% A)) / (n - 1)
  est_d <- sum(diag(Pc_ %*% O)) / (n - 1)
  data.frame(s2_rep = s2r, genetic = gen, estimation = est,
             estimation_share = est / (gen + est),
             estimation_share_diff = est_d / (gen_d + est_d),
             heritability_diff = gen_d / (gen_d + est_d),
             genetic_asv = gen_d, estimation_asv = est_d, row.names = NULL)
}))

## the Omega-only cost times the estimation share; the quantity the paper uses
## is wk_tab below
combined <- data.frame(
  s2_rep = S2R_GRID,
  unw_cost_on_Omega   = round(eff_tab$unw_vs_gls - 1, 4),
  estimation_share    = round(share_tab$estimation_share, 4),
  implied_cost        = round((eff_tab$unw_vs_gls - 1) *
                              share_tab$estimation_share, 5),
  implied_cost_on_SD  = round(sqrt(1 + (eff_tab$unw_vs_gls - 1) *
                              share_tab$estimation_share) - 1, 5))

## ---------- 5. the efficiency with the kernel ----------
## The second stage fits a genetic kernel, so the covariance a scheme gets
## wrong is not Omega but its share of s2_g A + Omega. The sandwich is taken
## with the true covariance s2_g A + Omega and each scheme's working one: the
## rotation s2_g A + Omega (generalized least squares, the minimum),
## inverse-variance weighting s2_g A + diag(Omega), the unweighted analysis
## s2_g A + c I with c the mean diagonal of Omega, since that scheme estimates
## its own scale. This quantity reproduces the simulated ratios marker by
## marker in scenario C (10d, 10e). Reported on the standard-deviation scale,
## mean and maximum over the simulation markers, beside what scenario B
## measured.
wk_tab <- do.call(rbind, lapply(S2R_GRID, function(s2r) {
  O  <- omega_at(s2r)
  Vt <- PARS$s2_g * A + O
  Wg <- solve(Vt)
  Wd <- solve(PARS$s2_g * A + diag(diag(O)))
  Wu <- solve(PARS$s2_g * A + mean(diag(O)) * diag(n))
  v  <- vapply(MARKERS, function(s) {
    z <- G[, s]
    c(g = var_beta(z, Vt, Wg), d = var_beta(z, Vt, Wd), u = var_beta(z, Vt, Wu))
  }, numeric(3))
  data.frame(s2_rep      = s2r,
             diag_SD     = mean(sqrt(v["d", ] / v["g", ])) - 1,
             diag_SD_max = max(sqrt(v["d", ] / v["g", ])) - 1,
             unw_SD      = mean(sqrt(v["u", ] / v["g", ])) - 1,
             unw_SD_max  = max(sqrt(v["u", ] / v["g", ])) - 1,
             row.names = NULL)
}))

## what scenario B measured on the same grid, read from its output
e09 <- new.env()
load(file.path(OUT, "09_simB.RData"), envir = e09)
b_meas <- e09$b_ratio[order(e09$b_ratio$s2_rep), c("s2_rep", "diag_over_rot", "unw_over_rot")]
rm(e09)
stopifnot(isTRUE(all.equal(b_meas$s2_rep, S2R_GRID, tolerance = 1e-6)))
bound_tab <- data.frame(
  s2_rep          = S2R_GRID,
  predicted_diag  = round(wk_tab$diag_SD, 5),
  measured_diag   = round(b_meas$diag_over_rot - 1, 5),
  predicted_unw   = round(wk_tab$unw_SD, 5),
  measured_unw    = round(b_meas$unw_over_rot - 1, 5),
  old_bound_diag  = round(sqrt(1 + (eff_tab$diag_vs_gls - 1) * share_tab$estimation_share) - 1, 5),
  old_bound_unw   = combined$implied_cost_on_SD)

## ---------- checks ----------
checks <- data.frame(
  Check = c("generalized least squares is the minimum at every level",
            "ignoring the correlation costs more as correlation rises",
            "markers lie within twice the chance share of the block space",
            "the estimation share is below one fifth throughout",
            "Omega is positive definite at every level"),
  Result = c(all(eff_tab$diag_vs_gls >= 1 - 1e-8) &&
               all(eff_tab$unw_vs_gls >= 1 - 1e-8),
             eff_tab$diag_vs_gls[length(S2R_GRID)] > eff_tab$diag_vs_gls[1],
             median(blk_share) < 2 * blk_expected,
             all(share_tab$estimation_share < 0.20),
             all(vapply(S2R_GRID, function(s)
               min(eigen(omega_at(s), symmetric = TRUE,
                         only.values = TRUE)$values) > 0, logical(1)))))

## ---------- report ----------
cat("\n=== S2 why the weighting schemes agree ===\n")

cat("\nexact variance of a marker coefficient relative to generalized least\n")
cat("squares, averaged over the", length(MARKERS), "simulation markers\n")
print(eff_tab, digits = 5, row.names = FALSE)
cat("   diag_vs_gls is the cost of ignoring the correlation, which is the only\n")
cat("   thing the rotation exploits that inverse-variance weighting cannot.\n")
cat("   unw_vs_gls is the cost of ignoring the heteroscedasticity as well.\n")

cat("\nare the markers aligned with the block structure?\n")
print(blk_tab, row.names = FALSE)
cat("\nwhere the correlation at the observed block variance comes from\n")
print(corr_source, digits = 5, row.names = FALSE)
cat("   unchanged with the replicate variance at zero: it is the fixed year effects\n")
cat("   estimated on an unbalanced design, not shared replicate effects\n")

cat("\nwhat share of an adjusted mean's variance can weighting act on?\n")
print(share_tab, digits = 5, row.names = FALSE)

cat("\nthe cost as first coded, (Omega-only cost - 1) times the estimation share\n")
print(combined, row.names = FALSE)
cat("   kept for the record; neither a first-order term nor a bound (see section 5)\n")
cat("\nthe efficiency with the kernel, standard-deviation scale, mean and maximum over the markers\n")
print(wk_tab, digits = 4, row.names = FALSE)
cat("\npredicted against what scenario B measured, on the standard-deviation scale\n")
print(bound_tab, row.names = FALSE)
cat(sprintf("   the measured values carry Monte Carlo error (%s replicates); the predictions carry none\n", format(NREP_B, big.mark = ",")))

cat("\nchecks\n"); print(checks, row.names = FALSE)

## saved before the checks are enforced, so a failed check does not lose the run
save(eff_tab, blk_tab, blk_share, share_tab, combined, wk_tab, bound_tab, checks,
     S2R_GRID, n_block_cols, blk_rank, blk_expected, corr_source,
     file = file.path(OUT, "S2_efficiency.RData"))

cat("\nsaved S2_efficiency.RData |",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

if (!all(checks$Result)) stop("efficiency checks failed")
