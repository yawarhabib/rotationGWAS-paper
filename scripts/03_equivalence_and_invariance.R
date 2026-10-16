################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##              Script 03: equivalence and invariance checks
##
##   Four checks on the rotation, all at the test marker:
##
##     1. ordinary least squares on the rotated data reproduces generalized
##        least squares on the unrotated means
##     2. the intercept must be rotated: a plain column of ones in its place
##        changes the marker effect
##     3. the three two-stage schemes against a single-stage fit to the plot
##        data; least squares on the means is unbiased under every weighting,
##        so an incorrect weighting shows in the standard error and in the
##        test statistic, not in the estimate
##     4. the result does not depend on which valid rotation is used
##
##   Output: 03_equivalence.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(lme4)
library(lme4qtl)
library(rrBLUP)

## Convergence of the lme4 fits: a fit that lme4 flags is refitted with the
## optimiser bobyqa, and the gradient scaled by the Hessian is reported. A
## boundary message (a variance component at zero) is not a failure.
converged   <- function(m) all(grepl("boundary", m@optinfo$conv$lme4$messages))
scaled_grad <- function(m) {
  d <- m@optinfo$derivs
  if (is.null(d)) NA_real_ else max(abs(solve(d$Hessian, d$gradient)))
}
BOBYQA <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

## A fit is also refitted when its scaled gradient exceeds GRAD_TOL.
GRAD_TOL <- 1e-3
fit_relmat <- function(form, data, kernel) {
  m <- relmatLmer(form, data = data, relmat = list(Id = kernel), REML = TRUE)
  refit <- !converged(m) || isTRUE(scaled_grad(m) > GRAD_TOL)
  if (refit)
    m <- relmatLmer(form, data = data, relmat = list(Id = kernel), REML = TRUE,
                    control = BOBYQA)
  list(m = m, refit = refit, ok = converged(m), grad = scaled_grad(m))
}

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "02_rotation.RData"))
t0 <- Sys.time()

snp  <- TEST_SNP
dose <- G[, snp]
mu   <- blues$BLUE

## The inverse of Omega, computed directly.
Omega_i <- solve(Omega)

## ---------- 1. rotated OLS versus GLS ----------
X_orig <- cbind(1, Q, dose)
b_gls  <- solve(t(X_orig) %*% Omega_i %*% X_orig) %*% (t(X_orig) %*% Omega_i %*% mu)

b_naive <- coef(lm(mu ~ Q + dose))

Z_rot_snp <- as.numeric(W %*% dose)
X_rot     <- cbind(ic_rot, Q_rot, Z_rot_snp)
b_rot     <- solve(crossprod(X_rot), crossprod(X_rot, y_rot))

gls_tab <- data.frame(
  Coefficient = c("Intercept", paste0("PC", 1:N_PC), "Marker"),
  Naive_OLS   = round(as.numeric(b_naive), 5),
  Rotated_OLS = round(as.numeric(b_rot), 5),
  GLS         = round(as.numeric(b_gls), 5))

dev_gls <- max(abs(as.numeric(b_rot) - as.numeric(b_gls)))

## ---------- 2. the intercept must be rotated ----------
snp_under <- function(intercept_vec) {
  X <- cbind(intercept_vec, Q_rot, Z_rot_snp)
  solve(crossprod(X), crossprod(X, y_rot))[ncol(X)]
}

intercept_tab <- data.frame(
  Intercept = c("GLS reference", "rotated, W times 1", "rotated, plain 1"),
  Marker_effect = round(c(b_gls[ncol(X_orig)],
                          snp_under(ic_rot),
                          snp_under(rep(1, n))), 5))
intercept_tab$Difference <- round(intercept_tab$Marker_effect -
                                    intercept_tab$Marker_effect[1], 5)

## ---------- 3. against a true single-stage fit ----------
## The kernel and the plot data must name the accessions identically.
stopifnot(identical(dimnames(A), list(rownames(A), rownames(A))),
          identical(levels(pheno$Id), rownames(A)))

dat <- pheno
dat$SNP <- dose[as.character(dat$Id)]
dat$PC1 <- Q[as.character(dat$Id), 1]
dat$PC2 <- Q[as.character(dat$Id), 2]
dat$PC3 <- Q[as.character(dat$Id), 3]

form_ss <- as.formula(paste(TRAIT,
  "~ 1 + SNP + PC1 + PC2 + PC3 + Year + (1|Year:Replication) + (1|Id:Year) + (1|Id)"))
## relmatLmer factorises the kernel, and A is singular, so it receives A + cJ,
## J the matrix of ones (script 00). At plot level the added term lies in the
## space of the intercept, so the marker effect, its standard error, the
## variance components and the REML criterion are those of A itself.
A_ss <- A + ONES_C * matrix(1, n, n)
f_ss    <- fit_relmat(form_ss, dat, A_ss)
fit_ss  <- f_ss$m
ss_coef <- summary(fit_ss)$coefficients["SNP", ]
b_ss    <- as.numeric(ss_coef["Estimate"])
se_ss   <- as.numeric(ss_coef["Std. Error"])

converged_ss <- f_ss$ok

## A transformed row is a combination of accessions, not a named line, so the
## rows get the labels r001, r002 and so on.
ids <- sprintf("r%03d", seq_len(n))

## The three schemes are fitted with the same function and the same kernel
## treatment, so they differ only in the transformation T. The transformed
## kernel T A T' is singular, so each receives T (A + cJ) T' (with_ones in
## script 00); the added term is the outer product of the transformed intercept
## and changes nothing reported here.
fit_two_stage <- function(Tm) {
  d <- data.frame(y   = as.numeric(Tm %*% mu),
                  ic  = as.numeric(Tm %*% rep(1, n)),
                  PC1 = (Tm %*% Q)[, 1], PC2 = (Tm %*% Q)[, 2], PC3 = (Tm %*% Q)[, 3],
                  SNP = as.numeric(Tm %*% dose),
                  Id  = factor(ids, levels = ids))
  Kn <- with_ones(A, Tm)
  dimnames(Kn) <- list(ids, ids)
  ff <- fit_relmat(y ~ 0 + ic + PC1 + PC2 + PC3 + SNP + (1|Id), d, Kn)
  f  <- ff$m
  cf <- summary(f)$coefficients["SNP", ]
  list(beta      = as.numeric(cf["Estimate"]),
       se        = as.numeric(cf["Std. Error"]),
       converged = ff$ok, refit = ff$refit, grad = ff$grad)
}

arm_unw  <- fit_two_stage(Tlist$unweighted)
arm_diag <- fit_two_stage(Tlist$diagonal)
arm_rot  <- fit_two_stage(Tlist$rotated)

max_grad <- max(c(f_ss$grad, arm_unw$grad, arm_diag$grad, arm_rot$grad), na.rm = TRUE)
n_refit  <- sum(c(f_ss$refit, arm_unw$refit, arm_diag$refit, arm_rot$refit))

## The single-stage fit is the reference. An incorrect weighting shows in the
## standard error and in the test statistic, not in the estimate.
b_arm  <- c(b_ss,  arm_unw$beta, arm_diag$beta, arm_rot$beta)
se_arm <- c(se_ss, arm_unw$se,   arm_diag$se,   arm_rot$se)

anchor_tab <- data.frame(
  Analysis = c("single-stage (plot level)", "two-stage, unweighted",
               "two-stage, inverse-variance", "two-stage, rotated"),
  Marker_effect = round(b_arm, 4),
  SE            = round(se_arm, 4),
  t             = round(b_arm / se_arm, 4),
  Distance_from_single_stage = round(abs(b_arm - b_ss), 4),
  SE_ratio                   = round(se_arm / se_ss, 4),
  t_ratio                    = round((b_arm / se_arm) / (b_ss / se_ss), 4),
  Converged = c(converged_ss, arm_unw$converged, arm_diag$converged, arm_rot$converged))

## ---------- 4. invariance to the choice of rotation ----------
## Any W with W'W equal to the inverse of Omega can be used. Four are compared:
## the symmetric square root used throughout, the eigenvector form, the
## inverse of the Cholesky factor, and the symmetric root premultiplied by a
## random orthogonal matrix.
sv <- svd(Omega)
set.seed(SEED_MISC)
Qo <- qr.Q(qr(matrix(rnorm(n * n), n)))

W_list <- list(
  symmetric = sv$u %*% diag(1 / sqrt(sv$d)) %*% t(sv$u),
  eigen     = diag(1 / sqrt(sv$d)) %*% t(sv$u),
  cholesky  = solve(t(chol(Omega))),
  random    = Qo %*% (sv$u %*% diag(1 / sqrt(sv$d)) %*% t(sv$u)))

## Each candidate is checked to be a rotation before the comparison.
validity <- sapply(W_list, function(Wm) max(abs(crossprod(Wm) - Omega_i)))
stopifnot(all(validity < 1e-8))

## mixed.solve tolerates a singular kernel, so A is passed as it is.
wald_under <- function(Wm) {
  y  <- as.numeric(Wm %*% mu)
  Xn <- cbind(as.numeric(Wm %*% rep(1, n)), Wm %*% Q)
  K  <- Wm %*% A %*% t(Wm)
  nf <- mixed.solve(y = y, X = Xn, Z = diag(n), K = K, return.Hinv = TRUE)
  X3 <- cbind(Xn, as.numeric(Wm %*% dose))
  p  <- ncol(X3); v2 <- n - p
  Wi <- solve(crossprod(X3, nf$Hinv %*% X3))
  b  <- Wi %*% crossprod(X3, nf$Hinv %*% matrix(y, n, 1))
  r  <- matrix(y, n, 1) - X3 %*% b
  s2 <- as.double(crossprod(r, nf$Hinv %*% r)) / v2
  se <- sqrt((s2 * Wi)[p, p])
  c(Vg = nf$Vu, Ve = nf$Ve, beta = b[p], se = se,
    lp = -log10(pbeta(v2 / (v2 + b[p]^2 / se^2), v2 / 2, 1 / 2)))
}

## The spread is computed on the unrounded values.
inv_raw <- t(sapply(W_list, wald_under))
inv_tab <- data.frame(Rotation = rownames(inv_raw), round(inv_raw, 8), row.names = NULL)

inv_spread <- data.frame(
  Quantity = c("Vg", "Ve", "beta", "se", "neglog10p"),
  Spread   = apply(inv_raw, 2, function(z) max(z) - min(z)),
  row.names = NULL)

dev_invariance <- max(inv_spread$Spread)

## ---------- report ----------
cat("\n=== 03 equivalence and invariance ===\n")
## the test marker on which every table below is based
cat("   anchored on", TEST_SNP, "| chromosome",
    map$Chromosome[match(TEST_SNP, map$SNP)], "position",
    format(map$Position[match(TEST_SNP, map$SNP)], big.mark = ","), "\n")
cat("\n1. rotated OLS versus GLS\n")
print(gls_tab, row.names = FALSE)
cat("   maximum difference, rotated against GLS:", signif(dev_gls, 3), "\n")

cat("\n2. the intercept must be rotated\n")
print(intercept_tab, row.names = FALSE)

cat("\n3. against a true single-stage fit\n")
print(anchor_tab, row.names = FALSE)
cat("   fits refitted with bobyqa:", n_refit,
    "| largest scaled gradient:", signif(max_grad, 3), "\n")

cat("\n4. invariance to the choice of rotation\n")
print(inv_tab, row.names = FALSE)
cat("\n   validity, max |W'W - Omega inverse| per candidate:\n")
print(data.frame(Rotation = names(validity), Deviation = signif(validity, 3)),
      row.names = FALSE)
cat("\n   spread across rotations\n")
print(inv_spread, row.names = FALSE)

checks <- data.frame(
  Check  = c("rotated OLS equals GLS",
             "all four candidates are valid rotations",
             "results invariant across rotations",
             "single-stage fit converged",
             "all two-stage arms converged",
             "scaled gradient small at every fit"),
  Result = c(dev_gls < 1e-8,
             all(validity < 1e-8),
             dev_invariance < 1e-6,
             converged_ss,
             all(c(arm_unw$converged, arm_diag$converged, arm_rot$converged)),
             max_grad < GRAD_TOL))

cat("\nchecks\n")
print(checks, row.names = FALSE)

if (!all(checks$Result)) stop("equivalence or invariance checks failed")

save(gls_tab, intercept_tab, anchor_tab, inv_tab, inv_spread, validity,
     checks, b_ss, se_ss, dat, dev_gls, dev_invariance, max_grad, n_refit,
     file = file.path(OUT, "03_equivalence.RData"))

cat("\nsaved 03_equivalence.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
