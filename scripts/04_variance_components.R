################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##               Script 04: variance components on the rotated scale
##
##   After rotation the residual variance has a known target. Var(W mu) is the
##   identity by construction, so a correctly specified second stage must return
##   a residual of one. This script fits the additive model and four extensions,
##   and asks in each case whether the residual is consistent with that target.
##
##     additive
##     additive + WW'         independent genetic effect, kernel = Omega inverse
##     additive + A#A         epistasis
##     additive + H           heterozygote structure
##     additive + A#A + WW'
##
##   Intervals for the residual come from profile likelihood on the REML
##   deviance, the 95% interval being the set of values within qchisq(0.95, 1)
##   of the minimum. A parametric bootstrap is run as a cross-check.
##
##   Output: 04_variance.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "02_rotation.RData"))
t0 <- Sys.time()

Xn <- cbind(ic_rot, Q_rot)
yv <- y_rot

## ---------- 1. kernels ----------
## The additive and the heterozygosity kernel are singular, centring having
## removed the direction of the ones vector. Each receives the matrix of ones,
## K + cJ with c = 1 (with_ones in script 00); the rotation transforms cJ into
## the outer product of the rotated intercept: the residual likelihood, the
## estimates and the deviance are the same for every c and the same as for the
## singular kernel. The kernels without the matrix of ones are kept for the
## rank table and for the bootstrap, which simulates from the fitted
## covariance.
## A kernel formed as a triple product is symmetrised, so that the Cholesky
## factorisation below does not depend on which triangle it reads.
sym <- function(M) (M + t(M)) / 2

## The epistatic and the heterozygosity kernel are scaled to the mean diagonal
## of A before rotation, so that their variance components are comparable with
## the additive one. The scale changes no deviance or AIC: the variance
## component absorbs it.
scale_to <- function(M, target) M * (target / mean(diag(M)))
mdA <- mean(diag(A))

K_add0 <- sym(K_rot)                     # the rotated additive kernel as it is, rank n - 1
K_add  <- with_ones(A, W)                # W(A + cJ)W', what the objective receives
K_ind  <- sym(WWt)                       # equals Omega inverse
K_epi <- sym(W %*% scale_to(A * A, mdA) %*% t(W))   # Hadamard square

## The heterozygosity kernel, het_kernel in script 00: one centred column per
## marker at which some but not all accessions are called heterozygous. Only
## called genotypes enter; an imputed cell takes the marker's mean and vanishes
## on centring. The kernel is a structure, not an estimate of a dominance
## component.
hk            <- het_kernel(G, MISS, mdA)
K_het0        <- hk$K
n_het_markers <- hk$m
n_het_calls   <- hk$n_het_calls
K_het_r0 <- sym(W %*% K_het0 %*% t(W))   # rotated, as it is, rank n - 1
K_het    <- with_ones(K_het0, W)         # with the matrix of ones, what the objective receives

## how distinct each structure is from the additive one, before rotation
offd <- function(M) M[upper.tri(M)]
cor_epi_A <- cor(offd(A), offd(scale_to(A * A, mdA)))
cor_het_A <- cor(offd(A), offd(K_het0))

## the rank table describes the kernels without the matrix of ones
kern_list <- list(additive = K_add0, epistasis = K_epi,
                  heterozygote = K_het_r0, independent = K_ind)
eig <- lapply(kern_list, function(M) eigen(M, symmetric = TRUE, only.values = TRUE)$values)

kernel_tab <- data.frame(
  Kernel         = c("additive", "epistasis A#A", "heterozygote", "independent WW'"),
  Effective_rank = sapply(eig, function(e) sum(e > 1e-6 * max(e))),
  Min_eigenvalue = signif(sapply(eig, min), 3),
  Max_eigenvalue = signif(sapply(eig, max), 3),
  row.names = NULL)

## ---------- 2. REML objective ----------
## The parameters are on the log scale, the last one the residual; y is an
## argument so that the bootstrap below can use the same function. With the
## matrix of ones on the singular kernels, H is positive definite for every
## residual, including zero. The determinants come from Cholesky factors,
## which fail if H is not positive definite.
## A residual below VE_BOUNDARY, three orders of magnitude under its target of
## one, counts as a boundary solution.
VE_BOUNDARY <- 1e-3

neg2 <- function(pars, kernels, y) {
  ve <- exp(pars[length(pars)])
  H <- diag(n) * ve
  for (k in seq_along(kernels)) H <- H + kernels[[k]] * exp(pars[k])
  val <- try({
    R     <- chol(H)
    Hi    <- chol2inv(R)
    XtHiX <- crossprod(Xn, Hi %*% Xn)
    Rx    <- chol(XtHiX)
    P     <- Hi - Hi %*% Xn %*% chol2inv(Rx) %*% crossprod(Xn, Hi)
    2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) +
      as.numeric(crossprod(y, P %*% y))
  }, silent = TRUE)
  if (inherits(val, "try-error") || !is.finite(val)) return(1e10)
  val
}

## Nelder-Mead can stop short of the optimum, so it is restarted from its own
## answer until the deviance stops improving.
optimise_reml <- function(start_log, kernels, y, max_restarts = 20) {
  par <- start_log; dev <- Inf; conv <- NA_integer_; used <- 0L
  for (i in seq_len(max_restarts)) {
    o <- optim(par, neg2, kernels = kernels, y = y, method = "Nelder-Mead",
               control = list(maxit = 5000, reltol = 1e-11))
    used <- i; conv <- o$convergence
    improved <- dev - o$value
    par <- o$par; dev <- o$value
    if (improved < 1e-9) break
  }
  list(par = par, dev = dev, conv = conv, restarts = used)
}

fit_model <- function(kernels, start) {
  o <- optimise_reml(log(start), kernels, yv)
  list(vc = exp(o$par), dev = o$dev, k = length(start),
       conv = o$conv, restarts = o$restarts)
}

models <- list(
  additive    = list(kern = list(K_add),               start = c(70, 1)),
  add_WWt     = list(kern = list(K_add, K_ind),        start = c(60, 5, 1)),
  add_epi     = list(kern = list(K_add, K_epi),        start = c(50, 20, 1)),
  add_het     = list(kern = list(K_add, K_het),        start = c(50, 20, 1)),
  add_epi_WWt = list(kern = list(K_add, K_epi, K_ind), start = c(40, 20, 5, 1))
)

fits <- lapply(models, function(m) fit_model(m$kern, m$start))

## AIC on the REML deviance, with the number of variance components as the
## penalty. All models have the same fixed effects, so their REML deviances
## are comparable. Each single extension is also given a likelihood-ratio test
## against the additive model below.
model_tab <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]; v <- f$vc
  data.frame(Model     = nm,
             Vg_add    = round(v[1], 3),
             Extra_1   = if (f$k >= 3) round(v[2], 3) else NA,
             Extra_2   = if (f$k >= 4) round(v[3], 3) else NA,
             Ve        = round(v[f$k], 4),
             deviance  = round(f$dev, 3),
             AIC       = round(f$dev + 2 * f$k, 3),
             Restarts  = f$restarts,
             Converged = f$conv == 0)
}))
best_model <- model_tab$Model[which.min(model_tab$AIC)]

## Each single extension reduces to the additive model when its added component
## is zero, a value on the boundary of the parameter space, so the
## likelihood-ratio statistic is referred to the equal mixture of chi-squared
## distributions on zero and one degree of freedom (Self and Liang, 1987; Stram
## and Lee, 1994).
lrt_tab <- do.call(rbind, lapply(c("add_WWt", "add_epi", "add_het"), function(nm) {
  lr <- fits$additive$dev - fits[[nm]]$dev
  data.frame(Model = nm, LR = round(lr, 3),
             p_mixture = signif(0.5 * pchisq(max(lr, 0), 1, lower.tail = FALSE), 3),
             dAIC = round((fits[[nm]]$dev + 2 * fits[[nm]]$k) -
                          (fits$additive$dev + 2 * fits$additive$k), 3),
             row.names = NULL)
}))

## ---------- 3. profile interval for the residual ----------
## The residual is fixed, the other components are re-optimised, and the limits
## are where the deviance increases by qchisq(0.95, 1) above its minimum. With
## one free component a one-dimensional optimiser is used.
LOGVG    <- log(c(1e-4, 1e5))     # bracket for the genetic variance
VE_UPPER <- 100                   # bracket for the upper end of the interval

profile_dev <- function(ve, kernels, start) {
  f <- function(p) neg2(c(p, log(ve)), kernels, yv)
  if (length(start) == 1) return(optimize(f, interval = LOGVG, tol = 1e-10)$objective)
  par <- log(start); dev <- Inf
  for (i in 1:10) {
    o <- optim(par, f, method = "Nelder-Mead",
               control = list(maxit = 5000, reltol = 1e-11))
    improved <- dev - o$value; par <- o$par; dev <- o$value
    if (improved < 1e-9) break
  }
  dev
}

## The optimum of the genetic variance is checked to lie inside its bracket.
vg_at_hat   <- optimize(function(p) neg2(c(p, log(fits$additive$vc[2])), list(K_add), yv),
                        interval = LOGVG, tol = 1e-10)$minimum
vg_interior <- vg_at_hat > LOGVG[1] + 1e-6 && vg_at_hat < LOGVG[2] - 1e-6

profile_interval <- function(kernels, start_g, ve_hat) {
  dev_min <- profile_dev(ve_hat, kernels, start_g)
  thresh  <- qchisq(0.95, 1)
  g <- function(ve) profile_dev(ve, kernels, start_g) - dev_min - thresh
  ## the lower limit is 0 when the deviance at a residual near zero is still
  ## inside the interval
  lo      <- if (g(1e-6) < 0) 0 else uniroot(g, c(1e-6, ve_hat))$root
  bounded <- g(VE_UPPER) > 0
  hi      <- if (bounded) uniroot(g, c(ve_hat, VE_UPPER))$root else NA_real_
  list(lower = lo, upper = hi, dev_min = dev_min, bounded = bounded)
}

ci_add <- profile_interval(list(K_add), fits$additive$vc[1], fits$additive$vc[2])
ci_wwt <- profile_interval(list(K_add, K_ind), fits$add_WWt$vc[1:2], fits$add_WWt$vc[3])

ci_tab <- data.frame(
  Model   = c("additive", "additive + WW'"),
  Ve      = round(c(fits$additive$vc[2], fits$add_WWt$vc[3]), 4),
  Lower   = round(c(ci_add$lower, ci_wwt$lower), 4),
  Upper   = round(c(ci_add$upper, ci_wwt$upper), 4),
  Bounded = c(ci_add$bounded, ci_wwt$bounded))
ci_tab$Covers_one <- ci_tab$Lower <= 1 & 1 <= ci_tab$Upper

## ---------- 3b. the Wald interval PROC MIXED would give ----------
## SAS PROC MIXED forms a Satterthwaite-Wald interval with nu = 2 Z^2 degrees
## of freedom, Z = Ve / se(Ve), and limits nu Ve / qchisq(1 - alpha/2, nu) and
## nu Ve / qchisq(alpha/2, nu). The standard error comes from the curvature of
## the profile deviance at the optimum: profile_dev returns -2 log L, so the
## observed information is half its second derivative.
## The second derivative is a central difference, computed at two step sizes;
## their agreement is one of the checks.
VE_CURV_STEP <- 0.05 * fits$additive$vc[2]

profile_curvature <- function(h) {
  v <- fits$additive$vc[2] + c(-h, 0, h)
  d <- vapply(v, function(x) profile_dev(x, list(K_add), fits$additive$vc[1]), 0)
  (d[1] - 2 * d[2] + d[3]) / h^2
}

d2_main <- profile_curvature(VE_CURV_STEP)
d2_half <- profile_curvature(VE_CURV_STEP / 2)
curv_agree <- abs(d2_main - d2_half) / abs(d2_main)

se_ve   <- sqrt(2 / d2_main)
nu_satt <- 2 * (fits$additive$vc[2] / se_ve)^2
wald_lo <- nu_satt * fits$additive$vc[2] / qchisq(0.975, nu_satt)
wald_hi <- nu_satt * fits$additive$vc[2] / qchisq(0.025, nu_satt)

wald_tab <- data.frame(
  Method  = c("profile likelihood", "Satterthwaite-Wald"),
  Lower   = round(c(ci_add$lower, wald_lo), 4),
  Upper   = round(c(ci_add$upper, wald_hi), 4),
  df      = c(NA_real_, round(nu_satt, 4)))
wald_tab$Covers_one <- wald_tab$Lower <= 1 & 1 <= wald_tab$Upper

## ---------- 4. deviance ladder around the additive model ----------
## the fitted value is added to the grid
ve_grid <- sort(unique(c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5,
                         round(fits$additive$vc[2], 4))))
ladder <- data.frame(
  Ve = ve_grid,
  Rise_above_minimum = round(sapply(ve_grid, function(v)
    profile_dev(v, list(K_add), fits$additive$vc[1])) - ci_add$dev_min, 4))
ladder$Inside_95 <- ladder$Rise_above_minimum <= qchisq(0.95, 1)

## ---------- 5. parametric bootstrap on the additive model ----------
set.seed(SEED_MISC + 1)
vg_hat <- fits$additive$vc[1]; ve_hat <- fits$additive$vc[2]
H_hat  <- sym(K_add0 * vg_hat + diag(n) * ve_hat)   # the fitted covariance, without the matrix of ones
Lh     <- t(chol(H_hat))
beta_hat <- solve(crossprod(Xn), crossprod(Xn, yv))
fit_mean <- as.numeric(Xn %*% beta_hat)

boot_ve   <- numeric(NREP_BOOT)
boot_conv <- integer(NREP_BOOT)
boot_dev  <- numeric(NREP_BOOT)
for (b in seq_len(NREP_BOOT)) {
  yv_b <- fit_mean + as.numeric(Lh %*% rnorm(n))
  o <- optimise_reml(log(c(vg_hat, ve_hat)), list(K_add), yv_b, max_restarts = 3)
  boot_ve[b]   <- exp(o$par[2])
  boot_conv[b] <- o$conv
  boot_dev[b]  <- o$dev
}

## With a genetic variance far larger than the residual, REML puts the residual
## on the boundary in a share of the resamples. Those draws are counted as
## boundary solutions, not failures; a failure is a draw that neither converged
## nor reached the boundary.
boot_boundary <- boot_ve < VE_BOUNDARY
boot_failed   <- boot_conv != 0 & !boot_boundary

boot_tab <- data.frame(
  Method   = "parametric bootstrap, additive model",
  Ve       = round(ve_hat, 4),
  Lower    = round(unname(quantile(boot_ve, 0.025)), 4),
  Upper    = round(unname(quantile(boot_ve, 0.975)), 4),
  Boundary = sum(boot_boundary),
  Failed   = sum(boot_failed))
boot_tab$Covers_one <- boot_tab$Lower <= 1 & 1 <= boot_tab$Upper

## ---------- checks ----------
checks <- data.frame(
  Check = c("all model fits converged",
            "genetic variance interior to its search bracket",
            "profile interval for the additive model bounded above",
            "profile interval for additive + WW' bounded above",
            "no bootstrap fit diverged",
            "every non-boundary bootstrap fit converged",
            "all kernels positive semidefinite",
            "profile curvature stable across the finite-difference step",
            "Satterthwaite degrees of freedom finite and positive"),
  Result = c(all(sapply(fits, function(f) f$conv) == 0),
             vg_interior,
             ci_add$bounded,
             ci_wwt$bounded,
             all(is.finite(boot_dev)) && min(boot_dev) > -1e6,
             !any(boot_failed),
             all(sapply(eig, min) > -1e-8 * max(sapply(eig, max))),
             is.finite(curv_agree) && curv_agree < 0.05,
             is.finite(nu_satt) && nu_satt > 0))

## ---------- report ----------
cat("\n=== 04 variance components ===\n")
cat("\nkernels\n");                          print(kernel_tab, row.names = FALSE)
cat("\nmodels\n");                           print(model_tab, row.names = FALSE)
cat("   Extra_1 and Extra_2 are the second and third kernels named in the model\n")
cat("   correlation with A before rotation: epistasis", round(cor_epi_A, 4),
    ", heterozygote", round(cor_het_A, 4), "
")
cat("   the heterozygote kernel spans", n_het_markers, "markers and",
    format(n_het_calls, big.mark = ","), "called heterozygotes; imputed cells contribute nothing\n")
cat("   best by AIC:", best_model, "\n")
cat("\nlikelihood-ratio tests of each single extension against the additive model\n")
print(lrt_tab, row.names = FALSE)
cat("\nprofile-likelihood intervals for the residual\n"); print(ci_tab, row.names = FALSE)
cat("\nthe interval PROC MIXED would give on the additive model\n")
print(wald_tab, row.names = FALSE)
cat("   se(Ve)", signif(se_ve, 4), "from a profile curvature stable to",
    signif(curv_agree, 2), "across a halved step\n")
cat("\ndeviance ladder, additive model\n");  print(ladder, row.names = FALSE)
cat("\nbootstrap\n");                        print(boot_tab, row.names = FALSE)
cat("\nchecks\n");                           print(checks, row.names = FALSE)

if (!all(checks$Result)) stop("variance-component checks failed")

save(kernel_tab, model_tab, ci_tab, wald_tab, ladder, boot_tab, checks, fits, boot_ve,
     nu_satt, se_ve, curv_agree, VE_CURV_STEP,
     boot_conv, boot_dev, boot_boundary, boot_failed,
     eig, best_model, lrt_tab, vg_interior,
     VE_BOUNDARY, ONES_C,
     ci_add, ci_wwt, K_add, K_add0, K_epi, K_het, K_het_r0, K_ind, K_het0,
     cor_epi_A, cor_het_A, n_het_markers, n_het_calls,
     file = file.path(OUT, "04_variance.RData"))

cat("\nsaved 04_variance.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
