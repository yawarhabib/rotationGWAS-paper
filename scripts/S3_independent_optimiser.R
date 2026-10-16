################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script S3: the multi-kernel likelihood against an independent
##                   implementation
##
##   Scenario D depends on ve_kernels, a REML likelihood coded directly in
##   script 07 and optimised by a restarted Nelder-Mead. This script fits the
##   same model with the propto structure of glmmTMB, an independent
##   implementation with a different parameterisation and a different
##   optimiser, and compares the estimates.
##
##   glmmTMB receives the rotated kernel directly and factorises it, which fails
##   for a singular kernel; each singular kernel is therefore passed as
##   K + cJ, J the matrix of ones (with_ones in script 00). Under rotation cJ
##   becomes the outer product of the rotated intercept, which lies in the
##   fixed-effect space, so the residual likelihood, the estimates and the
##   deviance are the same for every c and the same as for the singular kernel
##   (script S4). Both implementations receive the identical rotated kernel,
##   ve_kernels by way of the unrotated K + cJ, which it rotates itself.
##
##   Both solutions are scored on the same REML deviance, and a replicate
##   passes if the two agree or if the directly coded solution has the lower
##   deviance.
##
##   Output: S3_optimiser.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans); library(rrBLUP); library(glmmTMB)})

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
t0 <- Sys.time()

NREP_S3 <- 5          # replicates per structure; each costs a glmmTMB fit
S2_TRUE <- 25         # the withheld variance, mid-range of scenario D
TOL_VE  <- 0.05       # on a residual whose target is one
TOL_REL <- 0.05       # on components of order ten to a hundred

P <- list(pheno = pheno, id_idx = ID_IDX, yr_idx = YR_IDX, yr_rep = YR_REP,
          id_yr = ID_YR, pars = PARS, rootA = ROOT_A, rootE = ROOT_E,
          rootD = ROOT_D, Gsub = G_SUB, Q = Q, A = A, n = n, L = L_mat)

ids <- sprintf("r%03d", seq_len(n))
sym <- function(M) (M + t(M)) / 2
J   <- matrix(1, n, n)

## a kernel that centring has made singular: the smallest eigenvalue is zero
## to rounding
is_singular <- function(M) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  min(ev) < 1e-8 * max(ev)
}
## the kernel on the original scale as ve_kernels receives it, and its rotation
## as glmmTMB receives it; the two are the same matrix
unrotated <- function(M) if (is_singular(M)) M + ONES_C * J else M
rotated   <- function(M, W) {
  Mr <- if (is_singular(M)) with_ones(M, W) else sym(W %*% M %*% t(W))
  dimnames(Mr) <- list(ids, ids)
  Mr
}

## The second component under propto. glmmTMB reports a covariance matrix
## proportional to the supplied kernel, so dividing elementwise and averaging
## recovers the scalar.
propto_component <- function(fit, grp, K) {
  vb <- VarCorr(fit)$cond[[grp]]
  s  <- attr(vb, "stddev"); cr <- attr(vb, "correlation")
  mean((diag(s) %*% cr %*% diag(s)) / K)
}

## The common measure: the REML deviance of ve_kernels, evaluated at the
## variance components it is handed. Both solutions are scored on it.
make_dev <- function(y, Xn, Ka, Kb) function(va, vb, ve) {
  H <- diag(length(y)) * ve + Ka * va + Kb * vb
  val <- try({
    R  <- chol(H); Hi <- chol2inv(R)
    XtHiX <- crossprod(Xn, Hi %*% Xn); Rx <- chol(XtHiX)
    Pm <- Hi - Hi %*% Xn %*% chol2inv(Rx) %*% crossprod(Xn, Hi)
    2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) +
      as.numeric(crossprod(y, Pm %*% y))
  }, silent = TRUE)
  if (inherits(val, "try-error") || !is.finite(val)) NA_real_ else val
}

STRUCT <- data.frame(source = c("epi", "het"),
                     label  = c("epistasis A#A", "heterozygote"),
                     k      = 1:2, stringsAsFactors = FALSE)

cat("S3:", NREP_S3, "replicates x", nrow(STRUCT), "structures, two engines each\n")
cat("   singular kernels receive K + cJ with c =", ONES_C, "; additive:", is_singular(A),
    "epistatic:", is_singular(E_had), "heterozygote:", is_singular(D_dom), "\n")

## ---------- one replicate ----------
one <- function(r, i) {
  src <- STRUCT$source[i]
  set.seed(SEED_MISC + 5 + r * 10 + STRUCT$k[i])

  b  <- setNames(rep(0, length(MARKERS)), MARKERS)
  s1 <- stage1_fast(simulate_y(b,
                               s2_epi = if (src == "epi") S2_TRUE else 0,
                               s2_dom = if (src == "het") S2_TRUE else 0,
                               P = P), P = P)
  Tr <- make_transforms(s1$Omega)$rotated
  K2 <- if (src == "het") D_dom else E_had

  ## the same kernels, reached from the unrotated side so that ve_kernels keeps
  ## its rotated intercept
  mine <- ve_kernels(s1, Tr, list(add = unrotated(A), extra = unrotated(K2)), P = P)

  Ka <- rotated(A, Tr); Kb <- rotated(K2, Tr)
  Qr <- Tr %*% P$Q
  y  <- as.numeric(Tr %*% s1$mu)
  ic <- as.numeric(Tr %*% rep(1, n))
  Xn <- cbind(ic, Qr)
  dv <- make_dev(y, Xn, Ka, Kb)

  df <- data.frame(y = y, ic = ic, PC1 = Qr[, 1], PC2 = Qr[, 2], PC3 = Qr[, 3],
                   Id = factor(ids, levels = ids), g = factor(1), g2 = factor(1))

  ## Warnings of glmmTMB are recorded per fit, not treated as failures.
  warn <- character(0)
  gt <- withCallingHandlers(
    tryCatch(
      glmmTMB(y ~ 0 + ic + PC1 + PC2 + PC3 +
                propto(0 + Id | g, Ka) + propto(0 + Id | g2, Kb),
              data = df, REML = TRUE),
      error = function(e) NULL),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  Va_t <- if (is.null(gt)) NA_real_ else propto_component(gt, "g",  Ka)
  V2_t <- if (is.null(gt)) NA_real_ else propto_component(gt, "g2", Kb)
  Ve_t <- if (is.null(gt)) NA_real_ else sigma(gt)^2

  data.frame(
    rep = r, source = src, s2_true = S2_TRUE, c = ONES_C,
    Ve_mine = as.numeric(mine["Ve"]),    Ve_tmb = Ve_t,
    Va_mine = as.numeric(mine["V_add"]), Va_tmb = Va_t,
    V2_mine = as.numeric(mine["V_extra"]), V2_tmb = V2_t,
    dev_mine = as.numeric(mine["deviance"]),
    dev_tmb  = if (is.null(gt)) NA_real_ else dv(Va_t, V2_t, Ve_t),
    conv_mine = as.numeric(mine["converged"]),
    tmb_ok    = as.numeric(!is.null(gt)),
    tmb_warn = if (length(warn)) paste(unique(warn), collapse = " | ") else "",
    row.names = NULL)
}

S3 <- do.call(rbind, lapply(seq_len(NREP_S3), function(r)
        do.call(rbind, lapply(seq_len(nrow(STRUCT)), function(i) one(r, i)))))

## ---------- agreement ----------
S3$d_Ve <- abs(S3$Ve_mine - S3$Ve_tmb)
S3$d_Va <- abs(S3$Va_mine - S3$Va_tmb) / pmax(S3$Va_tmb, 1)
S3$d_V2 <- abs(S3$V2_mine - S3$V2_tmb) / pmax(S3$V2_tmb, 1)

S3$matched <- S3$d_Ve < TOL_VE & S3$d_Va < TOL_REL & S3$d_V2 < TOL_REL
## where they differ, a negative gap means the directly coded solution has the
## lower deviance
S3$dev_gap <- S3$dev_mine - S3$dev_tmb
S3$mine_better <- !S3$matched & S3$dev_gap < 0

agree <- data.frame(
  Quantity = c("residual, absolute", "additive component, relative",
               "second component, relative"),
  Max_matched = c(max(S3$d_Ve[S3$matched], na.rm = TRUE),
                  max(S3$d_Va[S3$matched], na.rm = TRUE),
                  max(S3$d_V2[S3$matched], na.rm = TRUE)),
  row.names = NULL)

## ---------- checks ----------
checks <- data.frame(
  Check = c("every glmmTMB fit returned a result",
            "every directly coded fit converged",
            "matched replicates agree on the residual",
            "matched replicates agree on both components",
            "where they differ, the directly coded solution has the lower deviance",
            "no missing values in the comparison"),
  Result = c(all(S3$tmb_ok == 1),
             all(S3$conv_mine == 1),
             max(S3$d_Ve[S3$matched], na.rm = TRUE) < TOL_VE,
             max(c(S3$d_Va[S3$matched], S3$d_V2[S3$matched]), na.rm = TRUE) < TOL_REL,
             all(S3$matched | S3$mine_better),
             !anyNA(S3[, c("Ve_mine", "Ve_tmb", "V2_mine", "V2_tmb", "dev_tmb")])))

## ---------- report ----------
cat("\n=== S3 the multi-kernel likelihood against glmmTMB ===\n")
cat("\nboth engines fitted to the identical rotated kernel, the singular kernels",
    "carrying the ones matrix with c =", ONES_C, "\n")
print(S3[, c("rep", "source", "Ve_mine", "Ve_tmb", "Va_mine", "Va_tmb",
             "V2_mine", "V2_tmb", "matched")], row.names = FALSE, digits = 5)

cat("\nreplicates where the two implementations matched:",
    sum(S3$matched), "of", nrow(S3), "\n")
print(agree, row.names = FALSE, digits = 3)
cat("   tolerances:", TOL_VE, "absolute on the residual,", TOL_REL,
    "relative on the components\n")

if (any(!S3$matched)) {
  cat("\nreplicates where they parted company, scored on the same deviance\n")
  print(S3[!S3$matched, c("rep", "source", "dev_mine", "dev_tmb", "dev_gap",
                          "mine_better")], row.names = FALSE, digits = 6)
  cat("   a negative gap means the directly coded solution is the better optimum.\n")
  w <- S3$tmb_warn[!S3$matched]
  cat("   glmmTMB reported on that fit:",
      if (all(w == "")) "nothing" else paste(w, collapse = " ; "), "\n")
  cat("   ve_kernels restarts until the deviance stops improving, which is the\n")
  cat("   policy this column is measuring.\n")
}

cat("\nglmmTMB fits returning a result:", sum(S3$tmb_ok), "of", nrow(S3), "\n")
cat("directly coded fits converged:", sum(S3$conv_mine), "of", nrow(S3), "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

## the warnings by replicate
warn_tab <- data.frame(
  rep = S3$rep, source = S3$source, matched = S3$matched,
  warning = ifelse(S3$tmb_warn == "", "none", substr(S3$tmb_warn, 1, 60)),
  row.names = NULL)

cat("\nglmmTMB warnings by replicate\n")
print(warn_tab, row.names = FALSE)

save(S3, agree, checks, warn_tab, NREP_S3, S2_TRUE, ONES_C, TOL_VE, TOL_REL,
     file = file.path(OUT, "S3_optimiser.RData"))

cat("\nsaved S3_optimiser.RData |",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

if (!all(checks$Result)) stop("S3 checks failed")
