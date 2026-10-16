################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##           Script S4: the matrix of ones in place of the nugget
##
##   Two kernels in this work are singular, the additive and the heterozygosity
##   kernel, each having lost one dimension to centring: K1 = 0. Scripts 03, 06
##   and S3 hand such a kernel to software that factorises it, lme4qtl and the
##   propto structure of glmmTMB. A nugget cI changes every eigenvalue and
##   moves the estimates; the matrix of ones, cJ with c any positive constant,
##   does not: (K + cJ)^-1 is a g-inverse of K, and J is removed from the
##   residual likelihood by the fixed intercept. The added term fills exactly
##   the dimension centring removed and lies in the column space of the
##   fixed-effect design. Under rotation W(K + cJ)W' = WKW' + c(W1)(W1)', and
##   W1 is the rotated intercept. Every error contrast annihilates it, so the
##   residual likelihood, the REML estimates, the deviance, the GLS estimates
##   and the standard errors of all fixed effects other than the intercept are
##   the same for every c and the same as for the singular kernel itself; only
##   the variance of the intercept estimate increases, by c times the kernel
##   variance. When the null space of K is spanned by 1, (K + cJ)^-1 = K^+ +
##   J / (c n^2), K^+ the Moore-Penrose inverse.
##
##   This is tested in three places.
##   1. The algebra on the observed kernels: the g-inverse identity, the closed
##      form of the inverse, the rotation of cJ, positive definiteness, and the
##      GLS identity on the observed trial: under T(A + cJ)T' the estimates and
##      their covariance equal those under TAT' except that the variance of the
##      intercept increases by exactly c times the kernel variance.
##   2. The multi-kernel comparison of S3 for c = 0.01, 1, 100 on the ten
##      simulated datasets of S3. Both implementations are given K + cJ; the
##      reference is the directly coded likelihood on the singular kernel,
##      which needs no factorisation. A control adds cJ after rotation, not
##      before, so that the added direction is 1 and not W1; that direction is
##      not in the fixed-effect space and c should then matter.
##   3. The scheme comparison of 06: its markers under the three
##      transformations with T(A + cJ)T' for the same three c, against the
##      nugget and against the singular kernel through the eigen route of
##      lme4qtl. The metrics of Table 2 and the anchor rows of Table 1 are
##      recomputed and set beside the stored ones.
##
##   Part 2 also evaluates both deviance functions at both solutions and along
##   the residual towards zero, and sets both against a third form of the same
##   deviance, through error contrasts, log|K'HK| + y'K(K'HK)^-1 K'y + log|X'X|
##   with K an orthonormal basis of the complement of X, which never meets the
##   singular direction. Where every kernel in the model is singular in the
##   same direction, the full-form evaluation on the singular kernel breaks
##   down as the residual approaches zero: log|H| and log|X'H^-1 X| diverge in
##   opposite directions and their sum, and the quadratic form, are lost to
##   rounding. With cJ the covariance matrix stays positive definite at a zero
##   residual and the deviance has a finite limit there.
##
##   Away from the boundary the K + cJ fits must reproduce the fits stored by
##   S3, and the c = 1 scheme of part 3 must reproduce the results stored by
##   06. The singular-kernel reference is fitted by a copy of the optimiser
##   that rotates the kernels as they are and floors the residual at 1e-6, so
##   this script does not depend on the ve_kernels of 07.
##
##   Output: S4_ones_matrix.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(lme4qtl); library(pbkrtest)
                  library(emmeans); library(rrBLUP); library(glmmTMB)})

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "02_rotation.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e03 <- new.env(); load(file.path(OUT, "03_equivalence.RData"), envir = e03)
e06 <- new.env(); load(file.path(OUT, "06_comparison.RData"), envir = e06)
eS3 <- new.env(); load(file.path(OUT, "S3_optimiser.RData"), envir = eS3)
sel <- e06$sel; ss <- e06$ss
t0 <- Sys.time()
stamp <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

C_SET   <- c(0.01, 1, 100)   # the kernels have a mean diagonal of about 1.8
NREP_S4 <- 5                 # the replicates of S3, same seeds
S2_TRUE <- 25
VE_GRID <- 10^seq(-1, -12, by = -1)
VE_SAFE <- 1e-3   # down to this residual the full-form evaluation on the singular kernel is
                  # accurate to TOL_EXACT
TOL_VE  <- 1e-3   # invariance on the residual, absolute: the optimiser noise
TOL_REL <- 1e-3   # on the components, relative
TOL_DEV <- 1e-3   # on the REML deviance between two optimisations
TOL_EXACT <- 1e-5 # between two evaluations of the same deviance at the same parameters
TOL_LP  <- 1e-4   # on -log10 p and the standard error, the precision reported
TOL_S3  <- 0.05   # the agreement criterion of S3 for the two implementations

sym <- function(M) (M + t(M)) / 2
ids <- sprintf("r%03d", seq_len(n))
J   <- matrix(1, n, n)
W   <- Tlist$rotated
w1  <- as.numeric(W %*% rep(1, n))

## ---------- 1. the algebra on the observed kernels ----------
pinv <- function(K) {
  e <- eigen(K, symmetric = TRUE); keep <- e$values > 1e-8 * max(e$values)
  e$vectors[, keep] %*% (t(e$vectors[, keep]) / e$values[keep])
}
KER <- list(additive = A, heterozygote = D_dom)
alg <- do.call(rbind, lapply(names(KER), function(nm) do.call(rbind, lapply(C_SET, function(cc) {
  K <- KER[[nm]]; Kc <- K + cc * J; Kci <- solve(Kc)
  Kr <- sym(W %*% Kc %*% t(W))
  data.frame(kernel = nm, c = cc,
             rank_K      = { ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values; sum(ev > 1e-8 * max(ev)) },
             K1          = max(abs(K %*% rep(1, n))),
             g_inverse   = max(abs(K %*% Kci %*% K - K)),
             closed_form = max(abs(Kci - (pinv(K) + J / (cc * n^2)))),
             rotation    = max(abs(Kr - (sym(W %*% K %*% t(W)) + cc * tcrossprod(w1)))),
             min_eig     = min(eigen(Kc, symmetric = TRUE, only.values = TRUE)$values),
             min_eig_rot = min(eigen(Kr, symmetric = TRUE, only.values = TRUE)$values),
             row.names = NULL)
}))))

## The GLS identity behind the invariance, on the observed trial. With the
## variance components of the rotated model at the test marker (script 03),
## (X'V_c^-1 X)^-1 = (X'V^-1 X)^-1 + U where U carries c times the kernel
## variance in the intercept cell and nothing else, and the GLS estimates are
## the same. Checked under all three transformations, whose intercept is T1.
vg <- e03$inv_tab$Vg[1]; ve <- e03$inv_tab$Ve[1]
gls <- function(Tm, cc) {
  X <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q, as.numeric(Tm %*% G[, TEST_SNP]))
  V <- vg * sym(Tm %*% (A + cc * J) %*% t(Tm)) + ve * diag(n)
  Vi <- chol2inv(chol(V)); C <- solve(crossprod(X, Vi %*% X))
  list(beta = as.numeric(C %*% crossprod(X, Vi %*% as.numeric(Tm %*% blues$BLUE))), C = C)
}
gls_tab <- do.call(rbind, lapply(names(Tlist), function(nm) do.call(rbind, lapply(C_SET, function(cc) {
  g0 <- gls(Tlist[[nm]], 0); g1 <- gls(Tlist[[nm]], cc); D <- g1$C - g0$C
  data.frame(scheme = nm, c = cc,
             d_beta       = max(abs(g1$beta - g0$beta)),
             d_cov        = max(abs(D[-1, -1])),          # every cell not involving the intercept
             d_cov_int    = max(abs(D[1, -1])),           # the intercept's covariances with the rest
             int_var_gain = D[1, 1], expected = cc * vg,
             se_marker    = sqrt(g1$C[5, 5]), row.names = NULL)
}))))
stamp("part 1 done")

## ---------- 2. the multi-kernel comparison of S3 under cJ ----------
P <- list(pheno = pheno, id_idx = ID_IDX, yr_idx = YR_IDX, yr_rep = YR_REP,
          id_yr = ID_YR, pars = PARS, rootA = ROOT_A, rootE = ROOT_E,
          rootD = ROOT_D, Gsub = G_SUB, Q = Q, A = A, n = n, L = L_mat)
STRUCT <- data.frame(source = c("epi", "het"), k = 1:2, stringsAsFactors = FALSE)

propto_component <- function(fit, grp, K) {
  vb <- VarCorr(fit)$cond[[grp]]
  s  <- attr(vb, "stddev"); cr <- attr(vb, "correlation")
  mean((diag(s) %*% cr %*% diag(s)) / K)
}
## the REML deviance as ve_kernels evaluates it, for two given rotated kernels
make_dev <- function(y, Xn, Ka, Kb) function(va, vb, ve) {
  H <- diag(length(y)) * ve + Ka * va + Kb * vb
  val <- try({
    R  <- chol(H); Hi <- chol2inv(R)
    XtHiX <- crossprod(Xn, Hi %*% Xn); Rx <- chol(XtHiX)
    Pm <- Hi - Hi %*% Xn %*% chol2inv(Rx) %*% crossprod(Xn, Hi)
    2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) + as.numeric(crossprod(y, Pm %*% y))
  }, silent = TRUE)
  if (inherits(val, "try-error") || !is.finite(val)) NA_real_ else val
}
## the same deviance through error contrasts, |H| |X'H^-1 X| = |K'HK| |X'X|
## and y'Py = y'K(K'HK)^-1 K'y for K an orthonormal basis of the complement of X
make_dev_ec <- function(y, Xn) {
  Kc  <- qr.Q(qr(Xn), complete = TRUE)[, -seq_len(ncol(Xn)), drop = FALSE]
  lxx <- as.numeric(determinant(crossprod(Xn), logarithm = TRUE)$modulus)
  z   <- crossprod(Kc, y)
  function(H) {
    R <- chol(crossprod(Kc, H %*% Kc))
    2 * sum(log(diag(R))) + as.numeric(crossprod(z, chol2inv(R) %*% z)) + lxx
  }
}

## The optimiser for the singular-kernel reference: kernels rotated as they
## are and the residual floored at 1e-6, kept here in full so that this script
## does not depend on the ve_kernels of 07, which adds the matrix of ones.
ve_kernels_old <- function(s1, Tm, kernels, P = NULL, start = NULL) {
  if (is.null(P)) P <- list(Q = Q, n = n)
  VE_FLOOR <- 1e-6
  sym <- function(M) (M + t(M)) / 2
  nn <- P$n
  y  <- as.numeric(Tm %*% s1$mu)
  Xn <- cbind(as.numeric(Tm %*% rep(1, nn)), Tm %*% P$Q)
  Ks <- lapply(kernels, function(M) sym(Tm %*% M %*% t(Tm)))
  nk <- length(Ks)
  f <- function(par) {
    ve <- exp(par[nk + 1])
    if (ve < VE_FLOOR) return(1e10)
    H <- diag(nn) * ve
    for (k in seq_len(nk)) H <- H + Ks[[k]] * exp(par[k])
    val <- try({
      R     <- chol(H)
      Hi    <- chol2inv(R)
      XtHiX <- crossprod(Xn, Hi %*% Xn)
      Rx    <- chol(XtHiX)
      Pm    <- Hi - Hi %*% Xn %*% chol2inv(Rx) %*% crossprod(Xn, Hi)
      2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) +
        as.numeric(crossprod(y, Pm %*% y))
    }, silent = TRUE)
    if (inherits(val, "try-error") || !is.finite(val)) return(1e10)
    val
  }
  if (is.null(start)) start <- c(50, rep(10, max(0, nk - 1)), 1)
  par <- log(start); dev <- Inf; conv <- NA_integer_
  for (i in 1:5) {
    o <- optim(par, f, method = "Nelder-Mead",
               control = list(maxit = 2000, reltol = 1e-10))
    improved <- dev - o$value; par <- o$par; dev <- o$value; conv <- o$convergence
    if (improved < 1e-9) break
  }
  vc <- exp(par)
  setNames(c(vc, as.numeric(conv == 0), dev),
           c(paste0("V_", names(kernels)), "Ve", "converged", "deviance"))
}

one <- function(r, i) {
  src <- STRUCT$source[i]
  set.seed(SEED_MISC + 5 + r * 10 + STRUCT$k[i])       # the datasets of S3
  b  <- setNames(rep(0, length(MARKERS)), MARKERS)
  s1 <- stage1_fast(simulate_y(b, s2_epi = if (src == "epi") S2_TRUE else 0,
                               s2_dom = if (src == "het") S2_TRUE else 0, P = P), P = P)
  Tr  <- make_transforms(s1$Omega)$rotated
  Ti1 <- as.numeric(solve(Tr) %*% rep(1, n))   # rotating its outer product gives J itself
  K2  <- if (src == "het") D_dom else E_had
  singular2 <- src == "het"

  raw <- ve_kernels_old(s1, Tr, list(add = A, extra = K2), P = P)

  y  <- as.numeric(Tr %*% s1$mu); ic <- as.numeric(Tr %*% rep(1, n)); Qr <- Tr %*% P$Q
  Xn <- cbind(ic, Qr)
  df <- data.frame(y = y, ic = ic, PC1 = Qr[, 1], PC2 = Qr[, 2], PC3 = Qr[, 3],
                   Id = factor(ids, levels = ids), g = factor(1), g2 = factor(1))
  Ka0 <- sym(Tr %*% A %*% t(Tr)); Kb0 <- sym(Tr %*% K2 %*% t(Tr))
  f0  <- make_dev(y, Xn, Ka0, Kb0)              # the singular-kernel deviance, full form
  fec <- make_dev_ec(y, Xn)                     # any kernel, through error contrasts

  tab <- do.call(rbind, lapply(C_SET, function(cc) {
    Ka_u <- A + cc * J
    Kb_u <- if (singular2) K2 + cc * J else K2
    mine <- ve_kernels_old(s1, Tr, list(add = Ka_u, extra = Kb_u), P = P)
    ctrl <- ve_kernels_old(s1, Tr, list(add = A + cc * tcrossprod(Ti1), extra = K2), P = P)

    Ka <- sym(Tr %*% Ka_u %*% t(Tr)); Kb <- sym(Tr %*% Kb_u %*% t(Tr))
    dimnames(Ka) <- dimnames(Kb) <- list(ids, ids)
    f1 <- make_dev(y, Xn, Ka, Kb)               # the K + cJ deviance, full form
    warn <- character(0)
    gt <- withCallingHandlers(
      tryCatch(glmmTMB(y ~ 0 + ic + PC1 + PC2 + PC3 +
                         propto(0 + Id | g, Ka) + propto(0 + Id | g2, Kb),
                       data = df, REML = TRUE), error = function(e) NULL),
      warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
    Va_t <- if (is.null(gt)) NA_real_ else propto_component(gt, "g",  Ka)
    V2_t <- if (is.null(gt)) NA_real_ else propto_component(gt, "g2", Kb)
    Ve_t <- if (is.null(gt)) NA_real_ else sigma(gt)^2
    stamp(sprintf("  rep %d %s c=%g  raw Ve %.5f  cJ Ve %.5f  tmb Ve %.5f  ctrl Ve %.5f",
                  r, src, cc, raw[["Ve"]], mine[["Ve"]], Ve_t, ctrl[["Ve"]]))
    data.frame(rep = r, source = src, c = cc,
               Ve_raw = raw[["Ve"]], Ve_mine = mine[["Ve"]], Ve_tmb = Ve_t, Ve_ctrl = ctrl[["Ve"]],
               Va_raw = raw[["V_add"]], Va_mine = mine[["V_add"]], Va_tmb = Va_t, Va_ctrl = ctrl[["V_add"]],
               V2_raw = raw[["V_extra"]], V2_mine = mine[["V_extra"]], V2_tmb = V2_t, V2_ctrl = ctrl[["V_extra"]],
               dev_raw = raw[["deviance"]], dev_mine = mine[["deviance"]],
               dev_tmb = if (is.null(gt)) NA_real_ else f1(Va_t, V2_t, Ve_t),
               ## both deviance functions at the other solution
               dev_cJ_at_raw = f1(raw[["V_add"]], raw[["V_extra"]], raw[["Ve"]]),
               dev_raw_at_cJ = f0(mine[["V_add"]], mine[["V_extra"]], mine[["Ve"]]),
               ## the error-contrast form at the cJ solution, evaluated on the
               ## singular kernel
               dev_ec_at_cJ  = fec(diag(n) * mine[["Ve"]] + Ka0 * mine[["V_add"]] + Kb0 * mine[["V_extra"]]),
               conv_mine = mine[["converged"]], conv_ctrl = ctrl[["converged"]],
               tmb_ok = as.numeric(!is.null(gt)),
               tmb_warn = if (length(warn)) paste(unique(warn), collapse = " | ") else "",
               row.names = NULL)
  }))

  ## all three evaluations along the residual, at the c = 1 solution
  m1  <- tab[tab$c == 1, ]
  Ka1 <- sym(Tr %*% (A + J) %*% t(Tr))
  Kb1 <- if (singular2) sym(Tr %*% (K2 + J) %*% t(Tr)) else Kb0
  f1  <- make_dev(y, Xn, Ka1, Kb1)
  H0  <- function(v) diag(n) * v + Ka0 * m1$Va_mine + Kb0 * m1$V2_mine
  H1  <- function(v) diag(n) * v + Ka1 * m1$Va_mine + Kb1 * m1$V2_mine
  path <- data.frame(rep = r, source = src, Ve = VE_GRID,
                     dev_raw    = sapply(VE_GRID, function(v) f0(m1$Va_mine, m1$V2_mine, v)),
                     dev_cJ     = sapply(VE_GRID, function(v) f1(m1$Va_mine, m1$V2_mine, v)),
                     dev_ec_raw = sapply(VE_GRID, function(v) fec(H0(v))),
                     dev_ec_cJ  = sapply(VE_GRID, function(v) fec(H1(v))),
                     row.names = NULL)
  list(tab = tab, path = path)
}

stamp("part 2:", NREP_S4, "replicates x", nrow(STRUCT), "structures x", length(C_SET), "values of c")
res  <- unlist(lapply(seq_len(NREP_S4), function(r)
          lapply(seq_len(nrow(STRUCT)), function(i) one(r, i))), recursive = FALSE)
S4   <- do.call(rbind, lapply(res, `[[`, "tab"))
path <- do.call(rbind, lapply(res, `[[`, "path"))

## the K + cJ fits at c = 1, with the floored optimiser, against the fits
## stored by S3, which uses the optimiser of 07; away from the boundary the two
## agree within the optimiser tolerance, so the boundary rows are left out
s3 <- merge(S4[S4$c == 1, c("rep", "source", "Ve_mine", "V2_mine")],
            eS3$S3[, c("rep", "source", "Ve_mine", "V2_mine")],
            by = c("rep", "source"), suffixes = c("", "_S3"))
s3 <- s3[s3$Ve_mine >= 1e-3 & s3$Ve_mine_S3 >= 1e-3, ]
d_S3 <- max(abs(s3$Ve_mine - s3$Ve_mine_S3), abs(s3$V2_mine - s3$V2_mine_S3) / pmax(s3$V2_mine_S3, 1))

## the directly coded implementation, cJ against the singular kernel
S4$boundary <- S4$Ve_raw < VE_SAFE | S4$Ve_mine < VE_SAFE
S4$d_Ve  <- abs(S4$Ve_mine - S4$Ve_raw)
S4$d_Va  <- abs(S4$Va_mine - S4$Va_raw) / pmax(S4$Va_raw, 1)
S4$d_V2  <- abs(S4$V2_mine - S4$V2_raw) / pmax(S4$V2_raw, 1)
S4$d_dev <- abs(S4$dev_mine - S4$dev_raw)
## the two full-form functions at the same parameters, the two solutions on the
## cJ function, and the cJ function against the error-contrast form
S4$d_fun <- abs(S4$dev_cJ_at_raw - S4$dev_raw)   # at the raw solution
S4$d_opt <- S4$dev_mine - S4$dev_cJ_at_raw       # negative or zero: cJ found at least as good a point
S4$d_ec  <- abs(S4$dev_mine - S4$dev_ec_at_cJ)   # at the cJ solution
## the two implementations on the same kernel, the S3 criterion
S4$matched <- abs(S4$Ve_tmb - S4$Ve_mine) < TOL_S3 &
              abs(S4$Va_tmb - S4$Va_mine) / pmax(S4$Va_mine, 1) < TOL_S3 &
              abs(S4$V2_tmb - S4$V2_mine) / pmax(S4$V2_mine, 1) < TOL_S3
S4$dev_gap <- S4$dev_mine - S4$dev_tmb
S4$mine_better <- !S4$matched & S4$dev_gap < 0
mm <- S4[S4$matched, ]
agree <- data.frame(
  Quantity = c("residual, absolute", "additive component, relative", "second component, relative"),
  Max_matched = c(max(abs(mm$Ve_tmb - mm$Ve_mine)),
                  max(abs(mm$Va_tmb - mm$Va_mine) / pmax(mm$Va_mine, 1)),
                  max(abs(mm$V2_tmb - mm$V2_mine) / pmax(mm$V2_mine, 1))), row.names = NULL)

## invariance to c, by replicate: the spread of each implementation across c
spread <- function(x, rel = FALSE) (max(x) - min(x)) / if (rel) max(1, max(x)) else 1
key <- interaction(S4$rep, S4$source)
inv <- do.call(rbind, lapply(split(S4, key), function(d) data.frame(
  rep = d$rep[1], source = d$source[1], all_matched = all(d$matched),
  mine_Ve = spread(d$Ve_mine), mine_Va = spread(d$Va_mine, TRUE), mine_V2 = spread(d$V2_mine, TRUE),
  tmb_Ve  = spread(d$Ve_tmb),  tmb_Va  = spread(d$Va_tmb, TRUE),  tmb_V2  = spread(d$V2_tmb, TRUE),
  ctrl_Ve = spread(d$Ve_ctrl), ctrl_Va = spread(d$Va_ctrl, TRUE), ctrl_V2 = spread(d$V2_ctrl, TRUE),
  row.names = NULL)))

## the breakdown of the full-form evaluation on the singular kernel towards a
## zero residual: the largest discrepancy from the error-contrast form, by
## structure and residual, beside the cJ full form and the invariance of the
## error-contrast form itself. With the epistatic structure only the additive
## kernel is singular and the epistatic kernel covers the direction; with the
## heterozygote structure both kernels are singular in the same direction.
path$d_raw <- abs(path$dev_raw - path$dev_ec_raw)
path$d_cJ  <- abs(path$dev_cJ  - path$dev_ec_cJ)
path$d_ec  <- abs(path$dev_ec_raw - path$dev_ec_cJ)
brk <- do.call(rbind, lapply(split(path, list(path$source, path$Ve)), function(d) data.frame(
  source = d$source[1], Ve = d$Ve[1], raw_vs_ec = max(d$d_raw), cJ_vs_ec = max(d$d_cJ),
  ec_raw_vs_ec_cJ = max(d$d_ec), row.names = NULL)))
brk <- brk[order(brk$source, -brk$Ve), ]
stamp("part 2 done")

## ---------- 3. the scheme comparison of 06 under cJ ----------
converged   <- function(m) length(m@optinfo$conv$lme4$messages) == 0
scaled_grad <- function(m) {
  d <- m@optinfo$derivs
  if (is.null(d)) NA_real_ else max(abs(solve(d$Hessian, d$gradient)))
}
BOBYQA <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
## the refit rule of 06: bobyqa when lme4 flags the fit or when the scaled
## gradient exceeds the tolerance
GRAD_TOL <- 1e-3
fit_relmat <- function(form, data, kernel) {
  m <- relmatLmer(form, data = data, relmat = list(Id = kernel), REML = TRUE)
  refit <- !converged(m) || isTRUE(scaled_grad(m) > GRAD_TOL)
  if (refit)
    m <- relmatLmer(form, data = data, relmat = list(Id = kernel), REML = TRUE, control = BOBYQA)
  list(m = m, refit = refit, ok = converged(m), grad = scaled_grad(m))
}

arm_fit <- function(Tm, Kn, label) {
  dimnames(Kn) <- list(ids, ids)
  base <- data.frame(y   = as.numeric(Tm %*% blues$BLUE),
                     ic  = as.numeric(Tm %*% rep(1, n)),
                     PC1 = (Tm %*% Q)[, 1], PC2 = (Tm %*% Q)[, 2], PC3 = (Tm %*% Q)[, 3],
                     Id  = factor(ids, levels = ids))
  f0 <- fit_relmat(y ~ 0 + ic + PC1 + PC2 + PC3 + (1|Id), base, Kn)
  m0 <- f0$m
  G_sel <- G[, sel, drop = FALSE]      # only these columns go to the cores
  res <- do.call(rbind, parallel::parLapply(cl, sel, function(s) {
    out <- try({
      d <- base; d$SNP <- as.numeric(Tm %*% G_sel[, s])
      ff <- fit_relmat(y ~ 0 + ic + PC1 + PC2 + PC3 + SNP + (1|Id), d, Kn)
      kr <- KRmodcomp(ff$m, m0)$stats
      cf <- summary(ff$m)$coefficients["SNP", ]
      data.frame(SNP = s, beta = cf[["Estimate"]], se = cf[["Std. Error"]],
                 df = kr$ddf, neglog10p = -log10(kr$p.value),
                 converged = ff$ok, refit = ff$refit, grad = ff$grad, row.names = NULL)
    }, silent = TRUE)
    if (inherits(out, "try-error"))
      data.frame(SNP = s, beta = NA_real_, se = NA_real_, df = NA_real_, neglog10p = NA_real_,
                 converged = FALSE, refit = NA, grad = NA_real_, row.names = NULL)
    else out
  }))
  cbind(kernel = label, res, conv_null = f0$ok)
}

## which route lme4qtl takes for the singular kernel: its Cholesky, or its
## eigen route when that fails
chol_route <- function(K) !inherits(try(Matrix::chol(K), silent = TRUE), "try-error")

stamp("part 3:", length(sel), "markers x 3 schemes x 5 kernels")
## the marker fits of each scheme run on several cores, as in 06; the number
## of cores changes the timing, not the result
ncore <- max(1, min(parallel::detectCores() - 1, 6))
cl <- parallel::makeCluster(ncore)
invisible(parallel::clusterEvalQ(cl, suppressMessages({library(lme4); library(lme4qtl); library(pbkrtest)})))
parallel::clusterExport(cl, c("fit_relmat", "converged", "scaled_grad", "BOBYQA", "GRAD_TOL"))
arms <- lapply(names(Tlist), function(nm) {
  Tm <- Tlist[[nm]]; K0 <- sym(Tm %*% A %*% t(Tm))
  ## the ones term is added after rotation, as with_ones in script 00 does, so
  ## that the c = 1 kernel is the one 06 fits
  t1 <- as.numeric(Tm %*% rep(1, n))
  kern <- c(list(nugget = regularise(K0), raw = K0),
            setNames(lapply(C_SET, function(cc) K0 + cc * tcrossprod(t1)),
                     paste0("ones_", C_SET)))
  out <- do.call(rbind, lapply(names(kern), function(k) {
    stamp("  ", nm, k); arm_fit(Tm, kern[[k]], k) }))
  cbind(scheme = nm, out, raw_chol = chol_route(K0))
})
arms <- do.call(rbind, arms)
parallel::stopCluster(cl)

## the c = 1 kernel against the results stored by 06, marker by marker
rep06 <- do.call(rbind, lapply(names(Tlist), function(nm) {
  a <- arms[arms$scheme == nm & arms$kernel == "ones_1", ]; b <- e06$ts[[nm]]
  a <- a[match(b$SNP, a$SNP), ]
  data.frame(scheme = nm, beta = max(abs(a$beta - b$beta)), se = max(abs(a$se - b$se)),
             df = max(abs(a$df - b$df)), neglog10p = max(abs(a$neglog10p - b$neglog10p)),
             row.names = NULL) }))

## the spread across c, and the distances to the singular kernel and the nugget
QTY <- c("beta", "se", "df", "neglog10p")
wide <- function(sch, k, q) { d <- arms[arms$scheme == sch & arms$kernel == k, ]; d[[q]][match(sel, d$SNP)] }
inv06 <- do.call(rbind, lapply(names(Tlist), function(sch) do.call(rbind, lapply(QTY, function(q) {
  M <- sapply(paste0("ones_", C_SET), function(k) wide(sch, k, q))
  data.frame(scheme = sch, quantity = q,
             spread_c       = max(apply(M, 1, function(v) max(v) - min(v))),
             ones_vs_raw    = max(abs(M[, "ones_1"] - wide(sch, "raw", q))),
             ones_vs_nugget = max(abs(M[, "ones_1"] - wide(sch, "nugget", q))),
             row.names = NULL)
}))))

## Table 2 and the anchor rows, recomputed under cJ (c = 1) beside the stored
## values of 06
metric_one <- function(two, quantity) {
  m <- merge(ss[, c("SNP", quantity)], two[, c("SNP", quantity)], by = "SNP", suffixes = c("_ss", "_ts"))
  a <- m[[paste0(quantity, "_ss")]]; b <- m[[paste0(quantity, "_ts")]]
  keep <- is.finite(a) & is.finite(b); a <- a[keep]; b <- b[keep]; d <- b - a
  se_bias <- sd(d) / sqrt(length(d))
  data.frame(Quantity = quantity, n = length(d), MSD = mean(d^2), Bias = mean(d),
             SE_bias = se_bias, t_bias = mean(d) / se_bias, Correlation = cor(b, a))
}
metrics_cJ <- do.call(rbind, lapply(names(Tlist), function(nm)
  cbind(Scheme = nm, do.call(rbind, lapply(c("beta", "se", "neglog10p"), function(q)
    metric_one(arms[arms$scheme == nm & arms$kernel == "ones_1", ], q))))))
num <- c("MSD", "Bias", "SE_bias", "t_bias", "Correlation")
metrics_cJ[, num] <- round(metrics_cJ[, num], 4)
metrics_cJ <- metrics_cJ[order(metrics_cJ$Quantity, metrics_cJ$Scheme), ]
m06 <- e06$metrics[order(e06$metrics$Quantity, e06$metrics$Scheme), ]
metrics_side <- data.frame(m06[, c("Scheme", "Quantity")],
                           Bias_06 = m06$Bias, Bias_cJ = metrics_cJ$Bias,
                           t_06 = m06$t_bias, t_cJ = metrics_cJ$t_bias,
                           MSD_06 = m06$MSD, MSD_cJ = metrics_cJ$MSD, row.names = NULL)
metrics_changed <- sum(metrics_side$Bias_06 != metrics_side$Bias_cJ | metrics_side$t_06 != metrics_side$t_cJ |
                       metrics_side$MSD_06 != metrics_side$MSD_cJ)

pick <- function(df, col) df[[col]][df$SNP == TEST_SNP]
anchor_cJ <- rbind(
  data.frame(Analysis = "single-stage (plot level)", beta = pick(ss, "beta"), se = pick(ss, "se"),
             neglog10p = pick(ss, "neglog10p")),
  do.call(rbind, lapply(names(Tlist), function(nm) {
    d <- arms[arms$scheme == nm & arms$kernel == "ones_1", ]
    data.frame(Analysis = paste("two-stage,", nm), beta = pick(d, "beta"), se = pick(d, "se"),
               neglog10p = pick(d, "neglog10p")) })))
anchor_cJ$Distance <- round(abs(anchor_cJ$beta - anchor_cJ$beta[1]), 4)
anchor_cJ$SE_ratio <- round(anchor_cJ$se / anchor_cJ$se[1], 4)
anchor_cJ[, 2:4]   <- round(anchor_cJ[, 2:4], 4)
anchor_changed <- sum(as.matrix(anchor_cJ[, -1]) != as.matrix(e06$anchor[, -1]))
stamp("part 3 done")

## ---------- checks ----------
off  <- S4[!S4$boundary, ]
safe <- path[path$Ve >= VE_SAFE, ]
## The breakdown is checked on the largest discrepancy at every residual from
## 1e-7 downwards, and on every dataset from 1e-8 downwards.
deep <- path[path$Ve <= 1e-7 & path$source == "het", ]
deep_broken <- all(tapply(deep$d_raw, deep$Ve, max) > 1) && min(deep$d_raw[deep$Ve <= 1e-8]) > 1
checks <- data.frame(
  Check = c("(K + cJ)^-1 is a g-inverse of K on both singular kernels",
            "(K + cJ)^-1 equals K^+ + J / (c n^2)",
            "W(K + cJ)W' equals WKW' + c(W1)(W1)'",
            "K + cJ and its rotation are positive definite for every c",
            "GLS estimates and every covariance but the intercept variance are invariant to c, all three schemes",
            "the intercept variance grows by exactly c times the kernel variance",
            "away from the boundary, the K + cJ fits reproduce the fits stored by S3",
            "every glmmTMB fit on K + cJ returned a result",
            "every directly coded fit converged",
            "directly coded estimates on K + cJ identical across c",
            "away from the boundary, the K + cJ fit reproduces the singular-kernel fit",
            "the two full-form functions agree at the same parameters where the residual is above 1e-3",
            "on the K + cJ likelihood its solution is at least as good as the singular-kernel solution",
            "the K + cJ full form equals the error-contrast form at every cJ solution",
            "the error-contrast deviance is identical for K and K + cJ along the whole path",
            "the K + cJ full form equals the error-contrast form down to a residual of 1e-12",
            "with both kernels singular the singular-kernel full form breaks down below a residual of 1e-7",
            "with one kernel positive definite the singular-kernel full form holds at every residual",
            "glmmTMB estimates are invariant to c where the engines agree",
            "the control, ones added after rotation, depends on c",
            "every lme4qtl fit converged",
            "the c = 1 arm reproduces the results stored by 06",
            "scheme comparison invariant to c at the reported precision",
            "scheme comparison on K + cJ equals the singular kernel through lme4qtl",
            "no missing values"),
  Result = c(max(alg$g_inverse) < 1e-8,
             max(alg$closed_form) < 1e-8,
             max(alg$rotation) < 1e-8,
             min(alg$min_eig, alg$min_eig_rot) > 0,
             max(gls_tab$d_beta, gls_tab$d_cov, gls_tab$d_cov_int) < 1e-8,
             max(abs(gls_tab$int_var_gain - gls_tab$expected) / gls_tab$expected) < 1e-8,
             d_S3 < 1e-3,
             all(S4$tmb_ok == 1),
             all(S4$conv_mine == 1),
             max(inv$mine_Ve, inv$mine_Va, inv$mine_V2) < 1e-8,
             max(off$d_Ve) < TOL_VE & max(off$d_Va, off$d_V2) < TOL_REL & max(off$d_dev) < TOL_DEV,
             max(S4$d_fun[S4$Ve_raw >= VE_SAFE]) < TOL_EXACT & max(safe$d_raw) < TOL_EXACT,
             max(S4$d_opt) < 1e-4,
             max(S4$d_ec) < TOL_EXACT,
             max(path$d_ec) < TOL_EXACT,
             max(path$d_cJ) < TOL_EXACT,
             deep_broken,
             max(path$d_raw[path$source == "epi"]) < TOL_EXACT,
             { ok <- inv[inv$all_matched, ]; max(ok$tmb_Ve) < TOL_VE & max(ok$tmb_Va, ok$tmb_V2) < TOL_REL },
             max(inv$ctrl_Ve) > TOL_VE,
             all(arms$converged) & all(arms$conv_null),
             max(rep06[, -1]) < 1e-6,
             max(inv06$spread_c[inv06$quantity %in% c("se", "neglog10p")]) < TOL_LP,
             max(inv06$ones_vs_raw[inv06$quantity %in% c("se", "neglog10p")]) < TOL_LP,
             !anyNA(S4[, c("Ve_mine", "Ve_tmb", "V2_mine", "V2_tmb", "dev_tmb", "dev_cJ_at_raw",
                           "dev_raw_at_cJ", "dev_ec_at_cJ")]) &
               !anyNA(arms$neglog10p) & !anyNA(path[, c("d_raw", "d_cJ", "d_ec")])))

## ---------- report ----------
cat("\n=== S4 the ones matrix in place of the nugget ===\n")
cat("\n1. algebra on the observed kernels\n"); print(alg, digits = 3, row.names = FALSE)
cat("\n   the GLS identity on the observed trial, test marker, Vg", signif(vg, 4), "Ve", signif(ve, 4), "\n")
print(gls_tab, digits = 4, row.names = FALSE)

cat("\n2. the multi-kernel comparison of S3, both engines on K + cJ\n")
cat("   K + cJ fits against those stored by S3, away from the boundary, largest difference:", signif(d_S3, 3), "\n")
print(S4[, c("rep", "source", "c", "Ve_raw", "Ve_mine", "Ve_tmb", "Ve_ctrl", "Va_raw", "Va_mine", "Va_tmb",
             "V2_raw", "V2_mine", "V2_tmb", "dev_raw", "dev_mine", "dev_tmb", "matched", "boundary")],
      digits = 6, row.names = FALSE)
cat("\n   directly coded engine, K + cJ against the singular kernel, maxima away from the boundary (",
    nrow(off), "of", nrow(S4), "fits )\n")
cat(sprintf("   residual %.2e absolute | additive %.2e relative | second %.2e relative | deviance %.2e\n",
            max(off$d_Ve), max(off$d_Va), max(off$d_V2), max(off$d_dev)))
cat("\n   the deviance functions at the same parameters\n")
print(S4[S4$c == 1, c("rep", "source", "Ve_raw", "dev_raw", "dev_cJ_at_raw", "d_fun", "dev_mine",
                      "dev_ec_at_cJ", "d_ec", "d_opt")], digits = 7, row.names = FALSE)
cat("   d_fun: |cJ full form - singular-kernel full form| at the singular-kernel solution\n")
cat("   d_ec : |cJ full form - error-contrast form| at the cJ solution\n")
cat("   d_opt: cJ solution minus singular-kernel solution, both on the cJ function; negative means cJ found the better point\n")
cat("\n   along the residual at the cJ solution: largest discrepancy from the error-contrast form, by structure\n")
print(brk, digits = 3, row.names = FALSE)
cat("\n   spread across c = 0.01, 1, 100 by replicate (absolute on the residual, relative on the components)\n")
print(inv, digits = 3, row.names = FALSE)
cat("\n   replicates where the two engines agree at every c:", sum(inv$all_matched), "of", nrow(inv), "\n")
print(agree, digits = 3, row.names = FALSE)
if (any(!S4$matched)) {
  cat("   where they parted company, scored on the same deviance\n")
  print(S4[!S4$matched, c("rep", "source", "c", "dev_mine", "dev_tmb", "dev_gap", "mine_better")],
        digits = 6, row.names = FALSE)
  cat("   glmmTMB reported:", unique(substr(S4$tmb_warn[!S4$matched], 1, 80)), "\n")
}

cat("\n3. the scheme comparison of 06 under T(A + cJ)T'\n")
cat("   lme4qtl route for the singular kernel, by scheme (TRUE = Cholesky succeeded):\n")
print(unique(arms[, c("scheme", "raw_chol")]), row.names = FALSE)
cat("   the c = 1 arm against the stored 06 results, largest differences\n"); print(rep06, digits = 3, row.names = FALSE)
cat("   maxima over the", length(sel), "markers\n"); print(inv06, digits = 3, row.names = FALSE)
cat("\n   Table 2 metrics, stored by 06 beside the c = 1 arm refitted here, both rounded as reported\n")
print(metrics_side, row.names = FALSE)
cat("   entries that differ at the reported precision:", metrics_changed, "of", 3 * nrow(metrics_side), "\n")
cat("\n   the test marker four ways, under cJ\n"); print(anchor_cJ, row.names = FALSE)
cat("   entries that differ from the stored anchor at the reported precision:", anchor_changed, "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

save(alg, gls_tab, S4, path, brk, inv, agree, d_S3, arms, rep06, inv06, metrics_cJ, metrics_side,
     metrics_changed, anchor_cJ, anchor_changed, checks, C_SET, NREP_S4, S2_TRUE, VE_GRID, VE_SAFE,
     TOL_VE, TOL_REL, TOL_DEV, TOL_EXACT, TOL_LP, TOL_S3,
     file = file.path(OUT, "S4_ones_matrix.RData"))
cat("\nsaved S4_ones_matrix.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
if (!all(checks$Result)) stop("S4 checks failed")
