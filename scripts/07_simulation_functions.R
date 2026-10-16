################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##                  Script 07: shared machinery for the simulations
##
##   Defines the functions and builds the objects used by scripts 08 to 11.
##   Runs no scenario itself.
##
##   Data are simulated on the real trial design, keeping the observed
##   accession-by-year-by-replicate incidence, the real relationship matrix and
##   the real principal components, and simulating only the response.
##
##   Omega is formed as L Var(beta) L', with L the contrast that averages the
##   fixed effects into the accession means; L depends only on the design and
##   is computed once. This is the unadjusted covariance; the analysis of the
##   observed data uses the Kenward-Roger adjusted form, which script 02 shows
##   to be larger by four parts in ten thousand on average.
##
##   The functions are saved in the results file that scripts 08 to 11 load.
##
##   Output: 07_sim_primitives.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(lme4)
library(lme4qtl)
library(emmeans)
library(rrBLUP)

load(file.path(OUT, "01_data.RData"))
t0 <- Sys.time()

## The kernel and the plot data must name the accessions identically.
stopifnot(identical(dimnames(A), list(rownames(A), rownames(A))),
          identical(levels(pheno$Id), rownames(A)))

converged <- function(m) length(m@optinfo$conv$lme4$messages) == 0

## ---------- 1. the averaging contrast, computed once ----------
form1  <- as.formula(paste(TRAIT, "~ Id + Year + (1|Year:Replication) + (1|Id:Year)"))
m_real <- lmer(form1, data = pheno, REML = TRUE)
L_mat  <- emmeans(m_real, ~ Id)@linfct
conv_real <- converged(m_real)

## ---------- 2. calibrate the generating model on the observed data ----------
## The fitting model carries the same terms as the generating model below,
## principal components included.
d0 <- pheno
d0$PC1 <- Q[as.character(d0$Id), 1]
d0$PC2 <- Q[as.character(d0$Id), 2]
d0$PC3 <- Q[as.character(d0$Id), 3]

form_ss <- as.formula(paste(TRAIT,
  "~ 1 + Year + PC1 + PC2 + PC3 + (1|Year:Replication) + (1|Id:Year) + (1|Id)"))
## relmatLmer factorises the kernel, and A is singular, so it receives A + cJ,
## J the matrix of ones (script 00); at plot level the added term lies in the
## space of the intercept and changes none of the calibrated values.
A_ss <- A + ONES_C * matrix(1, nrow(A), ncol(A))
m_ss <- relmatLmer(form_ss, data = d0, relmat = list(Id = A_ss), REML = TRUE)
conv_ss <- converged(m_ss)

vc <- as.data.frame(lme4::VarCorr(m_ss))
fx <- lme4::fixef(m_ss)

PARS <- list(
  mu      = as.numeric(fx[1]),
  s2_g    = vc$vcov[vc$grp == "Id"],
  s2_ge   = vc$vcov[vc$grp == "Id:Year"],
  s2_rep  = vc$vcov[vc$grp == "Year:Replication"],
  s2_e    = attr(lme4::VarCorr(m_ss), "sc")^2,
  year_fx = c(0, as.numeric(fx[grep("^Year", names(fx))])),
  pc_beta = as.numeric(fx[paste0("PC", 1:N_PC)])
)

## ---------- 3. marker selection ----------
## Causal markers are on separate chromosomes, and scripts 08 to 11 simulate
## one at a time; null markers are drawn from chromosomes carrying no causal
## marker. All are drawn above the minor allele frequency MAF_MIN, so that a
## simulated effect is carried by enough accessions to be estimated.
set.seed(SEED_MISC + 4)

poly   <- names(which(apply(G, 2, var) > 0))
pfreq  <- colMeans(G[, poly, drop = FALSE]) / 2
maf    <- pmin(pfreq, 1 - pfreq)
usable <- poly[maf > MAF_MIN]
chr_u  <- map$Chromosome[match(usable, map$SNP)]

chr_ok  <- names(which(table(chr_u) >= 5))
chr_cau <- sample(chr_ok, 6)
CAUSAL  <- unname(vapply(chr_cau, function(cc) sample(usable[chr_u == cc], 1), character(1)))
NULLS   <- sample(usable[!chr_u %in% chr_cau], 100)
MARKERS <- c(CAUSAL, NULLS)

chr_causal <- map$Chromosome[match(CAUSAL, map$SNP)]
chr_nulls  <- map$Chromosome[match(NULLS,  map$SNP)]

## ---------- 4. genetic structures and their matrix roots ----------
## A matrix root from the eigen decomposition, with negative eigenvalues set
## to zero, so that a singular kernel has one too.
root_of <- function(M) {
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)))
}

## Three structures: the additive relationship matrix, the additive-by-additive
## epistatic kernel as its Hadamard square, and the heterozygosity kernel. The
## two non-additive ones are scaled to the mean diagonal of A, as in script 04,
## so that a variance component of a given size adds the same total genetic
## variance whichever kernel carries it. Their correlation with A over the
## off-diagonal elements, reported in kernel_tab below, measures how much of
## each the additive kernel can absorb. Neither is a claim about a genetic mode
## of action: they are structures missing from the fitted model.
scale_to <- function(M, target) M * (target / mean(diag(M)))
mdA <- mean(diag(A))

E_had <- scale_to(A * A, mdA)

## The heterozygosity kernel, het_kernel in script 00, the same matrix as in
## script 04: only called genotypes enter, imputed cells contribute nothing.
hk            <- het_kernel(G, MISS, mdA)
D_dom         <- hk$K
n_het_markers <- hk$m
n_het_calls   <- hk$n_het_calls

## sommer::D.mat on the coding -1, 0, 1, computed only for its correlation with
## A and with the heterozygosity kernel.
D_som      <- sommer::D.mat(G - 1)
cor_dmat_A <- cor(D_som[upper.tri(A)], A[upper.tri(A)])
cor_dmat_H <- cor(D_som[upper.tri(A)], D_dom[upper.tri(A)])
rm(D_som)

## The Vitezica dominance coding, minus two q squared for AA, two p q for Aa
## and minus two p squared for aa, on called genotypes, a missing call taking
## the marker's mean coded value; only its correlation with A is kept.
p_vit <- colMeans(G) / 2
q_vit <- 1 - p_vit
V_cod <- matrix(NA_real_, nrow(G), ncol(G), dimnames = dimnames(G))
for (g in 0:2) {
  cell <- !MISS & G == g
  val  <- switch(g + 1, -2 * p_vit^2, 2 * p_vit * q_vit, -2 * q_vit^2)
  V_cod[cell] <- rep(val, each = nrow(G))[cell]
}
cm_v <- colMeans(V_cod, na.rm = TRUE)
ix_v <- which(is.na(V_cod), arr.ind = TRUE)
V_cod[ix_v] <- cm_v[ix_v[, 2]]
## centred before the cross-product, as the heterozygosity kernel is
V_cod <- scale(V_cod, center = TRUE, scale = FALSE)
V_vit <- scale_to(tcrossprod(V_cod) / ncol(V_cod), mdA)
cor_vit_A <- cor(A[upper.tri(A)], V_vit[upper.tri(V_vit)])
rm(V_cod, V_vit, cm_v, ix_v)

ROOT_A <- root_of(A)
ROOT_E <- root_of(E_had)
ROOT_D <- root_of(D_dom)

dev_root_A <- max(abs(tcrossprod(ROOT_A) - A))
dev_root_E <- max(abs(tcrossprod(ROOT_E) - E_had))
dev_root_D <- max(abs(tcrossprod(ROOT_D) - D_dom))

## how distinct the three structures are, measured on the off-diagonal elements
offd <- function(M) M[upper.tri(M)]
kernel_tab <- data.frame(
  Kernel       = c("additive", "epistasis A#A", "heterozygote"),
  Mean_diag    = round(c(mean(diag(A)), mean(diag(E_had)), mean(diag(D_dom))), 4),
  Eff_rank     = c(sum(eigen(A, symmetric = TRUE, only.values = TRUE)$values > 1e-6 * mdA),
                   sum(eigen(E_had, symmetric = TRUE, only.values = TRUE)$values > 1e-6 * mdA),
                   sum(eigen(D_dom, symmetric = TRUE, only.values = TRUE)$values > 1e-6 * mdA)),
  Cor_with_A   = round(c(1, cor(offd(A), offd(E_had)), cor(offd(A), offd(D_dom))), 4))

## ---------- 5. design index vectors ----------
YR_REP <- interaction(pheno$Year, pheno$Replication, drop = TRUE)
ID_YR  <- interaction(pheno$Id,   pheno$Year,        drop = TRUE)
ID_IDX <- as.integer(pheno$Id)
YR_IDX <- as.integer(pheno$Year)
G_SUB  <- G[, MARKERS, drop = FALSE]

## ---------- 6. simulate one plot-level response ----------
simulate_y <- function(beta_vec, s2_epi = 0, s2_dom = 0, s2_rep = NULL, P = NULL) {
  if (is.null(P)) P <- list(pheno = pheno, id_idx = ID_IDX, yr_idx = YR_IDX,
                            yr_rep = YR_REP, id_yr = ID_YR, pars = PARS,
                            rootA = ROOT_A, rootE = ROOT_E, rootD = ROOT_D,
                            Gsub = G_SUB, Q = Q, n = n)
  ## beta_vec must be named by marker
  stopifnot(!is.null(names(beta_vec)))

  p <- P$pars
  if (is.null(s2_rep)) s2_rep <- p$s2_rep
  nn <- P$n

  u   <- as.numeric(P$rootA %*% rnorm(nn)) * sqrt(p$s2_g)
  v   <- if (s2_epi > 0) as.numeric(P$rootE %*% rnorm(nn)) * sqrt(s2_epi) else rep(0, nn)
  w   <- if (s2_dom > 0) as.numeric(P$rootD %*% rnorm(nn)) * sqrt(s2_dom) else rep(0, nn)
  str <- as.numeric(P$Q %*% p$pc_beta)
  gsum <- u + v + w + str

  r_eff  <- rnorm(nlevels(P$yr_rep), 0, sqrt(s2_rep))
  ge_eff <- rnorm(nlevels(P$id_yr),  0, sqrt(p$s2_ge))

  mk <- names(beta_vec)[beta_vec != 0]
  snp <- if (length(mk)) as.numeric(P$Gsub[, mk, drop = FALSE] %*% beta_vec[mk]) else rep(0, nn)

  p$mu + p$year_fx[P$yr_idx] + gsum[P$id_idx] + snp[P$id_idx] +
    r_eff[as.integer(P$yr_rep)] + ge_eff[as.integer(P$id_yr)] +
    rnorm(nrow(P$pheno), 0, sqrt(p$s2_e))
}

## ---------- 7. first stage, fast path ----------
## Returns the adjusted means, Omega and four diagnostics: whether the fit is
## singular (the replicate-within-year variance is close to zero and often
## estimates to zero, which is not a failure), the gradient scaled by the
## Hessian, whether the fit was repeated, and whether lme4 gave any other
## message. A fit whose scaled gradient exceeds GRAD_TOL is refitted with the
## optimiser bobyqa.
stage1_fast <- function(y, P = NULL) {
  if (is.null(P)) P <- list(pheno = pheno, L = L_mat)
  GRAD_TOL <- 1e-3
  d <- P$pheno; d$ysim <- y

  form <- ysim ~ Id + Year + (1|Year:Replication) + (1|Id:Year)
  sg <- function(m) {
    dv <- m@optinfo$derivs
    if (is.null(dv)) NA_real_ else max(abs(solve(dv$Hessian, dv$gradient)))
  }

  m <- lmer(form, data = d, REML = TRUE)
  g <- sg(m)
  refit <- !is.na(g) && g > GRAD_TOL
  if (refit) {
    m <- lmer(form, data = d, REML = TRUE,
              control = lme4::lmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 1e5)))
    g <- sg(m)
  }

  list(mu       = as.numeric(P$L %*% lme4::fixef(m)),
       Omega    = as.matrix(P$L %*% as.matrix(vcov(m)) %*% t(P$L)),
       singular = lme4::isSingular(m),
       grad     = g,
       refit    = refit,
       ## any lme4 message other than the boundary one
       ok       = all(grepl("boundary", m@optinfo$conv$lme4$messages)))
}

## ---------- 8. scan the selected markers under one transformation ----------
## The scan of the selected markers, vectorised as in script 05: for a
## partitioned model the residual sum of squares after fitting the marker
## equals y'Py - b z'Py, with P removing the fixed effects under Hinv.
## mixed.solve tolerates a singular kernel, so the kernel is passed as it is.
fit_arm <- function(s1, Tm, P = NULL) {
  if (is.null(P)) P <- list(Q = Q, Gsub = G_SUB, A = A, n = n)
  nn <- P$n
  y  <- as.numeric(Tm %*% s1$mu)
  Xn <- cbind(as.numeric(Tm %*% rep(1, nn)), Tm %*% P$Q)
  Zr <- Tm %*% P$Gsub
  K  <- Tm %*% P$A %*% t(Tm)

  nf <- mixed.solve(y = y, X = Xn, Z = diag(nn), K = K, return.Hinv = TRUE)
  Hi <- nf$Hinv

  XtHiX <- crossprod(Xn, Hi %*% Xn)
  Pm    <- Hi - Hi %*% Xn %*% solve(XtHiX) %*% crossprod(Xn, Hi)

  PZ  <- Pm %*% Zr
  zPz <- colSums(Zr * PZ)
  zPy <- as.numeric(crossprod(PZ, y))
  yPy <- as.numeric(crossprod(y, Pm %*% y))

  p  <- ncol(Xn) + 1
  v2 <- nn - p
  ok <- zPz > 1e-8 * max(zPz)

  beta <- se <- lp <- lo <- hi <- rep(NA_real_, ncol(Zr))
  b  <- zPy[ok] / zPz[ok]
  s2 <- (yPy - b * zPy[ok]) / v2
  sb <- sqrt(s2 / zPz[ok])
  tq <- qt(0.975, v2)

  beta[ok] <- b
  se[ok]   <- sb
  lp[ok]   <- -log10(pbeta(v2 / (v2 + (b / sb)^2), v2 / 2, 1 / 2))
  lo[ok]   <- b - tq * sb
  hi[ok]   <- b + tq * sb

  data.frame(marker = colnames(Zr), beta = beta, se = se,
             neglog10p = lp, ci_lo = lo, ci_hi = hi,
             Vg = nf$Vu, Ve = nf$Ve, row.names = NULL)
}

## ---------- 9. rotated null model with an arbitrary set of kernels ----------
## Takes a named list of kernels on the original scale and rotates each, so the
## same function fits the additive model alone, the additive model plus the
## structure that generated the data, or the additive model plus a generic
## structure. diag(n) gives the independent genetic kernel, whose rotated form
## is W W', the inverse of Omega; that effect is identifiable when the second
## stage carries Omega, by rotation or by generalised least squares, while the
## unweighted analysis can estimate only its sum with the residual.
## A kernel that centring has made singular receives the matrix of ones before
## rotation, K + cJ with c = 1, as with_ones in script 00 does; the constant is
## set inside the function, which must be self-contained on the cores. The
## added term rotates to the outer product of the rotated intercept, so the
## residual likelihood, the estimates and the deviance are the same for every
## c. The determinants come from Cholesky factors, which fail if the
## covariance matrix is not positive definite.
ve_kernels <- function(s1, Tm, kernels, P = NULL, start = NULL) {
  if (is.null(P)) P <- list(Q = Q, n = n)
  ONES_C <- 1
  sym <- function(M) (M + t(M)) / 2

  nn <- P$n
  y  <- as.numeric(Tm %*% s1$mu)
  Xn <- cbind(as.numeric(Tm %*% rep(1, nn)), Tm %*% P$Q)
  Ks <- lapply(kernels, function(M) {
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < 1e-8 * max(ev)) M <- M + ONES_C * matrix(1, nn, nn)
    sym(Tm %*% M %*% t(Tm))
  })
  nk <- length(Ks)

  f <- function(par) {
    ve <- exp(par[nk + 1])
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

  ## Each chain is Nelder-Mead restarted until the deviance stops improving, as
  ## in 04. With two or more kernels the likelihood can have a local optimum at
  ## a zero residual beside a better interior one, so three further chains run
  ## from other starts; the best chain is polished by L-BFGS-B on the log
  ## scale, and the polish is kept only if it lowers the deviance. "converged"
  ## reports the Nelder-Mead status of the best chain.
  chain <- function(st) {
    par <- log(st); dev <- Inf; conv <- NA_integer_
    for (i in 1:5) {
      o <- optim(par, f, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-10))
      improved <- dev - o$value; par <- o$par; dev <- o$value; conv <- o$convergence
      if (improved < 1e-9) break
    }
    list(par = par, dev = dev, conv = conv)
  }
  starts <- list(if (is.null(start)) c(50, rep(10, max(0, nk - 1)), 1) else start)
  if (nk > 1) starts <- c(starts, list(c(50, rep(1, nk - 1), 1),
                                       c(50, rep(30, nk - 1), 0.3),
                                       c(20, rep(5, nk - 1), 3)))
  runs <- lapply(starts, chain)
  best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "dev"))]]
  pol <- try(optim(pmax(best$par, -14), f, method = "L-BFGS-B",
                   lower = rep(-14, nk + 1), upper = rep(8, nk + 1),
                   control = list(factr = 1e3, maxit = 2000)), silent = TRUE)
  if (!inherits(pol, "try-error") && is.finite(pol$value) && pol$value < best$dev)
    best <- list(par = pol$par, dev = pol$value, conv = best$conv)
  vc <- exp(best$par)
  setNames(c(vc, as.numeric(best$conv == 0), best$dev),
           c(paste0("V_", names(kernels)), "Ve", "converged", "deviance"))
}

## ---------- 10. checks ----------
checks <- data.frame(
  Check = c("first-stage model for the contrast converged",
            "calibration model converged",
            "all calibrated variance components positive",
            "causal markers on distinct chromosomes",
            "no null shares a chromosome with a causal marker",
            "every causal and null marker is above the MAF floor",
            "eigen roots reproduce their matrices",
            "averaging contrast has one row per accession",
            "non-additive kernels rescaled to the mean diagonal of A",
            "epistatic kernel carries structure distinct from A",
            "heterozygote kernel carries structure distinct from A"),
  Result = c(conv_real,
             conv_ss,
             all(c(PARS$s2_g, PARS$s2_ge, PARS$s2_rep, PARS$s2_e) > 0),
             length(unique(chr_causal)) == length(CAUSAL),
             !any(chr_nulls %in% chr_causal),
             all(MARKERS %in% usable),
             max(dev_root_A, dev_root_E, dev_root_D) < 1e-8,
             nrow(L_mat) == n,
             max(abs(c(mean(diag(E_had)), mean(diag(D_dom))) - mdA)) < 1e-8,
             abs(cor(offd(A), offd(E_had))) < 0.95,
             abs(cor(offd(A), offd(D_dom))) < 0.50))

## ---------- 11. report and save ----------
cal_tab <- data.frame(
  Component = c("mu", "s2_g", "s2_ge", "s2_rep", "s2_e",
                paste0("pc_beta", 1:N_PC)),
  Value = round(c(PARS$mu, PARS$s2_g, PARS$s2_ge, PARS$s2_rep, PARS$s2_e,
                  PARS$pc_beta), 4))

cat("\n=== 07 simulation primitives ===\n")
cat("\ncalibration, from a single-stage fit to the observed data\n")
print(cal_tab, row.names = FALSE)
cat("   year effects span", round(diff(range(PARS$year_fx)), 3), "days\n")

cat("\nmarkers: floor MAF >", MAF_MIN, "leaves", length(usable), "of", length(poly),
    "polymorphic\n")
print(data.frame(marker = CAUSAL, chromosome = chr_causal), row.names = FALSE)
cat("null markers:", length(NULLS),
    "| overlap with causal chromosomes:", sum(chr_nulls %in% chr_causal), "\n")

cat("\ngenetic kernels\n"); print(kernel_tab, row.names = FALSE)
cat("   the heterozygote kernel is the centred indicator of a called heterozygote over",
    n_het_markers, "markers and", format(n_het_calls, big.mark = ","), "calls\n")
cat("   sommer::D.mat on the documented coding: correlation with A", round(cor_dmat_A, 4),
    ", with the kernel used here", round(cor_dmat_H, 4), "\n")
cat("   Vitezica coding on called genotypes: correlation with A", round(cor_vit_A, 4), "\n")
cat("   eigen roots, max |R R' - M|:",
    signif(max(dev_root_A, dev_root_E, dev_root_D), 3), "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

if (!all(checks$Result)) stop("simulation-primitive checks failed")

save(L_mat, PARS, CAUSAL, NULLS, MARKERS, G_SUB,
     ROOT_A, ROOT_E, ROOT_D, E_had, D_dom, mdA, kernel_tab,
     n_het_markers, n_het_calls, cor_dmat_A, cor_dmat_H, cor_vit_A,
     YR_REP, ID_YR, ID_IDX, YR_IDX, cal_tab, checks, usable,
     chr_causal, chr_nulls, dev_root_A, dev_root_E, dev_root_D,
     simulate_y, stage1_fast, fit_arm, ve_kernels, root_of,
     file = file.path(OUT, "07_sim_primitives.RData"))

cat("\nsaved 07_sim_primitives.RData |",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
