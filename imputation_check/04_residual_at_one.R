################################################################################
##   Imputation check, step 4: what does holding the residual at one cost?
##
##   The rotated residual variance of the observed trial under the relationship
##   matrix in use and under the one built from the Beagle-imputed genotypes:
##   the deviance paid for fixing the residual at one, set against 2.71, the
##   one-sided 5% bound, and 3.84.
##
##   The residual likelihood of the rotated second-stage model is coded here
##   directly and checked against rrBLUP::mixed.solve before it is used.
################################################################################
source("scripts/00_config.R")
suppressMessages(library(rrBLUP))
HERE <- file.path(ROOT, "imputation_check")
load(file.path(OUT, "01_data.RData"))
load(file.path(HERE, "02_check_and_refit.RData"))          # Ab
e2 <- new.env(); load(file.path(OUT, "02_rotation.RData"), envir = e2)
Omega <- e2$Omega; blues <- e2$blues
sym <- function(M) (M + t(M)) / 2
psd <- function(M) { e <- eigen(sym(M), symmetric = TRUE); K <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors); dimnames(K) <- dimnames(M); K }
Tm <- make_transforms(Omega)$rotated
y  <- as.numeric(Tm %*% blues$BLUE); X <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q)

m <- rowMeans(MISS); p <- colMeans(G) / 2; w <- 2 * p * (1 - p)
As <- A / tcrossprod(1 - m); diag(As) <- diag(A) / (1 - m)
Zc <- sweep(G, 2, 2 * p); Zc[MISS] <- 0; Ob <- (!MISS) * 1
Apc <- tcrossprod(Zc) / (Ob %*% (t(Ob) * w)); dimnames(Apc) <- dimnames(A)
KS <- list("A in use" = A, "As" = psd(As), "Apc pairwise" = psd(Apc), "A from Beagle" = Ab)

## minus twice the residual log-likelihood of y = X b + u + e, Var = vg K + ve I
dev <- function(vg, ve, K) {
  V  <- vg * K + ve * diag(n); R <- chol(V); Vi <- chol2inv(R)
  XtViX <- crossprod(X, Vi %*% X); P <- Vi - Vi %*% X %*% solve(XtViX, crossprod(X, Vi))
  2 * sum(log(diag(R))) + as.numeric(determinant(XtViX, logarithm = TRUE)$modulus) + sum(y * (P %*% y))
}
out <- do.call(rbind, lapply(names(KS), function(k) {
  K  <- sym(Tm %*% KS[[k]] %*% t(Tm))
  ms <- mixed.solve(y = y, X = X, Z = diag(n), K = K)
  free <- optim(log(c(ms$Vu, max(ms$Ve, 1e-3))), function(p) dev(exp(p[1]), exp(p[2]), K), method = "BFGS",
                control = list(reltol = 1e-12, maxit = 500))
  one  <- optimize(function(lv) dev(exp(lv), 1, K), interval = log(c(1, 400)), tol = 1e-10)
  ## the profile interval: the residual values whose deviance is within 3.84
  ## of the minimum
  grid <- exp(seq(log(1e-4), log(8), length.out = 400))
  prof <- vapply(grid, function(ve) optimize(function(lv) dev(exp(lv), ve, K), interval = log(c(1, 400)), tol = 1e-8)$objective, numeric(1))
  inside <- grid[prof - free$value <= 3.84]
  data.frame(kernel = k, residual_mixed.solve = ms$Ve, residual_here = exp(free$par[2]), Vg_here = exp(free$par[1]),
             Vg_at_one = exp(one$minimum), cost_of_one = one$objective - free$value,
             interval_low = min(inside), interval_high = max(inside))
}))
options(width = 200)
cat("\nrotated second stage on the observed trial: the deviance paid for fixing the residual at one\n")
print(format(out, digits = 4), row.names = FALSE)
cat("\nbounds: 2.71 (one-sided 5%), 3.84 (two-sided 5%). The manuscript reports 0.49 and [0, 1.96] for the kernel in use.\n")
save(out, file = file.path(HERE, "04_residual_at_one.RData"))
