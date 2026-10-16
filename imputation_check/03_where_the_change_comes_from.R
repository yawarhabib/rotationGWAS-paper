################################################################################
##   Imputation check, step 3: where do the changes in the scan come from?
##
##   Imputing with Beagle changes two things at once: the relationship matrix,
##   and the genotype column of every marker that had missing calls. They are
##   separated here by crossing them, for the rotated scheme:
##     kernel in use  x  markers in use        (the main analysis)
##     kernel Beagle  x  markers in use
##     kernel in use  x  markers Beagle
##     kernel Beagle  x  markers Beagle
##   and the agreement of the scan is tabulated by how many calls a marker
##   lacks.
################################################################################
source("scripts/00_config.R")
suppressMessages(library(rrBLUP))
HERE <- file.path(ROOT, "imputation_check")
load(file.path(OUT, "01_data.RData"))
load(file.path(HERE, "02_check_and_refit.RData"))          # Ab, Gb
e2 <- new.env(); load(file.path(OUT, "02_rotation.RData"), envir = e2)
Omega <- e2$Omega; blues <- e2$blues
options(width = 200); sym <- function(M) (M + t(M)) / 2
Tm <- make_transforms(Omega)$rotated
scan <- function(K0, GG) {
  y  <- as.numeric(Tm %*% blues$BLUE); X0 <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q)
  Zt <- Tm %*% GG; K <- sym(Tm %*% K0 %*% t(Tm))
  nf <- mixed.solve(y = y, X = X0, Z = diag(n), K = K, return.Hinv = TRUE); Hi <- nf$Hinv
  P  <- Hi - Hi %*% X0 %*% solve(crossprod(X0, Hi %*% X0)) %*% crossprod(X0, Hi)
  PZ <- P %*% Zt; zPz <- colSums(Zt * PZ); zPy <- as.numeric(crossprod(PZ, y)); yPy <- as.numeric(crossprod(y, P %*% y))
  v2 <- n - 5; b <- zPy / zPz; s2 <- (yPy - b * zPy) / v2; se <- sqrt(s2 / zPz)
  list(Vg = nf$Vu, Ve = nf$Ve, lp = setNames(-log10(pbeta(v2 / (v2 + (b / se)^2), v2 / 2, 1 / 2)), colnames(GG)))
}
bonf <- -log10(ALPHA / ncol(G))
S <- list("kernel in use, markers in use" = scan(A, G), "kernel Beagle, markers in use" = scan(Ab, G),
          "kernel in use, markers Beagle" = scan(A, Gb), "kernel Beagle, markers Beagle" = scan(Ab, Gb))
base <- S[[1]]$lp
cat("\nrotated scheme, four combinations\n")
print(format(do.call(rbind, lapply(names(S), function(k) data.frame(combination = k, Vg = S[[k]]$Vg, residual = S[[k]]$Ve,
      above_threshold = sum(S[[k]]$lp > bonf), shared_with_manuscript = sum(S[[k]]$lp > bonf & base > bonf),
      cor_lp_with_manuscript = cor(S[[k]]$lp, base)))), digits = 4), row.names = FALSE)

mr <- colMeans(MISS); full <- S[[4]]$lp
bins <- cut(mr, c(-0.001, 0, 0.01, 0.03, 0.06, 0.11), labels = c("none", "up to 1%", "1 to 3%", "3 to 6%", "6 to 10%"))
cat("\nagreement of the scan (manuscript against Beagle), by the share of calls a marker lacks\n")
print(format(data.frame(missing_calls = levels(bins), markers = as.integer(table(bins)),
                        cor_lp = as.numeric(tapply(seq_along(mr), bins, function(i) cor(base[i], full[i]))),
                        mean_abs_change = as.numeric(tapply(abs(full - base), bins, mean))), digits = 3), row.names = FALSE)

sig <- union(names(which(base > bonf)), names(which(full > bonf)))
tab <- data.frame(marker = sub("glyma.Wm82.gnm4.", "", sig), missing = round(mr[sig], 3), lp_manuscript = round(base[sig], 2),
                  lp_kernel_only = round(S[[2]]$lp[sig], 2), lp_markers_only = round(S[[3]]$lp[sig], 2), lp_beagle = round(full[sig], 2),
                  maf_in_use = round(pmin(colMeans(G[, sig]) / 2, 1 - colMeans(G[, sig]) / 2), 3),
                  maf_beagle = round(pmin(colMeans(Gb[, sig]) / 2, 1 - colMeans(Gb[, sig]) / 2), 3))
cat("\nevery marker significant under either, rotated scheme, threshold", round(bonf, 2), "\n")
print(tab[order(-tab$lp_beagle), ], row.names = FALSE)
