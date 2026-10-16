################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##            Script 06b: the equivalence at every marker of the set
##
##   Script 03 shows at the test marker that least squares on the rotated data
##   reproduces generalized least squares on the unrotated means. This script
##   repeats the comparison at every marker of the comparison set of script 06
##   and reports the largest deviation, for every coefficient.
##
##   Output: equivalence_all_markers.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
load(file.path(OUT, "01_data.RData"))     # G, Q, A, TEST_SNP
load(file.path(OUT, "02_rotation.RData")) # Omega, W, blues, ic_rot, Q_rot, y_rot
load(file.path(OUT, "06_comparison.RData"))  # sel, the comparison set

mu      <- blues$BLUE
Omega_i <- solve(Omega)                   # computed directly, as in script 03
stopifnot(length(sel) > 0)
cat(sprintf("markers in the comparison set: %d\n", length(sel)))
cat(sprintf("coefficients per fit         : %d (intercept, %d PCs, marker)\n\n",
            2 + ncol(Q), ncol(Q)))

dev <- matrix(NA_real_, length(sel), 2 + ncol(Q),
              dimnames = list(sel, c("Intercept", paste0("PC", seq_len(ncol(Q))), "Marker")))
bgls <- brot <- numeric(length(sel))

for (i in seq_along(sel)) {
  dose <- G[, sel[i]]
  X_orig <- cbind(1, Q, dose)
  b_gls  <- solve(t(X_orig) %*% Omega_i %*% X_orig, t(X_orig) %*% Omega_i %*% mu)
  X_rot  <- cbind(ic_rot, Q_rot, as.numeric(W %*% dose))
  b_rot  <- solve(crossprod(X_rot), crossprod(X_rot, y_rot))
  dev[i, ] <- abs(as.numeric(b_rot) - as.numeric(b_gls))
  bgls[i]  <- b_gls[ncol(X_orig)]
  brot[i]  <- b_rot[ncol(X_rot)]
}

cat("maximum absolute deviation, rotated least squares against generalised least squares\n")
for (j in seq_len(ncol(dev)))
  cat(sprintf("   %-10s max %.3g   median %.3g\n", colnames(dev)[j],
              max(dev[, j]), median(dev[, j])))
cat(sprintf("\noverall maximum over all %d markers and all %d coefficients : %.3g\n",
            length(sel), ncol(dev), max(dev)))
cat(sprintf("largest marker-effect deviation                              : %.3g\n",
            max(abs(brot - bgls))))
cat(sprintf("relative to the marker effects themselves, which range %.3f to %.3f\n",
            min(bgls), max(bgls)))
cat(sprintf("the largest RELATIVE deviation is                            : %.3g\n",
            max(abs(brot - bgls) / pmax(abs(bgls), 1e-8))))
cat(sprintf("\nfor reference, script 03 reports %.3g at %s alone\n", 1.2e-12, TEST_SNP))

eq_all <- list(sel = sel, dev = dev, bgls = bgls, brot = brot,
               max_dev = max(dev), max_marker_dev = max(abs(brot - bgls)))
save(eq_all, file = file.path(OUT, "equivalence_all_markers.RData"))
cat("\nsaved equivalence_all_markers.RData\n")
