################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##         Script 06c: the equivalence with the variance components fixed
##
##   The single-stage analysis of the plot data and the rotated two-stage
##   analysis are compared with every variance component held at one common
##   set of values, so that only the algebra can differ. Both are solved from
##   the mixed model equations. The single-stage marginal covariance is
##     V = Zp A Zp' s2g + Zr Zr' s2r + Zv Zv' s2gt + I s2e
##   and the two-stage one is V2 = A s2g + Omega, with Omega the first-stage
##   covariance at the same s2r, s2gt and s2e. Omega enters without the
##   Kenward-Roger adjustment; the effect of the adjusted Omega is reported at
##   the end.
##
##   Output: equivalence_fixed_components.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
load(file.path(OUT, "01_data.RData"))        # pheno, A, Q, G, TEST_SNP
load(file.path(OUT, "02_rotation.RData"))    # Omega, Omega_naive, blues, varcomp, kr_factor
load(file.path(OUT, "06_comparison.RData"))  # sel

## ---------- components, one common set for both analyses ----------
## The first-stage components are taken unrounded from the stored fit of
## script 02.
suppressMessages(library(lme4))
vc   <- as.data.frame(VarCorr(model1))
s2r  <- vc$vcov[vc$grp == "Year:Replication"]
s2gt <- vc$vcov[vc$grp == "Id:Year"]
s2e  <- sigma(model1)^2
stopifnot(length(s2r) == 1, length(s2gt) == 1, is.finite(s2e))
S2G  <- 70.879                                # the rotated second-stage estimate of the genetic variance
cat(sprintf("components held fixed in BOTH analyses\n"))
cat(sprintf("   s2_g %.6f   s2_r %.10f   s2_gtau %.10f   s2_e %.10f\n", S2G, s2r, s2gt, s2e))
cat(sprintf("   (varcomp rounds these to 0.0543, 18.5397 and 20.5585; the rounded values give a\n"))
cat(sprintf("    difference of 1e-6 instead of 1e-12, which is what that rounding is worth)\n\n"))

## ---------- plot-level design ----------
ph  <- pheno
ph$Id <- factor(ph$Id, levels = rownames(A))
stopifnot(!any(is.na(ph$Id)))
Zp  <- model.matrix(~ Id - 1, ph);                colnames(Zp) <- levels(ph$Id)
Zr  <- model.matrix(~ interaction(Year, Replication, drop = TRUE) - 1, ph)
Zv  <- model.matrix(~ interaction(Id, Year, drop = TRUE) - 1, ph)
y   <- ph$Vegetation_period
Xyr <- model.matrix(~ Year, ph)                  # intercept plus year effects
Qp  <- Zp %*% Q                                  # principal components at plot level
n_p <- nrow(ph)
cat(sprintf("plot records %d | accessions %d | year-replicate levels %d | Id:Year levels %d\n\n",
            n_p, ncol(Zp), ncol(Zr), ncol(Zv)))

## the part of V that does not depend on the marker
V_base <- Zp %*% A %*% t(Zp) * S2G +
          tcrossprod(Zr) * s2r +
          tcrossprod(Zv) * s2gt +
          diag(n_p) * s2e
Vi <- chol2inv(chol(V_base))

## ---------- two-stage side ----------
eta <- blues$BLUE
V2  <- A * S2G + Omega_naive
V2i <- chol2inv(chol(V2))
V2k <- A * S2G + Omega                           # with the Kenward-Roger adjusted Omega
V2ki <- chol2inv(chol(V2k))

gls <- function(X, W, z) {
  XtW <- crossprod(X, W)
  C   <- solve(XtW %*% X)
  list(b = as.numeric(C %*% (XtW %*% z)), se = sqrt(diag(C)))
}

markers <- unique(c(TEST_SNP, sel))
res <- data.frame(SNP = markers, b_ss = NA_real_, b_2s = NA_real_,
                  se_ss = NA_real_, se_2s = NA_real_, b_2s_kr = NA_real_)
for (i in seq_along(markers)) {
  dose <- G[, markers[i]]
  ss <- gls(cbind(Xyr, Qp, Zp %*% dose), Vi,  y)
  ts <- gls(cbind(1, Q, dose),           V2i, eta)
  tk <- gls(cbind(1, Q, dose),           V2ki, eta)
  res$b_ss[i]  <- tail(ss$b, 1);  res$se_ss[i] <- tail(ss$se, 1)
  res$b_2s[i]  <- tail(ts$b, 1);  res$se_2s[i] <- tail(ts$se, 1)
  res$b_2s_kr[i] <- tail(tk$b, 1)
}

db  <- abs(res$b_2s - res$b_ss)
dse <- abs(res$se_2s - res$se_ss)
dkr <- abs(res$b_2s_kr - res$b_ss)
cat(sprintf("markers compared: %d\n\n", nrow(res)))
cat("SINGLE-STAGE PLOT LEVEL against ROTATED TWO-STAGE, all components fixed\n")
cat(sprintf("   marker effect   max |difference| %.3g   median %.3g\n", max(db), median(db)))
cat(sprintf("   standard error  max |difference| %.3g   median %.3g\n", max(dse), median(dse)))
cat(sprintf("   marker effects span %.3f to %.3f, so the largest RELATIVE difference is %.3g\n",
            min(res$b_ss), max(res$b_ss), max(db / pmax(abs(res$b_ss), 1e-8))))
cat(sprintf("\nwith the Kenward-Roger adjusted Omega instead, max |difference| %.3g\n", max(dkr)))
cat(sprintf("   (the adjustment inflates Omega by a factor of %.6f)\n", kr_factor))

k <- which(res$SNP == TEST_SNP)
cat(sprintf("\nat %s, the marker the manuscript anchors on\n", TEST_SNP))
cat(sprintf("   single-stage  %.10f  (se %.10f)\n", res$b_ss[k], res$se_ss[k]))
cat(sprintf("   two-stage     %.10f  (se %.10f)\n", res$b_2s[k], res$se_2s[k]))
cat(sprintf("   difference    %.3g   (the manuscript reports 0.012 with components estimated)\n",
            abs(res$b_2s[k] - res$b_ss[k])))

## the identity at a different s2g
S2G2 <- 35
Vi2  <- chol2inv(chol(Zp %*% A %*% t(Zp) * S2G2 + tcrossprod(Zr) * s2r +
                      tcrossprod(Zv) * s2gt + diag(n_p) * s2e))
V2i2 <- chol2inv(chol(A * S2G2 + Omega_naive))
dose <- G[, TEST_SNP]
a <- gls(cbind(Xyr, Qp, Zp %*% dose), Vi2, y); b <- gls(cbind(1, Q, dose), V2i2, eta)
cat(sprintf("\nrepeated at s2_g = %g, a value nothing was estimated at: difference %.3g\n",
            S2G2, abs(tail(a$b, 1) - tail(b$b, 1))))

eq_fixed <- res
save(eq_fixed, S2G, s2r, s2gt, s2e, file = file.path(OUT, "equivalence_fixed_components.RData"))
cat("\nsaved equivalence_fixed_components.RData\n")
