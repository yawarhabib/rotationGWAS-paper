################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##                 Script 02: first stage and the rotation
##
##   Fits the first-stage model to the plot data, extracts the adjusted
##   accession means and their full variance-covariance matrix Omega, and builds
##   the three weighting schemes as transformations of the second-stage inputs.
##
##   Omega is the Kenward-Roger adjusted variance-covariance matrix of the
##   estimated marginal means. The adjustment uses the package pbkrtest,
##   which emmeans calls itself; if pbkrtest is not installed, emmeans
##   returns the unadjusted matrix, so the script checks the adjustment
##   through its inflation factor, which must exceed one.
##
##   Output: 02_rotation.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(lme4)
library(emmeans)

load(file.path(OUT, "01_data.RData"))
t0 <- Sys.time()

## ---------- 1. first-stage model ----------
form1  <- as.formula(paste(TRAIT, "~ Id + Year + (1|Year:Replication) + (1|Id:Year)"))
model1 <- lmer(form1, data = pheno, REML = TRUE)

## The convergence of the first-stage model is recorded, since Omega and every
## later result depend on it.
conv_msg  <- model1@optinfo$conv$lme4$messages
converged <- length(conv_msg) == 0
if (!converged) cat("first-stage convergence messages:", paste(conv_msg, collapse = "; "), "\n")

vc1 <- as.data.frame(lme4::VarCorr(model1))
stopifnot(all(c("Id:Year", "Year:Replication") %in% vc1$grp))

varcomp <- data.frame(
  Component = c("Id:Year", "Year:Replication", "Residual"),
  Variance  = round(c(vc1$vcov[vc1$grp == "Id:Year"],
                      vc1$vcov[vc1$grp == "Year:Replication"],
                      attr(lme4::VarCorr(model1), "sc")^2), 4))

## ---------- 2. adjusted means and Omega ----------
emm      <- emmeans(model1, ~ Id, mode = "kenward-roger")
emm_summ <- summary(emm)

blues <- data.frame(Id = emm_summ$Id, BLUE = emm_summ$emmean)
Omega <- as.matrix(vcov(emm))

stopifnot(identical(as.character(blues$Id), rownames(A)))

## The inflation factor is the mean ratio of the adjusted to the unadjusted
## variances of the means; the unadjusted matrix is L Var(beta) L', with L the
## contrast that averages the fixed effects into the accession means.
L_contrast <- emm@linfct
Omega_naive <- as.matrix(L_contrast %*% as.matrix(vcov(model1)) %*% t(L_contrast))
kr_factor   <- mean(diag(Omega) / diag(Omega_naive))
## Without the adjustment the factor is exactly one, so a value of one means
## that the adjustment was not applied, and the script stops.
if (!(kr_factor > 1 + 1e-9))
  stop("Omega is not Kenward-Roger adjusted: inflation factor ", format(kr_factor, digits = 10),
       ". Is pbkrtest installed and enabled?")

## ---------- 3. the three transformations ----------
Tlist <- make_transforms(Omega)
W     <- Tlist$rotated

## ---------- 4. rotated second-stage inputs ----------
y_rot  <- as.numeric(W %*% blues$BLUE)
ic_rot <- as.numeric(W %*% rep(1, n))
Q_rot  <- W %*% Q
K_rot  <- W %*% A %*% t(W)
WWt    <- W %*% t(W)          # the inverse of Omega; checked in section 6

## ---------- 5. structure of Omega ----------
R    <- cov2cor(Omega)
offd <- abs(R[upper.tri(R)])
evO  <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
evR  <- eigen(R,     symmetric = TRUE, only.values = TRUE)$values

omega_tab <- data.frame(
  Quantity = c("BLUE minimum", "BLUE maximum",
               "diag(Omega) minimum", "diag(Omega) maximum",
               "variance ratio, max over min",
               "mean absolute correlation", "maximum absolute correlation",
               "pairs above 0.05",
               "condition number of Omega", "condition number of corr(Omega)",
               "Kenward-Roger inflation factor"),
  ## as text, so that each quantity prints with its own number of decimals
  Value = c(format(round(c(min(blues$BLUE), max(blues$BLUE),
                           min(diag(Omega)), max(diag(Omega)),
                           max(diag(Omega)) / min(diag(Omega))), 3), nsmall = 3),
            format(round(c(mean(offd), max(offd)), 4), nsmall = 4),
            format(sum(offd > 0.05)),
            format(round(c(max(evO) / min(evO), max(evR) / min(evR)), 3), nsmall = 3),
            format(round(kr_factor, 6), nsmall = 6)))

## ---------- 6. checks on the rotation ----------
dev_identity <- max(abs(W %*% Omega %*% t(W) - diag(n)))
dev_symmetry <- max(abs(W - t(W)))
dev_inverse  <- max(abs(WWt - solve(Omega)))
dev_corr     <- max(abs(Tlist$diagonal %*% Omega %*% t(Tlist$diagonal) - R))
ic_unw       <- as.numeric(Tlist$unweighted %*% rep(1, n))

rot_checks <- data.frame(
  Check = c("W Omega W' = I",
            "W symmetric",
            "W W' = Omega inverse",
            "unweighted intercept is constant",
            "rotated intercept is not constant",
            "diagonal transform gives corr(Omega)",
            "first-stage model converged"),
  MaxDeviation = c(dev_identity, dev_symmetry, dev_inverse, NA, NA, dev_corr, NA),
  Result = c(dev_identity < 1e-8,
             dev_symmetry < 1e-10,
             dev_inverse  < 1e-8,
             length(unique(round(ic_unw, 10))) == 1,
             length(unique(round(ic_rot, 8))) > 1,
             dev_corr < 1e-10,
             converged))

cat("\n=== 02 first stage and rotation ===\n")
cat("\nfirst-stage variance components\n")
print(varcomp, row.names = FALSE)
cat("\nstructure of Omega\n")
print(omega_tab, row.names = FALSE)
cat("\nrotation checks\n")
print(rot_checks, row.names = FALSE)

if (!all(rot_checks$Result)) stop("rotation checks failed")

save(model1, blues, Omega, Omega_naive, kr_factor, Tlist, W,
     y_rot, ic_rot, Q_rot, K_rot, WWt, varcomp, omega_tab, rot_checks,
     file = file.path(OUT, "02_rotation.RData"))

cat("\nsaved 02_rotation.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
