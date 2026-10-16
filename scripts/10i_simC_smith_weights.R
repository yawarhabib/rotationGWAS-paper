################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10i, scenario C: inverse-variance weights against Smith weights
##
##   The inverse-variance scheme weights each adjusted mean by the inverse of
##   its own variance, the inverse of diag(Omega). The variety-trial literature
##   (Smith, Cullis and Gilmour 2001; Moehring and Piepho 2009, method 4;
##   Damesa et al. 2017, 2019; Welham et al. 2010) uses the diagonal of the
##   inverse of Omega. This script compares the two on the designs of
##   scenario C.
##
##   Exact, without simulation: for each of the five panels and each selected
##   marker, the variance of the marker coefficient under a working covariance
##   is evaluated under the true covariance s2_g A + Omega, as in 10d, and set
##   against the rotation. Smith weights enter as the working residual matrix
##   S = [diag(Omega^-1)]^-1:
##     smith_fixed   s2_g A + S, the residual held at S
##     smith_scaled  s2_g A + c S, c = mean diag(Omega) / mean diag(S), the
##                   scale matched to that of Omega
##     smith_best    the same with c chosen to minimise the marker's variance
##
##   Output: 10i_simC_smith_weights.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))

YEARS_GRID <- c(4, 8, 16, 32, 64); NLEV <- 5
YRS_PER_GENO <- 3; N_REP_PLOT <- 3; CAUSAL1 <- CAUSAL[1]

make_panel <- function(n_years, lev) {                 # as in scripts 10 and 10d
  set.seed(SEED_C + 900000 + lev)
  gids <- rownames(A); ng <- length(gids)
  starts <- 1 + floor((seq_len(ng) - 1) * (n_years - YRS_PER_GENO + 1) / ng)
  rows <- do.call(rbind, lapply(seq_len(ng), function(i)
    expand.grid(Id = gids[i], Year = starts[i]:(starts[i] + YRS_PER_GENO - 1),
                Replication = seq_len(N_REP_PLOT), stringsAsFactors = FALSE)))
  ph <- data.frame(Year        = factor(rows$Year, levels = sort(unique(rows$Year))),
                   Replication = factor(rows$Replication),
                   Id          = factor(rows$Id, levels = gids))
  ph$tmp <- rnorm(nrow(ph))
  m <- lmer(tmp ~ Id + Year + (1|Year:Replication) + (1|Id:Year), data = ph, REML = TRUE)
  L <- emmeans(m, ~ Id, rg.limit = 2e4)@linfct
  list(pheno = ph, L = L,
       yr_rep = interaction(ph$Year, ph$Replication, drop = TRUE),
       id_yr  = interaction(ph$Id,   ph$Year,        drop = TRUE))
}

sandwich <- function(Xm, Sigma, Mw) {                  # variance of the last coefficient
  bread <- solve(crossprod(Xm, Mw %*% Xm))
  (bread %*% crossprod(Xm, Mw %*% Sigma %*% Mw %*% Xm) %*% bread)[ncol(Xm), ncol(Xm)]
}

res <- do.call(rbind, lapply(seq_len(NLEV), function(k) {
  Pc <- suppressMessages(make_panel(YEARS_GRID[k], k))
  d  <- Pc$pheno; N <- nrow(d)
  Zr <- model.matrix(~ 0 + Pc$yr_rep); Zge <- model.matrix(~ 0 + Pc$id_yr)
  X1 <- model.matrix(~ Id + Year, d)
  V1 <- PARS$s2_ge * tcrossprod(Zge) + PARS$s2_rep * tcrossprod(Zr) + PARS$s2_e * diag(N)
  O  <- Pc$L %*% solve(crossprod(X1, chol2inv(chol(V1)) %*% X1)) %*% t(Pc$L); O <- (O + t(O)) / 2
  G2 <- PARS$s2_g * A
  Sig <- G2 + O
  D  <- diag(diag(O)); S <- diag(1 / diag(solve(O)))
  cS <- mean(diag(O)) / mean(diag(S))
  W  <- list(rot = solve(Sig), inv = solve(G2 + D), smith_fixed = solve(G2 + S), smith_scaled = solve(G2 + cS * S))
  do.call(rbind, lapply(MARKERS, function(s) {
    Xm <- cbind(1, Q, G[, s])
    v  <- vapply(W, function(w) sandwich(Xm, Sig, w), numeric(1))
    best <- optimize(function(lc) sandwich(Xm, Sig, solve(G2 + exp(lc) * S)), c(-3, 3), tol = 1e-8)
    data.frame(n_years = YEARS_GRID[k], marker = s,
               inv = v["inv"] / v["rot"], smith_fixed = v["smith_fixed"] / v["rot"],
               smith_scaled = v["smith_scaled"] / v["rot"], smith_best = best$objective / v["rot"],
               c_best = exp(best$minimum), c_scaled = cS,
               diag_ratio = mean(diag(S)) / mean(diag(O)), row.names = NULL)
  }))
}))

summ <- do.call(rbind, lapply(split(res, res$n_years), function(d) {
  i <- d$marker == CAUSAL1
  data.frame(accessions_per_year = round(n * YRS_PER_GENO / d$n_years[1], 1),
             S_over_Omega_diag = d$diag_ratio[1],
             inv_median = median(d$inv), smith_fixed_median = median(d$smith_fixed),
             smith_scaled_median = median(d$smith_scaled), smith_best_median = median(d$smith_best),
             inv_causal = d$inv[i], smith_fixed_causal = d$smith_fixed[i],
             smith_scaled_causal = d$smith_scaled[i], smith_best_causal = d$smith_best[i],
             smith_scaled_better = mean(d$smith_scaled < d$inv))
}))
options(width = 220)
cat("\n=== 10i: variance of the marker coefficient relative to the rotation, scenario C designs, 106 markers ===\n")
cat("    true covariance s2_g A + Omega; working covariance s2_g A + diag(Omega) (inverse-variance) or s2_g A + c S, S = [diag(Omega^-1)]^-1 (Smith)\n\n")
print(format(summ, digits = 4), row.names = FALSE)
cat("\nS_over_Omega_diag: mean diagonal of S over that of Omega, how much smaller the Smith residual variances are\n")
cat("smith_scaled_better: share of the 106 markers at which the scale-matched Smith weights beat inverse-variance weights\n")
save(res, summ, file = file.path(OUT, "10i_simC_smith_weights.RData"))
