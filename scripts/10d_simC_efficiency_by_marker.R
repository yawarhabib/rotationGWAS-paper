################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10d: is the marker tested in scenario C a typical marker?
##
##   Scenario C simulates one causal marker at every level and in every
##   replicate. The efficiency of a scheme at a marker depends only on the
##   design and on the marker genotypes, so whether that marker is typical can
##   be answered exactly, without simulation.
##
##   For each of the five synthetic panels Omega is obtained from the design at
##   the calibrated variance components, as in script S2. For each selected
##   marker two ratios are computed, each a variance under a scheme divided by
##   the variance under the rotation:
##     omega_only  -- the sandwich variance under the true Omega with M = T'T;
##     with_kernel -- the same sandwich with the genetic variance included, the
##                    true covariance of the adjusted means being s2_g A + Omega
##                    and a scheme using s2_g A + D or s2_g A + mean(D) I; this
##                    is the ratio the simulated standard deviations estimate.
##
##   Output: 10d_simC_efficiency_by_marker.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e10 <- new.env(); load(file.path(OUT, "10_simC.RData"), envir = e10)

YEARS_GRID <- c(4, 8, 16, 32, 64); NLEV <- 5
YRS_PER_GENO <- 3; N_REP_PLOT <- 3; CAUSAL1 <- CAUSAL[1]
YEAR_SD <- sd(PARS$year_fx)

make_panel <- function(n_years, lev) {                 # as in script 10
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

sandwich <- function(z, Sigma, Mw) {                   # variance of the marker coefficient
  Xm    <- cbind(1, Q, z)
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
  Dg <- diag(diag(O)); Sig <- PARS$s2_g * A + O
  M  <- list(rot = solve(O), inv = diag(1 / diag(O)), unw = diag(n))
  Mk <- list(rot = solve(Sig), inv = solve(PARS$s2_g * A + Dg), unw = solve(PARS$s2_g * A + mean(diag(O)) * diag(n)))
  do.call(rbind, lapply(MARKERS, function(s) {
    z <- G[, s]
    a <- vapply(M,  function(m) sandwich(z, O,   m), numeric(1))
    b <- vapply(Mk, function(m) sandwich(z, Sig, m), numeric(1))
    data.frame(n_years = YEARS_GRID[k], marker = s, maf = min(mean(z) / 2, 1 - mean(z) / 2),
               omega_only_inv = a["inv"] / a["rot"], omega_only_unw = a["unw"] / a["rot"],
               with_kernel_inv = b["inv"] / b["rot"], with_kernel_unw = b["unw"] / b["rot"],
               v_rot = b["rot"], v_inv = b["inv"], v_unw = b["unw"],
               estimation_share = mean(diag(O)) / (PARS$s2_g * mean(diag(A)) + mean(diag(O))), row.names = NULL)
  }))
}))

options(width = 200)
q <- function(x) quantile(x, c(0, 0.25, 0.5, 0.75, 1))
tab <- do.call(rbind, lapply(split(res, res$n_years), function(d) {
  i <- d$marker == CAUSAL1
  data.frame(accessions_per_year = round(n * YRS_PER_GENO / d$n_years[1], 1),
             estimation_share = d$estimation_share[1],
             kernel_min = min(d$with_kernel_inv), kernel_median = median(d$with_kernel_inv),
             kernel_max = max(d$with_kernel_inv), kernel_causal = d$with_kernel_inv[i],
             causal_percentile = 100 * mean(d$with_kernel_inv <= d$with_kernel_inv[i]),
             omega_median = median(d$omega_only_inv), omega_causal = d$omega_only_inv[i],
             unw_kernel_median = median(d$with_kernel_unw), unw_kernel_causal = d$with_kernel_unw[i])
}))
cat("\n=== 10d: efficiency by marker on the scenario C designs, 106 markers ===\n")
cat("the marker simulated in scenario C:", CAUSAL1, "| MAF", round(res$maf[res$marker == CAUSAL1][1], 3),
    "| MAF of the 106:", paste(round(q(unique(res[, c("marker", "maf")])$maf), 3), collapse = " / "), "(min/q1/median/q3/max)\n")
cat("\nvariance under inverse-variance weighting relative to the rotation (kernel = with the genetic variance; omega = as defined in the Methods)\n")
print(format(tab, digits = 4), row.names = FALSE)

## against the simulation: script 10 reports the ratio of the standard
## deviations of the null estimates, averaged over the markers
sim <- e10$c_ratio
chk <- data.frame(n_years = YEARS_GRID,
                  exact_sd_ratio_inv = tapply(sqrt(res$with_kernel_inv), res$n_years, mean),
                  simulated_sd_ratio_inv = sim$diag_over_rot[match(YEARS_GRID, sim$n_years)],
                  exact_sd_ratio_unw = tapply(sqrt(res$with_kernel_unw), res$n_years, mean),
                  simulated_sd_ratio_unw = sim$unw_over_rot[match(YEARS_GRID, sim$n_years)])
cat("\ncross-check: exact ratio of standard deviations (mean over markers) against the simulated one\n")
print(format(chk, digits = 4), row.names = FALSE)
## Script 10 takes, within a replicate, the standard deviation of the null
## estimates over the markers and averages it over replicates; the exact
## quantity corresponding to it is the root of a ratio of sums over markers,
## not a mean of ratios.
pooled <- do.call(rbind, lapply(split(res, res$n_years), function(d) data.frame(
  n_years = d$n_years[1],
  mean_of_ratios_inv = mean(sqrt(d$v_inv / d$v_rot)), ratio_of_sums_inv = sqrt(sum(d$v_inv) / sum(d$v_rot)),
  mean_of_ratios_unw = mean(sqrt(d$v_unw / d$v_rot)), ratio_of_sums_unw = sqrt(sum(d$v_unw) / sum(d$v_rot)))))
pooled$simulated_inv <- sim$diag_over_rot[match(pooled$n_years, sim$n_years)]
pooled$simulated_unw <- sim$unw_over_rot[match(pooled$n_years, sim$n_years)]
cat("\nthe same cross-check with the exact quantity defined as script 10 defines the simulated one\n")
print(format(pooled, digits = 4), row.names = FALSE)
cat(sprintf("\nmarker-by-level cases in which the with-kernel ratio exceeds the Omega-only ratio: %d of %d (inverse-variance), %d of %d (unweighted)\n",
            sum(res$with_kernel_inv > res$omega_only_inv), nrow(res), sum(res$with_kernel_unw > res$omega_only_unw), nrow(res)))
cat(sprintf("\ncorrelation of the with-kernel ratio with MAF at the extreme level: %.3f\n",
            with(res[res$n_years == 64, ], cor(with_kernel_inv, maf))))
save(res, tab, chk, pooled, file = file.path(OUT, "10d_simC_efficiency_by_marker.RData"))
