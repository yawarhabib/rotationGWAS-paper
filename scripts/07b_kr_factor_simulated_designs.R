################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script 07b: the Kenward-Roger factor on the simulated designs
##
##   The observed data use the Kenward-Roger adjusted Omega, the simulations
##   the unadjusted one, L Var(beta) L'. On the observed trial the adjusted
##   variances are larger by four parts in ten thousand on average. This script
##   checks the ratio on the simulated designs. It simulates one dataset per
##   design and block variance at the calibrated variance components and
##   averages the ratio of the adjusted to the unadjusted variance of each
##   adjusted mean over the accessions, for
##     - the observed design (scenarios A and B), at the observed block
##       variance and at 40
##     - the synthetic rolling panels of scenario C at 119.2 and 7.5 accessions
##       per environment
##
##   Output: the printed table; nothing is saved
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))

ea <- eigen(A, symmetric = TRUE); ea$values[ea$values < 0] <- 0
Ahalf <- ea$vectors %*% diag(sqrt(ea$values))

kr_ratio <- function(d, s2r) {
  yr <- interaction(d$Year, d$Replication, drop = TRUE)
  ge <- interaction(d$Id, d$Year, drop = TRUE)
  u  <- as.numeric(Ahalf %*% rnorm(n)) * sqrt(PARS$s2_g); names(u) <- rownames(A)
  d$y <- u[as.character(d$Id)] + rnorm(nlevels(yr), 0, sqrt(s2r))[yr] +
         rnorm(nlevels(ge), 0, sqrt(PARS$s2_ge))[ge] + rnorm(nrow(d), 0, sqrt(PARS$s2_e))
  m <- suppressMessages(lmer(y ~ Id + Year + (1|Year:Replication) + (1|Id:Year),
                             data = d, REML = TRUE))
  L <- emmeans(m, ~ Id, rg.limit = 2e4)@linfct
  r <- diag(L %*% as.matrix(pbkrtest::vcovAdj(m)) %*% t(L)) /
       diag(L %*% as.matrix(vcov(m)) %*% t(L))
  c(mean = mean(r), max = max(r))
}

## observed design, scenarios A and B
d_obs <- droplevels(pheno[, c("Id", "Year", "Replication")])
d_obs$Id <- factor(d_obs$Id, levels = rownames(A))
obs_tab <- do.call(rbind, lapply(c(PARS$s2_rep, 40), function(s2r) {
  set.seed(1000 + round(s2r))
  data.frame(design = "observed", block_variance = s2r, t(kr_ratio(d_obs, s2r)))
}))

## synthetic rolling panels, scenario C
YEARS_GRID <- c(4, 8, 16, 32, 64); YRS_PER_GENO <- 3; N_REP_PLOT <- 3
make_panel <- function(n_years, lev) {                 # as in script 10
  set.seed(SEED_C + 900000 + lev)
  gids <- rownames(A); ng <- length(gids)
  starts <- 1 + floor((seq_len(ng) - 1) * (n_years - YRS_PER_GENO + 1) / ng)
  rows <- do.call(rbind, lapply(seq_len(ng), function(i)
    expand.grid(Id = gids[i], Year = starts[i]:(starts[i] + YRS_PER_GENO - 1),
                Replication = seq_len(N_REP_PLOT), stringsAsFactors = FALSE)))
  data.frame(Year = factor(rows$Year, levels = sort(unique(rows$Year))),
             Replication = factor(rows$Replication),
             Id = factor(rows$Id, levels = gids))
}
c_tab <- do.call(rbind, lapply(c(1, 5), function(k) {
  d <- make_panel(YEARS_GRID[k], k)
  set.seed(20261003 + k)
  data.frame(design = sprintf("scenario C, %.1f per environment", n * YRS_PER_GENO / YEARS_GRID[k]),
             block_variance = PARS$s2_rep, t(kr_ratio(d, PARS$s2_rep)))
}))

kr_tab <- rbind(obs_tab, c_tab)
cat("\n=== Kenward-Roger factor: adjusted over unadjusted variance of the adjusted means ===\n\n")
print(format(kr_tab, digits = 8), row.names = FALSE)
