################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10j, scenarios B and C: where the off-diagonal part of Omega
##   comes from
##
##   Computed on the five scenario C panels at the calibrated variance
##   components (as in 10d and 10i):
##     1. With few accessions per environment the estimated year effects carry
##        error shared by the accessions grown in that year, so the means of
##        accessions that share a year are strongly correlated (point 3: the
##        others are correlated too, positively or negatively by distance).
##     2. Raising the replicate-within-year variance (scenario B's lever) adds
##        to Omega only a covariance common to all adjusted means, which the
##        intercept absorbs: after centring the change is zero to rounding.
##     3. The correlation is not confined to accessions that share a year. The
##        year effects are linked through the accessions grown in neighbouring
##        years, so the correlation falls with the distance between the
##        accessions' years, from positive for nearby years to negative for
##        distant ones. The mean over the pairs that share no year is near zero
##        only because these signs cancel; its absolute value is large. dist_tab
##        gives the correlation by the distance between the accessions' first
##        years.
##
##   Output: 10j_simC_covariance_structure.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))

YEARS_GRID <- c(4, 8, 16, 32, 64); NLEV <- 5
YRS_PER_GENO <- 3; N_REP_PLOT <- 3

make_panel <- function(n_years, lev) {                 # as in scripts 10, 10d and 10i
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

Pc_ <- diag(n) - matrix(1 / n, n, n)
GAP_BREAKS <- c(-1, 2, 5, 10, 20, 40, 70)
GAP_LABELS <- c("0-2 (share a year)", "3-5", "6-10", "11-20", "21-40", "41-61")
res <- lapply(seq_len(NLEV), function(k) {
  Pn <- suppressMessages(make_panel(YEARS_GRID[k], k))
  d  <- Pn$pheno; N <- nrow(d)
  Zr <- model.matrix(~ 0 + Pn$yr_rep); Zge <- model.matrix(~ 0 + Pn$id_yr)
  X1 <- model.matrix(~ Id + Year, d)
  omega <- function(s2r) {
    V1 <- PARS$s2_ge * tcrossprod(Zge) + s2r * tcrossprod(Zr) + PARS$s2_e * diag(N)
    O  <- Pn$L %*% solve(crossprod(X1, chol2inv(chol(V1)) %*% X1)) %*% t(Pn$L); (O + t(O)) / 2
  }
  O   <- omega(PARS$s2_rep); O40 <- omega(40)
  R   <- cov2cor(O); up <- upper.tri(O)
  Y   <- (table(d$Id, d$Year) > 0) + 0
  sh  <- (Y %*% t(Y))[up]
  g_d <- PARS$s2_g * sum(diag(Pc_ %*% A)) / (n - 1); e_d <- sum(diag(Pc_ %*% O)) / (n - 1)
  ## distance between the accessions' first years; with windows of three
  ## consecutive years, a distance of two or less is the same as sharing at
  ## least one year
  st  <- tapply(as.integer(as.character(d$Year)), d$Id, min)   # Id level order, as in L and in Y
  gap <- abs(outer(st, st, "-"))[up]
  stopifnot(all((gap <= YRS_PER_GENO - 1) == (sh >= 1)))
  b   <- cut(gap, GAP_BREAKS, labels = GAP_LABELS)
  apy <- round(n * YRS_PER_GENO / YEARS_GRID[k], 1)
  dist <- data.frame(accessions_per_year = apy, distance_first_years = levels(b),
                     mean_cor = as.numeric(tapply(R[up], b, mean)),
                     mean_abs_cor = as.numeric(tapply(abs(R[up]), b, mean)),
                     pairs = as.numeric(table(b)))
  dist <- dist[dist$pairs > 0, ]
  row <- data.frame(n_years = YEARS_GRID[k], accessions_per_year = apy,
             mean_var = mean(diag(O)),
             estimation_share = mean(diag(O)) / (PARS$s2_g * mean(diag(A)) + mean(diag(O))),
             heritability_diff = g_d / (g_d + e_d),
             cov_sharing = mean(O[up][sh >= 1]), cov_not_sharing = mean(O[up][sh == 0]),
             cor_sharing = mean(R[up][sh >= 1]), cor_not_sharing = mean(R[up][sh == 0]),
             abs_cor_all = mean(abs(R[up])), abs_cor_not_sharing = mean(abs(R[up][sh == 0])),
             share_pairs_sharing = mean(sh >= 1),
             block40_change_max = max(abs(O40 - O)),
             block40_change_centred_max = max(abs(Pc_ %*% (O40 - O) %*% Pc_)))
  list(row = row, dist = dist)
})
cov_tab  <- do.call(rbind, lapply(res, `[[`, "row"))
dist_tab <- do.call(rbind, lapply(res, `[[`, "dist")); rownames(dist_tab) <- NULL

options(width = 200)
cat("\n=== 10j: covariance among adjusted means by shared years, and what the block variance adds ===\n\n")
print(format(cov_tab, digits = 4), row.names = FALSE)
cat("\ncor_sharing / cor_not_sharing: mean correlation of accession means that share at least one year / none\n")
cat("block40_change_centred_max: largest entry of the change in Omega from block variance 0.0375 to 40, after\n")
cat("   centring; zero to rounding means the change is common to all adjusted means\n")
cat("abs_cor_all: mean absolute correlation over all pairs (the measure reported in the Results for scenario C)\n")
cat("\n=== correlation by the distance between the accessions' first years ===\n\n")
print(format(dist_tab, digits = 3), row.names = FALSE)
save(cov_tab, dist_tab, file = file.path(OUT, "10j_simC_covariance_structure.RData"))
