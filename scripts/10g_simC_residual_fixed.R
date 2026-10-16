################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10g, scenario C: the genetic variance with the residual at one
##
##   On every dataset of scenario C, null and causal, the rotated second stage
##   is fitted twice with the marker among the fixed effects, once with the
##   residual free and once with it fixed at one, its value under a correctly
##   specified model, and the genetic variance of each fit is set beside the
##   single-stage estimate of the same dataset stored by script 10b.
##
##   Output: 10g_simC_residual_fixed.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans); library(parallel)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e10  <- new.env(); load(file.path(OUT, "10_simC.RData"), envir = e10)
e10b <- new.env(); load(file.path(OUT, "10b_simC_single_stage.RData"), envir = e10b)
t0 <- Sys.time()

NREP_G  <- as.integer(Sys.getenv("SS_NREP", NREP_C))
ncore   <- max(1, min(detectCores() - 1, NCORE_SIM))
OUTFILE <- file.path(OUT, if (NREP_G == NREP_C) "10g_simC_residual_fixed.RData" else "10g_simC_residual_fixed_smoke.RData")
YEARS_GRID <- e10$YEARS_GRID; NLEV <- length(YEARS_GRID)
YRS_PER_GENO <- 3; N_REP_PLOT <- 3; BETA_TEST <- e10$BETA_TEST; CAUSAL1 <- CAUSAL[1]
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
  ph$tmp <- NULL
  pars <- PARS; pars$year_fx <- rnorm(n_years, 0, YEAR_SD)
  list(pheno = ph, L = L,
       yr_rep = interaction(ph$Year, ph$Replication, drop = TRUE),
       id_yr  = interaction(ph$Id,   ph$Year,        drop = TRUE),
       id_idx = as.integer(ph$Id), yr_idx = as.integer(ph$Year),
       pars = pars, rootA = ROOT_A, rootE = ROOT_E, rootD = ROOT_D,
       Gsub = G_SUB, Q = Q, A = A, n = n)
}
PANELS <- suppressMessages(lapply(seq_len(NLEV), function(k) make_panel(YEARS_GRID[k], k)))
zC <- as.numeric(G_SUB[, CAUSAL1])

## The residual likelihood of y = X b + u + e, Var = vg K + ve I, in the
## spectral form: with S = I - X (X'X)^-1 X' and S K S = U diag(d) U', minus
## twice the residual log-likelihood is sum log(vg d_i + ve) + sum eta_i^2 /
## (vg d_i + ve) over the n - p non-null directions.
fit_two <- function(s1) {
  Tm <- make_transforms(s1$Omega)$rotated
  y  <- as.numeric(Tm %*% s1$mu)
  X  <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q, as.numeric(Tm %*% zC))
  K  <- Tm %*% A %*% t(Tm); K <- (K + t(K)) / 2
  S  <- diag(n) - X %*% solve(crossprod(X), t(X))
  ee <- eigen(S %*% K %*% S, symmetric = TRUE)
  keep <- seq_len(n - ncol(X))
  d   <- pmax(ee$values[keep], 0); eta <- as.numeric(crossprod(ee$vectors[, keep], y))
  dev <- function(vg, ve) sum(log(vg * d + ve)) + sum(eta^2 / (vg * d + ve))
  free <- optim(c(log(60), 0), function(p) dev(exp(p[1]), exp(p[2])), method = "BFGS",
                control = list(reltol = 1e-12, maxit = 500))
  one  <- optimize(function(lv) dev(exp(lv), 1), interval = log(c(0.5, 800)), tol = 1e-10)
  c(vg_free = exp(free$par[1]), ve_free = exp(free$par[2]), vg_one = exp(one$minimum),
    cost_of_one = one$objective - free$value)
}

rep_G <- function(r) {
  out <- try({
    rows <- list()
    for (lev in seq_len(NLEV)) {
      Pc <- PANELS[[lev]]; ny <- YEARS_GRID[lev]
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)
      for (kind in c("null", "causal")) {
        if (kind == "null") { set.seed(SEED_C + r * 10 + lev); bb <- b0 }
        else { set.seed(SEED_C + 100000 + r * 10 + lev); bb <- b0; bb[CAUSAL1] <- BETA_TEST }
        s1 <- stage1_fast(simulate_y(bb, P = Pc), P = Pc)
        rows[[length(rows) + 1]] <- data.frame(rep = r, n_years = ny, kind = kind, t(fit_two(s1)))
      }
    }
    do.call(rbind, rows)
  }, silent = TRUE)
  if (inherits(out, "try-error")) as.character(out) else out
}

cat("\nscenario C, rotated second stage with the residual free and fixed at one -", NREP_G, "replicates,", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, suppressMessages(library(lme4))))
clusterExport(cl, c("PANELS", "YEARS_GRID", "NLEV", "SEED_C", "MARKERS", "CAUSAL1", "BETA_TEST", "A", "Q", "n", "zC",
                    "fit_two", "simulate_y", "stage1_fast", "make_transforms"))
res <- parLapply(cl, seq_len(NREP_G), rep_G)
stopCluster(cl)
bad <- vapply(res, is.character, logical(1))
if (any(bad)) cat("failed replicates:", sum(bad), "| first message:", res[bad][[1]], "\n")
GG <- do.call(rbind, res[!bad]); rownames(GG) <- NULL

ss <- e10b$SS[e10b$SS$rep <= NREP_G, c("rep", "n_years", "kind", "s2_g")]
pr <- merge(GG, ss, by = c("rep", "n_years", "kind"))
sm <- function(d) data.frame(n = nrow(d),
  ss_mean = mean(d$s2_g), free_mean = mean(d$vg_free), one_mean = mean(d$vg_one),
  free_vs_ss_pct = 100 * (mean(d$vg_free) / mean(d$s2_g) - 1), one_vs_ss_pct = 100 * (mean(d$vg_one) / mean(d$s2_g) - 1),
  cor_free_ss = cor(d$vg_free, d$s2_g), cor_one_ss = cor(d$vg_one, d$s2_g),
  rmsd_free = sqrt(mean((d$vg_free - d$s2_g)^2)), rmsd_one = sqrt(mean((d$vg_one - d$s2_g)^2)),
  median_abs_pct_free = 100 * median(abs(d$vg_free / d$s2_g - 1)), median_abs_pct_one = 100 * median(abs(d$vg_one / d$s2_g - 1)),
  ve_free_median = median(d$ve_free), cost_above_2.71 = 100 * mean(d$cost_of_one > 2.71))
g_level <- do.call(rbind, lapply(split(pr, list(pr$n_years, pr$kind)), function(d) cbind(n_years = d$n_years[1], kind = d$kind[1], sm(d))))
g_level <- g_level[order(g_level$kind, g_level$n_years), ]; rownames(g_level) <- NULL

options(width = 230)
cat("\n=== 10g scenario C: genetic variance of the rotated second stage against single-stage ===\n")
cat("generating genetic variance:", round(PARS$s2_g, 2), "\n\n")
print(format(g_level, digits = 4), row.names = FALSE)
save(GG, g_level, file = OUTFILE)
cat("\nsaved", basename(OUTFILE), "|", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
