################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10f, scenario C: the two-stage schemes with the marker in the model
##
##   The two-stage schemes of script 10 hold the variance components at their
##   values under the model without the marker; the single-stage fit of script
##   10b estimates them with the marker in. This script removes that difference
##   and nothing else: on the same datasets each two-stage scheme is refitted by
##   REML with the marker among the fixed effects and tested by a Wald t on
##   n - 5 degrees of freedom. The results are set against the single-stage
##   analysis and against the fits of script 10 with the components held.
##
##   Both the causal and the null datasets of script 10 are used, the null ones
##   at the marker scenario C simulates, so that the standard deviation of the
##   estimate can be set against generalised least squares with every component
##   true (script 10e) and the size of the test is known.
##
##   Output: 10f_simC_full_fit_arms.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans); library(rrBLUP); library(parallel)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e10  <- new.env(); load(file.path(OUT, "10_simC.RData"), envir = e10)
e10b <- new.env(); load(file.path(OUT, "10b_simC_single_stage.RData"), envir = e10b)
e10e <- new.env(); load(file.path(OUT, "10e_simC_null_by_marker.RData"), envir = e10e)
t0 <- Sys.time()

NREP_F  <- as.integer(Sys.getenv("SS_NREP", NREP_C))
ncore   <- max(1, min(detectCores() - 1, NCORE_SIM))
OUTFILE <- file.path(OUT, if (NREP_F == NREP_C) "10f_simC_full_fit_arms.RData" else "10f_simC_full_fit_arms_smoke.RData")
YEARS_GRID <- e10$YEARS_GRID; NLEV <- length(YEARS_GRID)
YRS_PER_GENO <- 3; N_REP_PLOT <- 3; BETA_TEST <- e10$BETA_TEST; CAUSAL1 <- CAUSAL[1]
YEAR_SD <- sd(PARS$year_fx); DF_T <- n - 5
bonf_scan <- e10$bonf_scan
iC <- match(CAUSAL1, MARKERS)

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

## one scheme: the marker among the fixed effects, every component estimated
## by REML
full_arm <- function(s1, Tm) {
  y <- as.numeric(Tm %*% s1$mu)
  X <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q, as.numeric(Tm %*% zC))
  K <- Tm %*% A %*% t(Tm); K <- (K + t(K)) / 2
  f <- mixed.solve(y = y, X = X, Z = diag(n), K = K, SE = TRUE)
  b <- f$beta[5]; s <- f$beta.SE[5]
  data.frame(beta = b, se = s, neglog10p = -log10(2 * pt(-abs(b / s), DF_T)),
             ci_lo = b - qt(0.975, DF_T) * s, ci_hi = b + qt(0.975, DF_T) * s, Vg = f$Vu, Ve = f$Ve)
}

rep_F <- function(r) {
  out <- try({
    rows <- list()
    for (lev in seq_len(NLEV)) {
      Pc <- PANELS[[lev]]; ny <- YEARS_GRID[lev]
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)
      for (kind in c("null", "causal")) {
        if (kind == "null") { set.seed(SEED_C + r * 10 + lev); bb <- b0 }
        else { set.seed(SEED_C + 100000 + r * 10 + lev); bb <- b0; bb[CAUSAL1] <- BETA_TEST }
        s1 <- stage1_fast(simulate_y(bb, P = Pc), P = Pc)
        Ts <- make_transforms(s1$Omega)
        for (nm in names(Ts))
          rows[[length(rows) + 1]] <- cbind(rep = r, n_years = ny, kind = kind, arm = nm, full_arm(s1, Ts[[nm]]))
      }
    }
    do.call(rbind, rows)
  }, silent = TRUE)
  if (inherits(out, "try-error")) as.character(out) else out
}

cat("\nscenario C, two-stage arms with the marker in the model -", NREP_F, "replicates,", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, suppressMessages({library(lme4); library(rrBLUP)})))
clusterExport(cl, c("PANELS", "YEARS_GRID", "NLEV", "SEED_C", "MARKERS", "CAUSAL1", "BETA_TEST", "A", "Q", "n", "zC",
                    "DF_T", "full_arm", "simulate_y", "stage1_fast", "make_transforms"))
res <- parLapply(cl, seq_len(NREP_F), rep_F)
stopCluster(cl)
bad <- vapply(res, is.character, logical(1))
if (any(bad)) cat("failed replicates:", sum(bad), "| first message:", res[bad][[1]], "\n")
FF <- do.call(rbind, res[!bad]); rownames(FF) <- NULL

## ---------- causal datasets, paired with single-stage and P3D ----------
ss  <- e10b$SS[e10b$SS$kind == "causal" & e10b$SS$rep <= NREP_F, c("rep", "n_years", "beta", "se", "neglog10p")]
p3d <- e10$C_ca[e10$C_ca$rep <= NREP_F, c("rep", "n_years", "arm", "beta", "se", "neglog10p")]
ca  <- FF[FF$kind == "causal", ]
pr  <- merge(merge(ca, ss, by = c("rep", "n_years"), suffixes = c("", "_ss")), p3d, by = c("rep", "n_years", "arm"), suffixes = c("", "_p3d"))
one <- function(d) { h <- d$neglog10p > bonf_scan; hs <- d$neglog10p_ss > bonf_scan; hp <- d$neglog10p_p3d > bonf_scan
  a <- sum(h & !hs); b <- sum(!h & hs)
  data.frame(n = nrow(d), power_full = 100 * mean(h), power_p3d = 100 * mean(hp), power_ss = 100 * mean(hs),
             full_minus_ss = 100 * (mean(h) - mean(hs)), full_only = a, ss_only = b,
             p_exact = if (a + b > 0) binom.test(a, a + b)$p.value else NA,
             beta_full_minus_ss = mean(d$beta - d$beta_ss), se_of_that = sd(d$beta - d$beta_ss) / sqrt(nrow(d)),
             beta_p3d_minus_ss = mean(d$beta_p3d - d$beta_ss),
             se_ratio_full = mean(d$se / d$se_ss), se_ratio_p3d = mean(d$se_p3d / d$se_ss)) }
f_level  <- do.call(rbind, lapply(split(pr, list(pr$n_years, pr$arm)), function(d) cbind(n_years = d$n_years[1], arm = d$arm[1], one(d))))
f_level  <- f_level[order(f_level$arm, f_level$n_years), ]; rownames(f_level) <- NULL
f_pooled <- do.call(rbind, lapply(split(pr, pr$arm), function(d) cbind(arm = d$arm[1], one(d)))); rownames(f_pooled) <- NULL

## ---------- null datasets: spread and size ----------
nu  <- FF[FF$kind == "null", ]
mk  <- e10e$MK[e10e$MK$marker == iC & e10e$MK$rep <= NREP_F, ]
sdo <- aggregate(beta ~ n_years, mk[mk$arm == "oracle", ], sd);  names(sdo)[2] <- "sd_oracle"
sdp <- aggregate(beta ~ n_years, mk[mk$arm == "rotated", ], sd); names(sdp)[2] <- "sd_rot_p3d"
sdf <- aggregate(beta ~ n_years, nu[nu$arm == "rotated", ], sd); names(sdf)[2] <- "sd_rot_full"
f_null <- Reduce(function(a, b) merge(a, b, by = "n_years"), list(sdo, sdp, sdf))
f_null$p3d_over_oracle  <- f_null$sd_rot_p3d  / f_null$sd_oracle
f_null$full_over_oracle <- f_null$sd_rot_full / f_null$sd_oracle
f_size <- aggregate(cbind(type1_05 = neglog10p > -log10(0.05), type1_01 = neglog10p > 2) ~ arm, nu, mean)

options(width = 220)
cat("\n=== 10f scenario C: two-stage arms with the marker in the model ===\n")
cat("\npaired with single-stage, pooled over the five levels (power at", round(bonf_scan, 3), ", percent)\n")
print(format(f_pooled, digits = 4), row.names = FALSE)
cat("\nby level\n"); print(format(f_level, digits = 4), row.names = FALSE)
cat("\nnull datasets, the simulated marker: SD of the rotated estimate against generalised least squares with every component true\n")
print(format(f_null, digits = 5), row.names = FALSE)
cat(sprintf("\nsize of the full-fit test at the simulated marker, %s null datasets per arm\n",
            format(sum(FF$kind == "null" & FF$arm == FF$arm[1]), big.mark = ","))); print(format(f_size, digits = 4), row.names = FALSE)
save(FF, f_pooled, f_level, f_null, f_size, file = OUTFILE)
cat("\nsaved", basename(OUTFILE), "|", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
