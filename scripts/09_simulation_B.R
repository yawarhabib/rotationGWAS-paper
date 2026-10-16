################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##       Script 09, scenario B: does the size of the block effects matter?
##
##   The replicate-within-year variance is close to zero in this trial.
##   Increasing it introduces correlation among the adjusted means, since
##   accessions
##   sharing a replicate share that effect, and correlation is what the exact
##   rotation can exploit and inverse-variance weighting cannot. The first table
##   reports that correlation across the grid.
##
##   Two sub-scenarios at each level: a global null, from which type-I error and
##   the spread of the estimates are read; and a single causal marker, giving
##   bias, coverage and power.
##
##   Output: 09_simB.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(dplyr)          # for the summaries; tidyr for the pivot below
library(tidyr)
library(parallel)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
t0 <- Sys.time()

## from the observed value to about twice the plot error
S2R_GRID  <- c(PARS$s2_rep, 2, 5, 10, 20, 40)
NLEV      <- length(S2R_GRID)
BETA_TEST <- -7.5
CAUSAL1   <- CAUSAL[1]
ncore     <- max(1, min(detectCores() - 1, NCORE_SIM))

## the design's own threshold, and the one the real genome scan uses
bonf_sim  <- -log10(ALPHA / length(MARKERS))
bonf_scan <- -log10(ALPHA / ncol(G))

P <- list(pheno = pheno, id_idx = ID_IDX, yr_idx = YR_IDX, yr_rep = YR_REP,
          id_yr = ID_YR, pars = PARS, rootA = ROOT_A, rootE = ROOT_E,
          rootD = ROOT_D, Gsub = G_SUB, Q = Q, A = A, n = n, L = L_mat)

## ---------- one replicate, all levels ----------
## Each replicate runs inside try(); a failed replicate is counted and the
## others go on.
## Seeds follow the rule of script 00: the base of the scenario, plus 100000
## for the second sub-scenario, plus ten times the replicate number plus the
## level.
rep_B <- function(r) {
  out <- try({
    om <- nl <- ca <- list()

    for (lev in seq_len(NLEV)) {
      s2r <- S2R_GRID[lev]

      ## --- global null ---
      set.seed(SEED_B + r * 10 + lev)
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)
      s0 <- stage1_fast(simulate_y(b0, s2_rep = s2r, P = P), P = P)

      O <- s0$Omega; Rc <- cov2cor(O); off <- abs(Rc[upper.tri(Rc)])
      evR <- eigen(Rc, symmetric = TRUE, only.values = TRUE)$values
      om[[lev]] <- data.frame(rep = r, s2_rep = s2r,
                              max_corr = max(off), mean_corr = mean(off),
                              frac_above_05 = mean(off > 0.05),
                              cond_corr = max(evR) / min(evR),
                              var_ratio = max(diag(O)) / min(diag(O)),
                              s1_grad = s0$grad, s1_singular = s0$singular,
                              s1_refit = s0$refit, row.names = NULL)

      T0 <- make_transforms(O)
      for (nm in names(T0)) {
        f <- fit_arm(s0, T0[[nm]], P = P)
        f <- f[!is.na(f$neglog10p), ]
        nl[[length(nl) + 1]] <- data.frame(
          rep = r, s2_rep = s2r, arm = nm, n_null = nrow(f),
          n_05     = sum(f$neglog10p > -log10(0.05)),
          n_bonf   = sum(f$neglog10p > bonf_sim),
          n_scan   = sum(f$neglog10p > bonf_scan),
          sd_beta  = sd(f$beta),
          mean_se  = mean(f$se),
          Ve = f$Ve[1], row.names = NULL)
      }

      ## --- one causal marker ---
      set.seed(SEED_B + 100000 + r * 10 + lev)
      b1 <- b0; b1[CAUSAL1] <- BETA_TEST
      s1 <- stage1_fast(simulate_y(b1, s2_rep = s2r, P = P), P = P)
      T1 <- make_transforms(s1$Omega)
      for (nm in names(T1)) {
        fr <- fit_arm(s1, T1[[nm]], P = P)
        fr <- fr[fr$marker == CAUSAL1, ]
        stopifnot(nrow(fr) == 1)
        ca[[length(ca) + 1]] <- data.frame(
          rep = r, s2_rep = s2r, arm = nm, true_beta = BETA_TEST,
          beta = fr$beta, se = fr$se, ci_lo = fr$ci_lo, ci_hi = fr$ci_hi,
          neglog10p = fr$neglog10p,
          s1_grad = s1$grad, s1_singular = s1$singular, s1_refit = s1$refit,
          row.names = NULL)
      }
    }
    list(omega = do.call(rbind, om), null = do.call(rbind, nl),
         causal = do.call(rbind, ca))
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

cat("running scenario B -", NREP_B, "replicates,", NLEV, "levels, on", ncore, "cores\n")
cl <- makeCluster(ncore)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, { library(lme4); library(emmeans); library(rrBLUP) }))
clusterExport(cl, c("P", "SEED_B", "S2R_GRID", "NLEV", "BETA_TEST", "CAUSAL1",
                    "MARKERS", "bonf_sim", "bonf_scan",
                    "simulate_y", "stage1_fast", "fit_arm", "make_transforms"))
res <- parLapply(cl, seq_len(NREP_B), rep_B)
stopCluster(cl)

n_failed <- sum(vapply(res, is.null, logical(1)))
res <- res[!vapply(res, is.null, logical(1))]

B_om <- do.call(rbind, lapply(res, `[[`, "omega"))
B_nl <- do.call(rbind, lapply(res, `[[`, "null"))
B_ca <- do.call(rbind, lapply(res, `[[`, "causal"))

## ---------- summaries ----------
## the correlation among the adjusted means at each level of the block variance
b_omega <- B_om %>% group_by(s2_rep) %>%
  summarise(n = n(), max_corr = mean(max_corr), mean_corr = mean(mean_corr),
            frac_above_05 = mean(frac_above_05), cond_corr = mean(cond_corr),
            var_ratio = mean(var_ratio), .groups = "drop") %>% as.data.frame()

b_null <- B_nl %>% group_by(s2_rep, arm) %>%
  summarise(type1_05    = sum(n_05) / sum(n_null),
            type1_bonf  = sum(n_bonf) / sum(n_null),
            type1_scan  = sum(n_scan) / sum(n_null),
            SD_beta     = mean(sd_beta),
            mean_se     = mean(mean_se),
            mean_Ve     = mean(Ve), .groups = "drop") %>%
  arrange(s2_rep, arm) %>% as.data.frame()

## SD_beta is the spread of the null estimates across the markers within a
## replicate, averaged over replicates; markers of different allele frequency
## have different standard errors, and since every scheme sees the same markers
## that part cancels in the ratios below. mean_se is reported beside it.
b_ratio <- b_null %>% dplyr::select(s2_rep, arm, SD_beta) %>%
  pivot_wider(names_from = arm, values_from = SD_beta) %>%
  mutate(diag_over_rot = diagonal / rotated,
         unw_over_rot  = unweighted / rotated) %>% as.data.frame()

b_ratio_se <- b_null %>% dplyr::select(s2_rep, arm, mean_se) %>%
  pivot_wider(names_from = arm, values_from = mean_se) %>%
  mutate(diag_over_rot = diagonal / rotated,
         unw_over_rot  = unweighted / rotated) %>% as.data.frame()

b_causal <- B_ca %>% filter(!is.na(beta)) %>% group_by(s2_rep, arm) %>%
  summarise(n = n(),
            bias        = mean(beta - true_beta),
            RMSE        = sqrt(mean((beta - true_beta)^2)),
            coverage    = mean(ci_lo <= true_beta & true_beta <= ci_hi),
            power_106   = mean(neglog10p > bonf_sim),
            power_59206 = mean(neglog10p > bonf_scan), .groups = "drop") %>%
  arrange(s2_rep, arm) %>% as.data.frame()

## ---------- first-stage diagnostics ----------
## One fit per replicate and level in each sub-scenario; a stalled fit is one
## whose scaled gradient stays above 1e-3 after the refit inside stage1_fast.
s1_all <- rbind(B_om[, c("rep", "s2_rep", "s1_grad", "s1_singular", "s1_refit")],
                unique(B_ca[, c("rep", "s2_rep", "s1_grad", "s1_singular", "s1_refit")]))
n_stage1         <- nrow(s1_all)
n_stage1_stalled <- sum(s1_all$s1_grad > 1e-3, na.rm = TRUE)
n_stage1_refit   <- sum(s1_all$s1_refit)
max_s1_grad      <- max(s1_all$s1_grad, na.rm = TRUE)

sing_by_level <- s1_all %>% group_by(s2_rep) %>%
  summarise(pct_singular = round(100 * mean(s1_singular), 1), .groups = "drop") %>%
  as.data.frame()

## ---------- checks ----------
checks <- data.frame(
  Check = c("no replicate failed",
            "no first stage stalled after refitting",
            "the block variance raises the correlation among the means",
            "type-I error at 0.05 between 0.02 and 0.10 everywhere",
            "coverage between 0.90 and 0.99 everywhere",
            "every grid level represented in every table",
            "no missing values in the summaries"),
  Result = c(n_failed == 0,
             n_stage1_stalled == 0,
             b_omega$mean_corr[NLEV] > 2 * b_omega$mean_corr[1],
             all(b_null$type1_05 > 0.02 & b_null$type1_05 < 0.10),
             all(b_causal$coverage > 0.90 & b_causal$coverage < 0.99),
             nrow(b_omega) == NLEV && nrow(b_null) == 3 * NLEV &&
               nrow(b_causal) == 3 * NLEV,
             !anyNA(b_omega) && !anyNA(b_null) && !anyNA(b_causal)))

## ---------- report ----------
cat("\n=== 09 scenario B: does the size of the block effects matter? ===\n")
cat("\nthresholds:", round(bonf_sim, 3), "over", length(MARKERS), "markers here,",
    round(bonf_scan, 3), "over", ncol(G), "in the genome scan\n")

cat("\nwhat the block variance does to Omega\n")
print(b_omega, digits = 4, row.names = FALSE)
cat("   correlation is what the rotation exploits and inverse-variance weighting\n")
cat("   cannot. If it does not rise across the grid, nothing below can separate\n")
cat("   those two schemes.\n")

cat("\ntype-I error and spread under a global null\n")
print(b_null, digits = 4, row.names = FALSE)

cat("\nprecision ratios from the realised spread (above one favours the rotation)\n")
print(b_ratio, digits = 4, row.names = FALSE)
cat("\nthe same ratios from the model standard error, as a cross-check\n")
print(b_ratio_se, digits = 4, row.names = FALSE)

cat("\none causal marker at beta =", BETA_TEST, "\n")
print(b_causal, digits = 4, row.names = FALSE)

cat("\nfirst-stage fits:", n_stage1, "| refitted:", n_stage1_refit,
    "| still stalled:", n_stage1_stalled,
    "| largest scaled gradient:", signif(max_s1_grad, 3), "\n")
cat("singular fits by level, expected to fall as the variance rises\n")
print(sing_by_level, row.names = FALSE)
cat("replicates failed outright:", n_failed, "of", NREP_B, "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

## saved before the checks are enforced, so that a failing check leaves the
## run on disk
save(B_om, B_nl, B_ca, b_omega, b_null, b_ratio, b_ratio_se, b_causal,
     sing_by_level, checks, S2R_GRID, BETA_TEST, bonf_sim, bonf_scan,
     n_failed, n_stage1, n_stage1_stalled, n_stage1_refit, max_s1_grad,
     file = file.path(OUT, "09_simB.RData"))

cat("
saved 09_simB.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min
")

if (!all(checks$Result)) stop("scenario B checks failed")
