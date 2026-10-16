################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##     Script 10, scenario C: how many accessions share each environment?
##
##   Adjusted means become correlated when they depend on the same imprecisely
##   estimated environment effects. In this trial three years carry more than
##   a hundred accessions each, so the year effects are estimated tightly
##   and contribute little shared uncertainty to the means. The quantity that
##   matters is how many accessions each environment holds.
##
##   Fewer accessions per environment means the year effects are estimated
##   from less information: Omega increases relative to the genetic variance,
##   and the share of an adjusted mean's variance that weighting can act on,
##   the estimation share of script S2, increases. The tables report the
##   correlation structure of Omega, its magnitude and the estimation share.
##
##   The scenario builds its own designs: synthetic rolling panels that hold
##   the number of accessions, the years per accession, the replication and
##   the total number of plots constant and vary only the number of years, so
##   that the number of accessions per environment falls across the grid.
##
##   Output: 10_simC.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(lme4)
library(emmeans)
library(dplyr)
library(tidyr)
library(parallel)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
t0 <- Sys.time()

YEARS_GRID   <- c(4, 8, 16, 32, 64)
NLEV         <- length(YEARS_GRID)
YRS_PER_GENO <- 3
N_REP_PLOT   <- 3
BETA_TEST    <- -7.5
CAUSAL1      <- CAUSAL[1]
ncore        <- max(1, min(detectCores() - 1, NCORE_SIM))

YEAR_SD   <- sd(PARS$year_fx)
bonf_sim  <- -log10(ALPHA / length(MARKERS))
bonf_scan <- -log10(ALPHA / ncol(G))

## the reference the estimation share is measured against
GEN_VAR <- PARS$s2_g * mean(diag(A))

## ---------- build one synthetic rolling panel ----------
## Each accession is given a window of YRS_PER_GENO consecutive years, the
## windows rolling evenly across the calendar, so that the panel holds
## accessions, years per accession, replication and total plots constant while
## the number of years varies.
## A model on random values gives the averaging contrast L in lme4's own
## parameterization, the one stage1_fast multiplies.
make_panel <- function(n_years, lev) {
  set.seed(SEED_C + 900000 + lev)          # its own seed
  gids <- rownames(A)
  ng   <- length(gids)
  starts <- 1 + floor((seq_len(ng) - 1) * (n_years - YRS_PER_GENO + 1) / ng)

  rows <- do.call(rbind, lapply(seq_len(ng), function(i)
    expand.grid(Id = gids[i], Year = starts[i]:(starts[i] + YRS_PER_GENO - 1),
                Replication = seq_len(N_REP_PLOT), stringsAsFactors = FALSE)))

  ph <- data.frame(Year        = factor(rows$Year, levels = sort(unique(rows$Year))),
                   Replication = factor(rows$Replication),
                   Id          = factor(rows$Id, levels = gids))
  stopifnot(all(table(ph$Id) > 0), nlevels(ph$Year) == n_years)

  ph$tmp <- rnorm(nrow(ph))
  m <- lmer(tmp ~ Id + Year + (1|Year:Replication) + (1|Id:Year), data = ph, REML = TRUE)
  ## the reference grid, accessions by years, exceeds the default limit of
  ## emmeans at 64 years, so the limit is increased
  L <- emmeans(m, ~ Id, rg.limit = 2e4)@linfct
  ph$tmp <- NULL
  stopifnot(nrow(L) == n)

  pars <- PARS
  pars$year_fx <- rnorm(n_years, 0, YEAR_SD)

  list(pheno = ph, L = L,
       yr_rep = interaction(ph$Year, ph$Replication, drop = TRUE),
       id_yr  = interaction(ph$Id,   ph$Year,        drop = TRUE),
       id_idx = as.integer(ph$Id), yr_idx = as.integer(ph$Year),
       pars = pars, rootA = ROOT_A, rootE = ROOT_E, rootD = ROOT_D,
       Gsub = G_SUB, Q = Q, A = A, n = n)
}

PANELS <- lapply(seq_len(NLEV), function(k) make_panel(YEARS_GRID[k], k))
names(PANELS) <- as.character(YEARS_GRID)

panel_tab <- data.frame(
  years               = YEARS_GRID,
  plots               = sapply(PANELS, function(p) nrow(p$pheno)),
  accessions          = sapply(PANELS, function(p) nlevels(p$pheno$Id)),
  years_per_accession = YRS_PER_GENO,
  accessions_per_year = round(n * YRS_PER_GENO / YEARS_GRID, 1))

## ---------- one replicate, all levels ----------
## Each replicate runs inside try(); a failed replicate is counted and the
## others go on; a level whose Omega holds a missing value is recorded as
## skipped.
## Seeds follow the rule of script 00: the base of the scenario, plus 100000
## for the second sub-scenario, plus ten times the replicate number plus the
## level.
rep_C <- function(r) {
  out <- try({
    om <- nl <- ca <- list(); skipped <- integer(0)

    for (lev in seq_len(NLEV)) {
      ny <- YEARS_GRID[lev]
      Pc <- PANELS[[as.character(ny)]]

      ## --- global null ---
      set.seed(SEED_C + r * 10 + lev)
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)
      s0 <- stage1_fast(simulate_y(b0, P = Pc), P = Pc)
      if (anyNA(s0$Omega) || anyNA(s0$mu)) { skipped <- c(skipped, lev); next }

      O <- s0$Omega; Rc <- cov2cor(O); off <- abs(Rc[upper.tri(Rc)])
      evR <- eigen(Rc, symmetric = TRUE, only.values = TRUE)$values
      om[[length(om) + 1]] <- data.frame(
        rep = r, n_years = ny,
        accessions_per_year = round(n * YRS_PER_GENO / ny, 1),
        max_corr = max(off), mean_corr = mean(off),
        frac_above_05 = mean(off > 0.05), cond_corr = max(evR) / min(evR),
        ## the magnitude of Omega and the estimation share
        mean_diag = mean(diag(O)),
        var_ratio = max(diag(O)) / min(diag(O)),
        estimation_share = mean(diag(O)) / (GEN_VAR + mean(diag(O))),
        s1_grad = s0$grad, s1_singular = s0$singular, s1_refit = s0$refit,
        row.names = NULL)

      T0 <- make_transforms(O)
      for (nm in names(T0)) {
        f <- fit_arm(s0, T0[[nm]], P = Pc)
        f <- f[!is.na(f$neglog10p), ]
        nl[[length(nl) + 1]] <- data.frame(
          rep = r, n_years = ny, arm = nm, n_null = nrow(f),
          n_05   = sum(f$neglog10p > -log10(0.05)),
          n_bonf = sum(f$neglog10p > bonf_sim),
          n_scan = sum(f$neglog10p > bonf_scan),
          sd_beta = sd(f$beta), mean_se = mean(f$se),
          Ve = f$Ve[1], row.names = NULL)
      }

      ## --- one causal marker ---
      set.seed(SEED_C + 100000 + r * 10 + lev)
      b1 <- b0; b1[CAUSAL1] <- BETA_TEST
      s1 <- stage1_fast(simulate_y(b1, P = Pc), P = Pc)
      T1 <- make_transforms(s1$Omega)
      for (nm in names(T1)) {
        fr <- fit_arm(s1, T1[[nm]], P = Pc)
        fr <- fr[fr$marker == CAUSAL1, ]
        stopifnot(nrow(fr) == 1)
        ca[[length(ca) + 1]] <- data.frame(
          rep = r, n_years = ny, arm = nm, true_beta = BETA_TEST,
          beta = fr$beta, se = fr$se, ci_lo = fr$ci_lo, ci_hi = fr$ci_hi,
          neglog10p = fr$neglog10p,
          ## the causal sub-scenario fits its own first stage, so it carries
          ## its own diagnostics
          s1_grad = s1$grad, s1_singular = s1$singular, s1_refit = s1$refit,
          row.names = NULL)
      }
    }
    list(omega = do.call(rbind, om), null = do.call(rbind, nl),
         causal = do.call(rbind, ca), skipped = skipped)
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

cat("\nsynthetic rolling panels\n"); print(panel_tab, row.names = FALSE)
cat("\nrunning scenario C -", NREP_C, "replicates,", NLEV, "levels, on", ncore, "cores\n")

cl <- makeCluster(ncore)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, { library(lme4); library(emmeans); library(rrBLUP) }))
clusterExport(cl, c("PANELS", "SEED_C", "YEARS_GRID", "NLEV", "BETA_TEST",
                    "CAUSAL1", "MARKERS", "bonf_sim", "bonf_scan", "n",
                    "YRS_PER_GENO", "GEN_VAR",
                    "simulate_y", "stage1_fast", "fit_arm", "make_transforms"))
res <- parLapply(cl, seq_len(NREP_C), rep_C)
stopCluster(cl)

n_failed  <- sum(vapply(res, is.null, logical(1)))
res       <- res[!vapply(res, is.null, logical(1))]
n_skipped <- sum(vapply(res, function(x) length(x$skipped), integer(1)))

C_om <- do.call(rbind, lapply(res, `[[`, "omega"))
C_nl <- do.call(rbind, lapply(res, `[[`, "null"))
C_ca <- do.call(rbind, lapply(res, `[[`, "causal"))

## ---------- summaries ----------
## estimation_share is the share of an adjusted mean's variance that weighting
## can act on; script S2 reports it for the observed trial.
c_omega <- C_om %>% group_by(n_years, accessions_per_year) %>%
  summarise(n = n(), mean_corr = mean(mean_corr), max_corr = mean(max_corr),
            frac_above_05 = mean(frac_above_05), cond_corr = mean(cond_corr),
            mean_diag = mean(mean_diag), var_ratio = mean(var_ratio),
            estimation_share = mean(estimation_share), .groups = "drop") %>%
  arrange(n_years) %>% as.data.frame()

c_null <- C_nl %>% group_by(n_years, arm) %>%
  summarise(n = n(),
            type1_05   = sum(n_05) / sum(n_null),
            type1_bonf = sum(n_bonf) / sum(n_null),
            type1_scan = sum(n_scan) / sum(n_null),
            SD_beta = mean(sd_beta), mean_se = mean(mean_se),
            mean_Ve = mean(Ve), .groups = "drop") %>%
  arrange(n_years, arm) %>% as.data.frame()

c_ratio <- c_null %>% dplyr::select(n_years, arm, SD_beta) %>%
  pivot_wider(names_from = arm, values_from = SD_beta) %>%
  mutate(diag_over_rot = diagonal / rotated,
         unw_over_rot  = unweighted / rotated) %>% as.data.frame()

c_ratio_se <- c_null %>% dplyr::select(n_years, arm, mean_se) %>%
  pivot_wider(names_from = arm, values_from = mean_se) %>%
  mutate(diag_over_rot = diagonal / rotated,
         unw_over_rot  = unweighted / rotated) %>% as.data.frame()

c_causal <- C_ca %>% filter(!is.na(beta)) %>% group_by(n_years, arm) %>%
  summarise(n = n(),
            bias = mean(beta - true_beta),
            RMSE = sqrt(mean((beta - true_beta)^2)),
            coverage = mean(ci_lo <= true_beta & true_beta <= ci_hi),
            power_106   = mean(neglog10p > bonf_sim),
            power_59206 = mean(neglog10p > bonf_scan), .groups = "drop") %>%
  arrange(n_years, arm) %>% as.data.frame()

## ---------- first-stage diagnostics ----------
## Each level fits two first stages, one for the global null and one for the
## causal marker; both are counted.
s1_cols <- c("rep", "n_years", "s1_grad", "s1_singular", "s1_refit")
s1_all  <- rbind(data.frame(arm = "null",   unique(C_om[, s1_cols])),
                 data.frame(arm = "causal", unique(C_ca[, s1_cols])))
n_stage1_null   <- sum(s1_all$arm == "null")
n_stage1_causal <- sum(s1_all$arm == "causal")
n_stage1         <- nrow(s1_all)
n_stage1_stalled <- sum(s1_all$s1_grad > 1e-3, na.rm = TRUE)
n_stage1_refit   <- sum(s1_all$s1_refit)
max_s1_grad      <- max(s1_all$s1_grad, na.rm = TRUE)

## ---------- checks ----------
checks <- data.frame(
  Check = c("no replicate failed",
            "no level skipped for a missing value",
            "no first stage stalled after refitting",
            "the panels hold plots and accessions constant",
            "fewer accessions per environment raises the estimation share",
            "type-I error at 0.05 between 0.02 and 0.10 everywhere",
            "every grid level represented in every table",
            "no missing values in the summaries"),
  Result = c(n_failed == 0,
             n_skipped == 0,
             n_stage1_stalled == 0,
             length(unique(panel_tab$plots)) == 1 &&
               length(unique(panel_tab$accessions)) == 1,
             c_omega$estimation_share[NLEV] > c_omega$estimation_share[1],
             all(c_null$type1_05 > 0.02 & c_null$type1_05 < 0.10),
             nrow(c_omega) == NLEV && nrow(c_null) == 3 * NLEV &&
               nrow(c_causal) == 3 * NLEV,
             !anyNA(c_omega) && !anyNA(c_null) && !anyNA(c_causal)))

## ---------- report ----------
cat("\n=== 10 scenario C: how many accessions share each environment? ===\n")
cat("\nthresholds:", round(bonf_sim, 3), "over", length(MARKERS), "markers here,",
    round(bonf_scan, 3), "over", ncol(G), "in the genome scan\n")

cat("\nwhat environment size does to Omega\n")
print(c_omega, digits = 4, row.names = FALSE)
cat("   estimation_share is the fraction of an adjusted mean's variance that\n")
cat("   weighting can act on, the rest being genetic. Script S2 reports it for\n")
cat("   the observed trial, where the three schemes are indistinguishable, and\n")
cat("   predicts they separate as it rises; the genetic reference is", round(GEN_VAR, 2), "\n")

cat("\ntype-I error and spread under a global null\n")
print(c_null, digits = 4, row.names = FALSE)

cat("\nprecision ratios from the realised spread (above one favours the rotation)\n")
print(c_ratio, digits = 4, row.names = FALSE)
cat("\nthe same ratios from the model standard error, as a cross-check\n")
print(c_ratio_se, digits = 4, row.names = FALSE)

cat("\none causal marker at beta =", BETA_TEST, "\n")
print(c_causal, digits = 4, row.names = FALSE)

cat("\nfirst-stage fits:", n_stage1,
    "(null", n_stage1_null, "+ causal", n_stage1_causal, ")",
    "| refitted:", n_stage1_refit,
    "| still stalled:", n_stage1_stalled,
    "| largest scaled gradient:", signif(max_s1_grad, 3), "\n")
cat("levels skipped for a missing value:", n_skipped,
    "| replicates failed outright:", n_failed, "of", NREP_C, "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

## saved before the checks are enforced, so a failed check does not lose the run
save(C_om, C_nl, C_ca, c_omega, c_null, c_ratio, c_ratio_se, c_causal,
     panel_tab, checks, YEARS_GRID, BETA_TEST, bonf_sim, bonf_scan, GEN_VAR,
     n_failed, n_skipped, n_stage1, n_stage1_null, n_stage1_causal,
     n_stage1_stalled, n_stage1_refit, max_s1_grad,
     file = file.path(OUT, "10_simC.RData"))

cat("\nsaved 10_simC.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

if (!all(checks$Result)) stop("scenario C checks failed")
