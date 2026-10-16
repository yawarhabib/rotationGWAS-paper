################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##            Script 08, scenario A: is the estimator correct?
##
##   A1  No marker has any effect, so every rejection is a false positive and
##       the type-I error can be read directly. The rotated residual is also
##       recorded: with the fitted model matching the generating model its
##       expected value is one.
##
##   A2  One marker is given a known effect, one marker per data set. Bias,
##       precision, interval coverage and power follow. The pairing of marker
##       to effect size changes with the replicate, so each effect size is
##       averaged over all six causal markers.
##
##   Two Bonferroni thresholds are reported: the one of this design, which
##   tests the selected markers, and the one of the genome scan of script 05.
##
##   Output: 08_simA.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(dplyr)          # for the summaries
library(parallel)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
t0 <- Sys.time()

BETA_LADDER <- c(-7.5, 5.0, -3.0, 2.0, -1.0, 0.5)
NBETA <- length(BETA_LADDER)
ncore <- max(1, min(detectCores() - 1, NCORE_SIM))

## the design's own threshold, and the one the real genome scan uses
bonf_sim  <- -log10(ALPHA / length(MARKERS))
bonf_scan <- -log10(ALPHA / ncol(G))

## The objects the simulation functions need, packed for the cores.
P <- list(pheno = pheno, id_idx = ID_IDX, yr_idx = YR_IDX, yr_rep = YR_REP,
          id_yr = ID_YR, pars = PARS, rootA = ROOT_A, rootE = ROOT_E,
          rootD = ROOT_D, Gsub = G_SUB, Q = Q, A = A, n = n, L = L_mat)

## ---------- A1: global null ----------
## Each replicate runs inside try(); a failed replicate is counted and the
## others go on.
rep_A1 <- function(r) {
  out <- try({
    set.seed(SEED_A + r * 10)
    b  <- setNames(rep(0, length(MARKERS)), MARKERS)
    s1 <- stage1_fast(simulate_y(b, P = P), P = P)
    Ts <- make_transforms(s1$Omega)
    do.call(rbind, lapply(names(Ts), function(nm)
      cbind(rep = r, arm = nm, s1_grad = s1$grad, s1_singular = s1$singular,
            s1_refit = s1$refit,
            fit_arm(s1, Ts[[nm]], P = P))))
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

## ---------- A2: one causal marker per data set ----------
## The marker keeps its index k while the effect it carries changes with the
## replicate, so each effect size is evaluated at every marker.
rep_A2 <- function(r) {
  out <- try({
    do.call(rbind, lapply(seq_along(CAUSAL), function(k) {
      j  <- ((k - 1 + r) %% NBETA) + 1
      mk <- CAUSAL[k]; bt <- BETA_LADDER[j]
      set.seed(SEED_A + 100000 + r * 10 + k)

      b <- setNames(rep(0, length(MARKERS)), MARKERS); b[mk] <- bt
      s1 <- stage1_fast(simulate_y(b, P = P), P = P)
      Ts <- make_transforms(s1$Omega)

      do.call(rbind, lapply(names(Ts), function(nm) {
        f  <- fit_arm(s1, Ts[[nm]], P = P)
        fr <- f[f$marker == mk, ]
        stopifnot(nrow(fr) == 1)
        data.frame(rep = r, arm = nm, marker = mk, true_beta = bt,
                   s1_grad = s1$grad, s1_singular = s1$singular,
                   s1_refit = s1$refit,
                   beta = fr$beta, se = fr$se, ci_lo = fr$ci_lo, ci_hi = fr$ci_hi,
                   neglog10p = fr$neglog10p, Vg = fr$Vg, Ve = fr$Ve,
                   row.names = NULL)
      }))
    }))
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

run_par <- function(fun, nrep, label) {
  cat("running", label, "-", nrep, "replicates on", ncore, "cores ...\n")
  cl <- makeCluster(ncore)
  on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, { library(lme4); library(emmeans); library(rrBLUP) }))
  clusterExport(cl, c("P", "SEED_A", "MARKERS", "CAUSAL", "BETA_LADDER", "NBETA",
                      "simulate_y", "stage1_fast", "fit_arm", "make_transforms"),
                envir = environment())
  res <- parLapply(cl, seq_len(nrep), fun)
  list(data = do.call(rbind, res), failed = sum(vapply(res, is.null, logical(1))))
}

r1 <- run_par(rep_A1, NREP_A_NULL,   "A1 global null")
r2 <- run_par(rep_A2, NREP_A_CAUSAL, "A2 known effects")
A1 <- r1$data; A2 <- r2$data

## ---------- summaries ----------
a1_type1 <- A1 %>% filter(!is.na(neglog10p)) %>% group_by(arm) %>%
  summarise(n_tests       = n(),
            at_0.05       = mean(neglog10p > -log10(0.05)),
            at_0.01       = mean(neglog10p > -log10(0.01)),
            at_bonf_106   = mean(neglog10p > bonf_sim),
            at_bonf_59206 = mean(neglog10p > bonf_scan), .groups = "drop") %>%
  as.data.frame()

## Vg and Ve come from the null model and are repeated on every marker row, so
## one row per replicate and scheme is taken before summarising.
a1_ve <- A1 %>% distinct(rep, arm, .keep_all = TRUE) %>% group_by(arm) %>%
  summarise(n_reps = n(), mean_Vg = mean(Vg), mean_Ve = mean(Ve),
            median_Ve = median(Ve), q025 = quantile(Ve, .025),
            q975 = quantile(Ve, .975), .groups = "drop") %>%
  as.data.frame()

a2 <- A2 %>% filter(!is.na(beta)) %>% group_by(true_beta, arm) %>%
  summarise(n           = n(),
            bias        = mean(beta - true_beta),
            SD          = sd(beta),
            RMSE        = sqrt(mean((beta - true_beta)^2)),
            coverage    = mean(ci_lo <= true_beta & true_beta <= ci_hi),
            power_106   = mean(neglog10p > bonf_sim),
            power_59206 = mean(neglog10p > bonf_scan), .groups = "drop") %>%
  arrange(desc(abs(true_beta)), arm) %>% as.data.frame()

## ---------- calibration of the standard error ----------
## Within a cell (one effect size, one scheme) the estimates share one
## sampling distribution, so the model standard error is compared with the
## realised spread cell by cell and then pooled over the cells of a scheme.
## The three schemes analyse the same simulated data, so they are reported
## separately and not pooled with one another. The realised spread is itself
## estimated from the replicates, so individual cells vary.
ratio_cell <- A2 %>% filter(!is.na(beta)) %>% group_by(arm, true_beta) %>%
  summarise(n        = n(),
            mean_se  = mean(se),
            sd_beta  = sd(beta),
            ratio    = mean(se) / sd(beta),
            coverage = mean(ci_lo <= true_beta & true_beta <= ci_hi),
            .groups  = "drop") %>%
  arrange(arm, desc(abs(true_beta))) %>% as.data.frame()

ratio_arm <- ratio_cell %>% group_by(arm) %>%
  summarise(n_cells       = n(),
            mean_ratio    = mean(ratio),
            se_ratio      = sd(ratio) / sqrt(n()),
            t_vs_one      = (mean(ratio) - 1) / (sd(ratio) / sqrt(n())),
            mean_coverage = mean(coverage),
            .groups       = "drop") %>% as.data.frame()

## how often each effect was evaluated at each marker
cell_tab <- with(A2[A2$arm == "rotated", ], table(true_beta, marker))

## One first-stage fit per replicate in A1 and one per replicate and marker in
## A2; the distinct fits are counted, not the rows that repeat them.
s1_A1 <- unique(A1[, c("rep", "s1_grad", "s1_singular", "s1_refit")])
s1_A2 <- unique(A2[, c("rep", "marker", "s1_grad", "s1_singular", "s1_refit")])

n_stage1         <- nrow(s1_A1) + nrow(s1_A2)
n_stage1_stalled <- sum(s1_A1$s1_grad > 1e-3, na.rm = TRUE) +
                    sum(s1_A2$s1_grad > 1e-3, na.rm = TRUE)
pct_singular     <- 100 * (sum(s1_A1$s1_singular) + sum(s1_A2$s1_singular)) / n_stage1
max_s1_grad      <- max(c(s1_A1$s1_grad, s1_A2$s1_grad), na.rm = TRUE)
n_stage1_refit   <- sum(s1_A1$s1_refit) + sum(s1_A2$s1_refit)
maf_causal <- pmin(colMeans(G[, CAUSAL]) / 2, 1 - colMeans(G[, CAUSAL]) / 2)

## ---------- checks ----------
## The bounds on the type-I error and the coverage are wide: they check that
## the simulation is sound; the tables report the result.
checks <- data.frame(
  Check = c("no replicate failed in A1",
            "no replicate failed in A2",
            "no first stage stalled, scaled gradient below 1e-3",
            "type-I error at 0.05 between 0.02 and 0.10, every arm",
            "interval coverage between 0.90 and 0.99, every row",
            "every effect evaluated at every marker",
            "no missing values in the summaries",
            "standard error within 10 percent of the realised spread"),
  Result = c(r1$failed == 0,
             r2$failed == 0,
             n_stage1_stalled == 0,
             all(a1_type1$at_0.05 > 0.02 & a1_type1$at_0.05 < 0.10),
             all(a2$coverage > 0.90 & a2$coverage < 0.99),
             all(cell_tab > 0),
             !anyNA(a1_type1) && !anyNA(a2),
             all(abs(ratio_arm$mean_ratio - 1) < 0.10)))

## ---------- report ----------
cat("\n=== 08 scenario A: is the estimator correct? ===\n")
cat("\nthresholds: this design tests", length(MARKERS), "markers, Bonferroni",
    round(bonf_sim, 3), "\n")
cat("            the genome scan tests", ncol(G), "markers, Bonferroni",
    round(bonf_scan, 3), "\n")

cat("\ntype-I error under a global null\n")
print(a1_type1, digits = 4, row.names = FALSE)

cat("\nvariance components under a correct model (true s2_g =",
    round(PARS$s2_g, 2), ", target Ve = 1)\n")
print(a1_ve, digits = 4, row.names = FALSE)

cat("\nbias, precision, coverage and power at known effects\n")
print(a2, digits = 4, row.names = FALSE)
cat("   each effect is averaged over all six causal markers, whose minor-allele\n")
cat("   frequencies run from", round(min(maf_causal), 3), "to",
    round(max(maf_causal), 3), "\n")

cat("\ncalibration of the standard error, by cell\n")
print(ratio_cell, digits = 4, row.names = FALSE)
cat("\npooled within each arm, five degrees of freedom\n")
print(ratio_arm, digits = 4, row.names = FALSE)
cat("   ratio is the mean model standard error over the realised spread of the\n")
cat("   estimates. Below one means the interval is too short. The arms share\n")
cat("   their data, so they are not pooled with one another.\n")

cat("\nreplicates per effect and marker, rotated arm\n")
print(cell_tab)

cat("\nfirst-stage fits:", n_stage1,
    "| refitted with bobyqa:", n_stage1_refit,
    "| still stalled:", n_stage1_stalled,
    "| largest scaled gradient:", signif(max_s1_grad, 3), "\n")
cat("singular fits:", round(pct_singular, 1),
    "percent, expected since the replicate-within-year variance is",
    round(PARS$s2_rep, 4), "\n")
cat("replicates failed outright: A1", r1$failed, ", A2", r2$failed, "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

## saved before the checks are enforced, so that a failing check leaves the
## run on disk
save(A1, A2, a1_type1, a1_ve, a2, ratio_cell, ratio_arm, cell_tab, checks,
     BETA_LADDER,
     bonf_sim, bonf_scan, maf_causal,
     n_stage1, n_stage1_stalled, n_stage1_refit, pct_singular, max_s1_grad,
     file = file.path(OUT, "08_simA.RData"))

cat("
saved 08_simA.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min
")

if (!all(checks$Result)) stop("scenario A checks failed")
