################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10b, scenario C: the single-stage analysis of the simulated data
##
##   Script 10 compares the three two-stage schemes with one another. This
##   script adds the single-stage analysis of the same simulated data: the same
##   panels and the same seeds, so that every single-stage fit is paired with
##   the three two-stage results stored in 10_simC.RData. For each replicate
##   and design level two plot-level models are fitted with the marker as a
##   fixed effect and the additive relationship matrix on the accession effect:
##     causal -- the dataset simulated with the causal marker, which gives the
##               estimate, its standard error and the test used for power;
##     null   -- the dataset simulated under the global null, tested at the same
##               marker, which gives the size of that test.
##
##   The test is a Wald t on n - 5 degrees of freedom, the reference
##   distribution of the two-stage schemes, so that all four analyses are
##   judged against one distribution. The two-stage schemes hold the variance
##   components at their values under the model without the marker and
##   re-estimate the scale; the single-stage fit estimates every component
##   with the marker in the model. The null datasets are tested at one marker.
##
##   Three checks: the rotated estimate is recomputed for the first CHECK_REPS
##   replicates and must equal the stored one; with the true variance
##   components the single-stage and the rotated two-stage estimates must
##   coincide; and every fit records whether it reached a boundary.
##
##   Scripts 10f, 10g and 10h read its output.
##
##   Output: 10b_simC_single_stage.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(lme4qtl); library(emmeans); library(rrBLUP); library(parallel)})

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e10 <- new.env(); load(file.path(OUT, "10_simC.RData"), envir = e10)
t0 <- Sys.time()

## The environment variable SS_NREP can set a smaller number of replicates for
## a test run; the full run uses all.
NREP_SS    <- as.integer(Sys.getenv("SS_NREP", NREP_C))
CHECK_REPS <- min(20L, NREP_SS)
OUTFILE    <- file.path(OUT, if (NREP_SS == NREP_C) "10b_simC_single_stage.RData"
                              else "10b_simC_single_stage_smoke.RData")

YEARS_GRID   <- e10$YEARS_GRID
NLEV         <- length(YEARS_GRID)
YRS_PER_GENO <- 3
N_REP_PLOT   <- 3
BETA_TEST    <- e10$BETA_TEST
CAUSAL1      <- CAUSAL[1]
DF_T         <- n - 5
YEAR_SD      <- sd(PARS$year_fx)
GRAD_TOL     <- 1e-3
BOUNDARY     <- 1e-3            # a variance below this is counted as at the boundary
ncore        <- max(1, min(detectCores() - 1, NCORE_SIM))
bonf_sim     <- e10$bonf_sim
bonf_scan    <- e10$bonf_scan
stopifnot(DF_T == 154)

## ---------- the panels of script 10, same code and seeds ----------
make_panel <- function(n_years, lev) {
  set.seed(SEED_C + 900000 + lev)
  gids <- rownames(A)
  ng   <- length(gids)
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
  pars <- PARS
  pars$year_fx <- rnorm(n_years, 0, YEAR_SD)
  list(pheno = ph, L = L,
       yr_rep = interaction(ph$Year, ph$Replication, drop = TRUE),
       id_yr  = interaction(ph$Id,   ph$Year,        drop = TRUE),
       id_idx = as.integer(ph$Id), yr_idx = as.integer(ph$Year),
       pars = pars, rootA = ROOT_A, rootE = ROOT_E, rootD = ROOT_D,
       Gsub = G_SUB, Q = Q, A = A, n = n)
}
PANELS <- suppressMessages(lapply(seq_len(NLEV), function(k) make_panel(YEARS_GRID[k], k)))
A_ss   <- A + ONES_C * matrix(1, n, n)      # A + cJ, as in script 06

## ---------- one single-stage fit ----------
FORM <- ysim ~ 1 + SNP + PC1 + PC2 + PC3 + Year + (1|Year:Replication) + (1|Id:Year) + (1|Id)

fit_ss <- function(y, Pc) {
  d <- Pc$pheno; d$ysim <- y
  d$PC1 <- Q[Pc$id_idx, 1]; d$PC2 <- Q[Pc$id_idx, 2]; d$PC3 <- Q[Pc$id_idx, 3]
  d$SNP <- as.numeric(G_SUB[Pc$id_idx, CAUSAL1])
  sg <- function(m) {
    dv <- m@optinfo$derivs
    if (is.null(dv)) NA_real_ else max(abs(solve(dv$Hessian, dv$gradient)))
  }
  ## A singular fit (the replicate variance at zero) is recorded, not treated as
  ## a failure. A fit away from the boundary is refitted when it carries another
  ## message or its scaled gradient exceeds the tolerance.
  other <- function(m) { x <- unlist(m@optinfo$conv$lme4$messages); x[!grepl("boundary (singular)", x, fixed = TRUE)] }
  tm <- system.time({
    m <- relmatLmer(FORM, data = d, relmat = list(Id = A_ss), REML = TRUE)
    g <- sg(m)
    refit <- !lme4::isSingular(m) && (length(other(m)) > 0 || isTRUE(g > GRAD_TOL))
    if (refit) {
      m <- relmatLmer(FORM, data = d, relmat = list(Id = A_ss), REML = TRUE,
                      control = lme4::lmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 1e5)))
      g <- sg(m)
    }
  })["elapsed"]
  cf <- summary(m)$coefficients["SNP", ]
  vc <- as.data.frame(lme4::VarCorr(m)); v <- function(k) vc$vcov[vc$grp == k]
  b  <- unname(cf["Estimate"]); s <- unname(cf["Std. Error"]); tq <- qt(0.975, DF_T)
  data.frame(beta = b, se = s,
             neglog10p = -log10(2 * pt(-abs(b / s), DF_T)),
             ci_lo = b - tq * s, ci_hi = b + tq * s,
             s2_g = v("Id"), s2_ge = v("Id:Year"), s2_rep = v("Year:Replication"), s2_e = v("Residual"),
             singular = lme4::isSingular(m), grad = g, refit = refit,
             ok = length(other(m)) == 0, seconds = unname(tm))
}

## ---------- the check with the true variance components ----------
oracle <- function(y, Pc) {
  p  <- Pc$pars; N <- nrow(Pc$pheno)
  d  <- Pc$pheno
  d$PC1 <- Q[Pc$id_idx, 1]; d$PC2 <- Q[Pc$id_idx, 2]; d$PC3 <- Q[Pc$id_idx, 3]
  z  <- as.numeric(G_SUB[Pc$id_idx, CAUSAL1])
  Zr <- model.matrix(~ 0 + Pc$yr_rep); Zge <- model.matrix(~ 0 + Pc$id_yr)
  V1 <- p$s2_ge * tcrossprod(Zge) + p$s2_rep * tcrossprod(Zr) + p$s2_e * diag(N)
  ## single stage
  Vi  <- chol2inv(chol(V1 + p$s2_g * A[Pc$id_idx, Pc$id_idx]))
  X0  <- model.matrix(~ PC1 + PC2 + PC3 + Year, d)
  ViX <- Vi %*% X0
  Pm  <- Vi - ViX %*% solve(crossprod(X0, ViX), t(ViX))
  zPz <- as.numeric(crossprod(z, Pm %*% z))
  ss  <- c(as.numeric(crossprod(z, Pm %*% y)) / zPz, sqrt(1 / zPz))
  ## two stage, rotated: generalised least squares on the means with known
  ## components
  X1  <- model.matrix(~ Id + Year, d)
  V1i <- chol2inv(chol(V1))
  Cb  <- solve(crossprod(X1, V1i %*% X1))
  mu  <- as.numeric(Pc$L %*% (Cb %*% crossprod(X1, V1i %*% y)))
  Om  <- Pc$L %*% Cb %*% t(Pc$L); Om <- (Om + t(Om)) / 2
  V2i <- solve(p$s2_g * A + Om)
  X2  <- cbind(1, Q, G_SUB[, CAUSAL1])
  C2  <- solve(crossprod(X2, V2i %*% X2))
  ts  <- c(as.numeric(C2 %*% crossprod(X2, V2i %*% mu))[5], sqrt(C2[5, 5]))
  c(oracle_beta = abs(ss[1] - ts[1]), oracle_se = abs(ss[2] - ts[2]))
}

## ---------- one replicate, all levels ----------
rep_ss <- function(r) {
  out <- try({
    rows <- list()
    for (lev in seq_len(NLEV)) {
      Pc <- PANELS[[lev]]; ny <- YEARS_GRID[lev]
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)

      set.seed(SEED_C + r * 10 + lev)                    # the null dataset of script 10
      y0 <- simulate_y(b0, P = Pc)
      rows[[length(rows) + 1]] <- cbind(rep = r, n_years = ny, kind = "null", true_beta = 0,
                                        fit_ss(y0, Pc), same_data = NA_real_,
                                        oracle_beta = NA_real_, oracle_se = NA_real_)

      set.seed(SEED_C + 100000 + r * 10 + lev)           # the causal dataset of script 10
      b1 <- b0; b1[CAUSAL1] <- BETA_TEST
      y1 <- simulate_y(b1, P = Pc)
      sd_chk <- NA_real_; orc <- c(oracle_beta = NA_real_, oracle_se = NA_real_)
      if (r <= CHECK_REPS) {
        s1  <- stage1_fast(y1, P = Pc)
        fr  <- fit_arm(s1, make_transforms(s1$Omega)$rotated, P = Pc)
        st  <- e10$C_ca[e10$C_ca$rep == r & e10$C_ca$n_years == ny & e10$C_ca$arm == "rotated", ]
        sd_chk <- abs(fr$beta[fr$marker == CAUSAL1] - st$beta)
        orc <- oracle(y1, Pc)
      }
      rows[[length(rows) + 1]] <- cbind(rep = r, n_years = ny, kind = "causal", true_beta = BETA_TEST,
                                        fit_ss(y1, Pc), same_data = sd_chk,
                                        oracle_beta = orc["oracle_beta"], oracle_se = orc["oracle_se"])
    }
    do.call(rbind, rows)
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

cat("\nscenario C, single-stage arm -", NREP_SS, "replicates,", NLEV, "levels,", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, suppressMessages({library(lme4); library(lme4qtl); library(rrBLUP)})))
clusterExport(cl, c("PANELS", "YEARS_GRID", "NLEV", "SEED_C", "MARKERS", "CAUSAL1", "BETA_TEST", "A",
                    "A_ss", "Q", "G_SUB", "n", "DF_T", "GRAD_TOL", "CHECK_REPS", "FORM", "e10",
                    "fit_ss", "oracle", "simulate_y", "stage1_fast", "fit_arm", "make_transforms"))

## run in batches with the results saved after each, so that an interrupted run
## can be continued (run_in_blocks in script 00)
res <- run_in_blocks(cl, seq_len(NREP_SS), rep_ss, sub("\\.RData$", "_partial.RData", OUTFILE), ncore * 5, t0)
stopCluster(cl)

n_failed <- sum(vapply(res, is.null, logical(1)))
SS <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
rownames(SS) <- NULL
SS$accessions_per_year <- round(n * YRS_PER_GENO / SS$n_years, 1)

## ---------- summaries ----------
ca <- SS[SS$kind == "causal", ]; nu <- SS[SS$kind == "null", ]
two <- e10$C_ca[e10$C_ca$rep <= NREP_SS, ]

## power, error and coverage, the four analyses side by side
arm_row <- function(d, arm) data.frame(
  n_years = d$n_years[1], arm = arm, n = nrow(d),
  bias = mean(d$beta - BETA_TEST), RMSE = sqrt(mean((d$beta - BETA_TEST)^2)),
  mean_se = mean(d$se), coverage = mean(d$ci_lo <= BETA_TEST & BETA_TEST <= d$ci_hi),
  power_scan = 100 * mean(d$neglog10p > bonf_scan), power_sim = 100 * mean(d$neglog10p > bonf_sim))
ss_causal <- do.call(rbind, lapply(YEARS_GRID, function(ny) rbind(
  arm_row(ca[ca$n_years == ny, ], "single-stage"),
  do.call(rbind, lapply(c("rotated", "diagonal", "unweighted"), function(a)
    arm_row(two[two$n_years == ny & two$arm == a, ], a))))))

## each two-stage scheme against the single-stage fit of the same dataset
pair <- merge(two[, c("rep", "n_years", "arm", "beta", "se", "neglog10p")],
              ca[, c("rep", "n_years", "beta", "se", "neglog10p")],
              by = c("rep", "n_years"), suffixes = c("", "_ss"))
ss_pair <- do.call(rbind, lapply(split(pair, list(pair$n_years, pair$arm)), function(d) data.frame(
  n_years = d$n_years[1], arm = d$arm[1], n = nrow(d),
  abs_diff_beta = mean(abs(d$beta - d$beta_ss)),
  se_ratio      = mean(d$se / d$se_ss),
  se_diff       = mean(d$se - d$se_ss),
  t_se          = mean(d$se - d$se_ss) / (sd(d$se - d$se_ss) / sqrt(nrow(d))),
  pct_se_below  = 100 * mean(d$se < d$se_ss),
  lp_diff       = mean(d$neglog10p - d$neglog10p_ss),
  discordant_scan = sum((d$neglog10p > bonf_scan) != (d$neglog10p_ss > bonf_scan)))))
ss_pair <- ss_pair[order(ss_pair$n_years, ss_pair$arm), ]; rownames(ss_pair) <- NULL

## the size of the single-stage test; the two-stage rates come from script 10
ss_null <- do.call(rbind, lapply(YEARS_GRID, function(ny) {
  d <- nu[nu$n_years == ny, ]
  data.frame(n_years = ny, n = nrow(d), type1_05 = mean(d$neglog10p > -log10(0.05)),
             type1_01 = mean(d$neglog10p > 2), sd_beta = sd(d$beta), mean_se = mean(d$se))
}))
ss_null_all <- data.frame(n = nrow(nu), type1_05 = mean(nu$neglog10p > -log10(0.05)),
                          type1_01 = mean(nu$neglog10p > 2))

## boundaries and convergence of the single-stage fits
ss_fit <- do.call(rbind, lapply(split(SS, list(SS$n_years, SS$kind)), function(d) data.frame(
  n_years = d$n_years[1], kind = d$kind[1], n = nrow(d),
  rep_at_zero = 100 * mean(d$s2_rep < BOUNDARY), ge_at_zero = 100 * mean(d$s2_ge < BOUNDARY),
  g_at_zero = 100 * mean(d$s2_g < BOUNDARY), singular = 100 * mean(d$singular),
  refitted = 100 * mean(d$refit), not_converged = sum(!d$ok),
  mean_s2_rep = mean(d$s2_rep), mean_s2_g = mean(d$s2_g), mean_seconds = mean(d$seconds))))
ss_fit <- ss_fit[order(ss_fit$kind, ss_fit$n_years), ]; rownames(ss_fit) <- NULL

near <- c(sim  = sum(abs(ca$neglog10p - bonf_sim)  < 0.045),
          scan = sum(abs(ca$neglog10p - bonf_scan) < 0.045))

chk <- ca[!is.na(ca$same_data), ]
checks <- data.frame(
  Check = c("no replicate failed",
            "same data as script 10 (rotated estimate reproduced)",
            "equivalence with true components, marker estimate below 1e-8",
            "equivalence with true components, standard error below 1e-10",
            "no single-stage fit carries a convergence message other than a boundary",
            "type-I error at 0.05 between 0.02 and 0.10 at every level"),
  Result = c(n_failed == 0,
             nrow(chk) > 0 && max(chk$same_data) < 1e-10,
             max(chk$oracle_beta) < 1e-8,
             max(chk$oracle_se) < 1e-10,
             all(SS$ok),
             all(ss_null$type1_05 > 0.02 & ss_null$type1_05 < 0.10)))

## ---------- report ----------
options(width = 200)
cat("\n=== 10b scenario C: the single-stage arm ===\n")
cat("\nthresholds:", round(bonf_scan, 3), "genome scan,", round(bonf_sim, 3), "over the 106 simulation markers\n")
cat("\none causal marker at beta =", BETA_TEST, "- four arms on the same datasets\n")
print(format(ss_causal, digits = 4), row.names = FALSE)
cat("\neach two-stage arm against the single-stage fit of the same dataset\n")
print(format(ss_pair, digits = 4), row.names = FALSE)
cat("\nsize of the single-stage test under the global null, one marker per dataset\n")
print(format(ss_null, digits = 4), row.names = FALSE)
cat("pooled over levels:\n"); print(format(ss_null_all, digits = 4), row.names = FALSE)
cat("\ntype-I error at 0.05 of the two-stage arms, as stored by script 10 (all 106 markers)\n")
print(format(e10$c_null[, c("n_years", "arm", "type1_05")], digits = 4), row.names = FALSE)
cat("\nsingle-stage fits: boundaries (percent of fits), convergence, time\n")
print(format(ss_fit, digits = 3), row.names = FALSE)
s1 <- unique(e10$C_ca[e10$C_ca$rep <= NREP_SS, c("rep", "n_years", "s1_singular")])
cat("\nfor comparison, first-stage fits of the two-stage arms at a boundary, percent, by level:",
    paste(sprintf("%.1f", 100 * tapply(s1$s1_singular, s1$n_years, mean)), collapse = ", "), "\n")
cat(sprintf("\ncausal replicates with a single-stage -log10 p within 0.045 of a threshold: %d at %.3f, %d at %.3f, of %d\n",
            near["sim"], bonf_sim, near["scan"], bonf_scan, nrow(ca)))
cat(sprintf("checked on %d datasets: same data, largest difference %.1e; equivalence, estimate %.1e, standard error %.1e\n",
            nrow(chk), max(chk$same_data), max(chk$oracle_beta), max(chk$oracle_se)))
cat("\nchecks\n"); print(checks, row.names = FALSE)

save(SS, ss_causal, ss_pair, ss_null, ss_null_all, ss_fit, near, checks, n_failed,
     DF_T, BOUNDARY, bonf_sim, bonf_scan, YEARS_GRID, BETA_TEST, file = OUTFILE)
cat("\nsaved", basename(OUTFILE), "|", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
unlink(done_name(sub("\\.RData$", "_partial.RData", OUTFILE)))   # the raw replicate results, kept by run_in_blocks until this save
