################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10e, scenario C: the null datasets marker by marker, two references
##
##   Script 10 keeps only per-replicate summaries of its null datasets. This
##   script regenerates the same datasets from their seeds and stores the
##   estimates marker by marker, for three questions:
##     1. Does the exact efficiency of script 10d hold marker by marker? The
##        standard deviation of each marker's null estimate over the replicates
##        is compared with the exact value.
##     2. What does estimating Omega cost? Besides the three schemes, the
##        rotation is run with the true Omega, the second stage still
##        estimating its own components, and generalised least squares is run
##        with every component true.
##     3. What is the size of each two-stage scheme at the marker that
##        scenario C simulates, on the datasets the single-stage analysis of
##        script 10b is tested on?
##
##   The stored summaries of script 10 are reproduced as a check that the data
##   are the same.
##
##   Output: 10e_simC_null_by_marker.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans); library(rrBLUP); library(parallel)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e10 <- new.env(); load(file.path(OUT, "10_simC.RData"), envir = e10)
e10d <- new.env(); load(file.path(OUT, "10d_simC_efficiency_by_marker.RData"), envir = e10d)
t0 <- Sys.time()

NREP_E  <- as.integer(Sys.getenv("SS_NREP", NREP_C))
ncore   <- max(1, min(detectCores() - 1, NCORE_SIM))   # the number of cores, SIM_CORES
OUTFILE <- file.path(OUT, if (NREP_E == NREP_C) "10e_simC_null_by_marker.RData"
                           else "10e_simC_null_by_marker_smoke.RData")
YEARS_GRID <- e10$YEARS_GRID; NLEV <- length(YEARS_GRID)
YRS_PER_GENO <- 3; N_REP_PLOT <- 3; CAUSAL1 <- CAUSAL[1]
YEAR_SD <- sd(PARS$year_fx)
iC <- match(CAUSAL1, MARKERS)

## ---------- the panels of script 10, with the two references ----------
make_panel <- function(n_years, lev) {
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
  P <- list(pheno = ph, L = L,
            yr_rep = interaction(ph$Year, ph$Replication, drop = TRUE),
            id_yr  = interaction(ph$Id,   ph$Year,        drop = TRUE),
            id_idx = as.integer(ph$Id), yr_idx = as.integer(ph$Year),
            pars = pars, rootA = ROOT_A, rootE = ROOT_E, rootD = ROOT_D,
            Gsub = G_SUB, Q = Q, A = A, n = n)

  ## with the true components: the first stage as a fixed linear map, its Omega,
  ## and the plot-level generalised least squares projection for all markers
  N   <- nrow(ph)
  Zr  <- model.matrix(~ 0 + P$yr_rep); Zge <- model.matrix(~ 0 + P$id_yr)
  V1  <- PARS$s2_ge * tcrossprod(Zge) + PARS$s2_rep * tcrossprod(Zr) + PARS$s2_e * diag(N)
  X1  <- model.matrix(~ Id + Year, ph)
  V1i <- chol2inv(chol(V1))
  Cb  <- solve(crossprod(X1, V1i %*% X1))
  P$H1 <- L %*% Cb %*% crossprod(X1, V1i)                 # true-components adjusted means = H1 y
  Om   <- L %*% Cb %*% t(L); Om <- (Om + t(Om)) / 2
  P$T_true <- make_transforms(Om)$rotated
  d   <- ph; d$PC1 <- Q[P$id_idx, 1]; d$PC2 <- Q[P$id_idx, 2]; d$PC3 <- Q[P$id_idx, 3]
  X0  <- model.matrix(~ PC1 + PC2 + PC3 + Year, d)
  Vi  <- chol2inv(chol(V1 + PARS$s2_g * A[P$id_idx, P$id_idx]))
  ViX <- Vi %*% X0
  Pm  <- Vi - ViX %*% solve(crossprod(X0, ViX), t(ViX))
  Zm  <- G_SUB[P$id_idx, MARKERS, drop = FALSE]
  P$W   <- Pm %*% Zm
  P$zPz <- colSums(Zm * P$W)
  P
}
PANELS <- suppressMessages(lapply(seq_len(NLEV), function(k) make_panel(YEARS_GRID[k], k)))

rep_E <- function(r) {
  out <- try({
    mk <- list(); ds <- list()
    for (lev in seq_len(NLEV)) {
      Pc <- PANELS[[lev]]; ny <- YEARS_GRID[lev]
      set.seed(SEED_C + r * 10 + lev)                    # the null dataset of script 10
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)
      y  <- simulate_y(b0, P = Pc)
      s0 <- stage1_fast(y, P = Pc)
      if (anyNA(s0$Omega) || anyNA(s0$mu)) next
      T0 <- make_transforms(s0$Omega)
      fits <- lapply(T0, function(Tm) fit_arm(s0, Tm, P = Pc))
      fits$rot_true_omega <- fit_arm(list(mu = as.numeric(Pc$H1 %*% y)), Pc$T_true, P = Pc)
      for (nm in names(fits)) {
        f <- fits[[nm]]; stopifnot(identical(f$marker, MARKERS))
        mk[[length(mk) + 1]] <- data.frame(rep = r, n_years = ny, arm = nm, marker = seq_along(MARKERS),
                                           beta = f$beta, se = f$se, neglog10p = f$neglog10p)
        ds[[length(ds) + 1]] <- data.frame(rep = r, n_years = ny, arm = nm, Vg = f$Vg[1], Ve = f$Ve[1],
                                           sd_beta = sd(f$beta), mean_se = mean(f$se),
                                           n_05 = sum(f$neglog10p > -log10(0.05), na.rm = TRUE))
      }
      bo <- as.numeric(crossprod(Pc$W, y)) / Pc$zPz; so <- sqrt(1 / Pc$zPz)
      mk[[length(mk) + 1]] <- data.frame(rep = r, n_years = ny, arm = "oracle", marker = seq_along(MARKERS),
                                         beta = bo, se = so, neglog10p = -log10(2 * pnorm(-abs(bo / so))))
    }
    list(mk = do.call(rbind, mk), ds = do.call(rbind, ds))
  }, silent = TRUE)
  if (inherits(out, "try-error")) as.character(out) else out
}

cat("\nscenario C, null datasets by marker -", NREP_E, "replicates,", NLEV, "levels,", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, suppressMessages({library(lme4); library(rrBLUP)})))
clusterExport(cl, c("PANELS", "YEARS_GRID", "NLEV", "SEED_C", "MARKERS",
                    "simulate_y", "stage1_fast", "fit_arm", "make_transforms"))
## run in batches with the results saved after each, so that an interrupted run
## can be continued (run_in_blocks in script 00)
res <- run_in_blocks(cl, seq_len(NREP_E), rep_E, sub("\\.RData$", "_partial.RData", OUTFILE), ncore * 20, t0)
stopCluster(cl)
bad <- vapply(res, is.character, logical(1))
if (any(bad)) cat("failed replicates:", sum(bad), "| first message:", res[bad][[1]], "\n")
MK <- do.call(rbind, lapply(res[!bad], `[[`, "mk"))
DS <- do.call(rbind, lapply(res[!bad], `[[`, "ds"))

## ---------- 0. same data as script 10 ----------
st <- e10$C_nl[e10$C_nl$rep <= NREP_E, c("rep", "n_years", "arm", "sd_beta", "n_05")]
cmp <- merge(DS[DS$arm %in% c("rotated", "diagonal", "unweighted"), ], st, by = c("rep", "n_years", "arm"),
             suffixes = c("", "_stored"))
same <- c(rows = nrow(cmp), sd_beta = max(abs(cmp$sd_beta - cmp$sd_beta_stored)),
          n_05 = sum(cmp$n_05 != cmp$n_05_stored))

## ---------- 1 and 2. SD of each marker over replicates ----------
SD <- aggregate(beta ~ n_years + arm + marker, MK, sd)
wide <- reshape(SD, idvar = c("n_years", "marker"), timevar = "arm", direction = "wide")
names(wide) <- sub("^beta\\.", "sd_", names(wide))
ex <- e10d$res; ex$marker <- match(ex$marker, MARKERS)
wide <- merge(wide, ex[, c("n_years", "marker", "maf", "v_rot", "v_inv", "v_unw")], by = c("n_years", "marker"))
wide$exact_inv <- sqrt(wide$v_inv / wide$v_rot); wide$exact_unw <- sqrt(wide$v_unw / wide$v_rot)
wide$sim_inv   <- wide$sd_diagonal / wide$sd_rotated; wide$sim_unw <- wide$sd_unweighted / wide$sd_rotated

by_marker <- do.call(rbind, lapply(split(wide, wide$n_years), function(d) {
  fit <- lm(log(sim_inv) ~ log(exact_inv), d)
  data.frame(n_years = d$n_years[1],
             exact_min = min(d$exact_inv), exact_median = median(d$exact_inv), exact_max = max(d$exact_inv),
             sim_median = median(d$sim_inv), cor_exact_sim = cor(d$exact_inv, d$sim_inv),
             slope_log = unname(coef(fit)[2]), mean_sim_over_exact = mean(d$sim_inv / d$exact_inv),
             ## median absolute relative difference between the simulated
             ## and the exact factor, in percent
             medabs_rel_pct = 100 * median(abs(d$sim_inv / d$exact_inv - 1)))
}))
cost <- do.call(rbind, lapply(split(wide, wide$n_years), function(d) data.frame(
  n_years = d$n_years[1],
  oracle_check   = mean(d$sd_oracle / sqrt(d$v_rot)),          # the SD with every component true against its exact value
  rot_true_omega = mean(d$sd_rot_true_omega / d$sd_oracle),    # cost of estimating the second-stage components
  rot_est_omega  = mean(d$sd_rotated / d$sd_oracle),           # cost of estimating Omega as well
  omega_alone    = mean(d$sd_rotated / d$sd_rot_true_omega),
  inv_over_oracle = mean(d$sd_diagonal / d$sd_oracle), unw_over_oracle = mean(d$sd_unweighted / d$sd_oracle))))
tail64 <- wide[wide$n_years == max(YEARS_GRID), ]
tail64 <- tail64[order(-tail64$exact_inv), c("marker", "maf", "exact_inv", "sim_inv", "exact_unw", "sim_unw")][1:8, ]

## ---------- 3. size at the marker scenario C simulates ----------
sz <- MK[MK$marker == iC, ]
size_level <- aggregate(cbind(type1_05 = neglog10p > -log10(0.05), type1_01 = neglog10p > 2) ~ n_years + arm, sz, mean)
size_all   <- aggregate(cbind(type1_05 = neglog10p > -log10(0.05), type1_01 = neglog10p > 2) ~ arm, sz, mean)
size_all$n <- as.numeric(table(sz$arm)[size_all$arm])
size_106   <- aggregate(cbind(type1_05 = neglog10p > -log10(0.05)) ~ n_years + arm, MK, mean)

## the free scale of each scheme: the residual variance of its second-stage fit
ve <- do.call(rbind, lapply(split(DS, list(DS$n_years, DS$arm)), function(d) data.frame(
  n_years = d$n_years[1], arm = d$arm[1], mean = mean(d$Ve), q05 = quantile(d$Ve, 0.05),
  median = median(d$Ve), q95 = quantile(d$Ve, 0.95), at_zero = 100 * mean(d$Ve < 1e-3))))
ve <- ve[order(ve$arm, ve$n_years), ]; rownames(ve) <- NULL

options(width = 200)
cat("\n=== 10e scenario C: null datasets by marker ===\n")
cat(sprintf("\nsame data as script 10: %d dataset-by-arm rows, largest difference in the within-replicate SD %.1e, rejection counts differing %d\n",
            same["rows"], same["sd_beta"], same["n_05"]))
cat("\n1. inverse-variance over rotation, SD of the null estimate, marker by marker: exact against simulated\n")
print(format(by_marker, digits = 4), row.names = FALSE)
cat("\n   the eight markers with the largest exact ratio at the extreme level\n")
print(format(tail64, digits = 4), row.names = FALSE)
cat("\n2. what estimation costs: SD relative to generalised least squares with every component true (mean over markers)\n")
print(format(cost, digits = 5), row.names = FALSE)
cat("\n3. size at the marker scenario C simulates, by arm, pooled over levels\n")
print(format(size_all, digits = 4), row.names = FALSE)
cat("\n   by level\n"); print(format(size_level[order(size_level$arm, size_level$n_years), ], digits = 4), row.names = FALSE)
cat("\n   for reference, type-I error at 0.05 over all 106 markers\n")
print(format(reshape(size_106, idvar = "n_years", timevar = "arm", direction = "wide"), digits = 4), row.names = FALSE)
cat("\nresidual variance of the second-stage fit, the free scale of each arm\n")
print(format(ve, digits = 4), row.names = FALSE)

save(MK, DS, wide, by_marker, cost, tail64, size_all, size_level, size_106, ve, same, file = OUTFILE)
cat("\nsaved", basename(OUTFILE), "|", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
unlink(done_name(sub("\\.RData$", "_partial.RData", OUTFILE)))   # the raw replicate results, kept by run_in_blocks until this save
