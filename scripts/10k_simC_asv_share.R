################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10k, scenario C: the estimation share on the average semivariance
##
##   The estimation share defined with the average semivariance (ASV; Piepho
##   2019, Biometrical Journal 61, Eq. 11), trace(M P) / (n - 1) with P the
##   centring matrix, in place of the mean diagonal. The ASV is free of any
##   covariance common to all adjusted means, which the intercept absorbs.
##
##   Script 10 keeps only the mean diagonal of each replicate's estimated Omega.
##   This script refits the first stage of every null dataset of script 10 under
##   the same seed and the same panel, checks that the mean diagonal and the
##   share on it equal the stored values, and reports the ASV share.
##
##   Output: 10k_simC_asv_share.RData (read by F1 for Table 3 and Figure 2A)
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
suppressMessages({library(lme4); library(emmeans); library(parallel)})
load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
e10 <- new.env(); load(file.path(OUT, "10_simC.RData"), envir = e10)
t0 <- Sys.time()

YEARS_GRID   <- c(4, 8, 16, 32, 64)                     # as in script 10
NLEV         <- length(YEARS_GRID)
YRS_PER_GENO <- 3
N_REP_PLOT   <- 3
YEAR_SD      <- sd(PARS$year_fx)
ncore        <- max(1, min(detectCores() - 1, NCORE_SIM))
NREP_K       <- if (exists("NREP_SMOKE")) NREP_SMOKE else NREP_C

Pc_  <- diag(n) - matrix(1 / n, n, n)
asv  <- function(M) sum(diag(Pc_ %*% M)) / (n - 1)
GEN_VAR <- PARS$s2_g * mean(diag(A))                   # script 10's reference, for the check
GEN_ASV <- PARS$s2_g * asv(A)
stopifnot(abs(GEN_VAR - e10$GEN_VAR) < 1e-10)

## ---------- the panels, built exactly as in script 10 ----------
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
  stopifnot(all(table(ph$Id) > 0), nlevels(ph$Year) == n_years)
  ph$tmp <- rnorm(nrow(ph))
  m <- lmer(tmp ~ Id + Year + (1|Year:Replication) + (1|Id:Year), data = ph, REML = TRUE)
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

## ---------- the null datasets of script 10, first stage only ----------
rep_K <- function(r) {
  out <- try({
    do.call(rbind, lapply(seq_len(NLEV), function(lev) {
      ny <- YEARS_GRID[lev]; Pc <- PANELS[[as.character(ny)]]
      set.seed(SEED_C + r * 10 + lev)                  # the global-null seed of script 10
      b0 <- setNames(rep(0, length(MARKERS)), MARKERS)
      s0 <- stage1_fast(simulate_y(b0, P = Pc), P = Pc)
      if (anyNA(s0$Omega) || anyNA(s0$mu)) return(NULL)
      O <- s0$Omega; md <- mean(diag(O)); a <- asv(O)
      data.frame(rep = r, n_years = ny, mean_diag = md, asv_omega = a,
                 share_mean_diag = md / (GEN_VAR + md), share_asv = a / (GEN_ASV + a))
    }))
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

cat("refitting the first stage of", NREP_K, "x", NLEV, "null datasets on", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, { library(lme4); library(emmeans) }))
clusterExport(cl, c("PANELS", "SEED_C", "YEARS_GRID", "NLEV", "MARKERS", "GEN_VAR", "GEN_ASV",
                    "Pc_", "asv", "n", "simulate_y", "stage1_fast"))
res <- parLapply(cl, seq_len(NREP_K), rep_K)
stopCluster(cl)
K_om <- do.call(rbind, res[!vapply(res, is.null, logical(1))])

## ---------- check: the same datasets as script 10 ----------
m <- merge(K_om, e10$C_om[, c("rep", "n_years", "mean_diag", "estimation_share")], by = c("rep", "n_years"))
check_tab <- data.frame(
  rows_here = nrow(K_om), rows_script10 = sum(e10$C_om$rep <= NREP_K), rows_matched = nrow(m),
  max_diff_mean_diag = max(abs(m$mean_diag.x - m$mean_diag.y)),
  max_diff_old_share = max(abs(m$share_mean_diag - m$estimation_share)))
cat("\n=== check against the stored values of script 10 ===\n"); print(check_tab, row.names = FALSE)
stopifnot(check_tab$rows_matched == check_tab$rows_script10, check_tab$rows_here == check_tab$rows_script10,
          check_tab$max_diff_mean_diag < 1e-8, check_tab$max_diff_old_share < 1e-10)

## ---------- the estimation share by level ----------
k_share <- do.call(rbind, lapply(split(K_om, K_om$n_years), function(d) data.frame(
  n_years = d$n_years[1], accessions_per_year = round(n * YRS_PER_GENO / d$n_years[1], 1), n = nrow(d),
  share_mean_diag = mean(d$share_mean_diag), share_asv = mean(d$share_asv),
  difference = mean(d$share_asv) - mean(d$share_mean_diag))))
options(width = 200)
cat("\n=== estimation share, mean over replicates (Table 3 and Figure 2A take share_asv) ===\n")
print(format(k_share, digits = 4), row.names = FALSE)
cat(sprintf("\ngenetic ASV, s2_g x ASV(A) = %.4f x %.4f = %.3f\n", PARS$s2_g, asv(A), GEN_ASV))

save(K_om, k_share, check_tab, GEN_ASV, file = file.path(OUT, if (NREP_K == NREP_C) "10k_simC_asv_share.RData" else "10k_simC_asv_share_smoke.RData"))
cat(sprintf("\nsaved 10k_simC_asv_share | %.2f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
