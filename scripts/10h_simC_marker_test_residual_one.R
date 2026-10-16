################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10h, scenario C: the marker test with the rotated residual at one
##
##   On every dataset of scenario C, null and causal, the rotation is fitted
##   with the marker among the fixed effects, with the residual free and with
##   it fixed at one, and the marker effect, its standard error and the Wald
##   test on n - 5 degrees of freedom are kept for both. With the residual
##   free, its scale absorbs any error in the overall scale of the estimated
##   Omega; fixed at one, nothing does, so the size of the test is the
##   question. The free fit repeats the rotated scheme of script 10f through a
##   different route, which checks the likelihood coded here; the fixed fit is
##   set against the single-stage analysis of script 10b.
##
##   Output: 10h_simC_marker_test_residual_one.RData
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
e10f <- new.env(); load(file.path(OUT, "10f_simC_full_fit_arms.RData"), envir = e10f)
t0 <- Sys.time()

NREP_H  <- as.integer(Sys.getenv("SS_NREP", NREP_C))
ncore   <- max(1, min(detectCores() - 1, NCORE_SIM))
OUTFILE <- file.path(OUT, if (NREP_H == NREP_C) "10h_simC_marker_test_residual_one.RData" else "10h_simC_marker_test_residual_one_smoke.RData")
YEARS_GRID <- e10$YEARS_GRID; NLEV <- length(YEARS_GRID)
YRS_PER_GENO <- 3; N_REP_PLOT <- 3; BETA_TEST <- e10$BETA_TEST; CAUSAL1 <- CAUSAL[1]
YEAR_SD <- sd(PARS$year_fx); DF_T <- n - 5
bonf_scan <- e10$bonf_scan

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
  gls <- function(vg, ve) {                           # the marker effect and its standard error at given components
    Vi <- chol2inv(chol(vg * K + ve * diag(n))); XtViX <- crossprod(X, Vi %*% X); C <- solve(XtViX)
    c(beta = as.numeric(C %*% crossprod(X, Vi %*% y))[5], se = sqrt(C[5, 5]))
  }
  gf <- gls(exp(free$par[1]), exp(free$par[2])); g1 <- gls(exp(one$minimum), 1)
  c(vg_free = exp(free$par[1]), ve_free = exp(free$par[2]), vg_one = exp(one$minimum), cost_of_one = one$objective - free$value,
    beta_free = gf[["beta"]], se_free = gf[["se"]], beta_one = g1[["beta"]], se_one = g1[["se"]])
}

rep_H <- function(r) {
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

cat("\nscenario C, rotated second stage with the residual free and fixed at one -", NREP_H, "replicates,", ncore, "cores\n")
cl <- makeCluster(ncore)
invisible(clusterEvalQ(cl, suppressMessages(library(lme4))))
clusterExport(cl, c("PANELS", "YEARS_GRID", "NLEV", "SEED_C", "MARKERS", "CAUSAL1", "BETA_TEST", "A", "Q", "n", "zC",
                    "fit_two", "simulate_y", "stage1_fast", "make_transforms"))
res <- parLapply(cl, seq_len(NREP_H), rep_H)
stopCluster(cl)
bad <- vapply(res, is.character, logical(1))
if (any(bad)) cat("failed replicates:", sum(bad), "| first message:", res[bad][[1]], "\n")
HH <- do.call(rbind, res[!bad]); rownames(HH) <- NULL
pv <- function(b, s) 2 * pt(-abs(b / s), DF_T)

ss <- e10b$SS[e10b$SS$rep <= NREP_H, c("rep", "n_years", "kind", "beta", "se", "neglog10p")]; names(ss)[4:6] <- c("beta_ss", "se_ss", "nlp_ss")
ff <- e10f$FF[e10f$FF$arm == "rotated" & e10f$FF$rep <= NREP_H, c("rep", "n_years", "kind", "beta", "se")]; names(ff)[4:5] <- c("beta_10f", "se_10f")
pr <- merge(merge(HH, ss, by = c("rep", "n_years", "kind")), ff, by = c("rep", "n_years", "kind"))
pr$p_free <- pv(pr$beta_free, pr$se_free); pr$p_one <- pv(pr$beta_one, pr$se_one); pr$p_ss <- 10^(-pr$nlp_ss)

h_check <- c(max_abs_beta = max(abs(pr$beta_free - pr$beta_10f)), max_rel_se = max(abs(pr$se_free / pr$se_10f - 1)))

nu <- pr[pr$kind == "null", ]; ca <- pr[pr$kind == "causal", ]
sz <- function(d) data.frame(n = nrow(d),
  size05_one = mean(d$p_one < 0.05), size05_free = mean(d$p_free < 0.05), size05_ss = mean(d$p_ss < 0.05),
  size01_one = mean(d$p_one < 0.01), size01_free = mean(d$p_free < 0.01), size01_ss = mean(d$p_ss < 0.01),
  size05_one_z = mean(2 * pnorm(-abs(d$beta_one / d$se_one)) < 0.05),
  cal_one = sum(d$se_one^2) / sum(d$beta_one^2), cal_free = sum(d$se_free^2) / sum(d$beta_free^2), cal_ss = sum(d$se_ss^2) / sum(d$beta_ss^2),
  se_ratio_one = mean(d$se_one / d$se_ss), se_ratio_free = mean(d$se_free / d$se_ss))
h_size <- rbind(cbind(n_years = "all", sz(nu)), do.call(rbind, lapply(split(nu, nu$n_years), function(d) cbind(n_years = as.character(d$n_years[1]), sz(d)))))
pw <- function(d) { h1 <- -log10(d$p_one) > bonf_scan; hf <- -log10(d$p_free) > bonf_scan; hs <- d$nlp_ss > bonf_scan
  data.frame(n = nrow(d), power_one = 100 * mean(h1), power_free = 100 * mean(hf), power_ss = 100 * mean(hs),
    one_only_vs_ss = sum(h1 & !hs), ss_only = sum(!h1 & hs), one_only_vs_free = sum(h1 & !hf), free_only = sum(!h1 & hf),
    beta_one_minus_ss = mean(d$beta_one - d$beta_ss), se_of_that = sd(d$beta_one - d$beta_ss) / sqrt(nrow(d)),
    mean_abs_dev_one = mean(abs(d$beta_one - d$beta_ss)), mean_abs_dev_free = mean(abs(d$beta_free - d$beta_ss)),
    se_ratio_one = mean(d$se_one / d$se_ss), se_ratio_free = mean(d$se_free / d$se_ss)) }
h_power <- rbind(cbind(n_years = "all", pw(ca)), do.call(rbind, lapply(split(ca, ca$n_years), function(d) cbind(n_years = as.character(d$n_years[1]), pw(d)))))

options(width = 250)
cat("\n=== 10h scenario C: the rotation with the residual fixed at one, the marker test ===\n")
cat("\nthe free fit against the rotated arm of 10f: largest absolute difference in the estimate", signif(h_check[1], 3),
    ", largest relative difference in the SE", signif(h_check[2], 3), "\n")
cat("\nnull datasets: size (Wald t,", DF_T, "df), model variance over empirical variance, SE against single-stage\n")
print(format(h_size, digits = 4), row.names = FALSE)
cat("\ncausal datasets: power at", round(bonf_scan, 3), "(percent), paired counts, estimate and SE against single-stage\n")
print(format(h_power, digits = 4), row.names = FALSE)
save(HH, h_check, h_size, h_power, file = OUTFILE)
cat("\nsaved", basename(OUTFILE), "|", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
