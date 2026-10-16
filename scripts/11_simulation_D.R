################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script 11, scenario D: does the residual detect an incomplete
##                               genetic model, and what repairs it?
##
##   After rotation the residual variance has a known target of one, so an
##   estimate that differs from one is evidence that the fitted genetic model
##   is missing something. This scenario withholds a known amount of variance
##   and asks four questions of the residual:
##
##     Does it increase when variance is withheld?  fit: additive alone
##     Does restoring the structure that generated
##       the data return it?                        fit: additive + that source
##     Does a plausible but wrong kernel return it? fit: additive + other source
##     Does a generic kernel that knows nothing
##       about what was withheld return it?         fit: additive + WW'
##
##   The independent kernel WW' equals the inverse of Omega and is the rotated
##   image of an effect independent across accessions; it is identifiable only
##   after rotation, because on the original scale it and the residual both
##   multiply the identity and only their sum can be estimated.
##
##   Two structures are withheld in turn, both built in script 07 and described
##   by the correlation of their off-diagonal elements with those of the
##   additive matrix, read from script 07 below. Neither is a claim about a
##   genetic mode of action.
##
##   Everything is measured against the control row, in which nothing is
##   withheld, and not against one: adding a kernel the data do not need
##   lowers the residual, so each fitted model carries its own control offset.
##
##   Output: 11_simD.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(parallel)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "07_sim_primitives.RData"))
t0 <- Sys.time()

## ---------- 1. the design ----------
## Seven truth conditions: nothing withheld, then each structure at three
## levels. The variances are comparable across structures because script 07
## rescales both kernels to the mean diagonal of A.
TRUTH <- rbind(
  data.frame(source = "none", s2 = 0),
  data.frame(source = "epi",  s2 = c(10, 25, 50)),
  data.frame(source = "het",  s2 = c(10, 25, 50)))
TRUTH$k <- seq_len(nrow(TRUTH))          # index below 10, as the seed rule requires
NTRUTH  <- nrow(TRUTH)

## the correlation of each withheld structure with the additive matrix, from
## script 07
KERNEL_COR <- c(epi = kernel_tab$Cor_with_A[kernel_tab$Kernel == "epistasis A#A"],
                het = kernel_tab$Cor_with_A[kernel_tab$Kernel == "heterozygote"])
stopifnot(length(KERNEL_COR) == 2, !anyNA(KERNEL_COR))

ncore <- max(1, min(detectCores() - 1, NCORE_SIM))

P <- list(pheno = pheno, id_idx = ID_IDX, yr_idx = YR_IDX, yr_rep = YR_REP,
          id_yr = ID_YR, pars = PARS, rootA = ROOT_A, rootE = ROOT_E,
          rootD = ROOT_D, Gsub = G_SUB, Q = Q, A = A, n = n, L = L_mat)

## ---------- 2. one replicate ----------
## The four fitted models share one first stage per truth condition. Each
## replicate runs inside try(); a failed replicate is counted and the others
## go on. For the control row nothing is missing, so "matched" fits a
## superfluous epistatic kernel and "wrong" a superfluous heterozygosity one.
rep_D <- function(r) {
  out <- try({
    do.call(rbind, lapply(seq_len(NTRUTH), function(i) {
      src <- TRUTH$source[i]; s2 <- TRUTH$s2[i]
      set.seed(SEED_D + r * 10 + TRUTH$k[i])

      b  <- setNames(rep(0, length(MARKERS)), MARKERS)
      s1 <- stage1_fast(simulate_y(b,
                                   s2_epi = if (src == "epi") s2 else 0,
                                   s2_dom = if (src == "het") s2 else 0,
                                   P = P), P = P)
      Tr <- make_transforms(s1$Omega)$rotated

      ## the structure that generated the data, and the other one
      if (src == "het") {
        K_match <- D_dom; nm_match <- "het"
        K_wrong <- E_had; nm_wrong <- "epi"
      } else {
        K_match <- E_had; nm_match <- "epi"
        K_wrong <- D_dom; nm_wrong <- "het"
      }

      f_add <- ve_kernels(s1, Tr, list(add = A), P = P)
      f_mat <- ve_kernels(s1, Tr, setNames(list(A, K_match), c("add", nm_match)), P = P)
      f_wrg <- ve_kernels(s1, Tr, setNames(list(A, K_wrong), c("add", nm_wrong)), P = P)
      f_gen <- ve_kernels(s1, Tr, list(add = A, ind = diag(n)), P = P)

      data.frame(rep = r, source = src, s2_true = s2,
                 s1_grad = s1$grad, s1_singular = s1$singular, s1_refit = s1$refit,
                 Ve_additive   = as.numeric(f_add["Ve"]),
                 Ve_matched    = as.numeric(f_mat["Ve"]),
                 Ve_wrong      = as.numeric(f_wrg["Ve"]),
                 Ve_generic    = as.numeric(f_gen["Ve"]),
                 Vg_add_only   = as.numeric(f_add["V_add"]),
                 V_recovered   = as.numeric(f_mat[paste0("V_", nm_match)]),
                 V_wrong       = as.numeric(f_wrg[paste0("V_", nm_wrong)]),
                 V_independent = as.numeric(f_gen["V_ind"]),
                 conv_add = as.numeric(f_add["converged"]),
                 conv_mat = as.numeric(f_mat["converged"]),
                 conv_wrg = as.numeric(f_wrg["converged"]),
                 conv_gen = as.numeric(f_gen["converged"]),
                 row.names = NULL)
    }))
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

## ---------- 3. run ----------
cat("scenario D:", NREP_D, "replicates x", NTRUTH, "truth conditions x 4 models =",
    NREP_D * NTRUTH * 4, "variance-component fits\n")

cl <- makeCluster(ncore)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, { library(lme4); library(emmeans); library(rrBLUP) }))
clusterExport(cl, c("P", "TRUTH", "NTRUTH", "MARKERS", "n", "Q", "A",
                    "E_had", "D_dom", "SEED_D",
                    "simulate_y", "stage1_fast", "make_transforms", "ve_kernels"))
## run in batches with the results saved after each, so that an interrupted run
## can be continued (run_in_blocks in script 00)
res <- run_in_blocks(cl, seq_len(NREP_D), rep_D, file.path(OUT, "11_simD_partial.RData"), ncore * 5, t0)
stopCluster(cl)

n_failed <- sum(vapply(res, is.null, logical(1)))
D <- do.call(rbind, res[!vapply(res, is.null, logical(1))])

## ---------- 4. summarise ----------
## means with their standard errors over replicates, by truth condition
mse <- function(v) c(m = mean(v, na.rm = TRUE),
                     se = sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))

cells <- split(D, list(D$source, D$s2_true), drop = TRUE)
raw <- do.call(rbind, lapply(cells, function(x) {
  data.frame(source = x$source[1], s2_true = x$s2_true[1], n = nrow(x),
             Ve_add = mse(x$Ve_additive)["m"], se_add = mse(x$Ve_additive)["se"],
             Ve_mat = mse(x$Ve_matched)["m"],  se_mat = mse(x$Ve_matched)["se"],
             Ve_wrg = mse(x$Ve_wrong)["m"],    se_wrg = mse(x$Ve_wrong)["se"],
             Ve_gen = mse(x$Ve_generic)["m"],  se_gen = mse(x$Ve_generic)["se"],
             V_rec  = mean(x$V_recovered), V_rec_zero = mean(x$V_recovered < 1e-3),
             V_wrg  = mean(x$V_wrong),     V_wrg_zero = mean(x$V_wrong < 1e-3),
             V_ind  = mean(x$V_independent), V_ind_zero = mean(x$V_independent < 1e-3),
             row.names = NULL)
}))
raw <- merge(TRUTH[, c("source", "s2")], raw,
             by.x = c("source", "s2"), by.y = c("source", "s2_true"), all.x = TRUE)
raw <- raw[order(match(paste(raw$source, raw$s2), paste(TRUTH$source, TRUTH$s2))), ]

ctrl <- raw[raw$source == "none", ]

d_tab <- data.frame(
  Withheld = raw$source, Variance = raw$s2, n = raw$n,
  Ve_additive = round(raw$Ve_add, 4), se_add = round(raw$se_add, 4),
  Ve_matched  = round(raw$Ve_mat, 4), se_mat = round(raw$se_mat, 4),
  Ve_wrong    = round(raw$Ve_wrg, 4), se_wrg = round(raw$se_wrg, 4),
  Ve_generic  = round(raw$Ve_gen, 4), se_gen = round(raw$se_gen, 4),
  row.names = NULL)

## The excess is measured against the control row of the same fitted model.
exc_tab <- data.frame(
  Withheld = raw$source, Variance = raw$s2,
  excess_additive = round(raw$Ve_add - ctrl$Ve_add, 4),
  excess_matched  = round(raw$Ve_mat - ctrl$Ve_mat, 4),
  excess_wrong    = round(raw$Ve_wrg - ctrl$Ve_wrg, 4),
  excess_generic  = round(raw$Ve_gen - ctrl$Ve_gen, 4),
  row.names = NULL)
## repaired_* divides one difference by another, so it is reported for every
## cell but checked only where the additive excess reaches REPAIR_MIN_EXCESS.
REPAIR_MIN_EXCESS <- 0.5

ea <- exc_tab$excess_additive
exc_tab$repaired_matched <- ifelse(abs(ea) < 0.05, NA, round(1 - exc_tab$excess_matched / ea, 3))
exc_tab$repaired_wrong   <- ifelse(abs(ea) < 0.05, NA, round(1 - exc_tab$excess_wrong   / ea, 3))
exc_tab$repaired_generic <- ifelse(abs(ea) < 0.05, NA, round(1 - exc_tab$excess_generic / ea, 3))

## A component that is truly zero lies on the boundary in a large share of
## replicates, so its mean and its share at zero are both reported.
comp_tab <- data.frame(
  Withheld = raw$source, Variance = raw$s2,
  V_recovered = round(raw$V_rec, 3), pct_zero_rec = round(100 * raw$V_rec_zero, 1),
  V_wrong     = round(raw$V_wrg, 3), pct_zero_wrg = round(100 * raw$V_wrg_zero, 1),
  V_indep     = round(raw$V_ind, 3), pct_zero_ind = round(100 * raw$V_ind_zero, 1),
  row.names = NULL)

## ---------- 4b. what the real-data null could have hidden ----------
## Script 04 fits the additive model to the real panel and reports a residual
## below one with a profile interval that covers one. Scenario D runs on the
## same design with the same kernels, so its residuals are comparable with
## 04's, up to the Kenward-Roger adjustment of 04's Omega (four parts in ten
## thousand, script 02). Two rates are reported:
##   power_cal  share of replicates whose residual exceeds the 95th percentile
##              of this scenario's control cell, a test at the 5 percent level
##              of the null that nothing is withheld; its 5 percent in the
##              control row sets the threshold and is not a finding
##   cross_04   share whose residual exceeds the upper end of the interval of
##              04, a crossing rate against one sample's bound; it is not the
##              rejection rule of 04 and is calibrated to no nominal level
## The bound and the residual of 04 are read from its saved output, into an
## environment of their own because 04 also saves VE_BOUNDARY and a kernel_tab.
e04 <- new.env()
load(file.path(OUT, "04_variance.RData"), envir = e04)
UPPER_04 <- as.numeric(e04$ci_add$upper)
VE_04    <- as.numeric(e04$fits$additive$vc[2])
rm(e04)

ctl_ve  <- D$Ve_additive[D$source == "none"]
Q95_CTL <- as.numeric(quantile(ctl_ve, 0.95))

sens <- do.call(rbind, lapply(split(D, list(D$source, D$s2_true), drop = TRUE),
  function(x) data.frame(
    source = x$source[1], s2 = x$s2_true[1], n = nrow(x),
    mean_ve = round(mean(x$Ve_additive), 3),
    sd_ve   = round(sd(x$Ve_additive), 3),
    power_cal = round(100 * mean(x$Ve_additive > Q95_CTL), 1),
    cross_04  = round(100 * mean(x$Ve_additive > UPPER_04), 1),
    row.names = NULL)))
sens <- sens[order(match(paste(sens$source, sens$s2),
                         paste(TRUTH$source, TRUTH$s2))), ]

sens_tab <- data.frame(
  Withheld = sens$source, Variance = sens$s2, n = sens$n,
  Ve_mean = sens$mean_ve, Ve_sd = sens$sd_ve,
  power_cal = sens$power_cal, cross_04 = sens$cross_04, row.names = NULL)

## ---------- 5. diagnostics ----------
s1_all <- unique(D[, c("rep", "source", "s2_true", "s1_grad", "s1_singular", "s1_refit")])
n_stage1         <- nrow(s1_all)
n_stage1_stalled <- sum(s1_all$s1_grad > 1e-3, na.rm = TRUE)
n_stage1_refit   <- sum(s1_all$s1_refit)
max_s1_grad      <- max(s1_all$s1_grad, na.rm = TRUE)
pct_singular     <- 100 * mean(s1_all$s1_singular)

## A residual below VE_BOUNDARY is a boundary solution, not a failure, as in
## scripts 04 and 05; a fit the optimiser flags away from the boundary counts
## as a failure, and the range of the flagged residuals is printed.
VE_BOUNDARY <- 1e-3

fit_status <- function(ve, cv) {
  bnd <- ve < VE_BOUNDARY
  list(boundary = sum(bnd), failed = sum(cv == 0 & !bnd), n = length(ve))
}
st <- list(additive = fit_status(D$Ve_additive, D$conv_add),
           matched  = fit_status(D$Ve_matched,  D$conv_mat),
           wrong    = fit_status(D$Ve_wrong,    D$conv_wrg),
           generic  = fit_status(D$Ve_generic,  D$conv_gen))

fit_tab <- data.frame(
  Model    = names(st),
  n        = vapply(st, function(x) x$n, numeric(1)),
  boundary = vapply(st, function(x) x$boundary, numeric(1)),
  failed   = vapply(st, function(x) x$failed, numeric(1)),
  row.names = NULL)
fit_tab$pct_boundary <- round(100 * fit_tab$boundary / fit_tab$n, 1)

n_fit_failed   <- sum(fit_tab$failed)
n_fit_boundary <- sum(fit_tab$boundary)

## the range of the flagged residuals
all_ve   <- c(D$Ve_additive, D$Ve_matched, D$Ve_wrong, D$Ve_generic)
all_conv <- c(D$conv_add,    D$conv_mat,   D$conv_wrg,  D$conv_gen)
flag_max <- if (any(all_conv == 0)) max(all_ve[all_conv == 0]) else NA_real_
flag_min <- if (any(all_conv == 0)) min(all_ve[all_conv == 0]) else NA_real_
conv_below <- 100 * mean(all_ve[all_conv == 1] < VE_BOUNDARY)

## ---------- 6. checks ----------
## The checks test the shape of the response against the control row, not its
## level.
epi <- exc_tab[exc_tab$Withheld == "epi", ]
het <- exc_tab[exc_tab$Withheld == "het", ]
gate_rep <- exc_tab$excess_additive >= REPAIR_MIN_EXCESS

checks <- data.frame(
  Check = c("no replicate failed",
            "no first stage stalled after refitting",
            "every non-boundary variance-component fit converged",
            "every truth condition has its full complement of replicates",
            "the residual rises with the variance withheld, epistatic structure",
            "the residual rises with the variance withheld, heterozygote structure",
            "restoring the right structure removes most of the rise, where there is a rise",
            "no missing values in the residual table",
            "detection rises with the variance withheld, epistatic structure",
            "detection rises with the variance withheld, heterozygote structure"),
  Result = c(n_failed == 0,
             n_stage1_stalled == 0,
             n_fit_failed == 0,
             all(raw$n == NREP_D),
             all(diff(epi$excess_additive) > 0),
             all(diff(het$excess_additive) > 0),
             all(exc_tab$repaired_matched[gate_rep] > 0.5),
             !anyNA(d_tab[, c("Ve_additive", "Ve_matched", "Ve_wrong", "Ve_generic")]),
             all(diff(sens_tab$power_cal[sens_tab$Withheld == "epi"]) > 0),
             all(diff(sens_tab$power_cal[sens_tab$Withheld == "het"]) > 0)))

## ---------- 7. report ----------
cat("\n=== 11 scenario D: the residual as a specification test ===\n")
cat("\nthe two withheld structures, by how much the additive kernel can absorb:\n")
cat("   epi  the epistatic kernel, correlation", KERNEL_COR["epi"], "with the additive matrix\n")
cat("   het  the heterozygote kernel, correlation", KERNEL_COR["het"], "\n")

cat("\nresidual variance under four fitted models, with standard errors\n")
print(d_tab, row.names = FALSE)
cat("   additive: the structure is withheld from the model\n")
cat("   matched:  the structure that generated the data is restored\n")
cat("   wrong:    the other structure is fitted instead\n")
cat("   generic:  WW' is fitted, which knows nothing about what was withheld\n")
cat("   the row with nothing withheld is the control, and every fitted model\n")
cat("   carries its own control offset because adding an unneeded kernel lowers\n")
cat("   the residual\n")

cat("\nexcess residual over the control row of the same fitted model\n")
print(exc_tab, row.names = FALSE)
cat("   repaired_*: share of the additive model's excess that the repair removes\n")
cat("   the check on repaired_matched covers the", sum(gate_rep), "cells whose additive",
    "excess reaches", REPAIR_MIN_EXCESS, "\n")
below_gate <- exc_tab[!gate_rep & exc_tab$Withheld != "none", ]
if (nrow(below_gate) > 0)
  cat("   reported but not asserted on, the additive excess being below the gate:",
      paste0(below_gate$Withheld, " s2 = ", below_gate$Variance, " (",
             below_gate$excess_additive, ")", collapse = ", "), "\n")
cat("   the control row is blank by construction, its excess being zero\n")

cat("\nthe extra variance component, mean and share resting at zero\n")
print(comp_tab, row.names = FALSE)
cat("   a component that is truly zero rests on the boundary in a large share of\n")
cat("   replicates, so the mean and the share at zero say different things\n")

cat("\nwhat the real-data null could have hidden\n")
print(sens_tab, row.names = FALSE)
cat("   04 additive residual", round(VE_04, 4), ", profile upper bound", round(UPPER_04, 4),
    "| the control 95th percentile here", round(Q95_CTL, 4), "\n")
cat("   power_cal: share above the control 95th percentile here, a 5 percent",
    "test of the null; quote this column\n")
cat("              its 5 percent in the control row sets the threshold and is",
    "not a finding\n")
cat("   cross_04:  share above the bound 04 produced, a crossing rate against",
    "one sample;\n")
cat("              not the rule 04 would have rejected, and calibrated to no",
    "nominal level\n")
cat("   neither inverts each replicate own profile interval, which is the",
    "strict power of that rule\n")

cat("\nfirst-stage fits:", n_stage1, "| refitted:", n_stage1_refit,
    "| still stalled:", n_stage1_stalled,
    "| largest scaled gradient:", signif(max_s1_grad, 3), "\n")
cat("singular first stages:", round(pct_singular, 1),
    "percent, expected since the replicate-within-year variance is",
    round(PARS$s2_rep, 4), "\n")
cat("\nvariance-component fits by fitted model\n")
print(fit_tab, row.names = FALSE)
cat("   boundary means the residual came to rest below", VE_BOUNDARY,
    ", which is not a failure;\n")
cat("   failed means flagged by the optimiser away from the boundary\n")
if (any(all_conv == 0)) {
  cat("   flagged residuals span", signif(flag_min, 3), "to", signif(flag_max, 3),
      "over", sum(all_conv == 0), "of", length(all_ve),
      "fits, against a threshold of", VE_BOUNDARY, "\n")
} else {
  cat("   no fit was flagged by the optimiser in", length(all_ve), "fits\n")
}
cat("  ", round(conv_below, 1), "percent of converged fits lie below the threshold\n")
cat("replicates failed outright:", n_failed, "of", NREP_D, "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

## saved before the checks are enforced, so that a failing check leaves the
## run on disk
save(D, d_tab, exc_tab, comp_tab, raw, ctrl, TRUTH, checks, fit_tab,
     n_failed, n_stage1, n_stage1_stalled, n_stage1_refit, max_s1_grad,
     pct_singular, n_fit_failed, n_fit_boundary, VE_BOUNDARY, flag_max, flag_min,
     sens_tab, UPPER_04, VE_04, Q95_CTL, REPAIR_MIN_EXCESS, gate_rep, KERNEL_COR,
     file = file.path(OUT, "11_simD.RData"))

cat("\nsaved 11_simD.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
unlink(file.path(OUT, "11_simD_done.RData"))   # the raw replicate results, kept by run_in_blocks until this save

if (!all(checks$Result)) stop("scenario D checks failed")
