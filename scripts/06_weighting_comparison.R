################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##              Script 06: the three weighting schemes against
##                        a true single-stage analysis
##
##   Every scheme is compared with a single-stage fit to the plot data over a
##   set of markers spanning the significance range. Denominator degrees of
##   freedom come from Kenward-Roger, and the model without a marker is fitted
##   once per scheme and reused.
##
##   Least squares on the adjusted means is unbiased under every scheme, so an
##   incorrect weighting shows in the standard error and in the test, not in
##   the estimate.
##
##   Output: 06_comparison.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(lme4)
library(lme4qtl)
library(pbkrtest)

load(file.path(OUT, "01_data.RData"))
load(file.path(OUT, "02_rotation.RData"))
load(file.path(OUT, "05_scan.RData"))
t0 <- Sys.time()

sym <- function(M) (M + t(M)) / 2
## A boundary message (a variance component at zero) is not a failure.
converged <- function(m) all(grepl("boundary", m@optinfo$conv$lme4$messages))

## The gradient scaled by the Hessian, reported for every fit.
scaled_grad <- function(m) {
  d <- m@optinfo$derivs
  if (is.null(d)) NA_real_ else max(abs(solve(d$Hessian, d$gradient)))
}

BOBYQA <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

## A fit that lme4 flags, or whose scaled gradient exceeds GRAD_TOL, is refitted
## with the optimiser bobyqa.
GRAD_TOL <- 1e-3
fit_relmat <- function(form, data, kernel) {
  m <- relmatLmer(form, data = data, relmat = list(Id = kernel), REML = TRUE)
  refit <- !converged(m) || isTRUE(scaled_grad(m) > GRAD_TOL)
  if (refit)
    m <- relmatLmer(form, data = data, relmat = list(Id = kernel), REML = TRUE,
                    control = BOBYQA)
  list(m = m, refit = refit, ok = converged(m), grad = scaled_grad(m))
}

## The kernel and the plot data must name the accessions identically.
stopifnot(identical(dimnames(A), list(rownames(A), rownames(A))),
          identical(levels(pheno$Id), rownames(A)))

## ---------- 1. marker selection ----------
## Every marker significant under any scheme, plus markers drawn from each bin
## of -log10p below the threshold: N_PER_BIN per bin, and a bin holding fewer
## is taken whole.
N_PER_BIN <- 50
set.seed(SEED_MISC + 3)
strong <- Reduce(union, sig)

rr    <- scans$rotated$res
stopifnot(bonf > 5)              # the bin edges below assume it
valid <- !is.na(rr$neglog10p) & rr$neglog10p < bonf & !(rr$SNP %in% strong)
bins  <- cut(rr$neglog10p[valid], breaks = c(0, 1, 2, 3, 4, 5, bonf),
             include.lowest = TRUE)
pool  <- split(rr$SNP[valid], bins)
nulls <- unlist(lapply(pool, function(p)
  if (length(p) <= N_PER_BIN) p else sample(p, N_PER_BIN)), use.names = FALSE)

sel <- unique(c(strong, nulls))

## section 6 reports the test marker, so it must be in the set
stopifnot(TEST_SNP %in% sel)

sel_tab <- data.frame(
  Quantity = c("significant under any scheme", "nulls drawn", "total"),
  Value = c(length(strong), length(nulls), length(sel)))
bin_tab <- data.frame(Bin = names(pool), In_scan = sapply(pool, length),
                      Drawn = sapply(pool, function(p) min(length(p), N_PER_BIN)),
                      row.names = NULL)

## ---------- 2. plot-level frame ----------
dat <- pheno
dat$PC1 <- Q[as.character(dat$Id), 1]
dat$PC2 <- Q[as.character(dat$Id), 2]
dat$PC3 <- Q[as.character(dat$Id), 3]
for (s in sel) dat[[s]] <- G[as.character(dat$Id), s]

## ---------- 3. single-stage reference ----------
## Each fit runs inside try(), so that a failed fit is recorded in the table
## and the others go on.
form_null <- as.formula(paste(TRAIT,
  "~ 1 + PC1 + PC2 + PC3 + Year + (1|Year:Replication) + (1|Id:Year) + (1|Id)"))
## relmatLmer factorises the kernel, and A is singular, so the single-stage fits
## receive A + cJ, J the matrix of ones (script 00). At plot level the added
## term lies in the space of the intercept, so the marker effects, standard
## errors, variance components, REML criterion and Kenward-Roger tests are those
## of A itself.
A_ss  <- A + ONES_C * matrix(1, n, n)
f0_ss <- fit_relmat(form_null, dat, A_ss)
m0_ss <- f0_ss$m

cat("single-stage fits for", length(sel), "markers ...\n")
## The fits are spread over cores; the number of cores changes the timing, not
## the result.
ncore <- max(1, min(parallel::detectCores() - 1, 6))
cl <- parallel::makeCluster(ncore)
invisible(parallel::clusterEvalQ(cl, suppressMessages({library(lme4); library(lme4qtl); library(pbkrtest)})))
parallel::clusterExport(cl, c("dat", "A_ss", "m0_ss", "TRAIT", "fit_relmat", "converged",
                              "scaled_grad", "BOBYQA", "GRAD_TOL"))
ss <- do.call(rbind, parallel::parLapply(cl, sel, function(s) {
  out <- try({
    f <- as.formula(paste(TRAIT, "~ 1 +", s,
      "+ PC1 + PC2 + PC3 + Year + (1|Year:Replication) + (1|Id:Year) + (1|Id)"))
    ff <- fit_relmat(f, dat, A_ss)
    kr <- KRmodcomp(ff$m, m0_ss)$stats
    cf <- summary(ff$m)$coefficients[s, ]
    data.frame(SNP = s, beta = cf["Estimate"], se = cf["Std. Error"],
               df = kr$ddf, neglog10p = -log10(kr$p.value),
               converged = ff$ok, refit = ff$refit, grad = ff$grad,
               row.names = NULL)
  }, silent = TRUE)
  if (inherits(out, "try-error"))
    data.frame(SNP = s, beta = NA_real_, se = NA_real_, df = NA_real_,
               neglog10p = NA_real_, converged = FALSE, refit = NA,
               grad = NA_real_, row.names = NULL)
  else out
}))

## ---------- 4. the three two-stage schemes, same fitting function ----------
## The transformed kernel T A T' is singular, so each scheme receives
## T (A + cJ) T' (with_ones in script 00); the added term is the outer product
## of the transformed intercept and changes nothing reported here.
## A transformed row is a combination of accessions, not a named line, so the
## rows get the labels r001, r002 and so on.
ids <- sprintf("r%03d", 1:n)

two_stage_arm <- function(Tm) {
  Kn <- with_ones(A, Tm)
  dimnames(Kn) <- list(ids, ids)
  base <- data.frame(y   = as.numeric(Tm %*% blues$BLUE),
                     ic  = as.numeric(Tm %*% rep(1, n)),
                     PC1 = (Tm %*% Q)[, 1], PC2 = (Tm %*% Q)[, 2], PC3 = (Tm %*% Q)[, 3],
                     Id  = factor(ids, levels = ids))
  f0 <- fit_relmat(y ~ 0 + ic + PC1 + PC2 + PC3 + (1|Id), base, Kn)
  m0 <- f0$m
  G_sel <- G[, sel, drop = FALSE]      # only these columns go to the cores
  res <- do.call(rbind, parallel::parLapply(cl, sel, function(s) {
    out <- try({
      d <- base; d$SNP <- as.numeric(Tm %*% G_sel[, s])
      ff <- fit_relmat(y ~ 0 + ic + PC1 + PC2 + PC3 + SNP + (1|Id), d, Kn)
      kr <- KRmodcomp(ff$m, m0)$stats
      cf <- summary(ff$m)$coefficients["SNP", ]
      data.frame(SNP = s, beta = cf["Estimate"], se = cf["Std. Error"],
                 df = kr$ddf, neglog10p = -log10(kr$p.value),
                 converged = ff$ok, refit = ff$refit, grad = ff$grad,
                 row.names = NULL)
    }, silent = TRUE)
    if (inherits(out, "try-error"))
      data.frame(SNP = s, beta = NA_real_, se = NA_real_, df = NA_real_,
                 neglog10p = NA_real_, converged = FALSE, refit = NA,
                 grad = NA_real_, row.names = NULL)
    else out
  }))
  list(res = res, conv_null = f0$ok, refit_null = f0$refit, grad_null = f0$grad)
}

cat("two-stage arms ...\n")
arms <- lapply(Tlist, two_stage_arm)
ts   <- lapply(arms, `[[`, "res")
parallel::stopCluster(cl)

## ---------- 5. metrics against the single-stage reference ----------
## A marker is dropped where either analysis returned no finite value; the
## count used is reported with each metric.
metric_one <- function(two, quantity) {
  m <- merge(ss[, c("SNP", quantity)], two[, c("SNP", quantity)],
             by = "SNP", suffixes = c("_ss", "_ts"))
  a <- m[[paste0(quantity, "_ss")]]; b <- m[[paste0(quantity, "_ts")]]
  keep <- is.finite(a) & is.finite(b)
  a <- a[keep]; b <- b[keep]; d <- b - a
  ## Bias is the mean of the paired differences; its standard error over the
  ## markers gives the t statistic against zero.
  se_bias <- sd(d) / sqrt(length(d))
  data.frame(Quantity = quantity, n = length(d),
             MSD = mean(d^2), Bias = mean(d), SE_bias = se_bias,
             t_bias = mean(d) / se_bias, Correlation = cor(b, a))
}

metrics <- do.call(rbind, lapply(names(ts), function(nm)
  cbind(Scheme = nm, do.call(rbind,
    lapply(c("beta", "se", "neglog10p"), function(q) metric_one(ts[[nm]], q))))))
num <- c("MSD", "Bias", "SE_bias", "t_bias", "Correlation")
metrics[, num] <- round(metrics[, num], 4)
metrics <- metrics[order(metrics$Quantity, metrics$Scheme), ]

## ---------- 6. the test marker, four ways ----------
pick <- function(df, col) df[[col]][df$SNP == TEST_SNP]
anchor <- rbind(
  data.frame(Analysis = "single-stage (plot level)",
             beta = pick(ss, "beta"), se = pick(ss, "se"),
             neglog10p = pick(ss, "neglog10p")),
  do.call(rbind, lapply(names(ts), function(nm)
    data.frame(Analysis = paste("two-stage,", nm),
               beta = pick(ts[[nm]], "beta"), se = pick(ts[[nm]], "se"),
               neglog10p = pick(ts[[nm]], "neglog10p")))))
## distance and ratios from the unrounded values
anchor$Distance <- round(abs(anchor$beta - anchor$beta[1]), 4)
anchor$SE_ratio <- round(anchor$se / anchor$se[1], 4)
anchor[, 2:4]   <- round(anchor[, 2:4], 4)

## ---------- checks ----------
n_failed_ss <- sum(!ss$converged)
n_failed_ts <- sapply(ts, function(d) sum(!d$converged))
n_refit     <- sum(ss$refit, na.rm = TRUE) + f0_ss$refit +
               sum(sapply(ts,   function(d) sum(d$refit, na.rm = TRUE))) +
               sum(sapply(arms, function(a) a$refit_null))
## The scaled gradients of every fit: the marker fits under the three schemes,
## the single-stage marker fits, the single-stage null and the null of each
## scheme.
all_grad <- c(ss$grad,
              unlist(lapply(ts, function(d) d$grad)),
              f0_ss$grad,
              sapply(arms, function(a) a$grad_null))
max_grad <- max(all_grad, na.rm = TRUE)
n_fits   <- length(all_grad)
n_used   <- min(metrics$n)

checks <- data.frame(
  Check = c("test marker present in the selected set",
            "null single-stage model converged",
            "every single-stage marker fit converged",
            "null model converged in all three arms",
            "every two-stage marker fit converged",
            "all markers usable in every metric",
            "scaled gradient small at every fit"),
  Result = c(TEST_SNP %in% sel,
             f0_ss$ok,
             n_failed_ss == 0,
             all(sapply(arms, `[[`, "conv_null")),
             all(n_failed_ts == 0),
             n_used == length(sel),
             max_grad < GRAD_TOL))

## ---------- report ----------
cat("\n=== 06 weighting comparison ===\n")
cat("\nmarkers\n"); print(sel_tab, row.names = FALSE); print(bin_tab, row.names = FALSE)
cat("   fits refitted with bobyqa:", n_refit, "\n")
cat("   largest scaled gradient over all", n_fits, "fits:", signif(max_grad, 3), "\n")
cat("\nagreement with the single-stage analysis\n")
print(metrics, row.names = FALSE)
cat("   the beta rows are a consistency check: least squares on the adjusted\n")
cat("   means is unbiased under every scheme, so they cannot separate them.\n")
cat("   the se and neglog10p rows are where a difference must appear.\n")
cat("   t_bias tests the mean paired difference against zero over the markers.\n")
## the share of markers at which a scheme's standard error is below the
## single-stage one
se_below <- sapply(ts, function(d) { m <- merge(ss, d, by = "SNP"); round(100 * mean(m$se.y < m$se.x), 1) })
cat("\nshare of markers with a standard error below the single-stage one (percent)\n")
print(se_below)
cat("\nthe test marker\n"); print(anchor, row.names = FALSE)
cat("\nchecks\n");          print(checks, row.names = FALSE)

## saved before the checks are enforced, so that a failing check leaves the
## fits on disk
save(sel, sel_tab, bin_tab, N_PER_BIN, ss, ts, arms, metrics, anchor, checks, se_below,
     n_failed_ss, n_failed_ts, n_used, n_refit, max_grad, n_fits, all_grad, dat,
     file = file.path(OUT, "06_comparison.RData"))

cat("\nsaved 06_comparison.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

if (!all(checks$Result)) stop("weighting-comparison checks failed")
