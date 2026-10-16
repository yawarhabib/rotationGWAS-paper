################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 08b, scenario A: the rotated residual variances near zero
##
##   Reads the saved output of script 08 and summarises the rotated residual
##   variance under the global null. Scenario A fits the rotated model with
##   rrBLUP::mixed.solve (fit_arm in script 07), which searches the ratio of
##   residual to genetic variance on a linear scale and cannot resolve a
##   residual variance below about 0.004 at a genetic variance near 70;
##   estimates below 0.01 lie at that floor and are boundary solutions.
##
##   Output: 08b_simA_residual_near_zero.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
e <- new.env(); load(file.path(OUT, "08_simA.RData"), envir = e)
d <- e$A1[e$A1$arm == "rotated", c("rep", "Ve")]
v <- d$Ve[!duplicated(d$rep)]                        # one residual per replicate (global null)
low <- v[v < 0.01]
near_zero <- data.frame(replicates = length(v), pct_below_0.01 = 100 * mean(v < 0.01),
                        n_below_0.01 = length(low), smallest = min(low), largest_below_0.01 = max(low),
                        next_above = min(v[v >= 0.01]), below_1e3 = sum(v < 1e-3),
                        mean = mean(v), median = median(v), q025 = unname(quantile(v, 0.025)), q975 = unname(quantile(v, 0.975)))
print(format(near_zero, digits = 4), row.names = FALSE)
cat("\nthe", length(low), "estimates below 0.01 lie between", signif(min(low), 2), "and", signif(max(low), 2),
    "(the floor of rrBLUP::mixed.solve's search); the next value is", signif(min(v[v >= 0.01]), 2), "\n")
stopifnot(near_zero$replicates == 1000, round(near_zero$pct_below_0.01, 1) == 5.9,
          min(low) >= 0.004, max(low) <= 0.008, min(v[v >= 0.01]) > 0.01)
save(near_zero, low, file = file.path(OUT, "08b_simA_residual_near_zero.RData"))
