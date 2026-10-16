################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10l, scenario C: the distance of the P3D scan from the threshold
##   where its decision differs from the single-stage analysis
##
##   Markers within WINDOW units of -log10 p of the threshold are refitted with
##   the residual fixed at one (scripts 05b and 10h). For that rule to reach
##   every dataset in which the P3D scan and the single-stage analysis decide
##   differently, the scan's -log10 p must lie within the window there. This
##   script reads the stored results and measures it: causal datasets at the
##   genome-wide threshold (script 10, rotated scheme, against 10b), null
##   datasets at the 5% and 1% levels and the genome-wide threshold (the null
##   scans of 10e, which are the datasets of script 10, against 10b). The
##   distance is measured from the threshold of each decision.
##
##   Output: 10l_simC_refit_window.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
load(file.path(OUT, "07_sim_primitives.RData"))
e10  <- new.env(); load(file.path(OUT, "10_simC.RData"),            envir = e10)
e10b <- new.env(); load(file.path(OUT, "10b_simC_single_stage.RData"), envir = e10b)
e10e <- new.env(); load(file.path(OUT, "10e_simC_null_by_marker.RData"), envir = e10e)

WINDOW <- 1.5                                          # units of -log10 p on either side of the threshold
thr_gw <- e10$bonf_scan                                 # the genome-wide threshold of the scan
idx    <- match(CAUSAL[1], MARKERS)                     # 10e stores the markers by position

## the same datasets: the null scans of 10e reproduce the rotated residual of
## script 10 exactly
chk <- merge(e10$C_nl[e10$C_nl$arm == "rotated", c("rep", "n_years", "Ve")],
             e10e$DS[e10e$DS$arm == "rotated", c("rep", "n_years", "Ve")], by = c("rep", "n_years"))
N_PAIRS <- NREP_C * length(e10$YEARS_GRID)                 # replicates x design levels
stopifnot(nrow(chk) == N_PAIRS, max(abs(chk$Ve.x - chk$Ve.y)) == 0)

ca <- merge(e10$C_ca[e10$C_ca$arm == "rotated", c("rep", "n_years", "neglog10p")],
            setNames(e10b$SS[e10b$SS$kind == "causal", c("rep", "n_years", "neglog10p")], c("rep", "n_years", "ss")), by = c("rep", "n_years"))
nl <- merge(e10e$MK[e10e$MK$arm == "rotated" & e10e$MK$marker == idx, c("rep", "n_years", "neglog10p")],
            setNames(e10b$SS[e10b$SS$kind == "null", c("rep", "n_years", "neglog10p")], c("rep", "n_years", "ss")), by = c("rep", "n_years"))
stopifnot(nrow(ca) == N_PAIRS, nrow(nl) == N_PAIRS)

one <- function(m, thr, label) {
  dis <- m[(m$neglog10p > thr) != (m$ss > thr), ]
  data.frame(datasets = label, threshold = thr, n = nrow(m), discordant = nrow(dis),
             scan_distance_max = if (nrow(dis)) max(abs(dis$neglog10p - thr)) else NA_real_,
             within_window = if (nrow(dis)) sum(abs(dis$neglog10p - thr) <= WINDOW) else 0L,
             single_stage_distance_max = if (nrow(dis)) max(abs(dis$ss - thr)) else NA_real_,
             abs_difference_median = median(abs(m$neglog10p - m$ss)), abs_difference_q99 = unname(quantile(abs(m$neglog10p - m$ss), 0.99)),
             abs_difference_max = max(abs(m$neglog10p - m$ss)))
}
window_tab <- rbind(one(ca, thr_gw, "causal, genome-wide threshold"),
                    one(nl, -log10(0.05), "null, 5% level"), one(nl, -log10(0.01), "null, 1% level"), one(nl, thr_gw, "null, genome-wide threshold"))

## distinct datasets among the discordant decisions
d_ca <- ca[(ca$neglog10p > thr_gw) != (ca$ss > thr_gw), c("rep", "n_years")]
d_n5 <- nl[(nl$neglog10p > -log10(0.05)) != (nl$ss > -log10(0.05)), c("rep", "n_years")]
d_n1 <- nl[(nl$neglog10p > -log10(0.01)) != (nl$ss > -log10(0.01)), c("rep", "n_years")]
null_distinct <- nrow(unique(rbind(d_n5, d_n1)))
totals <- data.frame(discordant_decisions = sum(window_tab$discordant), discordant_datasets = nrow(d_ca) + null_distinct,
                     causal_datasets = nrow(d_ca), null_datasets = null_distinct, null_overlap_5_and_1 = nrow(d_n5) + nrow(d_n1) - null_distinct,
                     scan_distance_max_overall = max(window_tab$scan_distance_max, na.rm = TRUE), all_within_window = all(window_tab$within_window == window_tab$discordant))

## by design level, causal datasets
by_level <- do.call(rbind, lapply(split(ca, ca$n_years), function(x) data.frame(
  n_years = x$n_years[1], accessions_per_year = e10$panel_tab$accessions_per_year[match(x$n_years[1], e10$panel_tab$years)], n = nrow(x),
  abs_difference_median = median(abs(x$neglog10p - x$ss)), abs_difference_q99 = unname(quantile(abs(x$neglog10p - x$ss), 0.99)),
  discordant = sum((x$neglog10p > thr_gw) != (x$ss > thr_gw)))))

options(width = 200)
cat("=== 10l: the P3D scan against the single-stage analysis where the decisions differ (window", WINDOW, "units of -log10 p) ===\n\n")
print(format(window_tab, digits = 4), row.names = FALSE)
cat("\n"); print(format(totals, digits = 4), row.names = FALSE)
cat("\ncausal datasets by design level\n"); print(format(by_level, digits = 3), row.names = FALSE)
cat(sprintf("\nIn words: the scan and the single-stage analysis disagreed on the marker in %d datasets (%d causal at the genome-wide threshold, %d null at the 5%% and 1%% levels, none at the genome-wide threshold);\n",
            totals$discordant_datasets, totals$causal_datasets, totals$null_datasets))
cat(sprintf("in every one of them the scan's -log10 p was within %.2f of the threshold of that decision; a window of %g units covers %s.\n",
            totals$scan_distance_max_overall, WINDOW,
            if (isTRUE(totals$all_within_window)) "them all" else sprintf("%d of them", sum(window_tab$within_window))))
save(window_tab, totals, by_level, WINDOW, file = file.path(OUT, "10l_simC_refit_window.RData"))
## saved before the check is reported
cat(sprintf("\nCHECK every discordant dataset within %g units of its threshold: %s\n", WINDOW,
            if (isTRUE(totals$all_within_window)) "PASS" else "FAIL (widen the window in the text)"))
cat("\nsaved 10l_simC_refit_window\n")
