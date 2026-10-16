################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##   Script 10m, scenario C: the median absolute difference between the
##   predicted and the simulated standard-deviation factor, marker by marker
##
##   Reads the saved output of script 10e and computes, for each design level,
##   the median absolute relative difference between the predicted
##   standard-deviation factor of inverse-variance weighting (exact, from the
##   design) and the simulated one over the selected markers.
##
##   Output: 10m_simC_efficiency_medabs.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
e <- new.env(); load(file.path(OUT, "10e_simC_null_by_marker.RData"), envir = e)
w <- e$wide
medabs <- do.call(rbind, lapply(split(w, w$n_years), function(d) data.frame(
  n_years = d$n_years[1], markers = nrow(d),
  median_predicted = median(d$exact_inv), median_simulated = median(d$sim_inv), cor = cor(d$exact_inv, d$sim_inv),
  medabs_rel_pct = 100 * median(abs(d$sim_inv / d$exact_inv - 1)), medabs_factor = median(abs(d$sim_inv - d$exact_inv)),
  max_rel_pct = 100 * max(abs(d$sim_inv / d$exact_inv - 1)))))
options(width = 200)
cat("=== 10m: predicted against simulated standard-deviation factor of inverse-variance weighting, by design level (106 markers) ===\n\n")
print(format(medabs, digits = 4), row.names = FALSE)
i <- which(medabs$n_years == 64)
cat(sprintf("\nAt 7.5 accessions per environment (64 years): median predicted %.3f, median simulated %.3f, correlation %.3f, median absolute relative difference %.3f%% (the Results' 2.0%%)\n",
            medabs$median_predicted[i], medabs$median_simulated[i], medabs$cor[i], medabs$medabs_rel_pct[i]))
save(medabs, file = file.path(OUT, "10m_simC_efficiency_medabs.RData"))
## the quoted values, checked against this run
stopifnot(medabs$markers[i] == 106)
cat(sprintf("\nCHECK against the quoted values at 7.5 accessions: median absolute %.1f (quoted 2.0), predicted %.3f (quoted 1.157), correlation %.3f (quoted 0.957): %s\n",
            medabs$medabs_rel_pct[i], medabs$median_predicted[i], medabs$cor[i],
            if (round(medabs$medabs_rel_pct[i], 1) == 2.0 && round(medabs$median_predicted[i], 3) == 1.157 && round(medabs$cor[i], 3) == 0.957) "same" else "moved, re-quote"))
cat("\nsaved 10m_simC_efficiency_medabs\n")
