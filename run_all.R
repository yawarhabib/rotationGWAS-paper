################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##                     run_all.R: every script, in order
##
##   Runs the scripts of the analysis one after another, each in its own R
##   session, from the repository root, and writes each script's output to
##   results/logs/<script>.log. A script that stops with an error stops the run.
##
##   Usage, from the repository root:
##     Rscript run_all.R            every script
##     Rscript run_all.R 05 06      only the scripts whose names start so
##
##   The imputation check needs Java and the Beagle jar (see
##   imputation_check/README.md); without them the run stops before F2.
##
##   The whole run takes about a day on six cores; README.md gives the time of
##   each script.
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

ROOT <- normalizePath(getwd(), winslash = "/")
if (!dir.exists(file.path(ROOT, "scripts")) || !dir.exists(file.path(ROOT, "data")))
  stop("run from the repository root: the folders data/ and scripts/ were not found in ", ROOT)
dir.create(file.path(ROOT, "results", "logs"), showWarnings = FALSE, recursive = TRUE)

ORDER <- c(
  "scripts/01_data_preparation.R",
  "scripts/02_stage1_and_rotation.R",
  "scripts/03_equivalence_and_invariance.R",
  "scripts/04_variance_components.R",
  "scripts/05_genome_scan.R",
  "scripts/06_weighting_comparison.R",
  "scripts/05b_trial_refit_residual_one.R",
  "scripts/06b_equivalence_all_markers.R",
  "scripts/06c_equivalence_fixed_components.R",
  "scripts/07_simulation_functions.R",
  "scripts/07b_kr_factor_simulated_designs.R",
  "scripts/08_simulation_A.R",
  "scripts/08b_simA_residual_near_zero.R",
  "scripts/09_simulation_B.R",
  "scripts/10_simulation_C.R",
  "scripts/10b_simulation_C_single_stage.R",
  "scripts/10d_simC_efficiency_by_marker.R",
  "scripts/10e_simC_null_by_marker.R",
  "scripts/10f_simC_full_fit_arms.R",
  "scripts/10g_simC_residual_fixed.R",
  "scripts/10h_simC_marker_test_residual_one.R",
  "scripts/10i_simC_smith_weights.R",
  "scripts/10k_simC_asv_share.R",
  "scripts/10l_simC_refit_window.R",
  "scripts/10m_simC_efficiency_medabs.R",
  "scripts/11_simulation_D.R",
  "scripts/11b_simulation_D_one_sided_test.R",
  "scripts/S1_heterozygosity_sweep.R",
  "scripts/S2_weighting_efficiency.R",
  "scripts/S3_independent_optimiser.R",
  "scripts/S4_ones_matrix.R",
  "scripts/F1_tables_and_figures.R",
  "imputation_check/01_write_vcf.R",
  "BEAGLE",
  "imputation_check/02_check_and_refit.R",
  "imputation_check/03_where_the_change_comes_from.R",
  "imputation_check/04_residual_at_one.R",
  "scripts/F2_supporting_tables.R",
  "scripts/F3_supporting_figures.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  keep <- vapply(ORDER, function(s) any(startsWith(basename(s), args)), logical(1))
  ORDER <- ORDER[keep]
}

rscript <- file.path(R.home("bin"), "Rscript")
jar     <- file.path(ROOT, "imputation_check", "beagle.27Feb25.75f.jar")

run_one <- function(script) {
  log <- file.path(ROOT, "results", "logs", sub("\\.R$", ".log", basename(script)))
  cat(format(Sys.time(), "%H:%M:%S"), script, "\n")
  status <- system2(rscript, script, stdout = log, stderr = log)
  if (status != 0) stop("script ", script, " stopped with an error; see ", log)
}

run_beagle <- function() {
  if (!nzchar(Sys.which("java")) || !file.exists(jar))
    stop("the imputation check needs Java and ", basename(jar),
         " in imputation_check/ (see imputation_check/README.md)")
  for (set in c("full", "masked")) {
    cat(format(Sys.time(), "%H:%M:%S"), "Beagle on geno_", set, ".vcf\n", sep = "")
    status <- system2("java", c("-Xmx6g", "-jar", jar,
                                paste0("gt=imputation_check/geno_", set, ".vcf"),
                                paste0("out=imputation_check/imputed_", set),
                                "nthreads=2", "seed=20261002"),
                      stdout = file.path(ROOT, "imputation_check", paste0("beagle_", set, ".log")),
                      stderr = file.path(ROOT, "imputation_check", paste0("beagle_", set, ".log")))
    if (status != 0) stop("Beagle stopped with an error on the ", set, " set")
  }
}

t0 <- Sys.time()
for (step in ORDER) if (step == "BEAGLE") run_beagle() else run_one(step)
cat("done |", round(difftime(Sys.time(), t0, units = "hours"), 2), "hours\n")
