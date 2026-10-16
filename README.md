# Exact weighting in two-stage genome-wide association analysis: when does it matter?

The scripts and the data that reproduce the numbers, tables and figures of the
paper. Every script is run from the repository root and writes into
`results/`; the imputation check writes into its own folder.

## Citation

The DOI of this repository is [to be added when Zenodo issues it].

If you use these scripts or the data, please cite the paper

[the citation of the paper, to be added]

## Authors

Yawar Habib and Hans-Peter Piepho

## Contact

- Yawar Habib, habib.yawr.10@gmail.com, for the scripts and the data
- Hans-Peter Piepho, hans-peter.piepho@uni-hohenheim.de, the corresponding author, for the paper

## License

GPL (>= 2), for the scripts and the data. The full text is in `LICENSE`.

## Contents

- `data/` the trial phenotypes (year, replication, genotype, identifier in the
  marker data, vegetation period), the genotypes of the 160 genotyped lines at
  59,206 markers after quality control (PLINK `--recode A` format and the map;
  159 of the lines have phenotype records and enter the analysis), the three
  principal components, and, for the sensitivity analysis of script S1, the
  65,905 markers before the heterozygosity filter with their observed
  heterozygosity
- `scripts/` the analysis, in the order of the table below; `00_config.R`
  holds the paths, the constants and the functions the scripts share
- `imputation_check/` the sensitivity analysis to the imputation of missing
  calls, with its own README
- `run_all.R` runs everything in order

## Requirements

R 4.4.0 or later with the packages lme4, lme4qtl, emmeans, pbkrtest, rrBLUP,
sommer, glmmTMB, MASS, dplyr, tidyr, data.table, genio, ggplot2, patchwork,
officer and writexl. lme4qtl is installed from GitHub
(`remotes::install_github("variani/lme4qtl")`); the others from CRAN. The
imputation check needs Java and Beagle 5.5 (see `imputation_check/README.md`).

The results were produced with R 4.4.0 and lme4 1.1.36, lme4qtl 0.2.2,
emmeans 1.10.6, pbkrtest 0.5.3, rrBLUP 4.6.3, sommer 4.4.3, glmmTMB 1.1.10,
Matrix 1.7.1, MASS 7.3.64, dplyr 1.1.4, tidyr 1.3.1, data.table 1.18.2.1,
genio 1.1.2, ggplot2 4.0.3, patchwork 1.3.0, officer 0.6.7 and writexl 1.5.1,
with Java 20 and Beagle 5.5 (27Feb25.75f).

The simulations run on six cores by default; the environment variable
`SIM_CORES` sets another number. The number of cores changes the timing, not
the results.

## Running

From the repository root:

    Rscript run_all.R

or one script at a time, in the order below, for example
`Rscript scripts/01_data_preparation.R`. Each script writes its results file
into `results/` and prints its tables and checks; a script stops if one of its
checks fails.

| Script | What it does | Reads | Time |
|---|---|---|---|
| 01_data_preparation | phenotypes, genotypes, relationship matrix, principal components | data | under a minute |
| 02_stage1_and_rotation | first stage, Omega, the three weighting schemes | 01 | under a minute |
| 03_equivalence_and_invariance | equivalence and invariance checks at the test marker | 01, 02 | under a minute |
| 04_variance_components | the residual variance on the rotated scale, intervals, bootstrap | 01, 02 | 7 min |
| 05_genome_scan | the genome scan under each scheme | 01, 02, 04 | 3 min |
| 06_weighting_comparison | the three schemes against a single-stage analysis | 01, 02, 05 | 19 min |
| 05b_trial_refit_residual_one | the refit with the residual at one on the observed trial | 01, 02, 05, 06 | seconds |
| 06b_equivalence_all_markers | the equivalence at every marker of the comparison set | 01, 02, 06 | seconds |
| 06c_equivalence_fixed_components | the equivalence with the variance components fixed | 01, 02, 06 | seconds |
| 07_simulation_functions | the functions and objects of the simulations | 01 | seconds |
| 07b_kr_factor_simulated_designs | the Kenward-Roger factor on the simulated designs | 01, 07 | seconds |
| 08_simulation_A, 08b | scenario A: is the estimator correct? | 01, 07 | 31 min |
| 09_simulation_B | scenario B: the size of the block effects | 01, 07 | 67 min |
| 10_simulation_C | scenario C: accessions per environment | 01, 07 | 69 min |
| 10b | scenario C, the single-stage analysis of the same data | 01, 07, 10 | 4 h |
| 10d to 10m | scenario C, further analyses of the same data | 01, 07, 10, 10b, 10d to 10f | 4.6 h together |
| 11_simulation_D | scenario D: the residual as a specification test | 01, 04, 07 | 9 h |
| 11b | scenario D, the one-sided test | 01, 07, 11 | 30 min |
| S1 to S4 | the supporting analyses | 01 to 07, 09, S3 | 15 min together |
| F1 | the tables and figures of the paper | 01, 05, 06, 10, 10k, 11, 11b, S2 | seconds |
| imputation_check 01 to 04 | the imputation of missing calls with Beagle | 01, 02 | 20 min |
| F2, F3 | the supporting tables and figures | 01, 03, 05, 06, S1, F1, imputation_check | seconds |

Times are for six cores; the whole run takes about a day.

## Output

`results/` holds the `.RData` files with the objects the scripts save, the logs in `results/logs/`, the tables of the paper in
`results/tables/`, the figures in `results/figures/`, and the supporting
tables and figures as `Supporting_Information_Tables.xlsx` and
`Supporting_Information_Figures.docx`.
