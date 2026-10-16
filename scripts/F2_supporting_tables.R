################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script F2: the Supporting Information tables
##
##   Builds Tables S1 to S4 from the saved output of scripts 03, 05, 06 and S1
##   and of the imputation check, so that every value in the supporting workbook
##   is the value a script produced. Computes nothing new: each panel is a
##   selection, a rounding or a relabelling of results already checked.
##
##     S1a  the five fixed-effect coefficients under the three fits      (03)
##     S1b  the four rotations and their deviations from W'W = Omega^-1  (03)
##     S1c  the marker effect with the rotated and the plain intercept   (03)
##     S1d  the single-stage anchor and the three two-stage schemes      (03)
##     S2a  markers passing the threshold under any scheme               (05)
##     S2b  the comparison set of script 06, all quantities, all schemes (06)
##     S3   the heterozygosity-threshold sweep of the quality control    (S1)
##     S4   the missing calls imputed with Beagle 5.5 against the marker mean
##          (imputation_check, scripts 02 and 04)
##
##   Output: results/Supporting_Information_Tables.xlsx and, for inspection,
##           results/tables/S1a.csv ... S4.csv
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(writexl)
t0 <- Sys.time()

BUILD <- file.path(ROOT, "results")
TAB   <- file.path(BUILD, "tables")
XLSX  <- file.path(ROOT, "results/Supporting_Information_Tables.xlsx")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

e01 <- new.env(); load(file.path(OUT, "01_data.RData"),        envir = e01)
e03 <- new.env(); load(file.path(OUT, "03_equivalence.RData"), envir = e03)
e05 <- new.env(); load(file.path(OUT, "05_scan.RData"),        envir = e05)
e06 <- new.env(); load(file.path(OUT, "06_comparison.RData"),  envir = e06)
eS1 <- new.env(); load(file.path(OUT, "S1_het_sweep.RData"),   envir = eS1)
IMP <- file.path(ROOT, "imputation_check")
eI2 <- new.env(); load(file.path(IMP, "02_check_and_refit.RData"), envir = eI2)
eI4 <- new.env(); load(file.path(IMP, "04_residual_at_one.RData"), envir = eI4)

LAB <- c(unweighted = "unweighted", diagonal = "inverse-variance", rotated = "rotated")
sci <- function(x) {                       # 1.4e-15 -> "1.4x10^-15", for the deviation column
  m <- formatC(x, format = "e", digits = 1)
  mant <- sub("e.*$", "", m); ex <- as.integer(sub("^.*e", "", m))
  paste0(mant, "\u00d710^", ex)
}

## ---------- Table S1 ----------
S1a <- e03$gls_tab
names(S1a) <- c("Coefficient", "Least squares, unweighted", "Least squares, rotated",
                "GLS with weight Omega^-1")
S1a[, 2:4] <- round(S1a[, 2:4], 5)

S1b <- e03$inv_tab
S1b$Rotation <- c("Symmetric square root", "Lambda^(-1/2) S'", "Cholesky factor of Omega^-1",
                  "PW, P random orthonormal")
S1b <- data.frame(Rotation = S1b$Rotation,
                  "Genetic variance"  = round(S1b$Vg, 4),
                  "Residual variance" = round(S1b$Ve, 5),
                  "Marker effect"     = round(S1b$beta, 5),
                  "Standard error"    = round(S1b$se, 5),
                  "-log10p"           = round(S1b$lp, 4),
                  "Deviation from W'W = Omega^-1" = sci(e03$validity),
                  check.names = FALSE, stringsAsFactors = FALSE)

S1c <- e03$intercept_tab
S1c <- data.frame("Intercept column" = c("GLS reference, panel (a)", "Rotated column W1",
                                         "Plain column of ones"),
                  "Marker effect"       = round(S1c$Marker_effect, 5),
                  "Difference from GLS" = round(S1c$Difference, 5),
                  check.names = FALSE, stringsAsFactors = FALSE)

S1d <- e03$anchor_tab
S1d <- data.frame(Analysis = c("Single-stage, plot data", "Two-stage, unweighted",
                               "Two-stage, inverse-variance", "Two-stage, rotated"),
                  "Marker effect"  = round(S1d$Marker_effect, 4),
                  "Standard error" = round(S1d$SE, 4),
                  "t"              = round(S1d$t, 4),
                  "Difference from single stage" = round(S1d$Distance_from_single_stage, 4),
                  "SE ratio"       = round(S1d$SE_ratio, 4),
                  "t ratio"        = round(S1d$t_ratio, 4),
                  check.names = FALSE, stringsAsFactors = FALSE)

## ---------- Table S2 ----------
sig <- unique(e05$sig_tab$SNP)
k   <- match(sig, e05$sig_tab$SNP)
S2a <- data.frame(Marker = sig, Chromosome = e05$sig_tab$Chromosome[k],
                  Position = e05$sig_tab$Position[k], stringsAsFactors = FALSE)
for (nm in names(LAB)) {
  r <- e05$scans[[nm]]$res
  S2a[[paste0("-log10p, ", LAB[nm])]] <- round(r$neglog10p[match(sig, r$SNP)], 4)
}
S2a <- S2a[order(S2a$Position), ]

ss <- e06$ss; ts <- e06$ts; map <- e01$map
i  <- match(ss$SNP, map$SNP)
## five decimals, so that the medians of paired differences quoted in the text
## reproduce
S2b <- data.frame(Marker = ss$SNP, Chromosome = map$Chromosome[i], Position = map$Position[i],
                  "Significant under any scheme" = ifelse(ss$SNP %in% sig, "yes", "no"),
                  "Estimate, single-stage" = round(ss$beta, 5),
                  "SE, single-stage"       = round(ss$se, 5),
                  "-log10p, single-stage"  = round(ss$neglog10p, 5),
                  check.names = FALSE, stringsAsFactors = FALSE)
for (nm in names(LAB)) {
  t <- ts[[nm]]; j <- match(ss$SNP, t$SNP)
  S2b[[paste0("Estimate, ", LAB[nm])]] <- round(t$beta[j], 5)
  S2b[[paste0("SE, ", LAB[nm])]]       <- round(t$se[j], 5)
  S2b[[paste0("-log10p, ", LAB[nm])]]  <- round(t$neglog10p[j], 5)
}
S2b <- S2b[order(-S2b$"-log10p, single-stage"), ]

## ---------- Table S3 ----------
## The sweep refits the rotated model with its own optimiser, so its variance
## components at 0.10 differ slightly from Table 2; the caption gives the
## bound and a check below enforces it.
sw <- eS1$sweep_tab
S3 <- data.frame(Threshold = ifelse(sw$Threshold >= 1, "none", formatC(sw$Threshold, format = "f", digits = 2)),
                 "Markers retained"            = sw$Markers,
                 "Mean minor allele frequency" = round(sw$Mean_MAF, 4),
                 "Mean diagonal of A"          = round(sw$Mean_diag, 4),
                 "Rank of A"                   = sw$Eff_rank,
                 "Genetic variance"            = round(sw$Vg, 3),
                 "Residual variance"           = round(sw$Ve, 4),
                 "Genomic inflation factor"    = round(sw$Lambda, 4),
                 check.names = FALSE, stringsAsFactors = FALSE)

## ---------- Table S4 ----------
## The missing calls imputed with Beagle 5.5, set beside the marker mean used in
## the main text. Counts, variance components and windows are those of
## imputation_check/02; the residual test statistic is that of
## imputation_check/04, taken as zero when the estimate is at or below one.
it  <- eI2$tab; isg <- eI2$sig; io <- eI4$out
cnt <- function(sch, imp) it$above_threshold[it$fit == paste(sch, imp)]
win <- function(imp) { s <- isg[grepl(paste0(" ", imp, "$"), isg$fit), ]
  paste0(paste(unique(s$Chromosome), collapse = ", "), ": ", format(min(s$Position), big.mark = ","), "–", format(max(s$Position), big.mark = ",")) }
rot <- function(imp, col) it[[col]][it$fit == paste("rotated", imp)]
stat1 <- function(k) { r <- io[io$kernel == k, ]; if (r$residual_here <= 1) 0 else r$cost_of_one }
f4  <- function(x) formatC(x, format = "f", digits = 4)
f3  <- function(x) formatC(x, format = "f", digits = 3)
col4 <- function(imp, k) c(cnt("rotated", imp), cnt("diagonal", imp), cnt("unweighted", imp), win(imp),
                           f4(mean(diag(if (imp == "in use") e01$A else eI2$Ab))), f3(rot(imp, "Vg")),
                           f4(rot(imp, "Ve")), f4(stat1(k)))
S4 <- data.frame(
  Quantity = c("Markers above the threshold, rotated", "Markers above the threshold, inverse-variance",
               "Markers above the threshold, unweighted", "Chromosome and window of these markers",
               "Mean diagonal of A", "Genetic variance, rotated second stage",
               "Residual variance, rotated second stage", "One-sided test statistic for the residual variance"),
  "Marker mean (main text)" = col4("in use", "A in use"),
  "Beagle 5.5"              = col4("Beagle", "A from Beagle"),
  check.names = FALSE, stringsAsFactors = FALSE)

## ---------- captions ----------
## every number in a caption is taken from the panel it describes, at the
## precision the panel prints
f4   <- function(x) formatC(x, format = "f", digits = 4)
r10  <- S3$Threshold %in% c("none", "0.50", "0.30", "0.20", "0.10")
s3at <- function(col, thr) S3[[col]][S3$Threshold == thr]
cap_S3 <- paste0(
  "Table S3. Sensitivity to the heterozygosity threshold of the quality control. The marker set, the relationship matrix A and the rotated model were rebuilt at each threshold on observed heterozygosity, from no filter through 0.50, 0.30, 0.20, 0.10, 0.05 and 0.02; 0.10 is the threshold applied. From no filter to 0.10, the genomic inflation factor of the rotated scan remained between ",
  f4(min(S3$"Genomic inflation factor"[r10])), " and ", f4(max(S3$"Genomic inflation factor"[r10])),
  ", the rotated residual variance between ", f4(min(S3$"Residual variance"[r10])), " and ", f4(max(S3$"Residual variance"[r10])),
  ", and the rank of A at ", S3$"Rank of A"[1], ". The mean diagonal of A increased as the threshold was lowered, from ",
  f4(s3at("Mean diagonal of A", "none")), " with no filter to ", f4(s3at("Mean diagonal of A", "0.10")), " at 0.10 and ",
  f4(s3at("Mean diagonal of A", "0.02")), " at 0.02, where ", format(s3at("Markers retained", "0.02"), big.mark = ","),
  " markers remained and the residual variance was ", f4(s3at("Residual variance", "0.02")),
  ". The variance components are from the fit of script S1, which uses its own optimiser; at 0.10 they differ from those of Table 2 by less than 0.02.")

cap_S4 <- paste0(
  "Table S4. Sensitivity to the imputation of missing genotype calls. In the main text, the ",
  format(sum(e01$MISS), big.mark = ","), " missing calls (", formatC(100 * mean(e01$MISS), format = "f", digits = 2),
  "%) are imputed by the marker mean. Here they were imputed with Beagle 5.5 (Browning et al. 2018), without a genetic map or reference panel and with default settings, and the relationship matrix A, the genome scan under the three schemes and the rotated second stage were refitted on the observed trial. Markers above the threshold are those with -log10p above ",
  formatC(e05$bonf, format = "f", digits = 3), ". The variance components are from the rotated second stage without a marker. The one-sided test statistic is the increase in the REML deviance when the residual variance is fixed at one, taken as zero when the estimate is at or below one, and is referred to 2.71 (Materials and Methods). Browning, B. L., Y. Zhou, and S. R. Browning, 2018 A one-penny imputed genome from next-generation reference panels. American Journal of Human Genetics 103: 338–348. doi:10.1016/j.ajhg.2018.07.015")

captions <- c(
  "Supporting Information Tables", "",
  paste0("Table S1. Equivalence and invariance checks at the marker ", TEST_SNP, ", the marker with the smallest p-value in the genome scan. All quantities are computed on the observed trial."),
  "(S1a) The five fixed-effect coefficients under least squares without weighting, least squares on the rotated data, and generalised least squares with weight Omega^-1 on the unrotated data. The genetic kernel is excluded, so that the equivalence of (6) and (7) is exact.",
  "(S1b) The four rotations of Materials and Methods, with the genetic kernel included. The deviation is the largest absolute element of W'W - Omega^-1.",
  "(S1c) The marker effect with the rotated intercept column W1, and with a plain column of ones in its place.",
  "(S1d) A single-stage analysis of the plot data and the three two-stage schemes, all with the genetic kernel included. Ratios are to the single-stage reference.", "",
  "Table S2. The genome scan and the single-stage comparison on the observed trial.",
  paste0("(S2a) Markers passing the Bonferroni threshold (-log10p = ", formatC(e05$bonf, format = "f", digits = 3), ") under any of the three schemes, with -log10p under each scheme."),
  paste0("(S2b) The comparison set of ", nrow(S2b), " markers: for each marker, the estimate, standard error and -log10p from the single-stage analysis of the plot data and from each two-stage scheme. Rows are ordered by the single-stage -log10p."), "",
  cap_S3, "",
  cap_S4, "",
  "Source: 03_equivalence.RData, 05_scan.RData, 06_comparison.RData, S1_het_sweep.RData, and imputation_check/02_check_and_refit.RData and 04_residual_at_one.RData, written by script F2_supporting_tables.R.")

## ---------- write ----------
panels <- list(S1a = S1a, S1b = S1b, S1c = S1c, S1d = S1d, S2a = S2a, S2b = S2b, S3 = S3, S4 = S4)
for (nm in names(panels)) write.csv(panels[[nm]], file.path(TAB, paste0(nm, ".csv")), row.names = FALSE)
cap_sheet <- setNames(data.frame(captions[-1], stringsAsFactors = FALSE), captions[1])
write_xlsx(c(list(Captions = cap_sheet), panels), path = XLSX, col_names = TRUE, format_headers = TRUE)

## ---------- checks against the figures quoted in the text ----------
## the values the paper quotes, checked against the tables
d_lp <- function(nm) S2b[[paste0("-log10p, ", nm)]] - S2b$"-log10p, single-stage"
top  <- e05$scans$rotated$res
checks <- data.frame(
  Check = c("S1: the test marker has the smallest p-value of the rotated scan",
            "S1a: rotated OLS equals GLS on all five coefficients (to 1e-5)",
            "S1a: marker effect 6.93926 rotated, 6.44069 unweighted",
            "S1c: plain intercept 12.45060",
            "S1d: single-stage 6.2263 (1.0793); SE ratios 1.0024 / 1.0001 / 0.9848",
            "S2a: seven markers, all on chromosome 19, window 48,498,036 to 48,929,148",
            "S2a: 48,815,157 is 0.0142 above the threshold unweighted, 0.1619 below rotated",
            "S2b: 227 markers, 7 flagged significant",
            "S2b: 45 / 42 / 42% of -log10p differences positive (unweighted / rotated / inverse-variance); medians -0.0095 / -0.0210 / -0.0185",
            "S3: threshold 0.10 retains 59,206 markers; rank 158 throughout; lambda 0.9630-0.9743 and residual 0.5209-0.5762 down to 0.10",
            "S3: mean diagonal increases at every step, 1.7274 to 1.8999; 30,302 markers and residual 1.0363 at 0.02",
            "S3: variance components at 0.10 within 0.02 of Table 2",
            "S4: marker mean 6 / 6 / 7 above the threshold, as in S2a, window as in S2a; Beagle 9 / 9 / 11 on chromosome 19 within the same peak",
            "S4: rotated fit with the marker mean equals Table 2; Beagle residual 1.4516, statistic 0.8474 below 2.71; mean diagonal 1.8417 / 1.9249"),
  Result = c(top$SNP[which.max(top$neglog10p)] == TEST_SNP,
             all(abs(S1a[, 3] - S1a[, 4]) < 1e-5),
             S1a[5, 3] == 6.93926 && S1a[5, 2] == 6.44069,
             S1c[3, 2] == 12.4506,
             S1d[1, 2] == 6.2263 && S1d[1, 3] == 1.0793 && all(S1d$"SE ratio" == c(1, 0.9848, 1.0001, 1.0024)),
             nrow(S2a) == 7 && all(S2a$Chromosome == 19) && min(S2a$Position) == 48498036 && max(S2a$Position) == 48929148,
             abs((S2a$"-log10p, unweighted"[S2a$Position == 48815157] - e05$bonf) - 0.0142) < 1e-4 &&
               abs((e05$bonf - S2a$"-log10p, rotated"[S2a$Position == 48815157]) - 0.1619) < 1e-4,
             nrow(S2b) == 227 && sum(S2b$"Significant under any scheme" == "yes") == 7,
             round(100 * mean(d_lp("unweighted") > 0)) == 45 && round(100 * mean(d_lp("rotated") > 0)) == 42 &&
               round(100 * mean(d_lp("inverse-variance") > 0)) == 42 &&
               round(median(d_lp("unweighted")), 4) == -0.0095 && round(median(d_lp("rotated")), 4) == -0.0210 &&
               round(median(d_lp("inverse-variance")), 4) == -0.0185,
             s3at("Markers retained", "0.10") == 59206 && all(S3$"Rank of A" == 158) &&
               min(S3$"Genomic inflation factor"[r10]) == 0.9630 && max(S3$"Genomic inflation factor"[r10]) == 0.9743 &&
               s3at("Genomic inflation factor", "0.10") == round(e05$scan_tab$lambda_P3D[e05$scan_tab$Scheme == "rotated"], 4) &&
               min(S3$"Residual variance"[r10]) == 0.5209 && max(S3$"Residual variance"[r10]) == 0.5762,
             all(diff(S3$"Mean diagonal of A") > 0) && s3at("Mean diagonal of A", "none") == 1.7274 &&
               s3at("Mean diagonal of A", "0.02") == 1.8999 && s3at("Markers retained", "0.02") == 30302 &&
               s3at("Residual variance", "0.02") == 1.0363,
             abs(s3at("Genetic variance", "0.10") - e05$scans$rotated$Vg) < 0.02 &&
               abs(s3at("Residual variance", "0.10") - e05$scans$rotated$Ve) < 0.02,
             all(S4[1:3, 2] == c("6", "6", "7")) &&
               all(as.integer(S4[1:3, 2]) == sapply(c("rotated", "inverse-variance", "unweighted"), function(nm) sum(S2a[[paste0("-log10p, ", nm)]] > e05$bonf, na.rm = TRUE))) &&
               S4[4, 2] == paste0("19: ", format(min(S2a$Position), big.mark = ","), "–", format(max(S2a$Position), big.mark = ",")) &&
               all(S4[1:3, 3] == c("9", "9", "11")) && all(isg$Chromosome == 19) &&
               min(isg$Position) >= 48.49e6 && max(isg$Position) <= 48.94e6,
             abs(rot("in use", "Vg") - e05$scans$rotated$Vg) < 1e-3 && abs(rot("in use", "Ve") - e05$scans$rotated$Ve) < 1e-4 &&
               S4[7, 3] == "1.4516" && S4[8, 3] == "0.8474" && S4[8, 2] == "0.0000" && as.numeric(S4[8, 3]) < 2.71 &&
               S4[5, 2] == "1.8417" && S4[5, 3] == "1.9249"))
print(checks, row.names = FALSE)
if (!all(checks$Result)) stop("supporting-table checks failed")
cat("\nwritten", XLSX, "|", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
