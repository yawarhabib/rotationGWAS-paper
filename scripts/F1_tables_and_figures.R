################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script F1: the tables and figures of the paper
##
##   Builds every table and figure of the paper from the saved output of the
##   analysis scripts. Computes nothing new: each item is a selection or a
##   rearrangement of results already checked.
##
##   Tables are written as CSV and as a markdown fragment. Figures are drawn
##   with ggplot2 under one theme and written as PDF and as PNG from the same
##   plot object.
##
##   Output: results/tables/*.csv and tables.md, results/figures/*.pdf and *.png
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
t0 <- Sys.time()

BUILD <- file.path(ROOT, "results")
TAB   <- file.path(BUILD, "tables")
FIG   <- file.path(BUILD, "figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

e01 <- new.env(); load(file.path(OUT, "01_data.RData"),      envir = e01)
e05 <- new.env(); load(file.path(OUT, "05_scan.RData"),      envir = e05)
e06 <- new.env(); load(file.path(OUT, "06_comparison.RData"),envir = e06)
e10 <- new.env(); load(file.path(OUT, "10_simC.RData"),      envir = e10)
e11 <- new.env(); load(file.path(OUT, "11_simD.RData"),      envir = e11)
e11b <- new.env(); load(file.path(OUT, "11b_simD_one_sided.RData"), envir = e11b)
eS2 <- new.env(); load(file.path(OUT, "S2_efficiency.RData"), envir = eS2)

## colours distinguishable in greyscale as well as in colour
COL <- c(rotated = "#0072B2", diagonal = "#D55E00", unweighted = "#000000")   # Okabe-Ito
PCH <- c(rotated = 19,        diagonal = 17,        unweighted = 15)
LTY <- c(rotated = 1,         diagonal = 2,         unweighted = 3)
ARMS <- c("rotated", "diagonal", "unweighted")
LAB  <- c(rotated = "rotated", diagonal = "inverse-variance", unweighted = "unweighted")

## ---------- Table 1: the trial ----------
t1 <- e01$design_tab
names(t1) <- c("Year", "Accessions", "Plots", "Blocks analysed")

## ---------- Table 2: the three schemes on the observed trial ----------
sig_n <- vapply(ARMS, function(a) {
  r <- e05$scans[[a]]$res
  sum(!is.na(r$neglog10p) & r$neglog10p > e05$bonf)
}, numeric(1))

m <- e06$metrics
pick <- function(a, q, col) m[[col]][m$Scheme == a & m$Quantity == q]

t2 <- data.frame(
  Scheme            = unname(LAB[ARMS]),
  Vg                = vapply(ARMS, function(a) round(e05$scans[[a]]$Vg, 3), numeric(1)),
  Ve                = vapply(ARMS, function(a) round(e05$scans[[a]]$Ve, 4), numeric(1)),
  Markers_above     = unname(sig_n),
  Bias_neglog10p    = vapply(ARMS, function(a) round(pick(a, "neglog10p", "Bias"),   4), numeric(1)),
  t_neglog10p       = vapply(ARMS, function(a) round(pick(a, "neglog10p", "t_bias"), 2), numeric(1)),
  Bias_se           = vapply(ARMS, function(a) round(pick(a, "se", "Bias"),   4), numeric(1)),
  t_se              = vapply(ARMS, function(a) round(pick(a, "se", "t_bias"), 2), numeric(1)),
  row.names = NULL)

## column heads as the paper prints them
names(t2) <- c("Scheme", "sigma2_g", "sigma2_e", "Markers above threshold",
               "Mean difference in -log10p", "t for -log10p",
               "Mean difference in SE", "t for SE")

## ---------- Table 3: the regime map ----------
## both halves sorted by the number of accessions per environment
co <- e10$c_omega[order(-e10$c_omega$accessions_per_year), ]
## The estimation share on the average semivariance (ASV) comes from script
## 10k, matched by the number of years.
e10k <- new.env(); load(file.path(OUT, "10k_simC_asv_share.RData"), envir = e10k)
ik <- match(co$n_years, e10k$k_share$n_years); stopifnot(!anyNA(ik), all(e10k$k_share$n[ik] == co$n))
co$estimation_share <- e10k$k_share$share_asv[ik]
ca <- e10$c_causal
gp <- function(a, col) ca[[col]][ca$arm == a][order(ca$n_years[ca$arm == a])]

## power in percent, to one decimal
pct <- function(x) formatC(100 * x, format = "f", digits = 1)
t3 <- data.frame(
  Accessions_per_environment = co$accessions_per_year,
  Mean_correlation  = round(co$mean_corr, 4),
  Estimation_share  = round(co$estimation_share, 4),
  Power_rotated     = pct(gp("rotated",    "power_59206")),
  Power_inverse_var = pct(gp("diagonal",   "power_59206")),
  Power_unweighted  = pct(gp("unweighted", "power_59206")),
  RMSE_rotated      = round(gp("rotated",    "RMSE"), 3),
  RMSE_inverse_var  = round(gp("diagonal",   "RMSE"), 3),
  RMSE_unweighted   = round(gp("unweighted", "RMSE"), 3),
  row.names = NULL, stringsAsFactors = FALSE)
names(t3) <- c("Accessions per environment", "Mean absolute correlation", "Estimation share",
               "Power rotated", "Power inverse-variance", "Power unweighted",
               "RMSE rotated", "RMSE inverse-variance", "RMSE unweighted")

## ---------- Table 4: what the specification test can detect ----------
## A selection from scenario D: power_cal is the detection rate calibrated on
## the scenario's own control (script 11); the one-sided test's column comes
## from 11b, matched to the same condition. The two structures are labelled by
## their correlation with A, which script 11 saves.
st <- e11$sens_tab
stopifnot(all(st$Withheld %in% c("none", "epi", "het")))
lab_cor <- function(k) sprintf("correlation %.2f with A", e11$KERNEL_COR[[k]])
os <- e11b$one_sided_tab
os$key <- ifelse(grepl("^epi", os$structure), "epi", ifelse(grepl("^het", os$structure), "het", "none"))
io <- match(paste(st$Withheld, st$Variance), paste(os$key, os$withheld))
stopifnot(!anyNA(io), all(os$n[io] == st$n), isTRUE(all.equal(os$threshold_pct[io], st$power_cal)))
t4 <- data.frame(
  Withheld            = ifelse(st$Withheld == "none", "none",
                        ifelse(st$Withheld == "epi", lab_cor("epi"), lab_cor("het"))),
  Variance            = st$Variance,
  Residual_mean       = st$Ve_mean,
  Residual_sd         = st$Ve_sd,
  One_sided           = sprintf("%.1f (%.1f)", os$one_sided_pct[io], os$one_sided_mcse[io]),
  Calibrated          = st$power_cal,
  row.names = NULL, stringsAsFactors = FALSE)
names(t4) <- c("Structure", "Variance withheld", "Residual mean", "Residual SD",
               "One-sided test (%)", "Calibrated threshold (%)")

for (nm in c("t1", "t2", "t3", "t4"))
  write.csv(get(nm), file.path(TAB, paste0(nm, ".csv")), row.names = FALSE)

## ---------- the tables as a markdown fragment ----------
md_table <- function(d, caption) {
  hdr <- paste0("| ", paste(gsub("_", " ", names(d)), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  rows <- apply(d, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |"))
  c(caption, "", hdr, sep, rows, "")
}
cap1 <- paste("Table 1  The trial. Accessions entering and leaving between seasons.",
              "Plots number fewer than accessions times blocks because not every",
              "accession has a recorded value in every block.")
## on the observed trial the single-stage fit is a comparator, so these are
## mean differences, not biases
cap2 <- paste0("Table 2  The three weighting schemes on the observed trial. Differences and t ",
               "values are paired comparisons with a single-stage analysis over ",
               length(e06$sel), " markers.")

cap3 <- "Table 3  The regime map. Power, in percent, is measured at the genome-scan threshold."

cap4 <- paste("Table 4  Scenario D. The one-sided test is the proportion of replicates in which the increase in",
              "deviance from fixing the residual variance at one exceeds 2.71, the increase being taken as zero",
              "when the estimate is at or below one (Monte Carlo standard error in parentheses); the calibrated",
              "threshold is the proportion in which the residual variance exceeds the 95th percentile of the control.")

## written as one fragment and as two, one for the Methods and one for the
## Results
writeLines(c(md_table(t1, cap1), md_table(t2, cap2),
             md_table(t3, cap3), md_table(t4, cap4)),
           file.path(TAB, "tables.md"))
writeLines(md_table(t1, cap1), file.path(TAB, "tables_methods.md"))
writeLines(c(md_table(t2, cap2), md_table(t3, cap3), md_table(t4, cap4)),
           file.path(TAB, "tables_results.md"))

## ---------- what Figure 1 shows, as numbers ----------
## the share of differences above zero, the median, and the mean without each
## scheme's largest difference
pd <- lapply(ARMS, function(a) e06$ts[[a]]$neglog10p - e06$ss$neglog10p)
names(pd) <- ARMS

figS1_tab <- data.frame(
  Scheme           = unname(LAB[ARMS]),
  Mean             = round(vapply(pd, mean, numeric(1)), 4),
  Median           = round(vapply(pd, median, numeric(1)), 4),
  Percent_above_0  = round(100 * vapply(pd, function(d) mean(d > 0), numeric(1)), 0),
  Largest          = round(vapply(pd, max, numeric(1)), 4),
  Mean_less_largest = round(vapply(pd, function(d) mean(d[-which.max(d)]),
                                   numeric(1)), 4),
  row.names = NULL)

## the paired difference against the strength of the marker: a least-squares
## slope per scheme, and the means below and above a single-stage -log10 p of 4
x_ss <- e06$ss$neglog10p
fig1_dep <- do.call(rbind, lapply(ARMS, function(a) {
  f <- lm(pd[[a]] ~ x_ss); strong <- x_ss > 4
  data.frame(Scheme = unname(LAB[a]),
             Slope = round(coef(f)[2], 4), P_slope = round(summary(f)$coefficients[2, 4], 3),
             Mean_weak = round(mean(pd[[a]][!strong]), 4), N_weak = sum(!strong),
             Mean_strong = round(mean(pd[[a]][strong]), 4), N_strong = sum(strong),
             row.names = NULL) }))

## ---------- figures ----------
## All six figures are drawn with ggplot2 under one theme and written as PDF
## and as PNG from the same plot object.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
set.seed(20260916)                       # the jitter in Figure 1B is reproducible
GCOL <- setNames(COL[ARMS], LAB[ARMS]); GSHP <- setNames(PCH[ARMS], LAB[ARMS])
GLTY <- setNames(c("solid", "22", "12"), LAB[ARMS])
th <- theme_classic(base_size = 9) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.text = element_text(colour = "black"), legend.title = element_blank(),
        legend.background = element_blank(), legend.key.width = unit(1.4, "lines"),
        legend.text = element_text(size = 8), plot.tag = element_text(face = "bold", size = 11),
        strip.background = element_blank(), strip.text = element_text(size = 9),
        plot.margin = margin(4, 6, 2, 4))
scheme_scales <- function(keys = names(GCOL), aes = c("colour", "shape", "linetype")) {
  ## only the aesthetics a figure maps, so that ggplot2 does not warn about
  ## the others
  sc <- list(colour   = scale_colour_manual(values = GCOL[keys], breaks = keys),
             shape    = scale_shape_manual(values = GSHP[keys], breaks = keys),
             linetype = scale_linetype_manual(values = GLTY[keys], breaks = keys))
  unname(sc[aes])
}
legend_at <- function(x, y) theme(legend.position = "inside", legend.position.inside = c(x, y),
                                  legend.justification = c(x, y))
scheme_df <- function(a, ...) { d <- data.frame(...); d$scheme <- factor(LAB[a], levels = names(GCOL)); d }
save_fig <- function(name, p, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
LP_X <- expression("single-stage " * -log[10] * italic(p))
LP_D <- expression("paired difference in " * -log[10] * italic(p))

## ---------- Figure 1: agreement with the single-stage analysis ----------
## (A) the paired difference in -log10 p against the single-stage -log10 p, one
## point per marker and scheme, with a least-squares line per scheme. Every
## marker is shown; the axis covers the whole range. (B) the paired difference
## in the standard error, with the mean and two standard errors on either side.
fig1 <- function() {
  x <- e06$ss$neglog10p
  da <- do.call(rbind, lapply(ARMS, function(a) scheme_df(a, x = x, y = pd[[a]])))
  pA <- ggplot(da, aes(x, y, colour = scheme, shape = scheme)) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
    geom_vline(xintercept = e05$bonf, linetype = "22", colour = "grey40") +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
    scheme_scales(aes = c("colour", "shape")) + coord_cartesian(xlim = range(x)) +
    labs(x = LP_X, y = LP_D, tag = "A") + th +
    ## above the panel: with every marker shown, no corner of the panel is
    ## free of points
    theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))
  ds <- do.call(rbind, lapply(ARMS, function(a) scheme_df(a, d = e06$ts[[a]]$se - e06$ss$se)))
  sm <- do.call(rbind, lapply(split(ds, ds$scheme), function(g) {
    se <- sd(g$d) / sqrt(nrow(g))
    data.frame(scheme = g$scheme[1], m = mean(g$d), lo = mean(g$d) - 2 * se, hi = mean(g$d) + 2 * se) }))
  pB <- ggplot(ds, aes(scheme, d, colour = scheme, shape = scheme)) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
    geom_jitter(width = 0.14, height = 0, alpha = 0.45, size = 1.5) +
    geom_errorbar(data = sm, aes(x = scheme, ymin = lo, ymax = hi, colour = scheme), inherit.aes = FALSE,
                  width = 0.12, linewidth = 0.6) +
    geom_errorbar(data = sm, aes(x = scheme, ymin = m, ymax = m, colour = scheme), inherit.aes = FALSE,
                  width = 0.6, linewidth = 1.1) +
    scheme_scales(aes = c("colour", "shape")) + labs(x = NULL, y = "paired difference in standard error", tag = "B") +
    th + theme(legend.position = "none")
  pA + pB + plot_layout(ncol = 2)
}

## ---------- Figure 2: the regime map ----------
## Three panels in one row on the design axis: the quantity that drives the
## result, then the two results. Accessions per environment on a log scale,
## increasing to the right, with the number of years of the panel on top,
## since the two are the same grid read two ways.
fig2 <- function() {
  xb <- co$accessions_per_year; yrs <- co$n_years
  sx <- function() scale_x_log10(breaks = xb, labels = formatC(xb, format = "f", digits = 1),
                                 sec.axis = sec_axis(~ ., breaks = xb, labels = yrs, name = "years"),
                                 expand = expansion(mult = 0.06))
  obs_share <- eS2$share_tab$estimation_share_diff[1]   # the observed trial, from S2, on the ASV
  pA <- ggplot(co, aes(accessions_per_year, estimation_share)) +
    geom_hline(yintercept = obs_share, linetype = "dotted", colour = "grey45") +
    annotate("text", x = max(xb), y = obs_share, label = "observed trial", hjust = 1, vjust = -0.6,
             size = 2.6, colour = "grey30") +
    geom_line(colour = "grey15", linewidth = 0.6) +
    geom_point(shape = 21, fill = "white", colour = "grey15", size = 2, stroke = 0.7) +
    sx() + scale_y_continuous(limits = c(0, 0.36), expand = c(0, 0)) +
    labs(x = NULL, y = "estimation share", tag = "A") + th
  ca <- do.call(rbind, lapply(ARMS, function(a) scheme_df(a, x = xb, y = gp(a, "power_59206"))))
  pB <- ggplot(ca, aes(x, y, colour = scheme, shape = scheme, linetype = scheme)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) + scheme_scales() + sx() +
    scale_y_continuous(limits = c(0.4, 0.95)) +
    labs(x = "accessions per environment", y = "power", tag = "B") + th + legend_at(0.98, 0.04)
  ## the standard deviation of the estimates under the null relative to the
  ## rotation: the scheme effect alone, since the design effect cancels
  cr <- e10$c_ratio[match(co$n_years, e10$c_ratio$n_years), ]
  rat <- rbind(scheme_df("rotated", x = xb, y = 1),
               scheme_df("diagonal", x = xb, y = cr$diag_over_rot),
               scheme_df("unweighted", x = xb, y = cr$unw_over_rot))
  pC <- ggplot(rat, aes(x, y, colour = scheme, shape = scheme, linetype = scheme)) +
    geom_line(linewidth = 0.6) + geom_point(data = rat[rat$scheme != "rotated", ], size = 2) +
    scheme_scales() + sx() + scale_y_continuous(limits = c(0.98, 1.20)) +
    labs(x = NULL, y = "SD of estimates relative to rotation", tag = "C") + th +
    theme(legend.position = "none")
  ## narrower margins between the three panels
  tight <- theme(plot.margin = margin(4, 2, 2, 2))
  (pA + tight) + (pB + tight) + (pC + tight) + plot_layout(ncol = 3)
}

## ---------- Figure 3: the residual as a specification test ----------
## Excess residual over the control row of the same fitted model, so each model
## is read against its own baseline. Two panels for the two withheld
## structures, sharing a scale so the difference in detectability is visible.
fig3 <- function() {
  ex <- e11$exc_tab
  srcs <- c(epi = lab_cor("epi"), het = lab_cor("het"))
  kinds <- c(additive = "additive only", matched = "matched kernel", wrong = "wrong kernel", generic = "WW'")
  ## the legend prints WW with the transpose as a superscript, as in the caption
  KLAB <- setNames(expression("additive only", "matched kernel", "wrong kernel", bold(WW)^T), kinds)
  long <- do.call(rbind, lapply(names(srcs), function(s) {
    d <- rbind(ex[ex$Withheld == "none", ], ex[ex$Withheld == s, ])
    do.call(rbind, lapply(names(kinds), function(k)
      data.frame(panel = unname(srcs[s]), x = d$Variance, y = d[[paste0("excess_", k)]], model = unname(kinds[k]))))
  }))
  long$panel <- factor(long$panel, levels = srcs); long$model <- factor(long$model, levels = kinds)
  KCOL <- setNames(c("#D55E00", "#0072B2", "#CC79A7", "#009E73"), kinds)   # Okabe-Ito
  KSHP <- setNames(c(16, 17, 15, 18), kinds)
  ggplot(long, aes(x, y, colour = model, shape = model)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) +
    facet_wrap(~ panel) +
    scale_colour_manual(values = KCOL, labels = KLAB) + scale_shape_manual(values = KSHP, labels = KLAB) +
    labs(x = "variance withheld", y = "excess residual") + th + legend_at(0.02, 0.98)
}

## ---------- Figure S1: the distribution of the paired differences ----------
## The cumulative distribution of the paired differences in -log10 p. The share
## of a scheme's differences above zero is read where its curve crosses zero,
## the median where it crosses 0.5. This is the distribution behind figS1_tab.
## Every marker is shown; the axis covers the whole range.
figS1 <- function() {
  da <- do.call(rbind, lapply(ARMS, function(a) scheme_df(a, d = pd[[a]])))
  med <- do.call(rbind, lapply(ARMS, function(a) scheme_df(a, m = median(pd[[a]]))))
  ggplot(da, aes(d, colour = scheme)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, colour = "grey75", linetype = "22") +
    stat_ecdf(geom = "step", linewidth = 0.7, pad = TRUE) +
    geom_point(data = med, aes(m, 0.5, colour = scheme, shape = scheme), size = 2.4) +
    scheme_scales(aes = c("colour", "shape")) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = LP_D, y = "proportion of markers") + th + legend_at(0.98, 0.04)
}

## ---------- Figure S2: the detected region ----------
## Every marker on chromosome 19 within the window that the significant markers
## span, widened by 0.3 Mb on each side, under the three schemes, with the
## Bonferroni threshold.
figS2 <- function() {
  map <- e01$map; sig <- e05$sig_tab
  lo <- min(sig$Position) - 3e5; hi <- max(sig$Position) + 3e5
  w  <- which(map$Chromosome == 19 & map$Position >= lo & map$Position <= hi)
  snp <- map$SNP[w]; pos <- map$Position[w] / 1e6
  dd <- do.call(rbind, lapply(ARMS, function(a) {
    r <- e05$scans[[a]]$res; scheme_df(a, pos = pos, y = r$neglog10p[match(snp, r$SNP)]) }))
  ggplot(dd, aes(pos, y, colour = scheme, shape = scheme)) +
    geom_hline(yintercept = e05$bonf, linetype = "22", colour = "grey30") +
    geom_point(alpha = 0.75, size = 1.8) + scheme_scales(aes = c("colour", "shape")) +
    labs(x = "position on chromosome 19 (Mb)", y = expression(-log[10] * italic(p))) +
    th + legend_at(0.98, 0.72)   # the empty band right of the locus, clear of the threshold line
}

## ---------- Figure S3: standard errors across the genome ----------
## The ratio of each scheme's standard error to the rotated one over all
## markers, as a density.
figS3 <- function() {
  sr <- do.call(rbind, lapply(names(e05$se_ratio), function(a) scheme_df(a, r = e05$se_ratio[[a]])))
  ggplot(sr, aes(r, colour = scheme, linetype = scheme)) +
    geom_vline(xintercept = 1, colour = "grey35", linewidth = 0.4) +
    geom_density(adjust = 1.2, linewidth = 0.7, key_glyph = "path") +
    scheme_scales(unname(LAB[names(e05$se_ratio)]), aes = c("colour", "linetype")) +
    labs(x = "standard error relative to the rotated scheme", y = "density") + th + legend_at(0.02, 0.98)
}

save_fig("figure1_single_stage",       fig1(),  7.2, 3.6)
save_fig("figure2_regime_map",         fig2(),  7.5, 2.85)
save_fig("figure3_specification_test", fig3(),  7.0, 3.6)
save_fig("figureS1_paired_differences", figS1(), 5.0, 3.8)
save_fig("figureS2_locus",             figS2(), 6.4, 3.8)
save_fig("figureS3_se_ratio",          figS3(), 6.4, 3.8)

## ---------- checks ----------
files <- c(file.path(TAB, c("t1.csv", "t2.csv", "t3.csv", "t4.csv", "tables.md",
                            "tables_methods.md", "tables_results.md")),
           file.path(FIG, paste0(rep(c("figure1_single_stage", "figure2_regime_map",
                                       "figure3_specification_test", "figureS1_paired_differences",
                                       "figureS2_locus", "figureS3_se_ratio"), each = 2),
                                 c(".pdf", ".png"))))

checks <- data.frame(
  Check = c("every display file written",
            "no file is empty",
            "Table 1 has one row per year",
            "Table 2 has one row per scheme",
            "Table 3 has one row per grid level",
            "Table 4 has one row per truth condition",
            "no missing values in any table",
            "the significant-marker counts match the scan",
            "Figure 1A: level dependence computed on every marker of the comparison set"),
  Result = c(all(file.exists(files)),
             all(file.info(files)$size > 0),
             nrow(t1) == nlevels(e01$pheno$Year),
             nrow(t2) == length(ARMS),
             nrow(t3) == nrow(co),
             nrow(t4) == nrow(e11$sens_tab),
             !anyNA(t1) && !anyNA(t2) && !anyNA(t3) && !anyNA(t4),
             identical(as.numeric(t2[["Markers above threshold"]]), as.numeric(sig_n)),
             all(fig1_dep$N_weak + fig1_dep$N_strong == length(e06$sel))))

cat("\n=== F1 tables and figures ===\n")
cat("\nTable 1, the trial\n");            print(t1, row.names = FALSE)
cat("\nTable 2, the observed trial\n");   print(t2, row.names = FALSE)
cat("\nTable 3, the regime map\n");       print(t3, row.names = FALSE)
cat("\nTable 4, what the test can detect\n"); print(t4, row.names = FALSE)
cat("\nFigure 1A, the paired difference against the strength of the marker\n")
print(fig1_dep, row.names = FALSE)
cat("\nthe distribution behind Figure S1\n")
print(figS1_tab, row.names = FALSE)
cat("   a mean can be moved by one point; the share above zero and the median",
    "say where the body sits\n")
cat("\nfiles written:", length(files), "\n")
cat("   tables  ->", TAB, "\n")
cat("   figures ->", FIG, "\n")

cat("\nchecks\n"); print(checks, row.names = FALSE)

save(t1, t2, t3, t4, fig1_dep, figS1_tab, checks, files, file = file.path(OUT, "F1_display.RData"))

cat("\nsaved F1_display.RData |",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

if (!all(checks$Result)) stop("display-item checks failed")
