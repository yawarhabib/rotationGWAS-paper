################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##        Script F3: the Supporting Information figures
##
##   Assembles Figures S1 to S3, drawn by script F1, into one Word document with
##   their captions, so that the figure files and the captions are produced
##   together and every number a caption quotes is read from the saved results.
##
##     S1  cumulative distribution of the paired differences in -log10 p (06)
##     S2  the genome scan on chromosome 19 around the detected locus    (05)
##     S3  genome-wide standard errors relative to the rotated scheme    (05)
##
##   Output: results/Supporting_Information_Figures.docx
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(officer)
t0 <- Sys.time()

FIG    <- file.path(ROOT, "results/figures")
OUTDOC <- file.path(ROOT, "results/Supporting_Information_Figures.docx")

e05 <- new.env(); load(file.path(OUT, "05_scan.RData"),       envir = e05)
e06 <- new.env(); load(file.path(OUT, "06_comparison.RData"), envir = e06)
ARMS <- c("rotated", "diagonal", "unweighted")

## ---------- the numbers the captions quote ----------
ss <- e06$ss
pd <- lapply(ARMS, function(a) { t <- e06$ts[[a]]; t$neglog10p[match(ss$SNP, t$SNP)] - ss$neglog10p })
n_set    <- nrow(ss)                                              # markers behind S1, all of them drawn
sig_pos  <- e05$sig_tab$Position
window   <- c(min(sig_pos) - 3e5, max(sig_pos) + 3e5) / 1e6      # the S2 window, as drawn by F1
n_scan   <- length(e05$se_ratio$unweighted)                       # markers behind S3
fmt      <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

## ---------- text styles ----------
body  <- fp_text(font.family = "Times New Roman", font.size = 11)
bold  <- update(body, bold = TRUE)
ital  <- update(body, italic = TRUE)
sub   <- update(body, vertical.align = "subscript")
title <- update(body, font.size = 14, bold = TRUE)
lp    <- function() list(ftext("\u2212log", body), ftext("10", sub), ftext("p", ital))   # -log10 p
cap   <- function(label, ...) {
  pieces <- list(ftext(label, bold), ftext(" ", body))
  for (x in list(...)) pieces <- c(pieces, if (is.character(x)) list(ftext(x, body)) else x)
  do.call(fpar, pieces)
}

## ---------- captions ----------
cap_S1 <- cap("Figure S1",
  "Cumulative distribution of the paired differences in ", lp(),
  sprintf(" between each two-stage scheme and the single-stage analysis over the %d markers of Figure 1. ", n_set),
  "The share of a scheme's differences above zero is read where its curve crosses the vertical line at zero, ",
  "and its median where it crosses the horizontal line at 0.5, marked by the symbol of that scheme.")
cap_S2 <- cap("Figure S2",
  sprintf("The genome scan on chromosome 19 around the locus for vegetation period, from %.1f to %.1f Mb, under the three weighting schemes. ", window[1], window[2]),
  "The dashed line is the Bonferroni threshold at ", lp(), sprintf(" = %.3f. ", e05$bonf),
  sprintf("Every marker passing the threshold under any scheme lies between positions %s and %s.", fmt(min(sig_pos)), fmt(max(sig_pos))))
cap_S3 <- cap("Figure S3",
  "Density of the standard error of the marker effect under inverse-variance weighting and under the unweighted analysis, ",
  "each relative to the standard error under the rotation, over the ", fmt(n_scan), " markers of the genome scan. ",
  "The vertical line marks a ratio of one.")

## ---------- the document ----------
figs <- list(list("figureS1_paired_differences.png", 5.0, 3.8, cap_S1),
             list("figureS2_locus.png",           6.4, 3.8, cap_S2),
             list("figureS3_se_ratio.png",        6.4, 3.8, cap_S3))
doc <- read_docx()
doc <- body_add_fpar(doc, fpar(ftext("Supporting Information", title)))
doc <- body_add_fpar(doc, fpar(ftext("Exact weighting in two-stage genome-wide association analysis: when does it matter?", ital)))
doc <- body_add_fpar(doc, fpar(ftext("Y. Habib and H. P. Piepho", body)))
doc <- body_add_fpar(doc, fpar(ftext("Tables S1 to S3 are provided in the workbook Supporting_Information_Tables.xlsx.", body)))
for (f in figs) {
  doc <- body_add_par(doc, "")
  doc <- body_add_img(doc, src = file.path(FIG, f[[1]]), width = f[[2]], height = f[[3]])
  doc <- body_add_fpar(doc, f[[4]])
}
print(doc, target = OUTDOC)

## ---------- checks ----------
checks <- data.frame(
  Check  = c("every figure file exists", "S1 covers every marker of the comparison set of 06, with no missing difference",
             "S2 window covers the significant markers", "S3 marker count equals the scan",
             "document written"),
  Result = c(all(file.exists(file.path(FIG, sapply(figs, `[[`, 1)))),
             n_set == length(e06$sel) && all(lengths(pd) == n_set) && !anyNA(unlist(pd)),
             window[1] < min(sig_pos) / 1e6 && window[2] > max(sig_pos) / 1e6,
             n_scan == sum(!is.na(e05$scans$rotated$res$SE) & !is.na(e05$scans$unweighted$res$SE) & !is.na(e05$scans$diagonal$res$SE)),
             file.exists(OUTDOC)))
print(checks, row.names = FALSE)
cat(sprintf("\nS1 markers %d | S2 window %.3f-%.3f Mb | S3 markers %s\n", n_set, window[1], window[2], fmt(n_scan)))
if (!all(checks$Result)) stop("supporting-figure checks failed")
cat("written", OUTDOC, "|", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
