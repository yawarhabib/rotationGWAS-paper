################################################################################
##   Imputation check, step 2: is the imputation sound, and what does it change?
##
##   A. Observed calls must come back unchanged; only the gaps are filled.
##   B. Accuracy, measured on the 2% of observed calls that were hidden before a
##      second imputation, against what mean imputation gives on the same cells.
##   C. The relationship matrix rebuilt from the completed genotypes, set beside
##      the one in use.
##   D. The observed-trial fits under both: the rotated second stage, the other
##      two schemes, the single-stage fit, and the genome scan.
##
##   Reads the stored objects of the analysis; writes only into this folder.
################################################################################

source("scripts/00_config.R")
suppressMessages({library(rrBLUP); library(lme4); library(lme4qtl)})
HERE <- file.path(ROOT, "imputation_check")
load(file.path(HERE, "01_vcf_inputs.RData"))             # G (160 x 59206, NA = missing), hide, bim
G160 <- G; rm(G)
load(file.path(OUT, "01_data.RData"))                    # pheno, G (159, mean-imputed), A, Q, MISS, map, n
e2 <- new.env(); load(file.path(OUT, "02_rotation.RData"), envir = e2)
Omega <- e2$Omega; blues <- e2$blues
options(width = 220)
sym <- function(M) (M + t(M)) / 2

read_beagle <- function(f) {
  x  <- readLines(gzfile(f)); x <- x[!startsWith(x, "##")]
  hd <- strsplit(sub("^#", "", x[1]), "\t")[[1]]
  b  <- do.call(rbind, strsplit(x[-1], "\t", fixed = TRUE))
  gt <- b[, 10:ncol(b)]
  D  <- (substr(gt, 1, 1) == "1") + (substr(gt, 3, 3) == "1")
  dimnames(D) <- list(b[, 3], hd[10:length(hd)])
  t(D)                                                    # individuals x markers
}

## ---------- A. observed calls unchanged ----------
Bf <- read_beagle(file.path(HERE, "imputed_full.vcf.gz"))
stopifnot(identical(dimnames(Bf), list(rownames(G160), bim$SNP)), !anyNA(Bf))
obs <- !is.na(G160)
cat("\nA. the full imputation\n")
cat(sprintf("   observed calls changed by Beagle: %d of %d\n", sum(Bf[obs] != G160[obs]), sum(obs)))
cat(sprintf("   gaps filled: %d | filled as 0 / 1 / 2: %s\n", sum(!obs),
            paste(round(100 * prop.table(table(factor(Bf[!obs], 0:2))), 1), collapse = " / ")))
cat(sprintf("   observed calls are 0 / 1 / 2: %s\n", paste(round(100 * prop.table(table(factor(G160[obs], 0:2))), 1), collapse = " / ")))

## ---------- B. accuracy on hidden calls ----------
Bm <- read_beagle(file.path(HERE, "imputed_masked.vcf.gz"))
truth <- G160[hide]; got <- Bm[hide]
Gm <- G160; Gm[hide] <- NA
mm <- colMeans(Gm, na.rm = TRUE)
mean_fill <- mm[col(G160)[hide]]
acc <- data.frame(method = c("Beagle", "marker mean"),
                  exact_genotype = c(mean(got == truth), mean(round(mean_fill) == truth)),
                  dosage_r = c(cor(got, truth), cor(mean_fill, truth)),
                  dosage_rmse = c(sqrt(mean((got - truth)^2)), sqrt(mean((mean_fill - truth)^2))))
cat("\nB. accuracy on", length(hide), "hidden observed calls\n"); print(format(acc, digits = 4), row.names = FALSE)
byclass <- sapply(0:2, function(k) mean(got[truth == k] == k))
cat(sprintf("   Beagle, exact by true genotype 0 / 1 / 2: %s (n = %s)\n", paste(round(byclass, 3), collapse = " / "),
            paste(table(factor(truth, 0:2)), collapse = " / ")))
maf_all <- pmin(mm / 2, 1 - mm / 2)[col(G160)[hide]]
cat(sprintf("   Beagle, exact by minor allele frequency below 0.1 / 0.1 to 0.25 / above 0.25: %s\n",
            paste(round(tapply(got == truth, cut(maf_all, c(0, 0.1, 0.25, 0.5)), mean), 3), collapse = " / ")))

## ---------- C. the relationship matrix from the completed genotypes ----------
ids <- rownames(A)
Gb  <- Bf[ids, colnames(G)]                                # the 159 analysed accessions, as in script 01
stopifnot(identical(dimnames(Gb), dimnames(G)), all(Gb[!MISS] == G[!MISS]))
vanraden <- function(M) { p <- colMeans(M) / 2; K <- tcrossprod(sweep(M, 2, 2 * p)) / (2 * sum(p * (1 - p))); dimnames(K) <- list(rownames(M), rownames(M)); K }
stopifnot(max(abs(vanraden(G) - A)) < 1e-8)                # the function reproduces the matrix in use
Ab <- vanraden(Gb)
m  <- rowMeans(MISS)
cat("\nC. the relationship matrix\n")
cat(sprintf("   mean diagonal: in use %.4f | from Beagle %.4f\n", mean(diag(A)), mean(diag(Ab))))
cat(sprintf("   correlation between the two: diagonal %.4f, off-diagonal %.5f\n", cor(diag(A), diag(Ab)), cor(A[upper.tri(A)], Ab[upper.tri(A)])))
cat(sprintf("   correlation of the diagonal with the missing rate: in use %.3f | from Beagle %.3f\n", cor(diag(A), m), cor(diag(Ab), m)))
cat(sprintf("   smallest eigenvalue from Beagle %.1e (a Gram matrix, so not negative beyond rounding)\n",
            min(eigen(Ab, symmetric = TRUE, only.values = TRUE)$values)))
maf_b <- pmin(colMeans(Gb) / 2, 1 - colMeans(Gb) / 2)
cat(sprintf("   markers monomorphic after imputation: %d | minimum minor allele frequency %.4f\n", sum(maf_b == 0), min(maf_b)))

## ---------- D. the observed-trial fits ----------
Ts <- make_transforms(Omega)
scan <- function(Tm, K0, GG) {                             # as scan_arm in script 05
  y  <- as.numeric(Tm %*% blues$BLUE); X0 <- cbind(as.numeric(Tm %*% rep(1, n)), Tm %*% Q)
  Zt <- Tm %*% GG; K <- sym(Tm %*% K0 %*% t(Tm))
  nf <- mixed.solve(y = y, X = X0, Z = diag(n), K = K, return.Hinv = TRUE); Hi <- nf$Hinv
  P  <- Hi - Hi %*% X0 %*% solve(crossprod(X0, Hi %*% X0)) %*% crossprod(X0, Hi)
  PZ <- P %*% Zt; zPz <- colSums(Zt * PZ); zPy <- as.numeric(crossprod(PZ, y)); yPy <- as.numeric(crossprod(y, P %*% y))
  v2 <- n - 5; ok <- zPz > 1e-8 * max(zPz)
  b <- se <- lp <- rep(NA_real_, ncol(Zt))
  b[ok] <- zPy[ok] / zPz[ok]; s2 <- (yPy - b[ok] * zPy[ok]) / v2; se[ok] <- sqrt(s2 / zPz[ok])
  lp[ok] <- -log10(pbeta(v2 / (v2 + (b[ok] / se[ok])^2), v2 / 2, 1 / 2))
  list(Vg = nf$Vu, Ve = nf$Ve, LL = nf$LL, beta = setNames(b, colnames(GG)), se = setNames(se, colnames(GG)), lp = setNames(lp, colnames(GG)))
}
bonf <- -log10(ALPHA / ncol(G))
fits <- list()
for (nm in names(Ts)) {
  fits[[paste(nm, "in use")]] <- scan(Ts[[nm]], A,  G)
  fits[[paste(nm, "Beagle")]] <- scan(Ts[[nm]], Ab, Gb)
}
tab <- do.call(rbind, lapply(names(fits), function(k) { f <- fits[[k]]
  data.frame(fit = k, Vg = f$Vg, Ve = f$Ve, logLik = f$LL, above_threshold = sum(f$lp > bonf, na.rm = TRUE),
             test_beta = f$beta[TEST_SNP], test_se = f$se[TEST_SNP], test_lp = f$lp[TEST_SNP], row.names = NULL) }))
cat("\nD. second-stage fits and the genome scan, threshold", round(bonf, 3), "\n"); print(format(tab, digits = 5), row.names = FALSE)
## markers above the threshold under each fit, with their positions (Table S4)
sig <- do.call(rbind, lapply(names(fits), function(k) { s <- names(which(fits[[k]]$lp > bonf))
  if (!length(s)) return(NULL); i <- match(s, map$SNP)
  data.frame(fit = k, SNP = s, Chromosome = map$Chromosome[i], Position = map$Position[i], neglog10p = unname(fits[[k]]$lp[s]), row.names = NULL) }))

for (nm in names(Ts)) {
  a <- fits[[paste(nm, "in use")]]; b <- fits[[paste(nm, "Beagle")]]
  sa <- names(which(a$lp > bonf)); sb <- names(which(b$lp > bonf))
  cat(sprintf("   %-10s significant in both %d | only with the kernel in use %d | only with Beagle %d | correlation of -log10 p over all markers %.4f\n",
              nm, length(intersect(sa, sb)), length(setdiff(sa, sb)), length(setdiff(sb, sa)), cor(a$lp, b$lp, use = "complete.obs")))
}

d <- pheno; d$PC1 <- Q[as.character(d$Id), 1]; d$PC2 <- Q[as.character(d$Id), 2]; d$PC3 <- Q[as.character(d$Id), 3]
f0 <- as.formula(paste(TRAIT, "~ 1 + PC1 + PC2 + PC3 + Year + (1|Year:Replication) + (1|Id:Year) + (1|Id)"))
J  <- matrix(1, n, n)
ss <- function(K0) { K1 <- K0 + ONES_C * J; dimnames(K1) <- dimnames(A)
  fit <- relmatLmer(f0, data = d, relmat = list(Id = K1), REML = TRUE); vc <- as.data.frame(VarCorr(fit))
  c(Vg = vc$vcov[vc$grp == "Id"], Vge = vc$vcov[vc$grp == "Id:Year"], Ve = vc$vcov[vc$grp == "Residual"], REML = as.numeric(REMLcrit(fit))) }
s_use <- ss(A); s_bg <- ss(Ab)
gap <- data.frame(kernel = c("in use", "Beagle"),
                  two_stage_Vg = c(fits[["rotated in use"]]$Vg, fits[["rotated Beagle"]]$Vg),
                  single_stage_Vg = c(s_use["Vg"], s_bg["Vg"]),
                  rotated_residual = c(fits[["rotated in use"]]$Ve, fits[["rotated Beagle"]]$Ve),
                  single_stage_REML = c(s_use["REML"], s_bg["REML"]))
gap$gap_pct <- 100 * (gap$two_stage_Vg / gap$single_stage_Vg - 1)
cat("\n   the two numbers in question\n"); print(format(gap, digits = 5), row.names = FALSE)

save(acc, tab, gap, sig, Ab, Gb, file = file.path(HERE, "02_check_and_refit.RData"))
