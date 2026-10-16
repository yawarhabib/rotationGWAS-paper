################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##     Script 00: configuration. Paths, constants and the functions shared by
##     the other scripts. Computes nothing.
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################


## The scripts are run from the repository root, which is the working directory.
ROOT <- normalizePath(getwd(), winslash = "/")
if (!dir.exists(file.path(ROOT, "data")) || !dir.exists(file.path(ROOT, "scripts")))
  stop("run the scripts from the repository root: the folders data/ and scripts/ were not found in ", ROOT)

## ---------- data files ----------
pheno_file <- file.path(ROOT, "data/trial_phenotypes_vegetation_period.csv")

geno_raw <- file.path(ROOT, "data/genotypes_160_lines_59206_markers.raw")
geno_bim <- file.path(ROOT, "data/genotypes_160_lines_59206_markers.bim")
geno_pcs <- file.path(ROOT, "data/principal_components_3.csv")

## ---------- output ----------
OUT <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------- random number generator ----------
## Fixed so that the seeds give the same draws in every version of R.
RNGkind(kind        = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection")

## ---------- analysis constants ----------
TRAIT    <- "Vegetation_period"
## The marker with the smallest p-value in the genome scan of script 05;
## script 03 uses it before the scan.
TEST_SNP <- "glyma.Wm82.gnm4.Gm19_48709768"
N_PC     <- 3
MAF_MIN  <- 0.15
ALPHA    <- 0.05
NUGGET   <- 1e-4

## ---------- replicate counts ----------
NREP_A_NULL   <- 1000
NREP_A_CAUSAL <- 1000
NREP_B        <- 1000
NREP_C        <- 1000
NREP_D        <- 1000
NREP_BOOT     <- 1000
NREP_P3D      <- 1500

## ---------- the simulation loop ----------
## The simulation scripts run their replicates in batches on SIM_CORES parallel
## R processes (default 6) and save the results after each batch, so that an
## interrupted run can be continued. Each replicate sets its own seed.
NCORE_SIM <- as.integer(Sys.getenv("SIM_CORES", "6"))

run_in_blocks <- function(cl, idx, fun, partial_file, block_size, t0 = Sys.time()) {
  key  <- as.character(idx)
  res  <- setNames(vector("list", length(idx)), key)
  done <- setNames(logical(length(idx)), key)
  if (file.exists(partial_file)) {
    e <- new.env(); load(partial_file, envir = e)
    if (exists("idx_partial", e) && identical(e$idx_partial, idx)) {
      done[] <- e$done_partial[key]
      res[done] <- e$res_partial[key[done]]
      cat(sprintf("  resumed from %s: %d of %d replicates already done\n",
                  basename(partial_file), sum(done), length(idx)))
    } else cat(sprintf("  %s ignored: not a partial file of this run\n", basename(partial_file)))
  }
  todo <- idx[!done]
  if (length(todo)) {
    blocks <- split(todo, ceiling(seq_along(todo) / block_size))
    for (b in seq_along(blocks)) {
      out <- parallel::parLapply(cl, blocks[[b]], fun)
      kb  <- as.character(blocks[[b]])
      res[kb] <- out; done[kb] <- TRUE
      res_partial <- res; done_partial <- done; idx_partial <- idx
      save(res_partial, done_partial, idx_partial, file = partial_file)
      cat(sprintf("  block %d of %d done, replicates up to %d, %.1f min\n", b, length(blocks),
                  max(blocks[[b]]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
  done_file <- done_name(partial_file)
  if (file.exists(partial_file)) file.rename(partial_file, done_file)
  cat(sprintf("  all %d replicates done, %s kept as %s until the script's own save\n",
              length(idx), basename(partial_file), basename(done_file)))
  unname(res)
}
done_name <- function(partial_file) {
  if (grepl("_partial\\.RData$", partial_file)) sub("_partial\\.RData$", "_done.RData", partial_file)
  else sub("\\.RData$", "_done.RData", partial_file)
}

## ---------- seed bases, one per scenario ----------
## The seed of a replicate is the base of its scenario, plus ten times the
## replicate number, plus a small index; the scenarios never overlap.
SEED_A    <- 1000000     # scenario A, type-I error and known effects
SEED_B    <- 2000000     # scenario B, block-variance sweep
SEED_C    <- 3000000     # scenario C, environment-size sweep
SEED_D    <- 4000000     # scenario D, residual as specification test
SEED_MISC <- 5000000     # bootstrap, marker selection, P3D sample, invariance

## ---------- the three weighting schemes, as transformations ----------
make_transforms <- function(Omega) {
  n  <- nrow(Omega)
  sv <- svd(Omega)

  ## Omega must be positive definite:
  ## the rotation divides by the square roots of its eigenvalues.
  ev <- eigen((Omega + t(Omega)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0 || max(ev) / min(ev) > 1e10)
    stop("Omega is not safely positive definite: smallest eigenvalue ",
         signif(min(ev), 3), ", condition number ",
         signif(max(ev) / min(ev), 3))

  list(unweighted = diag(n),
       diagonal   = diag(1 / sqrt(diag(Omega))),
       rotated    = sv$u %*% diag(1 / sqrt(sv$d)) %*% t(sv$u))
}

## ---------- the matrix of ones for singular kernels ----------
## K + cJ, with J the matrix of ones and c > 0, is positive definite,
## and its inverse is a generalised inverse of K. Under the rotation cJ
## becomes the outer product of the rotated intercept, a column of the
## fixed-effect design matrix, so the residual likelihood is unchanged and
## only the variance of the intercept estimate depends on c.
## Used where the fitting software factorises the kernel (lme4qtl, glmmTMB).
ONES_C <- 1
with_ones <- function(M, Tm) {
  Mt <- Tm %*% M %*% t(Tm); Mt <- (Mt + t(Mt)) / 2
  Mt + ONES_C * tcrossprod(as.numeric(Tm %*% rep(1, nrow(M))))
}

## ---------- the heterozygosity kernel, from called genotypes only ----------
## H = C C' / m, with C the centred indicator of a heterozygous call over the
## m markers at which both classes occur, scaled to the mean diagonal `target`.
## A missing call (the indicator MISS of script 01) contributes nothing.
het_kernel <- function(G, MISS, target) {
  stopifnot(identical(dim(G), dim(MISS)), identical(dimnames(G), dimnames(MISS)),
            all(G[!MISS] %in% c(0, 1, 2)))
  I <- (G == 1) * 1
  I[MISS] <- NA
  n_het <- colSums(I, na.rm = TRUE)
  n_obs <- colSums(!MISS)
  keep  <- n_het > 0 & n_het < n_obs
  I  <- I[, keep, drop = FALSE]
  cm <- colMeans(I, na.rm = TRUE)
  ix <- which(is.na(I), arr.ind = TRUE)
  I[ix] <- cm[ix[, 2]]
  C <- scale(I, center = TRUE, scale = FALSE)
  H <- tcrossprod(C) / ncol(I)
  H <- H * (target / mean(diag(H)))
  dimnames(H) <- list(rownames(G), rownames(G))
  list(K = H, m = ncol(I), n_het_calls = sum(G[!MISS] == 1))
}

## A nugget regularisation of a singular kernel; script S4 compares it with the
## matrix of ones.
regularise <- function(M) {
  e <- max(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
  M + diag(nrow(M)) * NUGGET * e
}

cat("config loaded ->", OUT, "\n")
cat(R.version.string, "| SIM_CORES =", Sys.getenv("SIM_CORES", "default 6"), "|",
    paste(vapply(c("lme4", "lme4qtl", "emmeans", "Matrix", "rrBLUP", "glmmTMB", "sommer"),
                 function(p) paste(p, as.character(packageVersion(p))), ""), collapse = ", "), "\n")
