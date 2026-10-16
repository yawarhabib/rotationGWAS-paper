################################################################################
##        Rotated two-stage GWAS for unbalanced multi-environment trials
##
##                      Script 01: data preparation
##
##   Loads the phenotypes and the genotypes (PLINK --recode A format, missing
##   calls as NA), restricts both to the accessions common to them, imputes
##   missing calls at the marker mean, builds the additive relationship matrix
##   and the principal-component covariates, and puts every object on one
##   accession ordering that all later scripts assume.
##
##   Output: 01_data.RData
##
##                  Authors: Yawar Habib, Hans-Peter Piepho
##                                R >= 4.4.0
################################################################################

source("scripts/00_config.R")
library(dplyr)
library(data.table)

t0 <- Sys.time()

## ---------- 1. phenotypes ----------
Kh <- read.csv(pheno_file, header = TRUE)

pheno_all <- Kh %>%
  dplyr::select(Year, Replication, Id = ID_in_Geno, Genotype, all_of(TRAIT)) %>%
  filter(!is.na(.data[[TRAIT]]))

## ---------- 2. genotypes ----------
raw <- as.data.frame(fread(geno_raw))
bim <- as.data.frame(fread(geno_bim, header = FALSE,
                           col.names = c("chr_str", "SNP", "cm", "Position", "A1", "A2")))

mk <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G  <- as.matrix(raw[, mk])
colnames(G) <- sub("_[ACGT]$", "", mk)      # remove the allele suffix that PLINK appends
rownames(G) <- as.character(raw$IID)

n_genotyped <- nrow(G)

## The map comes from the .bim file and the dosages from the .raw file;
## the check below confirms that the two correspond marker for marker.
stopifnot(identical(colnames(G), bim$SNP),
          !any(duplicated(colnames(G))))

map <- data.frame(SNP        = bim$SNP,
                  Chromosome = as.integer(sub(".*Gm0?", "", bim$chr_str)),
                  Position   = bim$Position)

PCA <- read.csv(geno_pcs, header = TRUE)

## ---------- 3. restrict to accessions present in both ----------
common    <- intersect(as.character(pheno_all$Id), rownames(G))
geno_only <- setdiff(rownames(G), as.character(pheno_all$Id))

pheno <- pheno_all %>% filter(as.character(Id) %in% common)
G     <- G[as.character(common), , drop = FALSE]

## ---------- 4. heterozygosity, measured before imputation ----------
## Measured on the raw calls, so that imputed values are not counted as
## heterozygous. The bound of the check below is slightly above the
## quality-control threshold of 0.10, because the statistic is computed on
## the analysed accessions only, one fewer than the genotyped set.
het_ind <- rowMeans(G == 1, na.rm = TRUE)
het_mk  <- colMeans(G == 1, na.rm = TRUE)

## ---------- 5. impute missing calls at the marker mean ----------
n_missing <- sum(is.na(G))

## The indicator of missing calls is kept, because after imputation an imputed
## value and an observed call look the same, and the heterozygosity kernel of
## script 04 uses observed calls only.
MISS <- is.na(G)

marker_mean <- colMeans(G, na.rm = TRUE)
stopifnot(!anyNA(marker_mean))              # a marker missing in every accession

ix <- which(is.na(G), arr.ind = TRUE)
G[ix] <- marker_mean[ix[, "col"]]

## ---------- 6. one accession ordering everywhere ----------
## sommer::A.mat expects the coding -1, 0, 1, so one is subtracted from the
## dosages for the call; G itself stays coded 0, 1, 2. The check below confirms
## that A is VanRaden's matrix computed on every marker.
A <- sommer::A.mat(G - 1)

pheno <- pheno %>% mutate(Id          = factor(as.character(Id), levels = rownames(A)),
                          Year        = as.factor(Year),
                          Replication = as.factor(Replication))

A    <- A[levels(pheno$Id), levels(pheno$Id)]
G    <- G[levels(pheno$Id), , drop = FALSE]
MISS <- MISS[levels(pheno$Id), , drop = FALSE]
map <- map[match(colnames(G), map$SNP), ]

PCA$Taxa <- as.character(PCA$Taxa)
PCA <- PCA[match(rownames(A), PCA$Taxa), ]
Q   <- as.matrix(PCA[, paste0("PC", 1:N_PC)])
rownames(Q) <- rownames(A)

n <- nrow(A)

## ---------- 7. checks ----------
checks <- data.frame(
  Check = c("A rownames == colnames",
            "A ordering == genotype ordering",
            "A ordering == PC ordering",
            "A ordering == levels(pheno$Id)",
            "map aligned to genotype columns",
            "no NA in A, genotypes or PCs",
            "A symmetric",
            "every accession has phenotype data",
            "no duplicated plot records",
            "dosages within [0, 2]",
            "test marker present",
            "chromosomes parsed",
            "marker heterozygosity consistent with the upstream filter",
            "A is VanRaden's matrix on every marker",
            "missing-call mask aligned with the genotypes and complete",
            "every observed call is an integer dosage"),
  Result = c(identical(rownames(A), colnames(A)),
             identical(rownames(A), rownames(G)),
             identical(as.character(PCA$Taxa), rownames(A)),
             identical(rownames(A), levels(pheno$Id)),
             identical(as.character(map$SNP), colnames(G)),
             !anyNA(A) && !anyNA(G) && !anyNA(Q),
             isSymmetric(unname(A), tol = 1e-10),
             all(table(pheno$Id) > 0),
             !any(duplicated(pheno[, c("Id", "Year", "Replication")])),
             min(G) >= 0 && max(G) <= 2,
             TEST_SNP %in% colnames(G),
             !anyNA(map$Chromosome),
             max(het_mk) < 0.11,
             { p_chk <- colMeans(G) / 2
               max(abs(A - tcrossprod(sweep(G, 2, 2 * p_chk)) / (2 * sum(p_chk * (1 - p_chk))))) < 1e-8 },
             identical(dimnames(MISS), dimnames(G)) && sum(MISS) == n_missing,
             all(G[!MISS] %in% c(0, 1, 2))))

## ---------- 8. summary ----------
ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values

summary_tab <- data.frame(
  Quantity = c("accessions genotyped",
               "accessions carried forward",
               "accessions dropped, no phenotype records",
               "markers", "plots", "years",
               "missing calls imputed", "missing rate",
               "accession heterozygosity, median",
               "accession heterozygosity, maximum",
               "marker heterozygosity, maximum",
               "markers above the 0.10 filter threshold",
               "A mean diagonal", "A effective rank"),
  ## as text, so that counts and proportions print side by side
  Value = c(format(c(n_genotyped, n, length(geno_only), ncol(G), nrow(pheno),
                     nlevels(pheno$Year), n_missing), big.mark = ",", trim = TRUE),
            format(round(n_missing / (n * ncol(G)), 4), nsmall = 4),
            format(round(median(het_ind), 4), nsmall = 4),
            format(round(max(het_ind), 4), nsmall = 4),
            format(round(max(het_mk), 4), nsmall = 4),
            format(sum(het_mk > 0.10), big.mark = ",", trim = TRUE),
            format(round(mean(diag(A)), 3), nsmall = 3),
            format(sum(ev > 1e-6 * max(ev)), trim = TRUE)))

## ---------- the layout of the trial ----------
## For each year: the number of blocks, the entries in each block, whether
## every block held every entry, and the plots grown and the plots scored for
## the trait, counted over the whole layout. The tables that follow count the
## analysed plots only.
plots_grown  <- nrow(Kh)
plots_scored <- sum(!is.na(Kh[[TRAIT]]))

layout_tab <- do.call(rbind, lapply(sort(unique(Kh$Year)), function(y) {
  x <- Kh[Kh$Year == y, ]
  per <- tapply(x$ID_in_Geno, x$Replication, function(z) length(unique(z)))
  sc  <- tapply(!is.na(x[[TRAIT]]), x$Replication, sum)
  data.frame(Year = y, Blocks = length(per),
             Entries_per_block = paste(per, collapse = ", "),
             Complete = all(per == length(unique(x$ID_in_Geno))),
             Plots_grown = nrow(x),
             Plots_scored = sum(!is.na(x[[TRAIT]])),
             Scored_per_block = paste(sc, collapse = ", "),
             row.names = NULL)
}))

incidence <- table(pheno$Id, pheno$Year)
design_tab <- data.frame(
  Year        = levels(pheno$Year),
  Accessions  = as.integer(colSums(incidence > 0)),
  Plots       = as.integer(colSums(incidence)),
  Replicates  = as.integer(sapply(levels(pheno$Year),
                  function(y) length(unique(pheno$Replication[pheno$Year == y])))))

years_per_accession <- rowSums(incidence > 0)
rolling_tab <- data.frame(
  Quantity = c("years per accession, median",
               "accessions grown in all years",
               "accessions grown in four years or fewer"),
  Value = c(median(years_per_accession),
            sum(years_per_accession == nlevels(pheno$Year)),
            sum(years_per_accession <= 4)))

cat("\n=== 01 data preparation ===\n")
print(checks, row.names = FALSE)
cat("\n")
print(summary_tab, row.names = FALSE)
if (length(geno_only))
  cat("\ndropped for lack of phenotype records:", paste(geno_only, collapse = ", "), "\n")
cat("\n")
print(design_tab, row.names = FALSE)
cat("\n")
print(rolling_tab, row.names = FALSE)

cat("\nthe design as laid out, before filtering\n")
print(layout_tab, row.names = FALSE)
cat("   plots grown:", format(plots_grown, big.mark = ","),
    "| plots with a", TRAIT, "value:", format(plots_scored, big.mark = ","), "\n")
cat("   a block grown but not scored is part of the design and absent from",
    "the analysis\n")

if (!all(checks$Result)) stop("checks failed")

save(pheno, G, A, Q, map, n, TRAIT, TEST_SNP, MISS,
     n_missing, n_genotyped, geno_only, het_ind, het_mk,
     summary_tab, design_tab, rolling_tab, layout_tab,
     plots_grown, plots_scored,
     file = file.path(OUT, "01_data.RData"))

cat("\nsaved 01_data.RData |", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
