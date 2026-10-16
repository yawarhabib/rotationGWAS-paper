################################################################################
##   Imputation check, step 1: write the genotypes as VCF for Beagle
##
##   The analysis fills a missing marker call with the marker mean. This folder
##   tests what changes if the gaps are filled by an LD-based method. It reads
##   the inputs of the analysis and writes only into this folder.
##
##   Two files are written from the same 160 genotyped individuals and 59,206
##   markers that script 01 reads:
##     geno_full.vcf    every observed call, missing calls as ./.
##     geno_masked.vcf  the same with a random 2% of the OBSERVED calls hidden,
##                      so that the accuracy of the imputation can be measured
##                      on calls whose true value is known
##   The alternate allele is the allele PLINK counted in the .raw file, so a
##   dosage of 0, 1, 2 becomes 0/0, 0/1, 1/1 and maps back without ambiguity.
################################################################################

source("scripts/00_config.R")
library(data.table)
HERE <- file.path(ROOT, "imputation_check")

raw <- as.data.frame(fread(geno_raw))
bim <- as.data.frame(fread(geno_bim, header = FALSE, col.names = c("chr", "SNP", "cm", "pos", "A1", "A2")))
mk  <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G   <- as.matrix(raw[, mk]); rownames(G) <- as.character(raw$IID)
counted <- sub(".*_([ACGT])$", "\\1", mk)
snp     <- sub("_[ACGT]$", "", mk)
stopifnot(identical(snp, bim$SNP), all(G[!is.na(G)] %in% c(0, 1, 2)),
          all(counted == bim$A1 | counted == bim$A2))
other <- ifelse(counted == bim$A1, bim$A2, bim$A1)

## Beagle needs the markers of a chromosome in order of position
ord <- order(match(bim$chr, unique(bim$chr)), bim$pos)
stopifnot(identical(ord, seq_along(ord)))                       # already in order
dup <- duplicated(paste(bim$chr, bim$pos))
cat("individuals", nrow(G), "| markers", ncol(G), "| chromosomes", length(unique(bim$chr)),
    "| markers sharing a position with an earlier one", sum(dup), "\n")
cat("missing calls", sum(is.na(G)), "of", length(G), sprintf("(%.4f)\n", mean(is.na(G))))

write_vcf <- function(M, file) {
  gt <- matrix("./.", nrow(M), ncol(M))
  gt[which(M == 0)] <- "0/0"; gt[which(M == 1)] <- "0/1"; gt[which(M == 2)] <- "1/1"
  body <- data.table(bim$chr, bim$pos, bim$SNP, other, counted, ".", "PASS", ".", "GT", t(gt))
  con <- file(file, "wb")
  writeLines(c("##fileformat=VCFv4.2",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", rownames(M)), collapse = "\t")),
             con, sep = "\n")
  close(con)
  fwrite(body, file, append = TRUE, sep = "\t", col.names = FALSE, quote = FALSE, eol = "\n")
}

write_vcf(G, file.path(HERE, "geno_full.vcf"))

set.seed(SEED_MISC + 900)                 # its own seed
obs  <- which(!is.na(G))
hide <- sample(obs, round(0.02 * length(obs)))
Gm <- G; Gm[hide] <- NA
write_vcf(Gm, file.path(HERE, "geno_masked.vcf"))
cat("masked", length(hide), "observed calls; masked file has", sum(is.na(Gm)), "missing calls\n")

save(G, hide, bim, counted, other, file = file.path(HERE, "01_vcf_inputs.RData"))
