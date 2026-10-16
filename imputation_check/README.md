# Imputation check

A sensitivity analysis. The main analysis fills a missing marker call with the
marker mean; this folder asks what changes if the gaps are filled by an LD-based
method, Beagle 5.5. It reads the inputs and the stored results of the main
analysis and writes only into this folder. Script F2 reads its results for
Table S4.

## How to run

Every command is run from the repository root, after scripts 01 and 02 of the
main analysis.

1. `Rscript imputation_check/01_write_vcf.R` writes `geno_full.vcf` (the 160
   genotyped individuals and 59,206 markers, missing calls as `./.`) and
   `geno_masked.vcf` (the same with a random 2% of the observed calls hidden,
   so that the accuracy of the imputation can be measured on calls whose true
   value is known).
2. Beagle 5.5 (`beagle.27Feb25.75f.jar`, from
   https://faculty.washington.edu/browning/beagle/; Browning, Zhou and Browning
   2018, American Journal of Human Genetics 103: 338-348), run with Java, no
   genetic map and no reference panel, all other settings default. With the jar
   in this folder:

       java -Xmx6g -jar imputation_check/beagle.27Feb25.75f.jar gt=imputation_check/geno_full.vcf out=imputation_check/imputed_full nthreads=2 seed=20261002
       java -Xmx6g -jar imputation_check/beagle.27Feb25.75f.jar gt=imputation_check/geno_masked.vcf out=imputation_check/imputed_masked nthreads=2 seed=20261002

   Each run takes about 15 minutes on two threads.
3. `Rscript imputation_check/02_check_and_refit.R` checks the imputation
   (observed calls unchanged, accuracy on the hidden calls) and refits the
   observed-trial models with the relationship matrix built from the completed
   genotypes.
4. `Rscript imputation_check/03_where_the_change_comes_from.R` separates the
   effect of the relationship matrix from the effect of the marker columns.
5. `Rscript imputation_check/04_residual_at_one.R` computes the deviance paid
   for fixing the rotated residual variance at one under each relationship
   matrix.

The jar, the VCF files and the `.RData` results are not part of the
repository; steps 1 to 5 regenerate them.
