This package implements posterior predictive checks (PPCs) for population admixture models as proposed in [this paper](http://arxiv.org/abs/1407.0050v1).
A *discrepancy function* measures some characteristic of a collection of observed alleles
and their inferred population assignments.
A PPC compares values of a discrepancy function applied to the actual data to values of that discrepancy function applied to *replicated* datasets, where the population assignment variables are the same but the observed alleles are actually sampled from the trained admixture model.
If the PPC fails, and the observed discrepancy value is far from the values for replicated data that supposedly has the same latent population structure, we may consider rejecting the admixture model.

This package consists of java classes and optional R scripts that produce graphical output. The R scripts require the `ggplot2` and `dplyr` packages. (The R scripts are mostly very similar, differing only in variable names.)

# Data formats

All PPCs take three files as input. Currently only unphased diploid data is supported.

* A [plink](http://pngu.mgh.harvard.edu/~purcell/plink/) formatted binary pedigree file (.bed)
* A file containing population weights for each individual (the .Q file from [admixture](https://www.genetics.ucla.edu/software/admixture/))
* A file containing population weights for each SNP (the .P file from [admixture](https://www.genetics.ucla.edu/software/admixture/))

# Discrepancy functions

This package supports five discrepancy functions, which are fully described in the paper accompanying this package.
All discrepancies produce population-specific results, so for a model with 3 populations, the result will be a vector of three values.
The one exception is the LD discrepancy, which produces values for each population at several inter-SNP distances, or "lags".

## Inter-genome similarity

The java package will write results to standard output:

    ./ppc sim hapmap3-files/hapmap3.bed hapmap3-files/hapmap3.3.Q hapmap3-files/hapmap3.3.P > sim.txt

You can then create a plot with the following R script:

    Rscript R/sim.R sim.txt HapMap

The argument `sim.txt` is the name of the file we created, `HapMap` is a title to put on the plot. The R script is a simplification of plots included in our paper that contained multiple collections and varying numbers of populations.
If you modify the first few lines to include multiple collections they will appear in separate panels within the plot.

## Linkage disequilibrium between nearby SNPs

Most datasets used in admixture modeling contain SNPs that are relatively distant from each other, since nearby SNPs may be strongly correlated. This discrepancy function asks whether the SNPs are sufficiently widely spaced, or if there is still predictive power between SNPs that are adjacent *in the filtered dataset*.
This discrepancy function is distinct from the other four functions because rather than a single number for each population, the function returns a number for each population and "lag", or spacing in the filtered dataset between comparable SNPs.
The PPC takes two extra arguments, the distance between SNPs and the size of the total window of comparison.
For example, this command

    ./ppc lag-mi hapmap3-files/hapmap3.bed hapmap3-files/hapmap3.3.Q hapmap3-files/hapmap3.3.P 5 30 > lag-mi.txt

will compare the SNP at position 0 to the SNPs at positions 5, 10, 15, 20, 25, and 30, and the SNP at position 5 to SNPs 10, 15, 20, 25, 30, and 35, and so forth.
Mutual information is averaged over each distinct *lag*, for example lag 5 will include MI between SNP 0 and SNP 5, SNP 5 and SNP 10, SNP 10 and SNP 15, and so forth.
You can create a plot of mutual informations at each lag with this command

    Rscript R/lag-mi.R lag-mi.txt HapMap

## Association between SNPs and observed genome labels (F_ST)

Unlike the other discrepancies, F_ST requires a file with string tags (one per line) for each individual. 

	./ppc fst hapmap3-files/hapmap3.bed hapmap3-files/hapmap3.3.Q hapmap3-files/hapmap3.3.P hapmap3-files/hapmap3.pops > fst.txt

If the tags file contains more than one column, you can also optionally specify the column number of the tag on the command line after the tag file.

    Rscript R/fst.R fst.txt HapMap


## Posterior entropy over topic assignments

    ./ppc entropy hapmap3-files/hapmap3.bed hapmap3-files/hapmap3.3.Q hapmap3-files/hapmap3.3.P  > entropy.txt


    Rscript R/entropy.R entropy.txt HapMap

## Association between synthetic phenotypes and SNPs

    ./ppc gwas hapmap3-files/hapmap3.bed hapmap3-files/hapmap3.3.Q hapmap3-files/hapmap3.3.P > gwas.txt



    Rscript R/gwas.R gwas.txt HapMap