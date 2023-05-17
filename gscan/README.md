# Genome-wide selection scan for *Arabidopsis thaliana*
Scripts to scan positive and balancing selection in the genome of _A. thaliana_.  

## Subdirectory

1. **./prep**  
Comparing the A. lyrata and A. thaliana genome to estimate ancestral alleles.  

2. **./rehh**  
R scripts to scan the signature of positive selection. The extended haplotype homozygosity (EHH: Sabeti et al. 2002) and its integrated score are calculated using the rehh package (Gautier et al. 2017).  

3. **./BetaScan**  
R and Python script to scan the signature of balancing selection. Beta1 statistics are calculated using the Python script (Siewert & Voight 2017). The input files are created from the rehh format.  

## Usage
Proceed with the order of the three sub-directories above. See ```README.md``` in each subdirectory to run each analysis. The R script ```intersect_scan_gwas.R``` was finally used to perform the Fisher tests and figure preparation. Another directory ```../neiGWAS``` needs to run before the Fisher tests in this ```gscan``` directory.  

## References
- Sabeti, P. C., Reich, D. E., Higgins, J. M., Levine, H. Z., Richter, D. J., Schaffner, S. F. et al. (2002). Detecting recent positive selection in the human genome from haplotype structure. Nature, 419(6909), 832-837.  
- Gautier, M., Klassmann, A., Vitalis R. (2017). “rehh 2.0: a reimplementation of the R package rehh to detect positive selection from haplotype structure.” Molecular Ecology Resources, 17(1), 78-90.  

- Siewert, K. M., & Voight, B. F. (2017). Detecting long-term balancing selection using allele frequency correlation. Molecular Biology and Evolution, 34(11), 2996-3005.  
