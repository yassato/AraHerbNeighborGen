# prep

how to prepare the input files for the rehh package

1. use ```GetAraGenome.sh``` to get a multiple fasta format (.maf) of _A. thaliana_ and _A. lyrata_ genome sequences.
2. move to ```arabidopsis_thaliana_TAIR10_vs_arabidopsis_lyrata_v_1_0_lastz_net```
3. use ```maf2txt.py``` to compile ancestral (= _A. lyrata_) and derived (= _A. thaliana_) alleles in a csv format named "chrX.csv".
4. use ```reshape_chr_ann.R``` to reshape the csv files "chrX.csv" into input format of the rehh package

_Note: some output csv files are so large that should be kept outside the GitHub repository before pushing._
