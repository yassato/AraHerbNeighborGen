# Preparetion of input genotypes

## Usage
1. First, run subsetSNP.py to extract all SNPs and posiitions from .hdf5.  
2. Then, run reshapeSNP.R to prepare SNPs and their positions for target accessions.  

_Note also that output files are so large that should be saved outside the GitHub repository (e.g., ```../data/```)._

## Genotype data  
The full imputed SNP data is available at https://aragwas.1001genomes.org/#/download-center  
1. Access via "For other genotype versions, please refer to this _google drive page_" (accessed on 26 Oct 2022).
2. Download a compressed file named "full_imputed_SNP_MATRIX_2Jun2022.tar.gz", which includes "all_chromosomes_binary.hdf5" used in subsetSNP.py.

## References
- Matteo Togninalli, Ümit Seren, Dazhe Meng, Joffrey Fitz, Magnus Nordborg, Detlef Weigel, Karsten Borgwardt, Arthur Korte, and Dominik G. Grimm (2018) The AraGWAS Catalog: a curated and standardized Arabidopsis thaliana GWAS catalog; Nucleic Acids Research, gkx954, https://doi.org/10.1093/nar/gkx954

- The 1001 Genomes Consortium (2016) 1,135 genomes reveal the global pattern of polymorphism in Arabidopsis thaliana. Cell 166, no. 2: 481-491.https://doi.org/10.1016/j.cell.2016.05.063
