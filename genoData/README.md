# Preparetion of input genotypes

## Usage
1. First, run `subsetSNP.py` to extract all SNPs and posiitions from .hdf5.  
2. Then, run `reshapeSNP.R` to prepare SNPs and their positions for target accessions.  

_Note also that output files are so large that should be saved outside the GitHub repository (e.g., ```../data/```)._

## Genotype data  
The full imputed SNP data is available at https://aragwas.1001genomes.org/#/download-center  
1. Access via "For other genotype versions, please refer to this _google drive page_" (accessed on 26 Oct. 2022).
2. Download a compressed file named "full_imputed_SNP_MATRIX_2Jun2022.tar.gz", which includes "all_chromosomes_binary.hdf5" used in subsetSNP.py.

## Intermediate files
- positionsMAF5.rds  
Position file for GWAS.   

- positionsMAF5LD80.rds  
Position file for LASSO. Adjacent loci with $r^2$ > 0.8 are prunned.  

- sub_snpMAF5.rds  
SNP matrix for GWAS.  

- sub_snpMAF5LD80.rds  
SNP matrix for LASSO. Adjacent loci with $r^2$ > 0.8 are prunned.  

## Annotaiton files
The following annotation files for the *A. thaliana* genome are available within this directory.
- Araport11_genes.201606.transcript.rep_ERCC_Virus7457_GFP_GUS.txt.gz  
A gene description file for *A. thaliana*.  

- TAIR10_GFF3_genes.gff.gz  
Gene feature format (.gff) for *A. thaliana*, which is obtained TAIR data.  

- ulgSLIM.TAIR_211201  
GO slim annotation obtained from TAIR data.  

- AtGWASlocality_added.csv  
A file describing the locality of *A. thaliana* accessions. The locality available at the AraGWAS Catalog.   

- gwasIDlist.csv  
A file describing GWAS genotype ID of *A. thaliana* accessions.  


## References
- Matteo Togninalli, Uemit Seren, Dazhe Meng, Joffrey Fitz, Magnus Nordborg, Detlef Weigel, Karsten Borgwardt, Arthur Korte, and Dominik G. Grimm (2018) The AraGWAS Catalog: a curated and standardized Arabidopsis thaliana GWAS catalog; Nucleic Acids Research, gkx954, <https://doi.org/10.1093/nar/gkx954>

- The 1001 Genomes Consortium (2016) 1,135 genomes reveal the global pattern of polymorphism in Arabidopsis thaliana. Cell 166, no. 2: 481-491. <https://doi.org/10.1016/j.cell.2016.05.063>
