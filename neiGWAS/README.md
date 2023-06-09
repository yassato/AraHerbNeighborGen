# Neighbor GWAS of the insect herbivory and abundance
Variation partitioning, association mapping, and post GWAS analyses with the neighbor genotype effects.  

## Input files
- Genotype inputs can be prepared using scripts in "../genoData/" subdirectory.  
- Phenotype inputs can be prepared using scripts in "../insectData/" subdirectory.  

_Note: some intermediate files are so large that should be save outside the GitHub repository (e.g., ```../output/```)._

## Scripts
In the line of scripts, s = 1 and s = 2 corresponds to J = 4 (i.e., a spatial scale up to the nearest individuals) and J = 12 (up to the second nearest individuals), respectively. Scripts should be run in the following order.    

1. calcPVE_CHZ(or JPN).R  
R script to calculate the proportion of variation explained (PVE) by self or neighbor genotypes in the Zurich (CHZ) or Japan (JPN) site. The reference space can be changed using the ```scale``` argument in the script. This generates output files named like "CHZ/JPNoutSx.csv".  

1. NeighborCHZ(or JPN).R  
R script to perform association mapping with focal or neighbor genotypes in the Zurich (CHZ) or Japan (JPN) site. The reference space can be changed using the ```scale``` argument in the script. This generates output files named like "(CHZ/JPN)neiGWAS_(TRAIT_NAME)_Sx.csv".  

1. neiGWAS_genelist.py  
Python3 script to list candidate genes from the output of "NeighborCHZ(or JPN).R"  

1. neiGWAS_ManhattanPlot.R  
R script to depict Manhattan plots and other supplementary figures for the results of Neighbor GWAS. To complete a composite figure, run ```../gscan/``` before this script.     
