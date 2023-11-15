# glmnet regression using Python and R    

R and Python scripts to perform the LASSO regressions and cross-validation. Note that R version of the glmnet package does not work with millions of explanatory variables; thus, [glmnetPy](https://github.com/bbalasub1/glmnet_python) (Balakumar et al. 2016) is used instead.

## Input files
Genotypes and phenotypes should be prepared beforehand using scripts in ```../genoData``` and ```../insectData``` subdirectories. Some intermediate files are so large that should be saved outside the GitHub repository (e.g., ```../output/```).

## Scripts
In the line of scripts, s = 1 and s = 2 corresponds to J = 4 (a spatial scale up to the nearest individuals) and J = 12 (up to the second nearest individuals), respectively. Scripts should be run in the following order.    

1. glmnetPy_export4CHZ(or JPN).R  
R script to export input files for Python for the Zurich (CHZ) or Japan (JPN) site. Note that the output pickle files are as large as several GBs.

1. glmnetPy_CHZ(or JPN).py  
Python3 script to perform LASSO regressions and to save the results as a zipped pickle file.  

1. glmnetPy_coefExport.py  
Python3 script to extract estimated coefficients from the pickled results.    

1. CHZ(or JPN)2019CV.R  
R script to perform the cross-validation of LASSO models with the 2019 dataset for the Zurich (CHZ) or Japan (JPN) site.

1. glmnetPy_LASSOfigs.R  
R script to depict figures for LASSO results  

1. glmnetPy_genelist.py  
Python script to list up candidate genes for the LASSO results. The same usage as ```../neiGWAS/neiGWAS_genelist.py```.

1. MixtureEffectSize.R  
R script to estimate pairwise effect sizes of mixed planting, with additional plot functions for supplementary figures.  

1. SimulateEffectSize.R  
R script to simulate and visualize relationships between the herbivore damage and plant genotypic richness.  

1. LASSO_GO.R  
R script to perform a gene ontology (GO) enrichment analysis for the LASSO results  

## References
- Balakumar, B.J., Hastie, T., Friedman, J., Tibshirani, Simon, N.R. (2016) Glmnet for Python. http://hastie.su.domains/glmnet_python/  
