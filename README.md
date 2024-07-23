# AraHerbNeighborGen

Source codes and phenotype data for the analysis of **Ara**bidopsis **Herb**ivory with **Neighbor** **Gen**otypic effects  

## Subdirectory

-   **./genoData**\
    Scripts to prepare SNP data

-   **./insectData**\
    Data and scripts to compile phenotype data

-   **./neiGWAS**\
    Scripts to perform variation partitioning and association mapping, with additional scripts for post-GWAS analyses and visualization

-   **./gscan**\
    Scripts to scan signatures of positive and balancing selection

-   **./glmnetPy**\
    Scripts to perform LASSO regressions and estimate effect sizes, with additional scripts for post-LASSO analyses and plot functions

-   **./MixPlt**\
    Data and scripts for mixed planting experiments

-   **./ChoiceExp**\
    Data and scripts for laboratory choice experiments

-   **./simulation**\
    Numerical analysis and simulations for supplementary materials

-   **./output**\
    Directory to store intermediate files


## The other files
- AraHerbGWAS.Rproj  
A RStudio project to run the analyses above

- coord.R  
An additional source code to draw a Manhattan plot

- README.md  
This description file. Note that each subdirectory includes its own 'README.md' to run a subset of the analyses.  

## Usage

The workflow proceeds as follows 1-8. For each step, please read `README.md` in each subdirectory. Before running each script, make sure to set current working directory here at `usr/AraHerbGWAS/`.

1.  see `./genoData`, and prepare input genotypes
2.  see `./insectData`, and prepare input phenotypes
3.  see `./neiGWAS`, and perform PVE calculation and association mapping
4.  see `./gscan`, and perform genome-wide selection scan
5.  see `./glmnetPy`, and perform LASSO and estimate the effects of mixed planting
6.  see `./MixPlt`, and analyze data on mixed planting experiments
7.  see `./ChoiceExp`, and analyze data on laboratory choice experiments
8.  see `./simulation`, and perform additional tests of the model performance

Please also make `../data/`, `../output/`, and `../figs/` above this working space to save large files outside the git repository. Too large input, intermediate, and figure files should be placed there. A part of output files neccesary for the downstream analysis are stored within this repository named `./output`.

## External packages, libraries, and data

### Statistical analysis

1.  Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2020). vegan: Community Ecology Package. <https://CRAN.R-project.org/package=vegan>

2.  Hervé Perdry and Claire Dandine-Roulland (2020). gaston: Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models. <https://CRAN.R-project.org/package=gaston>

3.  Sato Y, Yamamoto E, Shimizu KK, Nagano AJ (2021) Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory. Heredity 126(4):597-614. <https://cran.r-project.org/package=rNeighborGWAS>

4.  Gautier M, Vitalis R (2012). rehh: An R package to detect footprints of selection in genome-wide SNP data from haplotype structure. Bioinformatics, 28(8) 1176-1177

5.  Gautier M, Klassmann A, Vitalis R (2017). rehh 2.0: a reimplementation of the R package rehh to detect positive selection from haplotype structure. Molecular Ecology Resources, 17(1) 78-90

6.  K. M. Siewert, B. F. Voight, Detecting long-term balancing selection using allele frequency correlation. Molecular Biology and Evolution. 34, 2996--3005 (2017). <https://github.com/ksiewert/BetaScan>

7.  Balakumar, B.J., Hastie, T., Friedman, J., Tibshirani, Simon, N.R. (2016). Glmnet for Python. <http://hastie.su.domains/glmnet_python/>; <https://github.com/bbalasub1/glmnet_python/> (accessed on 19 Sept. 2019)

8.  F. Supek, M. Bošnjak, N. Škunca, T. Šmuc, REVIGO summarizes and visualizes long lists of gene ontology terms. PLoS One. 6, e21800 (2011).

9.  S. Sayols, rrvgo: A bioconductor package to reduce and visualize gene ontology terms (2020; <https://ssayols.github.io/rrvgo>).

10. Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. <doi:10.18637/jss.v067.i01>.

11. Kuznetsova A, Brockhoff PB, Christensen RHB (2017). "lmerTest Package: Tests in Linear Mixed Effects Models." *Journal of Statistical Software*, *82*(13), 1-26. doi: 10.18637/jss.v082.i13 (URL: <https://doi.org/10.18637/jss.v082.i13>).

12. Russell V. Lenth (2021). emmeans: Estimated Marginal Means, aka Least-Squares Means. <https://CRAN.R-project.org/package=emmeans>

13. Bates D, Maechler M, Jagan M (2023). Matrix: Sparse and Dense Matrix Classes and Methods. <https://CRAN.R-project.org/package=Matrix>.


### Visualization

1.  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

2.  Thomas Lin Pedersen (2020). patchwork: The Composer of Plots. <https://CRAN.R-project.org/package=patchwork>

### Utilities

1.  Kevin Ushey, JJ Allaire and Yuan Tang (2021). reticulate: Interface to 'Python'. <https://CRAN.R-project.org/package=reticulate>

2.  Randy Lai (2020). arrangements: Fast Generators and Iterators for Permutations, Combinations, Integer Partitions and Compositions. <https://CRAN.R-project.org/package=arrangements>

3.  Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, <https://doi.org/10.21105/joss.01686>

4.  Hadley Wickham, Jim Hester and Winston Chang (2021). devtools: Tools to Make Developing R Packages Easier. <https://CRAN.R-project.org/package=devtools>

5. Yasuhiro Sato, & Atsushi J. Nagano. (2023). GOfisher (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.7901509  

### External data

1.  Togninalli M, Seren U, Meng D, Fitz J, Nordborg M, Weigel D, Borgwardt K, Korte A, Grimm DG. (2018). The AraGWAS catalog: A curated and standardized *Arabidopsis thaliana* GWAS catalog. Nucleic Acids Research 46:D1150--56. <https://aragwas.1001genomes.org/>

2.  Togninalli M, Seren U, Freudenthal JA, Monroe JG, Meng D, Nordborg M, Weigel D, Borgwardt K, Korte A, Grimm DG. (2020). AraPheno and the AraGWAS Catalog 2020: a major database update including RNA-Seq and knockout mutation data for *Arabidopsis thaliana*. Nucleic acids research, 48(D1), D1063-D1068. <https://arapheno.1001genomes.org/>

3.  T. Berardini, L. Reiser, E. Huala, TAIR functional annotation data (2021), <doi:10.5281/zenodo.7159104>.

4.  M. Carlson, GO.db: A set of annotation maps describing the entire Gene Ontology (2020), <https://doi.org/10.18129/B9.bioc.GO.db>).

5.  M. Carlson, org.At.tair.db: Genome wide annotation for Arabidopsis (2019), <https://doi.org/10.18129/B9.bioc.org.At.tair.db>).

6. Horton, M. W. et al. Genome-wide patterns of genetic variation in worldwide Arabidopsis thaliana accessions from the RegMap panel. Nature Genetics 44, 212–216 (2012). Data available at <http://bergelson.uchicago.edu/?page_id=790> (accessed on 9 Feb. 2017)


## Working environment and software version

### For association mapping, LASSO, and selection scan

Intel(R) Xeon(R) Gold 6126 CPU @ 2.60GHz $\times$ 24 cores  
Memory 384 GB with CentOS7  

**R v3.6.0**
- gaston: v1.5.5
- rNeighborGWAS: v1.2.4
- Matrix: v1.2-17
- reticulate: v1.15
- rehh: v3.1.2

**Python v3.6.8**
- Anaconda 3
- glmnet for Python: v1.0

### For the other analyses

Apple M2 Max $\times$ 12 CPU cores  
Memory 96GB with macOS Ventura  

**R v4.3.0**
- gaston: v1.5.9
- rNeighborGWAS: v1.2.4
- tidyverse: v2.0.0
- patchwork: v1.1.2
- arrangement: v1.1.9
- devtools: v2.4.5
- lme4: v1.1-34
- lmerTest: v3.1-3
- emmeans: v1.8.7
- vegan: v2.6-4
- GOfisher: v0.1.0
- GO.db: v3.17.0
- rrvgo: v1.12.0
- org.At.tair.db: v3.17.0
