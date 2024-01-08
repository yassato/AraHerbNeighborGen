# Compiling insect survey data as GWAS phenotypes
The original survey data are stored as .csv files named "AthInsectSurveyXXXX.csv". R scripts can be used to convert the original files as input phenotypes for GWAS.  

## Scripts
1. AthInsectData.R  
R script to convert the original files 

2. InsectRDA.R  
R script to perform an ordination analysis by RDA and then depict figures. "InsectRDAandNo.pdf" is its output figure. Another directory ```../neiGWAS``` needs to run beforehand to complete a composite figure.  

## Data files

1. AthInsectSurveyXXXX.csv  
csv files of the original survey data.  
Years: 2018; 2017; or 2019.  
Site: ZH, Zurich, Switzerland; or JP, Otsu, Japan.  

2. SurveyXXXX4GWAS_max.csv  
Output files from "AthInsectData.R". Phenotypes are summerized as the maximum number throughout the trial.  

## Header info for csv files
- IndivID: Individual plant ID  
- Name: Name of accession  
- Source: CS ID by the ABRC stock center      
- gwasID: GWAS ID linked for genotypes  
- Site: Zurich (ZH) or Japan (JP)  
- Year: Study year  
- Date: Study date of each survey; or starting date for "SurveyXXXX4GWAS_max.csv"
- Block: Experimental block ID
- position_X: x-axis position within the block  
- position_Y: y-axis position within the block  
- edge: Edge of the block (1) or not (0)  
- InitLeafLen: The length of the largest leaf (mm) at the start of the field experiment  
- Bolting: The presence (1) or absence (0) of bolting 2-wk after the start of the experiment  
- Score: Leaf damage score from 0 (almost intact) to 5 (almost complete loss of leaf area); available only for the JP site  
- Holes: No. of leaf holes made by flea beetles  
- Mines: No. of mining traces on leaves  
- Ps: No. of adult *Phyllotreta striolata* (incl. _P. undulata_ in Zurich)  
- Pa: No. of adult *Phyllotreta atra*
- Pp: No. of adult *Psylliodes punctifrons*  
- Lc: No. of adult *Listroderes costirostris*  
- Px: No. of adult *Plutella xylostella*  
- Pr_e: No. of eggs of *Pieris rapae*  
- Pr_l: No. of larval *Pieris rapae*  
- Ar: No. of larval *Athalia rosae*
- Tni: No. of larval *Trichoplusia ni*
-	Syrphinae_l: No. of larval hover flies  
- Cv: No. of parasites by *Cotesia vestalis*  
-	Fi: No. of *Frankliniella intonsa*  
-	Fo: No. of *Frankliniella	occidentalis*  
-	Bh: No. of *Bourletiella	hortensis*   	
- Er:	No. of *Eurydema 	rugosa*  
- Mp: No. of wingless *Myzus persicae*  
- Mp_w: No. of winged *Myzus persicae*	
- Le: No. of wingless *Lipaphis erysimi*  	
- Le_w: No. of winged *Lipaphis erysimi*	
- Bb: No. of wingless *Brevicoryne brassicae*  
- Bb_w: No. of winged *Brevicoryne brassicae*  	
- Mummy: No. of mummified aphids	

*Additional headers for rare species*  
- Wasps: No. of parasitoid wasps  
- crab_spider_l: No. of crab spiders inhabiting on the bottom of the pot (*Xysticus* spp.?)
- crab_spider_u: No. of crab spiders on the top of the plant (*Ebrechtella tricuspidata* or *Diaea subdola*?)
- sheetweb_spider_u: No. of sheetweb spiders (Linyphiidae sp.) inhabiting on the upper part of the plant  
- sheetweb_spider_l: No. of sheetweb spiders (Linyphiidae sp.) inhabiting on the lower part of the plant  
- Syrphinae_a: No. of hover fies visiting to the flowers of *A. thaliana*  
- red_mite: No. of red mites (not identified)  
- stink_bug_large: No. of large stink bugs (not identified)  
- stink_bug_smallblack: No. of small black stink bugs (not identified)  
- ant: No. of ants on the pot (not identified)  
- Mummy_black: No. of mummified aphis with black colors
- Al: No. of piggyback grasshoppers (*Atractomorpha lata*)  
- Cs: No. of seven-spotted ladybirds (*Coccinella septempunctata*)  
- leaf_hopper: No. of leaf hoppers (*Trialeurodes vaporariorum*?)
- lacewing_l: No. of larvae of lace wings  

*Indices calculated from R scripts*  
- divH: Shannon's diversity index, *H*  
- divHexp: exponential Shannon's diversity index  
- divD: Simpson's diversity index, *D*  
- richness: Species richness = Total no. of insect species  
- PxPr: Sum of Px and Pr   
- PsPa: Sum of Ps and Pa  
- chewer: No. of leaf chewing herbivores  
- sucker: No. of sap-sucking herbivores and internal feeders  
