################################
#export Japan data for glmnetPy#
################################

library(rNeighborGWAS)

#load data
geno_d = readRDS(file="./genoData/sub_snpMAF5LD80.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5LD80.rds")

pheno_d = read.csv("./insectData/SurveyCombined4GWAS_max.csv",header=TRUE)
pheno_d = subset(pheno_d, Site=="JP")
pheno_d$gwasID = paste0("X",pheno_d$gwasID)

#reshape pheno.data
pheno_d = subset(pheno_d, gwasID!="X7329")
pheno_d = subset(pheno_d, gwasID!="X3")
naID = which(is.na(pheno_d$InitLeafLen))
pheno_d = pheno_d[-naID,]

n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)

geno = geno_d[,as.character(pheno_d$gwasID)]
geno = t(geno)

smap = cbind(pheno_d$position_X,pheno_d$position_Y)

# change spatial scale sqrt(2)+0.01; sqrt(8)+0.01; sqrt(18)+0.01 for s = 1; s = 2; or s = 3; respectively
scale = sqrt(2)+0.01

rm(geno_d)
gc();gc()

g_nei = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)
X = as.matrix(model.matrix(~factor(Block)+scale(InitLeafLen)+Bolting+edge-1,data=pheno_d))

# export for python
library(reticulate)
geno = t(geno)
g_nei = t(g_nei)

cov = stats::model.matrix(~factor(Block)+scale(InitLeafLen)+Bolting+edge-1, data=pheno_d)
cov = r_to_py(cov)

self = r_to_py(geno)
nei = r_to_py(g_nei)
pheno = r_to_py(pheno_d)

# Note: Do not run through Rstudio. It may cause an error; Outputs are so large as several GB
py_save_object(self, filename="../data/glmnetPyJPNselfMAF5.pkl", pickle="pickle")
py_save_object(nei, filename="../data/glmnetPyJPNneiS1MAF5.pkl", pickle="pickle")
py_save_object(pheno, filename="../data/glmnetPyJPNphenoMAF5map.pkl", pickle="pickle")
py_save_object(cov, filename="../data/glmnetPyJPNcovMAF5map.pkl", pickle="pickle")
