#############################
#permutation for subset data#
#############################

i = commandArgs(trailingOnly=TRUE)
print(i)

library(rNeighborGWAS)

#load data
geno_d = readRDS(file="../data/sub_snpMAF5.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="../data/positionsMAF5.rds")

pheno_d = read.csv("./insectData/SurveyCombined4GWAS_max.csv",header=TRUE)
pheno_d = subset(pheno_d, Site=="ZH")
pheno_d = pheno_d[1:400,] # comment out this line to analyze full data
smap = cbind(pheno_d$position_X,pheno_d$position_Y)
perm_id = c(); for(j in 0:(nrow(pheno_d)/200 - 1)) perm_id = c(perm_id,sample(1:200,200)+j*200) # comment out this line to analyze w/o permutations
pheno_d = pheno_d[perm_id,] # comment out this line to analyze w/o permutations

smap = smap[-which(is.na(pheno_d$gwasID)),]
pheno_d = pheno_d[-which(is.na(pheno_d$gwasID)),]
pheno_d$gwasID = paste0("X",pheno_d$gwasID)

#reshape pheno.data
smap = smap[pheno_d$gwasID!="X7329",]
pheno_d = subset(pheno_d, gwasID!="X7329")

smap = smap[pheno_d$gwasID!="X3",]
pheno_d = subset(pheno_d, gwasID!="X3")

naID = which(is.na(pheno_d$InitLeafLen))
smap = smap[-naID,]
pheno_d = pheno_d[-naID,]

n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)

geno = geno_d[,as.character(pheno_d$gwasID)]
geno = t(geno)

# change spatial scale sqrt(2)+0.01; sqrt(8)+0.01; sqrt(18)+0.01 for s = 1(J=4); s = 2(J=12); or s = 3; respectively
scale = sqrt(8)+0.01

g_nei = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)
X = as.matrix(model.matrix(~factor(Block)+scale(InitLeafLen)+Bolting+edge-1,data=pheno_d))

res = nei_lmm(geno=geno,g_nei=g_nei,pheno=scale(log(pheno_d$Holes+1)),addcovar=X,n_core=10L)
res = data.frame(position,res)
colnames(res) = c("Chr","Position","MAF","beta_self","beta_nei","P_self","P_nei")
saveRDS(res,file=paste0("../output/neiGWAS_perm/partCHZ",i,".rds"),version=2)



