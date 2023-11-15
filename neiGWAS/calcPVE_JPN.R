####################
#PVE for Japan data#
####################

library(rNeighborGWAS)
library(gaston)

#load data
geno_d = readRDS(file="./genoData/sub_snpMAF5.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5.rds")

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

# set spatial map
smap = cbind(pheno_d$position_X,pheno_d$position_Y)

# change spatial scale sqrt(2)+0.01; sqrt(8)+0.01; sqrt(18)+0.01 for s = 1 (J=4); s = 2 (J=12); or s = 3; respectively
scale = sqrt(2)+0.01 

g_nei = nei_coval(geno,smap=smap,scale=scale,grouping=pheno_d$Block,n_core=1L)

q = ncol(geno)
K_self = tcrossprod(geno)
K_self = ((q-1)/2+K_self/2)/(q-1)
K_self = as.matrix(Matrix::nearPD(K_self,maxit=10^6)$mat)

K_nei = tcrossprod(g_nei)/(q-1)
K_nei = as.matrix(Matrix::nearPD(K_nei,maxit=10^6)$mat)

X = as.matrix(model.matrix(~factor(pheno_d$Block)+scale(pheno_d$InitLeafLen)+pheno_d$Bolting+pheno_d$edge))

resList = c()
Y = scale(pheno_d$Score)
res0 = lmm.aireml(Y=Y, X=X, K=list(K_self), verbose=TRUE)
res1 = lmm.aireml(Y=Y, X=X, K=list(K_self, K_nei), verbose=TRUE)
resList = rbind(resList,
                c(0,0,1,res0$logL0))
resList = rbind(resList,
                c(res0$tau,0,res0$sigma2,res0$logL))
resList = rbind(resList,
                c(res1$tau,res1$sigma2,res1$logL))

Y = scale(log(pheno_d$PxPr+1))
res0 = lmm.aireml(Y=Y, X=X, K=list(K_self), verbose=TRUE)
res1 = lmm.aireml(Y=Y, X=X, K=list(K_self, K_nei), verbose=TRUE)
resList = rbind(resList,
                c(0,0,1,res0$logL0))
resList = rbind(resList,
                c(res0$tau,0,res0$sigma2,res0$logL))
resList = rbind(resList,
                c(res1$tau,res1$sigma2,res1$logL))

Y = scale(log(pheno_d$sucker+1))
res0 = lmm.aireml(Y=Y, X=X, K=list(K_self), verbose=TRUE)
res1 = lmm.aireml(Y=Y, X=X, K=list(K_self, K_nei), verbose=TRUE)
resList = rbind(resList,
                c(0,0,1,res0$logL0))
resList = rbind(resList,
                c(res0$tau,0,res0$sigma2,res0$logL))
resList = rbind(resList,
                c(res1$tau,res1$sigma2,res1$logL))

Y = scale(log(pheno_d$chewer+1))
res0 = lmm.aireml(Y=Y, X=X, K=list(K_self), verbose=TRUE)
res1 = lmm.aireml(Y=Y, X=X, K=list(K_self, K_nei), verbose=TRUE)
resList = rbind(resList,
                c(0,0,1,res0$logL0))
resList = rbind(resList,
                c(res0$tau,0,res0$sigma2,res0$logL))
resList = rbind(resList,
                c(res1$tau,res1$sigma2,res1$logL))

Y = scale(log(pheno_d$richness+1))
res0 = lmm.aireml(Y=Y, X=X, K=list(K_self), verbose=TRUE)
res1 = lmm.aireml(Y=Y, X=X, K=list(K_self, K_nei), verbose=TRUE)
resList = rbind(resList,
                c(0,0,1,res0$logL0))
resList = rbind(resList,
                c(res0$tau,0,res0$sigma2,res0$logL))
resList = rbind(resList,
                c(res1$tau,res1$sigma2,res1$logL))

colnames(resList) = c("sigma_1","sigma_2","sigma_e","LL")

write.csv(resList,"../output/JPNoutS1.csv",row.names=FALSE)
