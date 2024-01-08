######################
#GWAS for Zurich data#
######################

library(rNeighborGWAS)

#load data
geno_d = readRDS(file="./genoData/sub_snpMAF5.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5.rds")

pheno_d = read.csv("./insectData/SurveyCombined4GWAS_max.csv",header=TRUE)
pheno_d = subset(pheno_d, Site=="ZH")
pheno_d = pheno_d[-which(is.na(pheno_d$gwasID)),]
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

# change spatial scale sqrt(2)+0.01; sqrt(8)+0.01; sqrt(18)+0.01 for s = 1(J=4); s = 2(J=12); or s = 3; respectively
scale = sqrt(2)+0.01

rm(geno_d)
gc();gc()

g_nei = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)
X = as.matrix(model.matrix(~factor(Block)+scale(InitLeafLen)+Bolting+edge-1,data=pheno_d))
pheno = scale(log(pheno_d$Holes+1))

q = n_marker
K_self = tcrossprod(geno)
K_self = ((q - 1)/2 + K_self/2)/(q - 1)
eiKs = eigen(K_self)

test_i = function(i) {
  X0 = cbind(X, geno[,i])
  X1 = cbind(X, geno[,i], g_nei[,i])
  
  res00 = gaston::lmm.diago(Y = pheno, X = X0, eigenK = eiKs, verbose = FALSE)
  LL00 = gaston::lmm.diago.profile.likelihood(tau = res00$tau, s2 = res00$sigma2, Y = pheno, X = X0, eigenK = eiKs)[1,1]
  
  res01 = gaston::lmm.diago(Y = pheno, X = X1, eigenK = eiKs, verbose = FALSE)
  LL01 = gaston::lmm.diago.profile.likelihood(tau = res01$tau, s2 = res01$sigma2, Y = pheno, X = X1, eigenK = eiKs)[1,1]
  
  return(pchisq(2*(LL01 - LL00),1,lower.tail = FALSE))
}

p1 = parallel::mcmapply(test_i, 1:dim(geno)[2],mc.cores=24L)
res = data.frame(position,p1)
write.csv(res,"../output/CHZneigGWAS_HolesS1_wo_K2.csv")

