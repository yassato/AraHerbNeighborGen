###################################
#permutation with genome shuffling# 
###################################

i = commandArgs(trailingOnly=TRUE)
print(i)

library(rNeighborGWAS)

#load data
geno_d = readRDS(file="./genoData/sub_snpMAF5.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5.rds")

rotate_i = function(i) {
  perm_chr = sample(1:5)
  perm_pos = c(which(position[,1]==perm_chr[1]),
               which(position[,1]==perm_chr[2]),
               which(position[,1]==perm_chr[3]),
               which(position[,1]==perm_chr[4]),
               which(position[,1]==perm_chr[5]))
  chr_rotate_i = geno_d[perm_pos,i]
  break_point = sample(length(perm_pos),1)
  chr_rotate_i = c(chr_rotate_i[(break_point+1):length(chr_rotate_i)], 
                   chr_rotate_i[1:(break_point)])
  return(chr_rotate_i)
}

geno_p = mapply(rotate_i,1:199)
colnames(geno_p) = colnames(geno_d)
geno_d = geno_p

rm(geno_p)
gc();gc()

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

rm(geno_d)
gc();gc()

# change spatial scale sqrt(2)+0.01; sqrt(8)+0.01; sqrt(18)+0.01 for s = 1 (J=4); s = 2 (J=12); or s = 3; respectively
scale = sqrt(2)+0.01
X = as.matrix(model.matrix(~factor(Block)+scale(InitLeafLen)+Bolting+edge-1,data=pheno_d))

res = calc_PVEnei(pheno=scale(log(pheno_d$Holes+1)),geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block,addcovar=X,n_core=8L)
saveRDS(res,file=paste0("./output/neiGWAS_rot/partCHZpve",i,".rds"),version=2)

g_nei = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)
res = nei_lmm(geno=geno,g_nei=g_nei,pheno=scale(log(pheno_d$Holes+1)),addcovar=X,n_core=8L)
res = data.frame(position,res)
colnames(res) = c("Chr","Position","MAF","beta_self","beta_nei","P_self","P_nei")
saveRDS(res,file=paste0("./output/neiGWAS_rot/partCHZ",i,".rds"),version=2)
