############################################
# Simulating effect size of mixed planting #
############################################

library(tidyverse)
library(arrangements)

#load SNP data
geno_d = readRDS("./genoData/sub_snpMAF5LD80.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

gwasid_list = read.csv("./genoData/gwasIDlist.csv",header=TRUE)

coef_path = "./output/"
f_name = "HolesS1CHZ_glmnetLassoMAF5_mean"
coef_d = read.csv(paste0(coef_path,f_name,".csv.gz"),header=TRUE)

self_which = which(coef_d$self_beta!=0)
nei_which = which(coef_d$nei_beta!=0)

self_beta = coef_d$self_beta[self_which]
nei_beta = coef_d$nei_beta[nei_which]

self_d = geno_d[self_which,]
nei_d = geno_d[nei_which,]


#define function 
res_pair_damage = function(selfID, neiID) {
  
  self_Dam = sum(self_beta*self_d[,selfID])
  
  line12 = nei_d[,selfID]*nei_d[,neiID]
  
  nei_Dam = sum(nei_beta*line12)
  
  # total_Dam = self_Dam + nei_Dam
  total_Dam = nei_Dam # nei eff only
  
  return(total_Dam)
}

n_geno = length(colnames(geno_d))

all_pair = permutations(colnames(geno_d),2,replace=TRUE)
poly_mat_elem = mapply(res_pair_damage, all_pair[,1], all_pair[,2])

poly_mat = matrix(poly_mat_elem, n_geno, n_geno, byrow=TRUE)
colnames(poly_mat) = colnames(geno_d)
rownames(poly_mat) = colnames(geno_d)

std_poly_mat = (poly_mat - mean(poly_mat))/sd(poly_mat)

# all genotypes
set.seed(1234)
eff_sim = c()
for(j in 1:9999) {
  n_geno = 8
  rand = sample(colnames(geno_d),n_geno)
  
  # 1 gennotype (= monoculture)
  e1 = mean(diag(std_poly_mat[rand,rand]))
  
  # n = 2
  sub_n = 2
  comb = combn(rand,sub_n)
  
  eff = c()
  for(i in 1:ncol(comb)) {
    eff = c(eff,mean(std_poly_mat[comb[,i],comb[,i]]))
  }
  e2 = mean(eff)
  
  # n = 4
  sub_n = 4
  comb = combn(rand,sub_n)
  
  eff = c()
  for(i in 1:ncol(comb)) {
    eff = c(eff,mean(std_poly_mat[comb[,i],comb[,i]]))
  }
  e4 = mean(eff)
  
  # 8 genotypes
  e8 = mean(std_poly_mat[rand,rand])
  
  eff_vec = c(e1,e2,e4,e8)
  eff_sim = rbind(eff_sim,eff_vec)
}

m = apply(eff_sim,2,mean)
er = apply(eff_sim,2,sd)

df1 = data.frame(x=c(1,2,4,8), y=m)
b <- ggplot(data=df1, aes(x=x, y=m)) + geom_line() + geom_point() +
  ylab("Simulated damage") + xlab("No. of genotypes") + xlim(1,8) +
  geom_errorbar(aes(ymax=m+er, ymin=m-er),width=0.0) + theme_test()

y1 = std_poly_mat[
  paste0("X",gwasid_list[gwasid_list$Name=="Bg-2",]$GWASid),
  paste0("X",gwasid_list[gwasid_list$Name=="Uod-1",]$GWASid)
]

y2 = std_poly_mat[
  paste0("X",gwasid_list[gwasid_list$Name=="Vastervik",]$GWASid),
  paste0("X",gwasid_list[gwasid_list$Name=="Jm-0",]$GWASid)
]

y3 = std_poly_mat[
  paste0("X",gwasid_list[gwasid_list$Name=="Bro1-6",]$GWASid),
  paste0("X",gwasid_list[gwasid_list$Name=="Bla-1",]$GWASid)
]

b = b + geom_text(data.frame(x=2.2,y=y3+0.1),mapping=aes(x=x,y=y),label="Bla-1 vs. Bro1-6",angle=0,size=3,hjust=0) +
  geom_point(data.frame(x=2,y=y3),mapping=aes(x=x,y=y),pch=1) +
  geom_text(data.frame(x=2.2,y=y2-0.1),mapping=aes(x=x,y=y),label="Vastervik vs. Jm-0",angle=0,size=3,hjust=0) +
  geom_point(data.frame(x=2,y=y2),mapping=aes(x=x,y=y),pch=1) +
  geom_text(data.frame(x=2.2,y=y1+0.1),mapping=aes(x=x,y=y),label="Bg-2 vs. Uod-1",angle=0,size=3,hjust=0) +
  geom_point(data.frame(x=2,y=y1),mapping=aes(x=x,y=y),pch=1)

# export for main Figure 4c
saveRDS(b,file="../figs/SimEffmain.rds",compress=TRUE,version=2)

#####################
# positive pairs only
relativeDam = read.csv(file=paste0(coef_path, f_name,"_relativeDam_nstd.csv"))
posIDs = subset(relativeDam,relativeDam_i>0)
posIDs$geno_i = paste0("X",posIDs$geno_i)
posIDs$geno_j = paste0("X",posIDs$geno_j)

allSet = intersect(unique(posIDs$geno_j),unique(posIDs$geno_j))

eff_sim = c()
for(j in 1:9999) {
  
  rand = c()
  while(length(rand)<8) {
    row_i = sample.int(nrow(posIDs),1)
    candidate = c(posIDs[row_i,3],posIDs[row_i,4])
    and = intersect(rand, candidate)
    if(length(and)==0) rand = c(rand,candidate)
  }
  
  # n_geno = 8
  # rand = sample(allSet,n_geno)
  
  # 1 gennotype (= monoculture)
  e1 = mean(diag(std_poly_mat[rand,rand]))
  
  # n = 2
  comb = combn(c(1,3,5,7),1)
  
  eff = c()
  for(i in comb) {
    eff = c(eff,mean(std_poly_mat[rand[i],rand[i+1]]))
  }
  e2 = mean(eff)
  
  # n = 4
  sub_n = 2
  comb = combn(c(1,3,5,7),sub_n)
  
  eff = c()
  for(i in 1:ncol(comb)) {
    IDs = c(rand[comb[,i]],rand[comb[,i]+1])
    eff = c(eff,mean(std_poly_mat[IDs,IDs]))
  }
  e4 = mean(eff)
  
  # 8 genotypes
  e8 = mean(std_poly_mat[rand,rand])
  
  eff_vec = c(e1,e2,e4,e8)
  eff_sim = rbind(eff_sim,eff_vec)
}

m = apply(eff_sim,2,mean)
er = apply(eff_sim,2,sd)

df2 = data.frame(x=c(1,2,4,8), y=m)
b = ggplot(data=df2, aes(x=x, y=m)) + geom_line() + geom_point() +
  ylab("Simulated damage") + xlab("No. of genotypes") + xlim(1,8) +
  geom_errorbar(aes(ymax=m+er, ymin=m-er),width=0.0) + theme_classic()

y1 = std_poly_mat[
  paste0("X",gwasid_list[gwasid_list$Name=="Bg-2",]$GWASid),
  paste0("X",gwasid_list[gwasid_list$Name=="Uod-1",]$GWASid)
]

y2 = std_poly_mat[
  paste0("X",gwasid_list[gwasid_list$Name=="Vastervik",]$GWASid),
  paste0("X",gwasid_list[gwasid_list$Name=="Jm-0",]$GWASid)
]

y3 = std_poly_mat[
  paste0("X",gwasid_list[gwasid_list$Name=="Bro1-6",]$GWASid),
  paste0("X",gwasid_list[gwasid_list$Name=="Bla-1",]$GWASid)
]

b = b + geom_text(data.frame(x=4,y=y3),mapping=aes(x=x,y=y),label="<- Bla-1 vs. Bro1-6",angle=0,size=3) +
  geom_point(data.frame(x=2,y=y3),mapping=aes(x=x,y=y),pch=1) +
  geom_text(data.frame(x=4,y=y2),mapping=aes(x=x,y=y),label="<- Vastervik vs. Jm-0",angle=0,size=3) +
  geom_point(data.frame(x=2,y=y2),mapping=aes(x=x,y=y),pch=1) +
  geom_text(data.frame(x=4,y=y1),mapping=aes(x=x,y=y),label="<- Bg-2 vs. Uod-1",angle=0,size=3) +
  geom_point(data.frame(x=2,y=y1),mapping=aes(x=x,y=y),pch=1)

# export for Figure S11f
saveRDS(b,file="../figs/SimEffsupp.rds",compress=TRUE,version=2)

