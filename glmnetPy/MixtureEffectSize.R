############################################
# Estimating effect size of mixed planting #
############################################

library(tidyverse)
library(patchwork)
library(vegan)

#load SNP data
geno_d = readRDS("./genoData/sub_snpMAF5LD80.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

gwasid_list = read.csv("./genoData/gwasIDlist.csv",header=TRUE)

coef_path = "./output/"
f_name = "HolesS1CHZ_glmnetLassoMAF5_mean"
coef_d = read.csv(paste0(coef_path,f_name,".csv.gz"), header=TRUE)

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
  
  total_Dam = self_Dam + nei_Dam
  # total_Dam = nei_Dam # nei eff only
  
  return(total_Dam)
}

library(arrangements)
n_geno = length(colnames(geno_d))

all_pair = permutations(colnames(geno_d),2,replace = TRUE)
poly_mat_elem = mapply(res_pair_damage, all_pair[,1], all_pair[,2])

poly_mat = matrix(poly_mat_elem, n_geno, n_geno, byrow=TRUE)
colnames(poly_mat) = colnames(geno_d)
rownames(poly_mat) = colnames(geno_d)

std_poly_mat = (poly_mat - mean(poly_mat))/sd(poly_mat)
#write.csv(std_poly_mat,file=paste0(f_name,"_stdNeiDamMat.csv"),quote=FALSE)

relativeDam = function(geno_i, geno_j) {
  genoi_relativeDam = std_poly_mat[geno_i, geno_i] - std_poly_mat[geno_i, geno_j]
  return(genoi_relativeDam)
}

relativeDam_elem = mapply(relativeDam, all_pair[,1], all_pair[,2])
relativeDamTable = data.frame(all_pair[,1], all_pair[,2], relativeDam_elem)
colnames(relativeDamTable) = c("geno_i", "geno_j", "relativeDam_i")

avg = aggregate(relativeDam_i~geno_i,relativeDamTable, mean)
temp_avg = data.frame(x=diag(std_poly_mat)[avg$geno_i],y=avg$relativeDam_i)
temp_avg = temp_avg[order(temp_avg$x),]
saveRDS(temp_avg,file="../output/temp_avg.rds",compress=TRUE)

relativeDamTable$geno_i = as.numeric(gsub("X", "", x=relativeDamTable$geno_i))
relativeDamTable$geno_j = as.numeric(gsub("X", "", x=relativeDamTable$geno_j))

# no. of positive-effect pairs / possible pairs
(sum(relativeDamTable$relativeDam_i>0)/2) / choose(199,2)

name_i = c()
for(i in relativeDamTable$geno_i) {
  if(i==1) { name_i = c(name_i,"gl1-1") 
  } else if(i==2) {
    name_i = c(name_i,"gl1-2")
  } else if(i==3) {
    name_i = c(name_i,"NA")
  } else {
    name_i = c(name_i, as.character(subset(gwasid_list, GWASid==i)$Name))
  }
}

name_j = c()
for(i in relativeDamTable$geno_j) {
  if(i==1) { name_j = c(name_j,"gl1-1") 
  } else if(i==2) {
    name_j = c(name_j,"gl1-2")
  } else if(i==3) {
    name_j = c(name_j,"NA")
  } else {
    name_j = c(name_j, as.character(subset(gwasid_list, GWASid==i)$Name))
  }
}

relativeDamTable = data.frame(name_i, name_j, relativeDamTable)

# write.csv(relativeDamTable, file=paste0(coef_path, f_name,"_relativeDam.csv"), row.names=FALSE, quote=FALSE)


# histogram of effect size
damhist = ggplot(relativeDamTable,aes(x=relativeDam_i)) + geom_histogram(binwidth=0.05) +
  theme_classic() + ylab("No. of genotype pairs") + xlab("Estimated effect size") +
  geom_vline(xintercept=0, lty=2, colour=grey(0.5,0.5)) +
  geom_text(data.frame(x=0.075,y=1000),mapping=aes(x=x,y=y),label="Bla-1 vs. Bro1-6",angle=90,size=3) +
  geom_text(data.frame(x=0.25,y=750),mapping=aes(x=x,y=y),label="Vastervik vs. Jm-0",angle=90,size=3) +
  geom_text(data.frame(x=0.8,y=500),mapping=aes(x=x,y=y),label="Bg-2 vs. Uod-1",angle=90,size=3)

saveRDS(damhist,file="../figs/EffectSizeDamage.rds",compress=TRUE,version=2)


# eff vs. kinship
K = (crossprod(geno_d)+1)/nrow(geno_d)
cor.test(K[upper.tri(K)],std_poly_mat[upper.tri(std_poly_mat)])
mantel(std_poly_mat,K,permutations=999)

# eff vs. geo
geo = read.csv("./genoData/AtGWASlocality_added.csv",header=TRUE) # locality info available via AraPheno database
rownames(std_poly_mat)
gdis = c()
for(i in 1:nrow(std_poly_mat)) {
  Xid = paste0("X",geo$GWAS_id)
  target = which(Xid == rownames(std_poly_mat)[i])
  gdis = rbind(gdis,geo[target,c("Lat","Long")])
}

gmat = as.matrix(dist(gdis,diag=TRUE, upper=TRUE))
cor.test(gmat[upper.tri(gmat)],std_poly_mat[upper.tri(std_poly_mat)])
mantel(std_poly_mat,gmat,permutations=999)
mantel.partial(std_poly_mat,gmat,K,permutations=999)
mantel.partial(std_poly_mat,K,gmat,permutations=999)

# for Figure S13c
eff_img = std_poly_mat %>% 
  as.data.frame() %>%
  rownames_to_column("f_id") %>%
  pivot_longer(-c(f_id), names_to = "samples", values_to = "pred.") %>%
  ggplot(aes(x=samples, y=f_id, fill=pred.)) + 
  geom_raster() +
  scale_fill_viridis_c() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ylab("accession i") + xlab("accession j")

# for Figure S13d
genp = ggplot(NULL,aes(x=K[upper.tri(K)],y=std_poly_mat[upper.tri(std_poly_mat)])) + geom_point(alpha=0.1) +
  theme_classic() + ylab("Esitimated effect size of mixing") + xlab("Genetic distance")

# for Figure S13e
geop = ggplot(NULL,aes(x=gmat[upper.tri(gmat)],y=std_poly_mat[upper.tri(std_poly_mat)])) + geom_point(alpha=0.1) +
  theme_classic() + ylab("Estimated effect size of mixing") + xlab("Geographical distance")

# load Figure S13f
simp = readRDS(file="../figs/SimEffsupp.rds")

temp_avg = readRDS("../output/temp_avg.rds")
avgp = ggplot(temp_avg,mapping=aes(x=x,y=y)) +
  geom_point() + theme_classic() + 
  ylab("Avg. estimated effect size") + xlab("Est. damage under monoculture") +
  geom_smooth(method=lm,se=TRUE) +
  geom_text(data.frame(x=1.5,y=-0.4),mapping=aes(x=x,y=y),label=paste0("r = ",round(cor(temp_avg$x, temp_avg$y),2),"***"),size=4)

cor.test(temp_avg$x, temp_avg$y)
cor.test(temp_avg[-c(195:199),]$x, temp_avg[-c(195:199),]$y)

saveRDS(avgp,file="../figs/avgp.rds")

# composite for Figure S13c-f
eff_biplot = (eff_img + ggtitle(substitute(paste(bold("c"))))) | (genp + ggtitle(substitute(paste(bold("d"))))) | (geop + ggtitle(substitute(paste(bold("e"))))) | (simp + ggtitle(substitute(paste(bold("f")))))
saveRDS(eff_biplot,file="../figs/EffectSizeBiplot.rds",compress=TRUE,version=2)
