################################
#analysis of permutation output#
################################
library(tidyverse)
library(patchwork)
source("coord.R")

ggMan = function(f,type="nei") {
  d = read.csv(f,header=TRUE)
  d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
  chr_rep = table(d$Chr)
  cols = c(rep(rgb(1,0,0,0.5), chr_rep[1]),
           rep(rgb(0,1,0,0.5), chr_rep[2]),
           rep(rgb(0,0,1,0.5), chr_rep[3]), 
           rep(rgb(0,0,0,0.5), chr_rep[4]),
           rep(rgb(1,0,1,0.5), chr_rep[5]))
  x = coord(d$Chr,d$Position)
  if(type=="nei") {
    y = -log10(d$P_nei)
  } else {
    y = -log10(d$P_self)
  }
  man = ggplot(NULL,aes(x=x$coord,y=y)) + geom_point(colour=cols) + theme_classic() + scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
    ylab(expression(-log[10]*(italic(p)))) + geom_hline(yintercept=-log10(0.05/nrow(d)),lty=2,colour="black")
  return(man)
}

###################
# 99-times permutations on the subset data
org = readRDS("../output/neiGWAS_perm/J4/partCHZ_org_.rds") # set J4 or J12 
bonf1 = -log10(0.05/nrow(org))

perm_p = c()
for(i in 1:99) {
  res = readRDS(paste0("../output/neiGWAS_perm/J4/partCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}
# saveRDS(perm_p,"./output/partCHZ_perm_p.rds")

perm_th1p = quantile(-log10(perm_p),0.95)

# output examples
# > quantile(-log10(perm_p),0.95)
# 95% 
# 6.738267 
# > quantile(-log10(perm_p),0.90)
# 90% 
# 6.574153 

perm_p = c()
for(i in 1:99) {
  res = readRDS(paste0("../output/neiGWAS_rot/J4/partCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}
# saveRDS(perm_p,"./output/partCHZ_rot_p.rds")

perm_th1r = quantile(-log10(perm_p),0.95)

# output examples
# > quantile(-log10(perm_p),0.95)
# 95% 
# 6.802742 
# > quantile(-log10(perm_p),0.90)
# 90% 
# 6.689917 

# Figure S6b
man1a = ggMan(f="./output/partCHZ_org.csv.gz")
man1a = man1a + geom_hline(yintercept=perm_th1p,lty=2,lwd=1.25,colour="blue") + 
  labs(subtitle=substitute(paste(bold("c"),"  Subset data on herbivore damage (Leaf holes)")))

man1b = ggMan(f="./output/partCHZ_org.csv.gz")
man1b = man1b + geom_hline(yintercept=perm_th1r,lty=2,lwd=1.25,colour="red") +
  labs(subtitle=substitute(paste(bold("f"),"  Subset data on herbivore damage (Leaf holes)")))


########################
# 20-times permutation on the full data
perm_p = c()
for(i in 1:20) {
  res = readRDS(paste0("../output/neiGWAS_perm/full/partCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}
# saveRDS(perm_p,"./output/fullCHZ_perm_p.rds")

perm_th2p = max(-log10(perm_p))

perm_p = c()
for(i in 1:20) {
  res = readRDS(paste0("../output/neiGWAS_rot/J4/fullCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}
# saveRDS(perm_p,"./output/fullCHZ_rot_p.rds")

perm_th2r = max(-log10(perm_p))


# Figure S6a
man2a = ggMan(f="./output/CHZneiGWAS_HolesS1.csv.gz")
man2a = man2a + geom_hline(yintercept=perm_th2p,lty=2,lwd=1.25,colour="blue") + 
  labs(subtitle=substitute(paste(bold("b"),"  Full data on herbivore damage (Leaf holes)")))

man2b = ggMan(f="./output/CHZneiGWAS_HolesS1.csv.gz")
man2b = man2b + geom_hline(yintercept=perm_th2r,lty=2,lwd=1.25,colour="red") +
  labs(subtitle=substitute(paste(bold("e"),"  Full data on herbivore damage (Leaf holes)")))


###############
# PVE permutations
perm_pve1 = c()
for(i in 1:99) {
  res = readRDS(paste0("../output/neiGWAS_perm/pve/fullCHZpve",i,".rds"))
  netPVE = (res$PVEself[2]+res$PVEnei[2]) - res$PVEself[1]
  perm_pve1 = c(perm_pve1,netPVE)
}
# saveRDS(perm_pve1,"./output/fullCHZ_perm_pve.rds")

perm_pve2 = c()
for(i in 1:99) {
  res = readRDS(paste0("../output/neiGWAS_rot/pve/fullCHZpve",i,".rds"))
  netPVE = (res$PVEself[2]+res$PVEnei[2]) - res$PVEself[1]
  perm_pve2 = c(perm_pve2,netPVE)
}
# saveRDS(perm_pve2,"./output/fullCHZ_rot_pve.rds")

# Figure S6c
pve_hist1 = ggplot(NULL,aes(x=perm_pve1)) + theme_classic() +
  geom_histogram(fill=rgb(0,0,1,0.25)) +  
  geom_vline(xintercept=(0.51330843-0.452626269),lwd=1.25,colour=grey(0,0.75)) + # PVE for holes at J = 4
  geom_vline(xintercept = quantile(perm_pve1,0.95),lty=2,lwd=1.25,colour="blue") +
  xlab("PVE by neighbor genotype effects") + ylab("No. of iterations") +
  labs(subtitle=substitute(paste(bold("a"),"  Full data on herbivore damage (Leaf holes)")))

pve_hist2 = ggplot(NULL,aes(x=perm_pve2)) + theme_classic() +
  geom_histogram(fill=rgb(1,0,0,0.25)) +  
  geom_vline(xintercept=(0.51330843-0.452626269),lwd=1.25,colour=grey(0,0.75)) + # PVE for holes at J = 4
  geom_vline(xintercept = quantile(perm_pve2,0.95),lty=2,lwd=1.25,colour="red") +
  xlab("PVE by neighbor genotype effects") + ylab("No. of iterations") +
  labs(subtitle=substitute(paste(bold("d"),"  Full data on herbivore damage (Leaf holes)")))


man = (pve_hist1 | man2a | man1a) / (pve_hist2 | man2b | man1b)
ggsave(man,filename="../figs/permJ4v2.jpg",width=12,height=6,dpi=300)

