#############################
#permutation for subset data#
#############################
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

perm_th1 = quantile(-log10(perm_p),0.95)

# output examples
# > quantile(-log10(perm_p),0.95)
# 95% 
# 6.738267 
# > quantile(-log10(perm_p),0.90)
# 90% 
# 6.574153 

# Figure S6b
man1 = ggMan(f="./output/partCHZ_org.csv.gz")
man1 = man1 + geom_hline(yintercept=perm_th1,lty=2,lwd=1.25,colour="grey") + labs(subtitle=substitute(paste(bold("c"),"  Subset data on herbivore damage (Leaf holes)")))

########################
# 20-times permutation on the full data
perm_p = c()
for(i in 1:20) {
  res = readRDS(paste0("../output/neiGWAS_perm/full/partCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}
# saveRDS(perm_p,"./output/fullCHZ_perm_p.rds")

perm_th2 = max(-log10(perm_p))
perm_th2

# Figure S6a
man2 = ggMan(f="./output/CHZneiGWAS_HolesS1.csv.gz")
man2 = man2 + geom_hline(yintercept=perm_th2,lty=2,lwd=1.25,colour="grey") + labs(subtitle=substitute(paste(bold("b"),"  Full data on herbivore damage (Leaf holes)")))

###############
# PVE permutations
perm_pve = c()
for(i in 1:99) {
  res = readRDS(paste0("../output/neiGWAS_perm/pve/fullCHZpve",i,".rds"))
  netPVE = (res$PVEself[2]+res$PVEnei[2]) - res$PVEself[1]
  perm_pve = c(perm_pve,netPVE)
}

# saveRDS(perm_pve,"./output/fullCHZ_perm_pve.rds")

# Figure S6c
pve_hist = ggplot(NULL,aes(x=perm_pve)) + geom_histogram() + theme_classic() + 
  geom_vline(xintercept=(0.51330843-0.452626269)) + # PVE for holes at J = 4
  geom_vline(xintercept = quantile(perm_pve,0.95),lty=2,lwd=1.25,colour="grey") +
  xlab("PVE by neighbor genotype effects") + ylab("No. of iterations") +
  labs(subtitle=substitute(paste(bold("a"),"  Full data on herbivore damage (Leaf holes)")))

man = pve_hist | man2 | man1
ggsave(man,filename="../figs/permJ4.jpg",width=12,height=3,dpi=300)

