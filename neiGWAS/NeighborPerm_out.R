#############################
#permutation for subset data#
#############################
source("coord.R")

# 99-times permutations on the subset data
org = readRDS("../output/neiGWAS_perm/J4/partCHZ_org_.rds") # set J4 or J4 
bonf1 = -log10(0.05/nrow(org))

perm_p = c()
for(i in 1:99) {
  res = readRDS(paste0("../output/neiGWAS_perm/J4/partCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}

perm_th1 = quantile(-log10(perm_p),0.95)

d = readRDS("../output/neiGWAS_perm/J4/partCHZ_org_.rds")
d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
chr_rep = cumsum(table(d$Chr))
cols = c(rgb(1,0,0, 2*d$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*d$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*d$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*d$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*d$MAF[(chr_rep[4]+1):(chr_rep[5])]))
x = coord(d$Chr,d$Position)
y = -log10(d$P_nei)
man = ggplot(NULL,aes(x=x,y=y)) + geom_point(colour=cols) + theme_classic() + ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) +
  ylab(expression(-log[10](p))) + xlab("Chromosomes") + geom_hline(yintercept=-log10(0.05/nrow(d)),lty=2,colour="black") + 
  geom_hline(yintercept=perm_th1,lty=2,colour="grey") + labs(subtitle="Leaf holes (subset data)")

# Figure S10
ggsave(man,filename="../figs/permJ4sub.jpg",width=6,height=3,dpi=300)


# 20-times permutation on the full data
perm_p = c()
for(i in 1:20) {
  res = readRDS(paste0("../output/neiGWAS_perm/full/partCHZ",i,".rds"))
  perm_p = c(perm_p,min(res$P_nei))
}

perm_th2 = max(-log10(perm_p))
perm_th2

