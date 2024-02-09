####################
# GWAS simulations #
####################

library(tidyverse)
library(patchwork)
library(rNeighborGWAS)

set.seed(1234)
pves = c(); slopes = c()
for(i in 1:30) {
  g_bin = readRDS("./simulation/Ara250kRegMap_chr12_MAF10.rds")
  perm = sample(colnames(g_bin),9*36,replace=FALSE)
  g_self = g_bin[,perm]
  g_self[g_self==0] = -1
  smap = data.frame(X=rep(1:9,each=36), Y=rep(1:36,9))
  g_self = t(g_self)
  
  pheno = nei_simu(geno=g_self,smap=smap,scale=sqrt(2),n_causal=50,pveB=0.6,pve=0.8)
  pheno = pheno[,1]
  
  scales = c(sqrt(2),sqrt(8),sqrt(18),sqrt(32),sqrt(50),sqrt(72),sqrt(128),sqrt(162))
  
  pve = calc_PVEnei(pheno=pheno,geno=g_self,smap=smap,scale=scales,
                    response="quantitative",n_core=8L)
  
  slope = c()
  for (j in scales) {
    g_nei = nei_coval(geno=g_self,smap=smap,scale=j,n_core=8L)
    res = nei_lmm(geno=g_self,g_nei=g_nei,pheno=pheno,response="quantitative",n_core=8L)
    
    p = res$p_nei
    x = -log(ppoints(length(p)),10)
    y = -log(sort(p,decreasing=FALSE),10)
    plot(x,y)
    abline(a=0,b=1)
    abline(a=coef(lm(y~x))[1],b=coef(lm(y~x))[2],lty=2)
    
    slope = c(slope, coef(lm(y~x))[2])
  }
  slopes = rbind(slopes,slope)
  pves = append(pves,list(pve))
}

saveRDS(slopes,file="./simulation/slopes.rds")
saveRDS(pves,file="./simulation/pve.rds")

slopes = readRDS("./simulation/slopes.rds")
pves = readRDS("./simulation/pve.rds")

scales = c(sqrt(2),sqrt(8),sqrt(18),sqrt(32),sqrt(50),sqrt(72),sqrt(128),sqrt(162))

# Figure S15
pdf(file="../figs/slopePVEsimu.pdf",width=8,height=4)
par(mfcol=(c(1,2)))
plot(round(scales,2),slopes[1,],type="b",ylim=c(0.5,1.5),las=1,
     ylab="slope of QQ-plot",xlab="spatial distance from a focal plant",main="a                                           ")
for(i in 2:30) { points(round(scales,2),slopes[i,],type="b") }
abline(a=1,b=0,lwd=2,col="blue",lty=2)
abline(v=sqrt(2),lwd=2,col="red",lty=2)

net_pves = c()
for(i in 1:30) {
  net_pves = rbind(net_pves,((pves[[i]]$PVEself + pves[[i]]$PVEnei) - pves[[i]]$PVEself[1]))
}

plot(c(0,round(scales,2)),net_pves[1,],type="b",ylim=c(0,0.3),las=1,
     ylab="PVE - h2",xlab="spatial distance from a focal plant",main="b                                           ")
for(i in 2:30) { points(c(0,round(scales,2)),net_pves[i,],type="b") }
abline(a=1,b=0,lwd=2,col="blue",lty=2)
abline(v=sqrt(2),lwd=2,col="red",lty=2)

dev.off()


# MAF and correl, r; to be done
g_bin = readRDS("./simulation/Ara250kRegMap_chr12_MAF10.rds")
perm = sample(colnames(g_bin),9*36,replace=FALSE)
g_self = g_bin[,perm]
af = apply(g_self,1,sum)/dim(g_self)[2]
g_self[g_self==0] = -1
smap = data.frame(X=rep(1:9,each=36), Y=rep(1:36,9))
g_self = t(g_self)
g_nei = nei_coval(geno=g_self,smap=smap,scale=sqrt(2),n_core=8L)

scales = c(sqrt(2),sqrt(8),sqrt(18),sqrt(32),sqrt(50),sqrt(72),sqrt(128),sqrt(162))

r_mat = c()
for(j in scales) {
  g_nei = nei_coval(geno=g_self,smap=smap,scale=j,n_core=8L)
  
  r_vec = c()
  for(i in 1:dim(g_self)[2]) {
    r_vec = c(r_vec, cor(g_self[,i],g_nei[,i])) 
  }
  r_mat = cbind(r_mat,r_vec)
}

colnames(r_mat) = scales

pr1 = ggplot(NULL,aes(x=af,y=r_mat[,1])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[1]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr2 = ggplot(NULL,aes(x=af,y=r_mat[,2])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[2]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr3 = ggplot(NULL,aes(x=af,y=r_mat[,3])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[3]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr4 = ggplot(NULL,aes(x=af,y=r_mat[,4])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[4]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr5 = ggplot(NULL,aes(x=af,y=r_mat[,5])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[5]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr6 = ggplot(NULL,aes(x=af,y=r_mat[,6])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[6]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr7 = ggplot(NULL,aes(x=af,y=r_mat[,7])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[7]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

pr8 = ggplot(NULL,aes(x=af,y=r_mat[,8])) + geom_point(alpha=0.1) + theme_classic() + ylim(-1,1) + 
  labs(subtitle=paste0("spatial distance = ",round(as.numeric(colnames(r_mat)[8]),2))) + 
  xlab("allele frequency") + ylab("Pearson's correlation coef.") + geom_smooth(method="gam")

# Figure S14
pr_all = (pr1 | pr2 | pr3 | pr4) / (pr5 | pr6 | pr7 | pr8)
ggsave(pr_all,filename="../figs/cor_af.jpg",width=10,height=5,dpi=300)

# test w/o K2
coef = c()
for(i in 1:30) {
  g_bin = readRDS("./simulation/Ara250kRegMap_chr12_MAF10.rds")
  perm = sample(colnames(g_bin),9*36,replace=FALSE)
  g_self = g_bin[,perm]
  g_self[g_self==0] = -1
  smap = data.frame(X=rep(1:9,each=36), Y=rep(1:36,9))
  
  g_self = t(g_self)
  
  pheno = nei_simu(geno=g_self,smap=smap,scale=sqrt(2),n_causal=50,pveB=0.6,pve=0.8)
  pheno = pheno[,1]
  
  g_nei = nei_coval(geno=g_self,smap=smap,scale=sqrt(2),n_core=8L)
  q = ncol(g_self)
  K_self = tcrossprod(g_self)
  K_self = ((q - 1)/2 + K_self/2)/(q - 1)
  X = matrix(1, nrow = length(pheno))
  eiKs = eigen(K_self)
  
  test_i = function(i) {
    X0 = cbind(X, g_self[,i])
    X1 = cbind(X, g_self[,i], g_nei[,i])
    
    res00 = gaston::lmm.diago(Y = pheno, X = X0, eigenK = eiKs, verbose = FALSE)
    LL00 = gaston::lmm.diago.profile.likelihood(tau = res00$tau, s2 = res00$sigma2, Y = pheno, X = X0, eigenK = eiKs)[1,1]
    
    res01 = gaston::lmm.diago(Y = pheno, X = X1, eigenK = eiKs, verbose = FALSE)
    LL01 = gaston::lmm.diago.profile.likelihood(tau = res01$tau, s2 = res01$sigma2, Y = pheno, X = X1, eigenK = eiKs)[1,1]
    
    return(pchisq(2*(LL01 - LL00),1,lower.tail = FALSE))
  }
  
  # w/o K2
  p1 = parallel::mcmapply(test_i, 1:dim(g_nei)[2],mc.cores=8L)
  p = p1
  x = -log(ppoints(length(p)),10)
  y = -log(sort(p,decreasing=FALSE),10)
  slope1 = coef(lm(y~x))[2]
  
  # w/ K2
  res = nei_lmm(geno=g_self,g_nei=g_nei,pheno=pheno,response="quantitative",n_core=8L)
  p2 = res$p_nei
  
  p = p2
  x = -log(ppoints(length(p)),10)
  y = -log(sort(p,decreasing=FALSE),10)
  slope2 = coef(lm(y~x))[2]
  
  coef = rbind(coef,c(slope1,slope2))
}

saveRDS(coef,file="./simulation/pK.rds")

## copy & paste from source data of Fig. S9 to draw the same figure
# d = read.table(pipe("pbpaste"),header=TRUE)
# p1 = 10^(-d$y1); p2 = 10^(-d$y2)

x1 = -log(ppoints(length(p1)),10)
y1 = -log(sort(p1,decreasing=FALSE),10)
slope1 = coef(lm(y1~x1))
x2 = -log(ppoints(length(p2)),10)
y2 = -log(sort(p2,decreasing=FALSE),10)
slope2 = coef(lm(y2~x2))
qqKsim = ggplot(NULL,aes(x=x1,y=y1)) + geom_point(alpha=0.2,color="red") + geom_abline(slope=slope1[2],intercept=slope1[1],color="red") +
  geom_point(aes(x=x2,y=y2),alpha=0.2,color="blue") + geom_abline(slope=slope2[2],intercept=slope2[1],color="blue") +
  geom_abline(slope=1,intercept=0,lty=2) + labs(subtitle="Simulated data") +
  theme_classic() + xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p))) +
  geom_label(data.frame(x=0.25,y=6),mapping=aes(x=x,y=y),label=paste0("w/o K2: slope = ",round(slope1[2],3)),size=3,color=rgb(1,0.0,0,1),hjust=0) +
  geom_label(data.frame(x=0.25,y=5.5),mapping=aes(x=x,y=y),label=paste0("w/ K2: slope = ",round(slope2[2],3)),size=3,color=rgb(0,0,1,1),hjust=0)

coef = readRDS("./simulation/pK.rds")
colnames(coef) = c("K1","K1_and_K2")
coef = as.data.frame(coef)
d = gather(coef,key="K",value="slope")

bp = ggplot(d,aes(x=K,y=slope)) + geom_boxplot() + theme_classic() + labs(subtitle="Simulation") +
  ylab("slope of QQ-plot") + xlab("K") + geom_jitter() + geom_hline(yintercept=1,lty=2)

dK2 = read.csv("./output/CHZneigGWAS_HolesS1_wo_K2.csv.gz",header=TRUE)
p = dK2$p1
x = -log(ppoints(length(p)),10)
y = -log(sort(p,decreasing=FALSE),10)
qqK2 = ggplot(NULL,aes(x=x,y=y)) + geom_point(alpha=0.2) + geom_abline(slope=1,intercept=0,lty=2) + labs(subtitle="Herbivore damage (Leaf holes)") +
  theme_classic() + xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

# Figure S16
bpK2 = (bp | qqKsim | qqK2) + plot_annotation(tag_level="a")
ggsave(bpK2,filename="../figs/bp_qqK2.jpg",width=12,height=4,dpi=300)

