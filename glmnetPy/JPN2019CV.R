##########################################
# cross validation using 2019 Japan data #
##########################################

# set function and load library
scaling = function(x){
  return((x - mean(x))/sd(x))
}

library(rNeighborGWAS)

#load data
geno_d = readRDS(file="./genoData/sub_snpMAF5LD80.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5LD80.rds")

pheno_d = read.csv("./insectData/Survey20194GWAS_max.csv",header=TRUE)
pheno_d = subset(pheno_d, Site=="JP")
pheno_d$gwasID = paste0("X",pheno_d$gwasID)

#reshape pheno.data
pheno_d = subset(pheno_d, gwasID!="X7329")
pheno_d = subset(pheno_d, gwasID!="X3")

n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)

geno = geno_d[,as.character(pheno_d$gwasID)]
geno = t(geno)

smap = cbind(pheno_d$position_X,pheno_d$position_Y)

# change spatial scale sqrt(2)+0.01; sqrt(8)+0.01; sqrt(18)+0.01 for s = 1; s = 2; or s = 3; respectively
scale = sqrt(2)+0.01
g_nei = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)

scale = sqrt(8)+0.01
g_nei2 = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)

X = as.matrix(model.matrix(~pheno_d$Bolting+scale(pheno_d$InitLeafLen)+factor(pheno_d$Block)+pheno_d$edge-1))

rm(geno_d)
gc();gc()

# standardize genotype values for glmnet comparisons
geno = mapply(function(x){ scaling(geno[,x]) },1:n_marker)
g_nei = mapply(function(x){ scaling(g_nei[,x]) },1:n_marker)
g_nei2 = mapply(function(x){ scaling(g_nei2[,x]) },1:n_marker)

geno = t(geno)
g_nei = t(g_nei)
g_nei2 = t(g_nei2)


######################
# set the target trait
Y = scale(log(pheno_d$richness+1)) # change trait names when needed

# beta i to all
y_obs = lm(Y~scale(InitLeafLen)+Bolting+factor(Block)+edge,data=pheno_d)$residuals
y_obs = as.numeric(y_obs)

corlist = c()

# input file must match the target trait
beta = read.csv("./output/richnessJPNself_glmnetLassoMAF5_coefALL.csv.gz",header=TRUE)
beta = beta[-c(1:19),]
r_self = c()
for(i in 1:100) {
  y_pred = beta[,i+1]%*%geno
  y_pred = lm(as.numeric(y_pred)~scale(pheno_d$InitLeafLen)+pheno_d$Bolting+factor(pheno_d$Block)+pheno_d$edge)$residuals
  y_pred = as.numeric(y_pred)
  
  r_self = c(r_self, cor(y_pred, y_obs, method="spearman"))
}
corlist = cbind(corlist,r_self)

beta = read.csv("./output/richnessS1JPN_glmnetLassoMAF5_coefALL.csv.gz",header=TRUE)
beta = beta[-c(1:19),]
nsnp = nrow(beta)/2
r_both = c()
for(i in 1:100) {
  y_pred = beta[(1:nsnp),i+1]%*%geno + beta[(nsnp+1):(nsnp*2),i+1]%*%g_nei
  y_pred = lm(as.numeric(y_pred)~scale(pheno_d$InitLeafLen)+pheno_d$Bolting+factor(pheno_d$Block)+pheno_d$edge)$residuals
  y_pred = as.numeric(y_pred)
  r_both = c(r_both, cor(y_pred, y_obs, method="spearman"))
}
corlist = cbind(corlist,r_both)

beta = read.csv("./output/richnessS2JPN_glmnetLassoMAF5_coefALL.csv.gz",header=TRUE)
beta = beta[-c(1:19),]
nsnp = nrow(beta)/2
r_both = c()
for(i in 1:100) {
  y_pred = beta[(1:nsnp),i+1]%*%geno + beta[(nsnp+1):(nsnp*2),i+1]%*%g_nei2
  y_pred = lm(as.numeric(y_pred)~scale(pheno_d$InitLeafLen)+pheno_d$Bolting+factor(pheno_d$Block)+pheno_d$edge)$residuals
  y_pred = as.numeric(y_pred)
  r_both = c(r_both, cor(y_pred, y_obs, method="spearman"))
}
corlist = cbind(corlist,r_both)

plot(corlist[,2], main="no. of richness", ylab="rho", xlab="lambda rank",ylim=c(0,0.425))
points(corlist[,3],pch=2)
points(corlist[,1],pch=20)

colnames(corlist) = c("richness_self","richness_bothS1","richness_bothS2")
write.csv(corlist,"../output/richnessJPNrhoLASSO.csv",row.names=FALSE)

# summary files
fn = system("ls ../output/*JPNrhoLASSO.csv", intern=TRUE)
rholist = c()
for(i in fn) {
  f = read.csv(i,header=TRUE)
  rholist = cbind(rholist,as.matrix(f))
  print(dim(f))
}

rholist = as.data.frame(rholist)
write.csv(rholist,"../output/SummaryJPNrhoLASSO.csv",row.names=FALSE)

fn = system("ls ../output/*JPN*DFlambda.csv", intern=TRUE)
dflist = c()
for(i in fn) {
  f = read.csv(i,header=TRUE)
  dflist = cbind(dflist,as.matrix(f))
  print(dim(f))
}

dflist = as.data.frame(dflist)
end = regexpr("_glmnet",fn)
header = paste0(rep(substr(fn,10+1,end-1),times=1, each=2), rep(c("_df","_lambda"),12))
colnames(dflist) = header
write.csv(dflist,"../output/SummaryJPNlambdaLASSO.csv",row.names=FALSE)

