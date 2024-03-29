#############################
# Figures for LASSO results #
#############################

# load plot library and functions
library(tidyverse)
library(patchwork)
source("coord.R")

#load data
geno_d = readRDS("./genoData/sub_snpMAF5LD80.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5LD80.rds")
gwasid_list = read.csv("./genoData/gwasIDlist.csv",header=TRUE)

dz = read.csv("./output/SummaryCHZrhoLASSO.csv",header=TRUE)
dj = read.csv("./output/SummaryJPNrhoLASSO.csv",header=TRUE)
which(dz$Holes_bothS1>dz$Holes_self)
which(dj$Score_bothS1>dj$Score_self)


chr_rep = table(position[,1])
cols = c(rep(rgb(1,0,0,0.5), chr_rep[1]),
         rep(rgb(0,1,0,0.5), chr_rep[2]),
         rep(rgb(0,0,1,0.5), chr_rep[3]), 
         rep(rgb(0,0,0,0.5), chr_rep[4]),
         rep(rgb(1,0,1,0.5), chr_rep[5]))

coef_path = "./output/"
f_name = "HolesS1CHZ_glmnetLassoMAF5_coefALL"
coef_d = read.csv(paste0(coef_path,f_name,".csv.gz"), header=TRUE)
coef_d = coef_d[-c(1:19),]
self_coef = coef_d[1:nrow(position),]
self_coef = self_coef[,-1]
nei_coef = coef_d[(nrow(position)+1):nrow(coef_d),]
nei_coef = nei_coef[,-1]
selfEff = apply(self_coef[,30:50],1,mean)
x = coord(position[,1],position[,2])
man1 = ggplot(NULL,aes(x=x$coord,y=selfEff)) + geom_point(colour=cols) + theme_classic() + 
  scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
  ylab(expression(hat(italic(beta))[1])) + ggtitle(substitute(paste(bold("a"))))

neiEff = apply(nei_coef[,30:50],1,mean)
man2 = ggplot(NULL,aes(x=x$coord,y=neiEff)) + geom_point(colour=cols) + theme_classic() + 
  scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
  ylab(expression(hat(italic(beta))[2]*" (J = 4)")) + ggtitle(substitute(paste(bold("b"))))

man = man1 + man2
saveRDS(man,file="../figs/HolesS1CHZ_glmnetLassoMAF5_mean.rds")

out = data.frame(position,selfEff,neiEff)
colnames(out)[4:5] = c("self_beta","nei_beta")
write.csv(out,"../output/HolesS1CHZ_glmnetLassoMAF5_mean.csv",row.names=FALSE)

biplot = readRDS("../figs/EffectSizeBiplot.rds")

# composite Figure S11
lasso_p = (man / biplot) + plot_layout(heights=c(1,1.5))
ggsave(lasso_p,filename="../figs/LASSOeffect.jpg",width=12,height=7,dpi=600)

###################
# figS: CV correlations
dz = read.csv("./output/SummaryCHZrhoLASSO.csv",header=TRUE)
lz = read.csv("./output/SummaryCHZlambdaLASSO.csv",header=TRUE)

# holes
maxy = max(c(dz$Holes_bothS1,dz$Holes_bothS2,dz$Holes_bothS3,dz$Holes_self),na.rm=TRUE)
p1 = ggplot(NULL,aes(x=log(lz$HolesS1CHZ_lambda),y=dz$Holes_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lz$HolesS2CHZ_lambda),y=dz$Holes_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lz$HolesCHZself_lambda),y=dz$Holes_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("a"),"  Zurich")),subtitle="Herbivore damage (Leaf holes)") + 
  geom_vline(xintercept = log(lz$HolesS2CHZ_lambda)[30],lty=2) +
  geom_vline(xintercept = log(lz$HolesS2CHZ_lambda)[50],lty=2) +
  geom_line(aes(x=log(lz$HolesS2CHZ_lambda)[30:50],y=dz$Holes_bothS1[30:50]),lwd=2) +
  geom_point(aes(x=log(lz$HolesS1CHZ_lambda)[39],y=dz$Holes_bothS1[39]),pch=1,size=4)

# chewer
maxy = max(c(dz$chewer_bothS1,dz$chewer_bothS2,dz$chewer_bothS3,dz$chewer_self),na.rm=TRUE)
p2 = ggplot(NULL,aes(x=log(lz$chewerS1CHZ_lambda),y=dz$chewer_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lz$HolesS2CHZ_lambda),y=dz$chewer_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lz$chewerCHZself_lambda),y=dz$chewer_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("b"))),subtitle="Individual no. of external feeders")

# sucker
maxy = max(c(dz$sucker_bothS1,dz$sucker_bothS2,dz$sucker_bothS3,dz$sucker_self),na.rm=TRUE)
p3 = ggplot(NULL,aes(x=log(lz$suckerS1CHZ_lambda),y=dz$sucker_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lz$suckerS2CHZ_lambda),y=dz$sucker_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lz$suckerCHZself_lambda),y=dz$sucker_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("c"))),subtitle="Individual no. of internal feeders")

# richness
maxy = max(c(dz$richness_bothS1,dz$richness_bothS2,dz$richness_bothS3,dz$richness_self),na.rm=TRUE)
p4 = ggplot(NULL,aes(x=log(lz$richnessS1CHZ_lambda),y=dz$richness_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lz$richnessS2CHZ_lambda),y=dz$richness_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lz$richnessCHZself_lambda),y=,dz$richness_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("d"))),subtitle="Total no. of insect species")

# figS: correlation at Japan
dj = read.csv("./output/SummaryJPNrhoLASSO.csv",header=TRUE)
lj = read.csv("./output/SummaryJPNlambdaLASSO.csv",header=TRUE)

# score
maxy = max(c(dj$Score_bothS1,dj$Score_bothS2,dj$Score_self),na.rm=TRUE)
p5 = ggplot(NULL,aes(x=log(lj$ScoreS1JPN_lambda),y=dj$Score_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lj$ScoreS2JPN_lambda),y=dj$Score_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lj$ScoreJPNself_lambda),y=dj$Score_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("e"),"  Otsu")),subtitle="Herbivore damage (Leaf area loss)")

# chewer
maxy = max(c(dj$chewer_bothS1,dj$chewer_bothS2,dj$chewer_self),na.rm=TRUE)
p6 = ggplot(NULL,aes(x=log(lj$chewerS1JPN_lambda),y=dj$chewer_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lj$chewerS2JPN_lambda),y=dj$chewer_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lj$chewerJPNself_lambda),y=dj$chewer_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("f"))),subtitle="Individual no. of external feeders")

# sucker
maxy = max(c(dj$sucker_bothS1,dj$sucker_bothS2,dj$sucker_self),na.rm=TRUE)
p7 = ggplot(NULL,aes(x=log(lj$suckerS1JPN_lambda),y=dj$sucker_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lj$suckerS2JPN_lambda),y=dj$sucker_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lj$suckerJPNself_lambda),y=dj$sucker_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("g"))),subtitle="Individual no. of internal feeders")

# richness
maxy = max(c(dj$richness_bothS1,dj$richness_bothS2,dj$richness_self),na.rm=TRUE)
p8 = ggplot(NULL,aes(x=log(lj$richnessS1JPN_lambda),y=dj$richness_bothS1)) + geom_line(lwd=1) + 
  geom_line(aes(x=log(lj$richnessS2JPN_lambda),y=,dj$richness_bothS2),colour="skyblue",lwd=0.8) +
  geom_line(aes(x=log(lj$richnessJPNself_lambda),y=dj$richness_self),colour="grey",lwd=1,lty=2) + 
  theme_classic() + ylab(expression(rho)) + xlab(expression(log(lambda))) + ylim(0,maxy) + 
  labs(title=substitute(paste(bold("h"))),subtitle="Total no. of insect species")

# Figure S10
rho_p = (p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8)
ggsave(rho_p,filename="../figs/CV2019rho.pdf",height=8,width=12)

