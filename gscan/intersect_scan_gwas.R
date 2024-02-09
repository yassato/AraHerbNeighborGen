#########################################################
# summarize and visualize the results of selection scan #
#########################################################

library(tidyverse)
library(patchwork)

# load and merge scan results
chr1_hh = read.csv("./output/scan_ihs200_chr1.csv")
chr2_hh = read.csv("./output/scan_ihs200_chr2.csv")
chr3_hh = read.csv("./output/scan_ihs200_chr3.csv")
chr4_hh = read.csv("./output/scan_ihs200_chr4.csv")
chr5_hh = read.csv("./output/scan_ihs200_chr5.csv")

ehh_all = rbind(chr1_hh,chr2_hh,chr3_hh,chr4_hh,chr5_hh)

ehh_all = na.omit(ehh_all)
ehh_out = ehh_all[ehh_all$IHS>quantile(ehh_all$IHS,0.95),]


chr1_beta = read.table("./output/BetaScanChr1out.txt", header=TRUE)
chr2_beta = read.table("./output/BetaScanChr2out.txt", header=TRUE)
chr3_beta = read.table("./output/BetaScanChr3out.txt", header=TRUE)
chr4_beta = read.table("./output/BetaScanChr4out.txt", header=TRUE)
chr5_beta = read.table("./output/BetaScanChr5out.txt", header=TRUE)
Chr = c(rep(1,nrow(chr1_beta)), rep(2,nrow(chr2_beta)), rep(3,nrow(chr3_beta)), rep(4,nrow(chr4_beta)), rep(5,nrow(chr5_beta)))

beta_all = rbind(chr1_beta,chr2_beta,chr3_beta,chr4_beta,chr5_beta)
beta_all = data.frame(Chr,beta_all)
beta_out = beta_all[beta_all$Beta1>quantile(beta_all$Beta1,0.95,na.rm=TRUE),]

# for all files
fn = system("ls ./output/*neiGWAS*.csv.gz", intern=TRUE)

ehh_SNP = paste0(ehh_out$CHR,"-",ehh_out$POSITION)
beta_SNP = paste0(beta_out$Chr,"-",beta_out$Position)

res = c()
for(i in fn) {
  gwas_out = read.csv(i, header=TRUE)
  gwas_out = gwas_out[gwas_out$P_nei<quantile(gwas_out$P_nei,0.001),]
  gwas_pos = gwas_out[gwas_out$beta_nei>0,]
  gwas_neg = gwas_out[gwas_out$beta_nei<0,]
  
  gwas_SNPpos = paste0(gwas_pos$Chr,"-",gwas_pos$Pos)
  gwas_SNPneg = paste0(gwas_neg$Chr,"-",gwas_neg$Pos)
  
  
  print(i)
  print(length(intersect(ehh_SNP,gwas_SNPpos)))
  print(length(intersect(beta_SNP,gwas_SNPpos)))
  print("----")
  print(length(intersect(ehh_SNP,gwas_SNPneg)))
  print(length(intersect(beta_SNP,gwas_SNPneg)))
  
  trait = i
  a = length(intersect(ehh_SNP,gwas_SNPpos))
  b = length(intersect(beta_SNP,gwas_SNPpos))
  c = length(intersect(ehh_SNP,gwas_SNPneg))
  d = length(intersect(beta_SNP,gwas_SNPneg))
  
  t_mat = matrix(c(a,b,c,d),2,2,byrow=TRUE)
  
  # "less" means H_0 "odd ratio (a*d)/(b*c) = 1; "H_1 "odd ratio (a*d)/(b*c) < 1" 
  # set "greater" when testing the opposite hypothesis
  p = fisher.test(t_mat,alternative="less")$p.value 
  
  res = rbind(res, c(trait,a,b,c,d,p))
}

res = as.data.frame(res)
colnames(res) = c("trait","beta2Pos_EHH","beta2Pos_Beta","beta2Neg_EHH","beta2Neg_Beta","pvalue")
res$beta2Pos_EHH = as.numeric(res$beta2Pos_EHH)
res$beta2Pos_Beta = as.numeric(res$beta2Pos_Beta)
res$beta2Neg_EHH = as.numeric(res$beta2Neg_EHH)
res$beta2Neg_Beta = as.numeric(res$beta2Neg_Beta)
res$pvalue = as.numeric(res$pvalue)

res$trait[res$pvalue < 0.05]
res$pvalue[res$pvalue < 0.05]


################################
# number and frequency at J = 12
s2 = list(); b2 = list()
fn = system("ls ./output/*_*S2.csv.gz", intern=TRUE)
for(i in fn[c(1,2,4,3,5,6,8,7)]) {
  gwas_out = read.csv(i, header=TRUE)
  gwas_out = gwas_out[gwas_out$P_nei<quantile(gwas_out$P_nei,0.001),]
  gwas_out$MAF[gwas_out$MAF>0.5] = 1 - gwas_out$MAF[gwas_out$MAF>0.5]
  print(wilcox.test(gwas_out$MAF~(gwas_out$beta_nei>0))$p.value)
  p = ggplot(gwas_out,aes(x=beta_nei)) + geom_histogram() + theme_classic() + 
    ylab("No. of SNPs (J = 12)") + xlab(expression(hat(italic(beta))[2]))
  
  b = ggplot(gwas_out,aes(x=(beta_nei>0),y=MAF,group=(beta_nei>0))) + ylim(0,0.6) +
    geom_boxplot()+geom_point()+theme_classic()+scale_x_discrete(labels=c("<0", ">0")) + 
    ylab("No. of SNPs (J = 12)") + xlab(expression(hat(italic(beta))[2])) +
    geom_text(x=1,y=0.6,label=paste0("(",sum(gwas_out$beta_nei<0),")"),size=4) +
    geom_text(x=2,y=0.6,label=paste0("(",sum(gwas_out$beta_nei>0),")"),size=4)
  
  s2 = append(s2,list(p))
  b2 = append(b2,list(b))
}

# for Figure S8a-h
beta2maf_hist = ((s2[[1]]+labs(title=substitute(paste(bold("a"),"  coef. at Zurich")), subtitle="Herbivore damage (Leaf holes)")) | 
  (s2[[2]]+labs(title=substitute(paste(bold("b"))),subtitle="Individual no. of external feeders")) |
  (s2[[3]]+labs(title=substitute(paste(bold("c"))),subtitle="Individual no. of internal feeders")) | 
  (s2[[4]]+labs(title=substitute(paste(bold("d"))),subtitle="Total no. of insect species"))) /
  ((s2[[5]]+labs(title=substitute(paste(bold("e"),"  coef. at Otsu")), subtitle="Herbivore damage (Leaf area loss)")) | 
  (s2[[6]]+labs(title=substitute(paste(bold("f"))),subtitle="Individual no. of external feeders")) |
  (s2[[7]]+labs(title=substitute(paste(bold("g"))),subtitle="Individual no. of internal feeders")) | 
  (s2[[8]]+labs(title=substitute(paste(bold("h"))),subtitle="Total no. of insect species"))) /
  ((b2[[1]]+labs(title=substitute(paste(bold("i"),"  MAF at Zurich")), subtitle="Herbivore damage (Leaf holes)")) |
  (b2[[2]]+labs(title=substitute(paste(bold("j"))),subtitle="Individual no. of external feeders")) |
  (b2[[3]]+labs(title=substitute(paste(bold("k"))),subtitle="Individual no. of internal feeders")) |
  (b2[[4]]+labs(title=substitute(paste(bold("l"))),subtitle="Total no. of insect species"))) /
  ((b2[[5]]+labs(title=substitute(paste(bold("m"),"  MAF at Otsu")), subtitle="Herbivore damage (Leaf area loss)")) |
  (b2[[6]]+labs(title=substitute(paste(bold("n"))),subtitle="Individual no. of external feeders")) |
  (b2[[7]]+labs(title=substitute(paste(bold("o"))),subtitle="Individual no. of internal feeders")) |
  (b2[[8]]+labs(title=substitute(paste(bold("p"))),subtitle="Total no. of insect species")))

ggsave(beta2maf_hist, filename="../figs/mafhistJ12.pdf",width=12,height=12)


########################
# selection scan figures

# for Figure S9a
ihs_p = ggplot(ehh_all,aes(x=IHS)) + geom_histogram() + 
  ylab("No. of SNPs")+xlab("iHS")+theme_classic()+ggtitle(substitute(paste(bold("a"),"  Positive selection"))) +
  geom_vline(xintercept=quantile(ehh_all$IHS,0.95),lty=2)

# for Figure S9b
beta_p = ggplot(beta_all,aes(x=Beta1)) + geom_histogram() + 
  ylab("No. of SNPs")+xlab("BetaScan")+theme_classic()+ggtitle(substitute(paste(bold("b"),"  Balancing selection"))) +
  geom_vline(xintercept=quantile(beta_all$Beta1,0.95),lty=2)

scan_hist = (ihs_p | beta_p)


# repanel (J=4)
repanel = rbind(res[1,], res[3,], res[7,], res[5,],
                res[9,], res[11,], res[15,], res[13,])
repanel$trait = factor(repanel$trait,levels=repanel$trait)

# barplots
bar = c()
for(i in repanel$trait) {
  sub = repanel[repanel$trait==i,2:5]
  df = data.frame(trait=rep(i,4),value=as.numeric(sub),group=colnames(sub))
  bar = rbind(bar,df)
}
bar$trait = factor(bar$trait,levels=repanel$trait)

which(repanel$pvalue[1:4]<0.05)

h1 = ggplot(bar[1:16,],aes(y=value,x=trait,fill=group))+geom_bar(stat="identity", position="dodge") + labs(title=substitute(paste(bold("c"),"  Scan at Zurich (J = 4)"))) +
  theme_classic() + ylab("No. of SNPs") + xlab("") + scale_x_discrete(labels=c("Herbivore damage (Leaf holes)","Individual no. of external feeders","Individual no. of internal feeders","Total no. of insect species")) +
  scale_fill_manual(values=c("red", "blue", "pink","skyblue"),name="Category",
                    breaks=c("beta2Pos_EHH", "beta2Pos_Beta", "beta2Neg_EHH", "beta2Neg_Beta"),
                    labels=c(expression(iHS:italic(beta)[2]>0), expression(BETA:italic(beta)[2]>0), expression(iHS:italic(beta)[2]<0), expression(BETA:italic(beta)[2]<0))) +
  geom_text(data.frame(x=1,y=17,group=NA),mapping=aes(x=x,y=y),label="*",size=8) +
  geom_text(data.frame(x=2,y=27,group=NA),mapping=aes(x=x,y=y),label="**",size=8) + 
  geom_text(data.frame(x=4,y=10,group=NA),mapping=aes(x=x,y=y),label="*",size=8) + theme(legend.position=c(0.9,0.9))


which(repanel$pvalue[5:8]<0.001)

h2 = ggplot(bar[17:32,],aes(y=value,x=trait,fill=group))+geom_bar(stat="identity", position="dodge") + labs(title=substitute(paste(bold("d"),"  Scan at Otsu (J = 4)"))) +
  theme_classic() + ylab("No. of SNPs") + xlab("") + scale_x_discrete(labels=c("Herbivore damage (Leaf area loss)","Individual no. of external feeders","Individual no. of internal feeders","Total no. of insect species")) +
  scale_fill_manual(values=c("red", "blue", "pink","skyblue"),name="Category",
                    breaks=c("beta2Pos_EHH", "beta2Pos_Beta", "beta2Neg_EHH", "beta2Neg_Beta"),
                    labels=c(expression(iHS:italic(beta)[2]>0), expression(BETA:italic(beta)[2]>0), expression(iHS:italic(beta)[2]<0), expression(BETA:italic(beta)[2]<0))) +
  geom_text(data.frame(x=1,y=27,group=NA),mapping=aes(x=x,y=y),label="***",size=8) + theme(legend.position="none")

# for Figure S9c-d
scan_bar1 = (h1 | h2)


# repanel (J=12)
repanel = rbind(res[2,], res[4,], res[8,], res[6,],
                res[10,], res[12,], res[16,], res[14,])
repanel$trait = factor(repanel$trait,levels=repanel$trait)

# barplots
bar = c()
for(i in repanel$trait) {
  sub = repanel[repanel$trait==i,2:5]
  df = data.frame(trait=rep(i,4),value=as.numeric(sub),group=colnames(sub))
  bar = rbind(bar,df)
}
bar$trait = factor(bar$trait,levels=repanel$trait)

which(repanel$pvalue[1:4]<0.05)

h1 = ggplot(bar[1:16,],aes(y=value,x=trait,fill=group))+geom_bar(stat="identity", position="dodge") + labs(title=substitute(paste(bold("e"),"  Scan at Zurich (J = 12)"))) +
  theme_classic() + ylab("No. of SNPs") + xlab("") + scale_x_discrete(labels=c("Herbivore damage (Leaf holes)","Individual no. of external feeders","Individual no. of internal feeders","Total no. of insect species")) +
  scale_fill_manual(values=c("red", "blue", "pink","skyblue"),name="Category",
                    breaks=c("beta2Pos_EHH", "beta2Pos_Beta", "beta2Neg_EHH", "beta2Neg_Beta"),
                    labels=c(expression(iHS:italic(beta)[2]>0), expression(BETA:italic(beta)[2]>0), expression(iHS:italic(beta)[2]<0), expression(BETA:italic(beta)[2]<0))) +
  geom_text(data.frame(x=1,y=38,group=NA),mapping=aes(x=x,y=y),label="**",size=8) + theme(legend.position=c(0.9,0.9))

which(repanel$pvalue[5:8]<0.001)

h2 = ggplot(bar[17:32,],aes(y=value,x=trait,fill=group))+geom_bar(stat="identity", position="dodge") + labs(title=substitute(paste(bold("f"),"  Scan at Otsu (J = 12)"))) +
  theme_classic() + ylab("No. of SNPs") + xlab("") + scale_x_discrete(labels=c("Herbivore damage (Leaf area loss)","Individual no. of external feeders","Individual no. of internal feeders","Total no. of insect species")) +
  scale_fill_manual(values=c("red", "blue", "pink","skyblue"),name="Category",
                    breaks=c("beta2Pos_EHH", "beta2Pos_Beta", "beta2Neg_EHH", "beta2Neg_Beta"),
                    labels=c(expression(iHS:italic(beta)[2]>0), expression(BETA:italic(beta)[2]>0), expression(iHS:italic(beta)[2]<0), expression(BETA:italic(beta)[2]<0))) +
  geom_text(data.frame(x=4,y=35,group=NA),mapping=aes(x=x,y=y),label="***",size=8) + theme(legend.position="none")

# for Figure S9e-f
scan_bar2 = (h1 | h2)

# Figure S9
scan_p = (scan_hist / scan_bar1 / scan_bar2) + plot_layout(heights=c(1,1,1))
ggsave(scan_p, filename="../figs/SelectionScan.pdf",width=18,height=9)
