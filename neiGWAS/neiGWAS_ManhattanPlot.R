##########################
# neiGWAS Manhattan plot #
##########################

# load plot library and functions
library(ggplot2)
library(patchwork)

coord = function(chr, pos) {
  if(length(pos)!=length(chr)) stop("chr and pos length differ")
  chr <- as.factor(chr)
  coord <- 0
  M <- 0
  for (i in 1:nlevels(chr)) {
    w <- (chr == levels(chr)[i])
    pos.c <- pos[w]
    coord[w] <- M + pos.c
    mx <- max(pos.c)
    M <- M + mx
  }
  coord <- coord/M
  return(coord)
}

#############
# PVE barplot
res = c()
for(i in 1:2) {
  tmp = read.csv(paste0("../output/CHZoutS",i,".csv"),header=TRUE)
  tmp$PVE = ((tmp$sigma_1+tmp$sigma_2) / (tmp$sigma_1+tmp$sigma_2+tmp$sigma_e))
  res = rbind(res,tmp)
}

# LL-tests
# J=0
pchisq(2*(res$LL[(res$sigma_2==0)&(res$sigma_1!=0)] - res$LL[(res$sigma_2==0)&(res$sigma_1==0)]),1,lower.tail=FALSE)
# J=4
pchisq(2*(res$LL[(res$sigma_2!=0)&(res$sigma_1!=0)] - res$LL[(res$sigma_2==0)&(res$sigma_1!=0)]),1,lower.tail=FALSE)[1:5]
# J=12
pchisq(2*(res$LL[(res$sigma_2!=0)&(res$sigma_1!=0)] - res$LL[(res$sigma_2==0)&(res$sigma_1!=0)]),1,lower.tail=FALSE)[6:10]

pveZH = c(res$PVE[c(2,15*c(0:1)+3)],
          res$PVE[c(5,15*c(0:1)+6)],
          res$PVE[c(8,15*c(0:1)+9)],
          res$PVE[c(11,15*c(0:1)+12)],
          res$PVE[c(14,15*c(0:1)+15)]
)

res = c()
for(i in 1:2) {
  tmp = read.csv(paste0("../output/JPNoutS",i,".csv"),header=TRUE)
  tmp$PVE = ((tmp$sigma_1+tmp$sigma_2) / (tmp$sigma_1+tmp$sigma_2+tmp$sigma_e))
  res = rbind(res,tmp)
}

# LL-tests
# J=0
pchisq(2*(res$LL[(res$sigma_2==0)&(res$sigma_1!=0)] - res$LL[(res$sigma_2==0)&(res$sigma_1==0)]),1,lower.tail=FALSE)
# J=4
pchisq(2*(res$LL[(res$sigma_2!=0)&(res$sigma_1!=0)] - res$LL[(res$sigma_2==0)&(res$sigma_1!=0)]),1,lower.tail=FALSE)[1:5]
# J=12
pchisq(2*(res$LL[(res$sigma_2!=0)&(res$sigma_1!=0)] - res$LL[(res$sigma_2==0)&(res$sigma_1!=0)]),1,lower.tail=FALSE)[6:10]

pveJP = c(res$PVE[c(2,15*c(0:1)+3)],
          res$PVE[c(5,15*c(0:1)+6)],
          res$PVE[c(8,15*c(0:1)+9)],
          res$PVE[c(11,15*c(0:1)+12)],
          res$PVE[c(14,15*c(0:1)+15)]
)

# PVE barplot (all)
# Zurich
bar1 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[1:3])) + geom_col() +
  ylim(0,1) + ylab("PVE") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=c(1.75,2.75),y=c(0.6,0.75)),mapping=aes(x=x,y=y),label="***",size=6) +
  labs(title="(A) Zurich",subtitle="Leaf holes")

bar2 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[10:12])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=c(1.75,2.75),y=c(0.25,0.25)),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="External feeders")

bar3 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[7:9])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  labs(subtitle="Internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[13:15])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=c(1.75,2.75),y=c(0.25,0.25)),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Species number")

barz = bar1 | bar2 | bar3 | bar4

# Japan
bar1 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[1:3])) + geom_col() +
  ylim(0,1) + ylab("PVE") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=1.75,y=0.30),mapping=aes(x=x,y=y),label="**",size=6) +
  geom_text(data.frame(x=2.75,y=0.30),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(title="(B) Otsu",subtitle="Leaf area loss")

bar2 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[10:12])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  labs(subtitle="External feeders")

bar3 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[7:9])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=1.75,y=0.2),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[13:15])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=1.75,y=0.25),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Species number")

barj = bar1 | bar2 | bar3 | bar4

bar_pve = barz / barj
ggsave(bar_pve,filename="../figs/PVEsupp.pdf",width=7,height=5,bg="transparent")


# PVE barplot (s=1)
# Zurich
bar1 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[1:2])) + geom_col() +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.55),mapping=aes(x=x,y=y),label="***",size=8) +
  labs(subtitle="Leaf holes")

bar2 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[10:11])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=8) +
  labs(subtitle="External feeders")

bar3 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[7:8])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("") + theme_classic() +
  labs(subtitle="Internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[13:14])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=8) +
  labs(subtitle="Species number")

barz = bar1 | bar2 | bar3 | bar4

# Japan
bar1 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[1:2])) + geom_col() +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.30),mapping=aes(x=x,y=y),label="**",size=8) +
  labs(subtitle="Leaf area loss")

bar2 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[10:11])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("") + theme_classic() +
  labs(subtitle="External feeders")

bar3 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[7:8])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.2),mapping=aes(x=x,y=y),label="*",size=8) +
  labs(subtitle="Internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[13:14])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=8) +
  labs(subtitle="Species number")

barj = bar1 | bar2 | bar3 | bar4

bar_pve = barz / barj
ggsave(bar_pve,filename="../figs/PVEmain.pdf",width=7,height=6,bg="transparent")


###################
# Manhattan plots for the main figure
ggMan = function(f) {
  d = read.csv(f,header=TRUE)
  d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
  chr_rep = cumsum(table(d$Chr))
  cols = c(rgb(1,0,0, 2*d$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*d$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*d$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*d$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*d$MAF[(chr_rep[4]+1):(chr_rep[5])]))
  x = coord(d$Chr,d$Position)
  man = ggplot(NULL,aes(x=x,y=-log10(d$P_nei))) + geom_point(colour=cols) + theme_classic() + ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) +
    ylab("") + xlab("Chromosomes") + geom_hline(yintercept=-log10(0.05/nrow(d)),lty=2,colour="black") + geom_hline(yintercept=-log10(quantile(d$P_nei,0.001)),lty=2,colour="grey")
  return(man)
}

man1 = ggMan(f="../output/CHZneiGWAS_HolesS1.csv")
man2 = ggMan(f="../output/CHZneiGWAS_chewerS1.csv")
man3 = ggMan(f="../output/CHZneiGWAS_suckerS1.csv")
man4 = ggMan(f="../output/CHZneiGWAS_richnessS1.csv")

man1 = man1 + labs(title="(A) Zurich",subtitle="Leaf holes")
man2 = man2 + labs(subtitle="External feeders")
man3 = man3 + labs(subtitle="Internal feeders")
man4 = man4 + labs(subtitle="Species number")
lab = ggplot(data.frame(l = "Association score", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")
man = lab + (man1 / man2 / man3 / man4) + plot_layout(widths = c(1, 30))
saveRDS(man,file="../figs/CHZ_GWAS_s1.rds")

man1 = ggMan(f="../output/JPNneiGWAS_ScoreS1.csv")
man2 = ggMan(f="../output/JPNneiGWAS_chewerS1.csv")
man3 = ggMan(f="../output/JPNneiGWAS_suckerS1.csv")
man4 = ggMan(f="../output/JPNneiGWAS_richnessS1.csv")

man1 = man1 + labs(title="(B) Otsu",subtitle="Leaf area loss")
man2 = man2 + labs(subtitle="External feeders")
man3 = man3 + labs(subtitle="Internal feeders")
man4 = man4 + labs(subtitle="Species number")
lab = ggplot(data.frame(l = "Association score", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")
man = lab + (man1 / man2 / man3 / man4) + plot_layout(widths = c(1, 30))
saveRDS(man,file="../figs/JPN_GWAS_s1.rds")


#########################
# Manhattan plots for supplements
ggMan = function(f,type="nei") {
  d = read.csv(f,header=TRUE)
  d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
  chr_rep = cumsum(table(d$Chr))
  cols = c(rgb(1,0,0, 2*d$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*d$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*d$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*d$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*d$MAF[(chr_rep[4]+1):(chr_rep[5])]))
  x = coord(d$Chr,d$Position)
  if(type=="nei") {
    y = -log10(d$P_nei)
  } else {
    y = -log10(d$P_self)
  }
  man = ggplot(NULL,aes(x=x,y=y)) + geom_point(colour=cols) + theme_classic() + ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) +
    ylab(expression(-log[10](p))) + xlab("Chromosomes") + geom_hline(yintercept=-log10(0.05/nrow(d)),lty=2,colour="black")
  return(man)
}

ggQQ = function(f,type="nei") {
  d = read.csv(f,header=TRUE)
  d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
  chr_rep = cumsum(table(d$Chr))
  cols = c(rgb(1,0,0, 2*d$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*d$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*d$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*d$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*d$MAF[(chr_rep[4]+1):(chr_rep[5])]))
  if(type=="nei") {
    p = d$P_nei
    cols = cols[order(d$P_nei,decreasing=FALSE)]
  } else {
    p = d$P_self
    cols = cols[order(d$P_self,decreasing=FALSE)]
  }

  o = -log(sort(p,decreasing=FALSE),10)
  e = -log(ppoints(length(p)),10)

  qq = ggplot(data=NULL, mapping=aes(x=e,y=o))+
    geom_point(colour=cols)+
    geom_abline(intercept=0,slope=1,linetype="dashed")+
    theme_classic()+
    xlab(expression("Expected "*-log[10](p)))+ylab(expression("Observed "*-log[10](p)))
  return(qq)
}

# Zurich at J = 0
man1 = ggMan(f="../output/CHZneiGWAS_HolesS1.csv",type="self")
man2 = ggMan(f="../output/CHZneiGWAS_chewerS1.csv",type="self")
man3 = ggMan(f="../output/CHZneiGWAS_suckerS1.csv",type="self")
man4 = ggMan(f="../output/CHZneiGWAS_richnessS1.csv",type="self")
man1 = man1 + labs(title="(A) J = 0",subtitle="Leaf holes")
man2 = man2 + labs(subtitle="External feeders")
man3 = man3 + labs(subtitle="Internal feeders")
man4 = man4 + labs(subtitle="Species number")

qq1 = ggQQ(f="../output/CHZneiGWAS_HolesS1.csv",type="self")
qq2 = ggQQ(f="../output/CHZneiGWAS_chewerS1.csv",type="self")
qq3 = ggQQ(f="../output/CHZneiGWAS_suckerS1.csv",type="self")
qq4 = ggQQ(f="../output/CHZneiGWAS_richnessS1.csv",type="self")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/CHZ_GWAS_s0.rds")


# Otsu at J = 0
man1 = ggMan(f="../output/JPNneiGWAS_ScoreS1.csv",type="self")
man2 = ggMan(f="../output/JPNneiGWAS_chewerS1.csv",type="self")
man3 = ggMan(f="../output/JPNneiGWAS_suckerS1.csv",type="self")
man4 = ggMan(f="../output/JPNneiGWAS_richnessS1.csv",type="self")
man1 = man1 + labs(title="(A) J = 0",subtitle="Leaf area loss")
man2 = man2 + labs(subtitle="External feeders")
man3 = man3 + labs(subtitle="Internal feeders")
man4 = man4 + labs(subtitle="Species number")

qq1 = ggQQ(f="../output/JPNneiGWAS_ScoreS1.csv",type="self")
qq2 = ggQQ(f="../output/JPNneiGWAS_chewerS1.csv",type="self")
qq3 = ggQQ(f="../output/JPNneiGWAS_suckerS1.csv",type="self")
qq4 = ggQQ(f="../output/JPNneiGWAS_richnessS1.csv",type="self")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/JPN_GWAS_s0.rds")


# Zurich at J = 4
qq1 = ggQQ(f="../output/CHZneiGWAS_HolesS1.csv",type="nei")
qq2 = ggQQ(f="../output/CHZneiGWAS_chewerS1.csv",type="nei")
qq3 = ggQQ(f="../output/CHZneiGWAS_suckerS1.csv",type="nei")
qq4 = ggQQ(f="../output/CHZneiGWAS_richnessS1.csv",type="nei")
qq1 = qq1 + labs(title="(B) J = 4",subtitle="Leaf holes")
qq2 = qq2 + labs(subtitle="External feeders")
qq3 = qq3 + labs(subtitle="Internal feeders")
qq4 = qq4 + labs(subtitle="Species number")

qq = qq1 / qq2 /qq3 /qq4
saveRDS(qq,file="../figs/CHZ_QQ_s1.rds")

# Otsu at J = 4
qq1 = ggQQ(f="../output/JPNneiGWAS_ScoreS1.csv",type="nei")
qq2 = ggQQ(f="../output/JPNneiGWAS_chewerS1.csv",type="nei")
qq3 = ggQQ(f="../output/JPNneiGWAS_suckerS1.csv",type="nei")
qq4 = ggQQ(f="../output/JPNneiGWAS_richnessS1.csv",type="nei")
qq1 = qq1 + labs(title="(B) J = 4",subtitle="Leaf area loss")
qq2 = qq2 + labs(subtitle="External feeders")
qq3 = qq3 + labs(subtitle="Internal feeders")
qq4 = qq4 + labs(subtitle="Species number")

qq = qq1 / qq2 /qq3 /qq4
saveRDS(qq,file="../figs/JPN_QQ_s1.rds")


# Zurich at J = 12
man1 = ggMan(f="../output/CHZneiGWAS_HolesS2.csv",type="nei")
man2 = ggMan(f="../output/CHZneiGWAS_chewerS2.csv",type="nei")
man3 = ggMan(f="../output/CHZneiGWAS_suckerS2.csv",type="nei")
man4 = ggMan(f="../output/CHZneiGWAS_richnessS2.csv",type="nei")
man1 = man1 + labs(title="(C) J = 12",subtitle="Leaf holes")
man2 = man2 + labs(subtitle="External feeders")
man3 = man3 + labs(subtitle="Internal feeders")
man4 = man4 + labs(subtitle="Species number")

qq1 = ggQQ(f="../output/CHZneiGWAS_HolesS2.csv",type="nei")
qq2 = ggQQ(f="../output/CHZneiGWAS_chewerS2.csv",type="nei")
qq3 = ggQQ(f="../output/CHZneiGWAS_suckerS2.csv",type="nei")
qq4 = ggQQ(f="../output/CHZneiGWAS_richnessS2.csv",type="nei")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/CHZ_GWAS_s2.rds")


# Otsu at J = 12
man1 = ggMan(f="../output/JPNneiGWAS_ScoreS2.csv",type="nei")
man2 = ggMan(f="../output/JPNneiGWAS_chewerS2.csv",type="nei")
man3 = ggMan(f="../output/JPNneiGWAS_suckerS2.csv",type="nei")
man4 = ggMan(f="../output/JPNneiGWAS_richnessS2.csv",type="nei")
man1 = man1 + labs(title="(C) J = 12",subtitle="Leaf area loss")
man2 = man2 + labs(subtitle="External feeders")
man3 = man3 + labs(subtitle="Internal feeders")
man4 = man4 + labs(subtitle="Species number")

qq1 = ggQQ(f="../output/JPNneiGWAS_ScoreS2.csv",type="nei")
qq2 = ggQQ(f="../output/JPNneiGWAS_chewerS2.csv",type="nei")
qq3 = ggQQ(f="../output/JPNneiGWAS_suckerS2.csv",type="nei")
qq4 = ggQQ(f="../output/JPNneiGWAS_richnessS2.csv",type="nei")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/JPN_GWAS_s2.rds")


# main fig
zman = readRDS(file="../figs/CHZ_GWAS_s1.rds")
jman = readRDS(file="../figs/JPN_GWAS_s1.rds")
man = (zman) | (jman)

scan = readRDS(file="../figs/SelectionSNPno.rds")

p = (man / scan) + plot_layout(heights=c(3,2))

ggsave(p,filename="../figs/GWASgscan.jpg",width=12,height=9,dpi=600,bg="transparent")

# suppl Zurich
a = readRDS(file="../figs/CHZ_GWAS_s0.rds")
b = readRDS(file="../figs/CHZ_QQ_s1.rds")
c = readRDS(file="../figs/CHZ_GWAS_s2.rds")
abc = (a | b | c) + plot_layout(widths=c(2,1,2))
ggsave(abc,filename="../figs/ZurichGWASall.jpg",width=16,height=8,dpi=300,bg="transparent")

# suppl Otsu
a = readRDS(file="../figs/JPN_GWAS_s0.rds")
b = readRDS(file="../figs/JPN_QQ_s1.rds")
c = readRDS(file="../figs/JPN_GWAS_s2.rds")
abc = (a | b | c) + plot_layout(widths=c(2,1,2))
ggsave(abc,filename="../figs/JapanGWASall.jpg",width=16,height=8,dpi=300,bg="transparent")

