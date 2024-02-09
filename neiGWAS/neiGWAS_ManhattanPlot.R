##########################
# neiGWAS Manhattan plot #
##########################

# load plot library and functions
library(tidyverse)
library(patchwork)
source("coord.R")

###################
# Manhattan plots for the main figure
ggMan = function(f) {
  d = read.csv(f,header=TRUE)
  d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
  chr_rep = table(d$Chr)
  cols = c(rep(rgb(1,0,0,0.5), chr_rep[1]),
           rep(rgb(0,1,0,0.5), chr_rep[2]),
           rep(rgb(0,0,1,0.5), chr_rep[3]), 
           rep(rgb(0,0,0,0.5), chr_rep[4]),
           rep(rgb(1,0,1,0.5), chr_rep[5]))
  x = coord(d$Chr,d$Position)
  man = ggplot(NULL,aes(x=x$coord,y=-log10(d$P_nei))) + geom_point(colour=cols) + theme_classic() + scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
    ylab(expression(-log[10]*(italic(p)))) + geom_hline(yintercept=-log10(0.05/nrow(d)),lty=2,colour="black")
  return(man)
}

ggHist = function(f) {
  gwas_out = read.csv(f, header=TRUE)
  gwas_out = gwas_out[gwas_out$P_nei<quantile(gwas_out$P_nei,0.001),]
  p = ggplot(gwas_out,aes(x=beta_nei))+geom_histogram()+theme_classic()+ylab("No. of SNPs")+xlab(expression(hat(italic(beta))[2]))
  return(p)
}

ggBox = function(f) {
  gwas_out = read.csv(f, header=TRUE)
  gwas_out = gwas_out[gwas_out$P_nei<quantile(gwas_out$P_nei,0.001),]
  gwas_out$MAF[gwas_out$MAF>0.5] = 1 - gwas_out$MAF[gwas_out$MAF>0.5]
  b = ggplot(gwas_out,aes(x=(beta_nei>0),y=MAF,group=(beta_nei>0))) +
    geom_boxplot(outlier.shape=NA) + geom_point() + theme_classic() + ylim(0,0.6) +
    scale_x_discrete(labels=c("<0", ">0")) + ylab("MAF") + xlab(expression(hat(italic(beta))[2])) +
    geom_text(x=1,y=0.58,label=paste0("(",sum(gwas_out$beta_nei<0),")"),size=3) +
    geom_text(x=2,y=0.58,label=paste0("(",sum(gwas_out$beta_nei>0),")"),size=3)
  return(b)
}

# Zurich
man1 = ggMan(f="./output/CHZneiGWAS_HolesS1.csv.gz")
man2 = ggMan(f="./output/CHZneiGWAS_chewerS1.csv.gz")
man3 = ggMan(f="./output/CHZneiGWAS_suckerS1.csv.gz")
man4 = ggMan(f="./output/CHZneiGWAS_richnessS1.csv.gz")

hist1 = ggHist(f="./output/CHZneiGWAS_HolesS1.csv.gz")
hist2 = ggHist(f="./output/CHZneiGWAS_chewerS1.csv.gz")
hist3 = ggHist(f="./output/CHZneiGWAS_suckerS1.csv.gz")
hist4 = ggHist(f="./output/CHZneiGWAS_richnessS1.csv.gz")

box1 = ggBox(f="./output/CHZneiGWAS_HolesS1.csv.gz")
box2 = ggBox(f="./output/CHZneiGWAS_chewerS1.csv.gz")
box3 = ggBox(f="./output/CHZneiGWAS_suckerS1.csv.gz")
box4 = ggBox(f="./output/CHZneiGWAS_richnessS1.csv.gz")

man1 = man1 + labs(title=substitute(paste(bold("a"),"  Zurich")),subtitle="Herbivore damage (Leaf holes)") 
man2 = man2 + labs(title=substitute(paste(bold("b"))),subtitle="Individual no. of external feeders")
man3 = man3 + labs(title=substitute(paste(bold("c"))),subtitle="Individual no. of internal feeders")
man4 = man4 + labs(title=substitute(paste(bold("d"))),subtitle="Total no. of insect species")

hist1 = hist1 + labs(title=substitute(paste(bold("e"))),subtitle="Herbivore damage (Leaf holes)")
hist2 = hist2 + labs(title=substitute(paste(bold("f"))),subtitle="Individual no. of external feeders")
hist3 = hist3 + labs(title=substitute(paste(bold("g"))),subtitle="Individual no. of internal feeders")
hist4 = hist4 + labs(title=substitute(paste(bold("h"))),subtitle="Total no. of insect species")

man1 = man1 + hist1 + box1 + plot_layout(widths = c(3,1.5,0.5))
man2 = man2 + hist2 + box2 + plot_layout(widths = c(3,1.5,0.5))
man3 = man3 + hist3 + box3 + plot_layout(widths = c(3,1.5,0.5))
man4 = man4 + hist4 + box4 + plot_layout(widths = c(3,1.5,0.5))

man = (man1 / man2 / man3 / man4)
saveRDS(man,file="../figs/CHZ_GWAS_s1.rds")

nman1 = plot_spacer() + hist1 + box1 + plot_layout(widths = c(3,1.5,0.5))
nman2 = plot_spacer() + hist2 + box2 + plot_layout(widths = c(3,1.5,0.5))
nman3 = plot_spacer() + hist3 + box3 + plot_layout(widths = c(3,1.5,0.5))
nman4 = plot_spacer() + hist4 + box4 + plot_layout(widths = c(3,1.5,0.5))
nman = (nman1 / nman2 / nman3 / nman4)
saveRDS(nman,file="../figs/CHZ_GWAS_s1_blank.rds")


# Otsu
man1 = ggMan(f="./output/JPNneiGWAS_ScoreS1.csv.gz")
man2 = ggMan(f="./output/JPNneiGWAS_chewerS1.csv.gz")
man3 = ggMan(f="./output/JPNneiGWAS_suckerS1.csv.gz")
man4 = ggMan(f="./output/JPNneiGWAS_richnessS1.csv.gz")

hist1 = ggHist(f="./output/JPNneiGWAS_ScoreS1.csv.gz")
hist2 = ggHist(f="./output/JPNneiGWAS_chewerS1.csv.gz")
hist3 = ggHist(f="./output/JPNneiGWAS_suckerS1.csv.gz")
hist4 = ggHist(f="./output/JPNneiGWAS_richnessS1.csv.gz")

box1 = ggBox(f="./output/JPNneiGWAS_ScoreS1.csv.gz")
box2 = ggBox(f="./output/JPNneiGWAS_chewerS1.csv.gz")
box3 = ggBox(f="./output/JPNneiGWAS_suckerS1.csv.gz")
box4 = ggBox(f="./output/JPNneiGWAS_richnessS1.csv.gz")

man1 = man1 + labs(title=substitute(paste(bold("i"),"  Otsu")),subtitle="Herbivore damage (Leaf area loss)")
man2 = man2 + labs(title=substitute(paste(bold("j"))),subtitle="Individual no. of external feeders")
man3 = man3 + labs(title=substitute(paste(bold("k"))),subtitle="Individual no. of internal feeders")
man4 = man4 + labs(title=substitute(paste(bold("l"))),subtitle="Total no. of insect species")

hist1 = hist1 + labs(title=substitute(paste(bold("m"))),subtitle="Herbivore damage (Leaf area loss)")
hist2 = hist2 + labs(title=substitute(paste(bold("n"))),subtitle="Individual no. of external feeders")
hist3 = hist3 + labs(title=substitute(paste(bold("o"))),subtitle="Individual no. of internal feeders")
hist4 = hist4 + labs(title=substitute(paste(bold("p"))),subtitle="Total no. of insect species")

man1 = man1 + hist1 + box1 + plot_layout(widths = c(3,1.5,0.5))
man2 = man2 + hist2 + box2 + plot_layout(widths = c(3,1.5,0.5))
man3 = man3 + hist3 + box3 + plot_layout(widths = c(3,1.5,0.5))
man4 = man4 + hist4 + box4 + plot_layout(widths = c(3,1.5,0.5))

man = (man1 / man2 / man3 / man4)
saveRDS(man,file="../figs/JPN_GWAS_s1.rds")

nman1 = plot_spacer() + hist1 + box1 + plot_layout(widths = c(3,1.5,0.5))
nman2 = plot_spacer() + hist2 + box2 + plot_layout(widths = c(3,1.5,0.5))
nman3 = plot_spacer() + hist3 + box3 + plot_layout(widths = c(3,1.5,0.5))
nman4 = plot_spacer() + hist4 + box4 + plot_layout(widths = c(3,1.5,0.5))
nman = (nman1 / nman2 / nman3 / nman4)
saveRDS(nman,file="../figs/JPN_GWAS_s1_blank.rds")

#########################
# Manhattan plots for supplements
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

ggQQ = function(f,type="nei") {
  d = read.csv(f,header=TRUE)
  d$MAF[d$MAF>0.5] = 1 - d$MAF[d$MAF>0.5]
  chr_rep = table(d$Chr)
  cols = c(rep(rgb(1,0,0,0.5), chr_rep[1]),
           rep(rgb(0,1,0,0.5), chr_rep[2]),
           rep(rgb(0,0,1,0.5), chr_rep[3]), 
           rep(rgb(0,0,0,0.5), chr_rep[4]),
           rep(rgb(1,0,1,0.5), chr_rep[5]))
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
man1 = ggMan(f="./output/CHZneiGWAS_HolesS1.csv.gz",type="self")
man2 = ggMan(f="./output/CHZneiGWAS_chewerS1.csv.gz",type="self")
man3 = ggMan(f="./output/CHZneiGWAS_suckerS1.csv.gz",type="self")
man4 = ggMan(f="./output/CHZneiGWAS_richnessS1.csv.gz",type="self")
man1 = man1 + labs(title=substitute(paste(bold("a"),"  J = 0")),subtitle="Herbivore damage (Leaf holes)")
man2 = man2 + labs(title=substitute(paste(bold("b"))),subtitle="Individual no. of external feeders")
man3 = man3 + labs(title=substitute(paste(bold("c"))),subtitle="Individual no. of internal feeders")
man4 = man4 + labs(title=substitute(paste(bold("d"))),subtitle="Total no. of insect species")

qq1 = ggQQ(f="./output/CHZneiGWAS_HolesS1.csv.gz",type="self")
qq2 = ggQQ(f="./output/CHZneiGWAS_chewerS1.csv.gz",type="self")
qq3 = ggQQ(f="./output/CHZneiGWAS_suckerS1.csv.gz",type="self")
qq4 = ggQQ(f="./output/CHZneiGWAS_richnessS1.csv.gz",type="self")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/CHZ_GWAS_s0.rds")


# Otsu at J = 0
man1 = ggMan(f="./output/JPNneiGWAS_ScoreS1.csv.gz",type="self")
man2 = ggMan(f="./output/JPNneiGWAS_chewerS1.csv.gz",type="self")
man3 = ggMan(f="./output/JPNneiGWAS_suckerS1.csv.gz",type="self")
man4 = ggMan(f="./output/JPNneiGWAS_richnessS1.csv.gz",type="self")
man1 = man1 + labs(title=substitute(paste(bold("a"),"  J = 0")),subtitle="Herbivore damage (Leaf area loss)")
man2 = man2 + labs(title=substitute(paste(bold("b"))),subtitle="Individual no. of external feeders")
man3 = man3 + labs(title=substitute(paste(bold("c"))),subtitle="Individual no. of internal feeders")
man4 = man4 + labs(title=substitute(paste(bold("d"))),subtitle="Total no. of insect species")

qq1 = ggQQ(f="./output/JPNneiGWAS_ScoreS1.csv.gz",type="self")
qq2 = ggQQ(f="./output/JPNneiGWAS_chewerS1.csv.gz",type="self")
qq3 = ggQQ(f="./output/JPNneiGWAS_suckerS1.csv.gz",type="self")
qq4 = ggQQ(f="./output/JPNneiGWAS_richnessS1.csv.gz",type="self")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/JPN_GWAS_s0.rds")


# Zurich at J = 4
qq1 = ggQQ(f="./output/CHZneiGWAS_HolesS1.csv.gz",type="nei")
qq2 = ggQQ(f="./output/CHZneiGWAS_chewerS1.csv.gz",type="nei")
qq3 = ggQQ(f="./output/CHZneiGWAS_suckerS1.csv.gz",type="nei")
qq4 = ggQQ(f="./output/CHZneiGWAS_richnessS1.csv.gz",type="nei")
qq1 = qq1 + labs(title=substitute(paste(bold("e"),"  J = 4")),subtitle="Herbivore damage (Leaf holes)")
qq2 = qq2 + labs(title=substitute(paste(bold("f"))),subtitle="Individual no. of external feeders")
qq3 = qq3 + labs(title=substitute(paste(bold("g"))),subtitle="Individual no. of internal feeders")
qq4 = qq4 + labs(title=substitute(paste(bold("h"))),subtitle="Total no. of insect species")

qq = qq1 / qq2 /qq3 /qq4
saveRDS(qq,file="../figs/CHZ_QQ_s1.rds")

# Otsu at J = 4
qq1 = ggQQ(f="./output/JPNneiGWAS_ScoreS1.csv.gz",type="nei")
qq2 = ggQQ(f="./output/JPNneiGWAS_chewerS1.csv.gz",type="nei")
qq3 = ggQQ(f="./output/JPNneiGWAS_suckerS1.csv.gz",type="nei")
qq4 = ggQQ(f="./output/JPNneiGWAS_richnessS1.csv.gz",type="nei")
qq1 = qq1 + labs(title=substitute(paste(bold("e"),"  J = 4")),subtitle="Herbivore damage (Leaf area loss)")
qq2 = qq2 + labs(title=substitute(paste(bold("f"))),subtitle="Individual no. of external feeders")
qq3 = qq3 + labs(title=substitute(paste(bold("g"))),subtitle="Individual no. of internal feeders")
qq4 = qq4 + labs(title=substitute(paste(bold("h"))),subtitle="Total no. of insect species")

qq = qq1 / qq2 /qq3 /qq4
saveRDS(qq,file="../figs/JPN_QQ_s1.rds")


# Zurich at J = 12
man1 = ggMan(f="./output/CHZneiGWAS_HolesS2.csv.gz",type="nei")
man2 = ggMan(f="./output/CHZneiGWAS_chewerS2.csv.gz",type="nei")
man3 = ggMan(f="./output/CHZneiGWAS_suckerS2.csv.gz",type="nei")
man4 = ggMan(f="./output/CHZneiGWAS_richnessS2.csv.gz",type="nei")
man1 = man1 + labs(title=substitute(paste(bold("i"),"  J = 12")),subtitle="Herbivore damage (Leaf holes)")
man2 = man2 + labs(title=substitute(paste(bold("j"))),subtitle="Individual no. of external feeders")
man3 = man3 + labs(title=substitute(paste(bold("k"))),subtitle="Individual no. of internal feeders")
man4 = man4 + labs(title=substitute(paste(bold("l"))),subtitle="Total no. of insect species")

qq1 = ggQQ(f="./output/CHZneiGWAS_HolesS2.csv.gz",type="nei")
qq2 = ggQQ(f="./output/CHZneiGWAS_chewerS2.csv.gz",type="nei")
qq3 = ggQQ(f="./output/CHZneiGWAS_suckerS2.csv.gz",type="nei")
qq4 = ggQQ(f="./output/CHZneiGWAS_richnessS2.csv.gz",type="nei")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/CHZ_GWAS_s2.rds")


# Otsu at J = 12
man1 = ggMan(f="./output/JPNneiGWAS_ScoreS2.csv.gz",type="nei")
man2 = ggMan(f="./output/JPNneiGWAS_chewerS2.csv.gz",type="nei")
man3 = ggMan(f="./output/JPNneiGWAS_suckerS2.csv.gz",type="nei")
man4 = ggMan(f="./output/JPNneiGWAS_richnessS2.csv.gz",type="nei")
man1 = man1 + labs(title=substitute(paste(bold("i"),"  J = 12")),subtitle="Herbivore damage (Leaf area loss)")
man2 = man2 + labs(title=substitute(paste(bold("j"))),subtitle="Individual no. of external feeders")
man3 = man3 + labs(title=substitute(paste(bold("k"))),subtitle="Individual no. of internal feeders")
man4 = man4 + labs(title=substitute(paste(bold("l"))),subtitle="Total no. of insect species")

qq1 = ggQQ(f="./output/JPNneiGWAS_ScoreS2.csv.gz",type="nei")
qq2 = ggQQ(f="./output/JPNneiGWAS_chewerS2.csv.gz",type="nei")
qq3 = ggQQ(f="./output/JPNneiGWAS_suckerS2.csv.gz",type="nei")
qq4 = ggQQ(f="./output/JPNneiGWAS_richnessS2.csv.gz",type="nei")

qman1 = (man1 + qq1) + plot_layout(widths=c(3,1))
qman2 = (man2 + qq2) + plot_layout(widths=c(3,1))
qman3 = (man3 + qq3) + plot_layout(widths=c(3,1))
qman4 = (man4 + qq4) + plot_layout(widths=c(3,1))

qman = qman1 /qman2 / qman3 / qman4
saveRDS(qman,file="../figs/JPN_GWAS_s2.rds")


# main Figure 3
zman = readRDS(file="../figs/CHZ_GWAS_s1.rds")
jman = readRDS(file="../figs/JPN_GWAS_s1.rds")
man = (zman) | (jman)
ggsave(man,filename="../figs/GWASmain.jpg",width=18,height=8,dpi=600,bg="transparent")

znman = readRDS(file="../figs/CHZ_GWAS_s1_blank.rds")
jnman = readRDS(file="../figs/JPN_GWAS_s1_blank.rds")
nman = (znman) | (jnman)
ggsave(nman,filename="../figs/GWASmain_blank.pdf",width=18,height=8,bg="transparent")

# suppl Zurich (Figure S4)
a = readRDS(file="../figs/CHZ_GWAS_s0.rds")
b = readRDS(file="../figs/CHZ_QQ_s1.rds")
c = readRDS(file="../figs/CHZ_GWAS_s2.rds")
abc = (a | b | c) + plot_layout(widths=c(2,0.75,2))
ggsave(abc,filename="../figs/ZurichGWASall.jpg",width=17,height=8,dpi=300,bg="transparent")

# suppl Otsu (Figure S5)
a = readRDS(file="../figs/JPN_GWAS_s0.rds")
b = readRDS(file="../figs/JPN_QQ_s1.rds")
c = readRDS(file="../figs/JPN_GWAS_s2.rds")
abc = (a | b | c) + plot_layout(widths=c(2,0.75,2))
ggsave(abc,filename="../figs/JapanGWASall.jpg",width=17,height=8,dpi=300,bg="transparent")

