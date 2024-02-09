#####################################
# Ordination analysis and bar plots #
#####################################

# load library
library(vegan)
library(tidyverse)
library(patchwork)

# load data on 2017 and 2018
d = read.csv("./insectData/SurveyCombined4GWAS_max.csv")
d2017 = subset(d, Year=="2017")
com2017 = d2017[,c("Pr","Mp","Bb","Le","Mummy","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Fi","Fo","Bh","Er")]

d2018 = subset(d, Year=="2018")
com2018 = d2018[,c("Pr","Mp","Bb","Le","Mummy","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Fi","Fo","Bh","Er")]

com = rbind(com2017,com2018)
d_env = rbind(d2017[,1:6],d2018[,1:6])

res_rda = rda(log(com+1)~factor(d_env$Site)+factor(d_env$Year),scale=TRUE)
# res_rda = rda(log(com+1),scale=TRUE) # w/o constraint factors
varp = res_rda$CA$eig^2 / sum(res_rda$CA$eig^2) # contributions calc. from the eigen vector
anova(res_rda)
plot(res_rda)

df = data.frame(res_rda$CCA$wa,d_env)
df$Year = factor(df$Year)

# for Figure S1a
pca_p = ggplot(data=df,aes(x=RDA1,y=RDA2,shape=Year,color=Site)) + geom_point(alpha=0.25) + 
  theme_classic() + geom_vline(xintercept=0,colour="grey",lty=2) + geom_hline(yintercept=0,colour="grey",lty=2) +
  xlab(paste0("RDA1 (",round(varp[1],2)*100,"%)")) + ylab(paste0("RDA2 (",round(varp[2],2)*100,"%)")) +
  scale_color_manual(values=c("red","blue"),labels=c("Otsu","Zurich")) + xlim(NA,0.08) + theme(legend.position=c(0.9,0.7))
  
# for Figure S1b-g
df = sort(apply(com2017[1601:3200,],2,sum), decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b1 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + 
  geom_bar(stat="identity",fill=c("green","blue","blue","green","blue","darkgray","green","green","green","blue")) + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Zurich 2017") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90)) + 
  geom_label(data.frame(x=8,y=1000),mapping=aes(x=x,y=y),label="External feeder",size=2,fill="lightgreen",hjust=0) +
  geom_label(data.frame(x=8,y=200),mapping=aes(x=x,y=y),label="Internal feeder",size=2,fill="skyblue",hjust=0) + 
  geom_label(data.frame(x=8,y=40),mapping=aes(x=x,y=y),label="Carnivore",size=2,fill="gray",hjust=0)

df = sort(apply(com2017[1:1600,],2,sum), decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b2 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + 
  geom_bar(stat="identity",fill=c("blue","green","blue","green","green","blue","green","blue","blue","green")) + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Otsu 2017") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

df = sort(apply(com2018[1601:3200,],2,sum),decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b3 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + 
  geom_bar(stat="identity",fill=c("blue","blue","green","green","blue","darkgray","green","darkgrey","green","blue")) + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Zurich 2018") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

df = sort(apply(com2018[1:1600,],2,sum),decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b4 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + 
  geom_bar(stat="identity",fill=c("green","blue","green","green","green","blue","darkgray","blue","green","green")) + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Otsu 2018") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

# load 2019 data
d2019 = read.csv("./insectData/Survey20194GWAS_max.csv")

dz2019 = subset(d2019, Site=="ZH")
df = sort(apply(dz2019[,c("Pr","Mp","Bb","Le","Mummy","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Fi","Fo","Bh","Er")],2,sum),decreasing = T)[1:10]
df = data.frame(names=names(df),no=df)
b5 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + 
  geom_bar(stat="identity",fill=c("green","green","blue","blue","green","blue","blue","green","green","green")) + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Zurich 2019") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

dj2019 = subset(d2019, Site=="JP")
df = sort(apply(dj2019[,c("Pr","Mp","Bb","Le","Mummy","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Fi","Fo","Bh","Er")],2,sum),decreasing = T)[1:10]
df = data.frame(names=names(df),no=df)
b6 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + 
  geom_bar(stat="identity",fill=c("blue","green","green","green","blue","green","blue","blue","blue","darkgray")) + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Otsu 2019") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

# Figure S1
p = (pca_p) / (b1+b2) / (b3+b4) / (b5+b6) + plot_annotation(tag_levels = "a") + plot_layout(heights=c(2,1,1,1))
ggsave(p,filename="../figs/InsectRDAandNo.pdf",width=8,height=10)


# trait values
h1 = ggplot(data=subset(d, Site=="ZH"),aes(x=Holes)) + geom_histogram() +
  theme_classic() + labs(title=substitute(paste(bold("a"),"  Zurich")),subtitle="Herbivore damage (Leaf holes)") + #scale_x_log10() +
  xlab("No. of leaf holes per plant") + ylab("No. of plants")
h2 = ggplot(data=subset(d, Site=="ZH"),aes(x=chewer)) + geom_histogram(binwidth=1) +
  theme_classic() + labs(title=substitute(paste(bold("b"))),subtitle="Individual no. of external feeders") + #scale_x_log10() +
  xlab("Individual no. of external feeders per plant") + ylab("No. of plants")
h3 = ggplot(data=subset(d, Site=="ZH"),aes(x=sucker)) + geom_histogram(binwidth=1) +
  theme_classic() + labs(title=substitute(paste(bold("c"))),subtitle="Individual no. of internal feeders") + #scale_x_log10() +
  xlab("Individual no. of internal feeders per plant") + ylab("No. of plants")
h4 = ggplot(data=subset(d, Site=="ZH"),aes(x=richness)) + geom_histogram(binwidth=1) +
  theme_classic() + labs(title=substitute(paste(bold("d"))),subtitle="Total no. of insect species") + #scale_x_log10() +
  xlab("Total no. of insect species per plant") + ylab("No. of plants")

h5 = ggplot(data=subset(d, Site=="ZH"),aes(x=divHexp)) + geom_histogram() + ggtitle(substitute(paste(bold("a"),"  Zurich"))) +
  theme_classic() + xlab("Exponential Shannon diversity") + ylab("No. of plants")
h6 = ggplot(data=subset(d, Site=="ZH"),aes(x=divD)) + geom_histogram() +
  theme_classic() + xlab("Simpson diversity") + ylab("No. of plants")


h7 = ggplot(data=subset(d, Site=="JP"),aes(x=Score)) + geom_histogram(binwidth=1) +
  theme_classic() + ggtitle(substitute(paste(bold("i"),"  Otsu"))) + labs(subtitle="Herbivore damage (Leaf area loss)") +
  xlab("Leaf area loss (score) per plant") + ylab("No. of plants")
h8 = ggplot(data=subset(d, Site=="JP"),aes(x=chewer)) + geom_histogram(binwidth=1) + # scale_x_log10() + 
  labs(title=substitute(paste(bold("j"))),subtitle="Individual no. of external feeders") + theme_classic() + 
  xlab("Individual no. of external feeders per plant") + ylab("No. of plants")
h9 = ggplot(data=subset(d, Site=="JP"),aes(x=sucker)) + geom_histogram(binwidth=1) + # scale_x_log10() + 
  labs(title=substitute(paste(bold("k"))),subtitle="Individual no. of internal feeders") + theme_classic() + 
  xlab("Individual no. of internal feeders per plant") + ylab("No. of plants")
h10 = ggplot(data=subset(d, Site=="JP"),aes(x=richness)) + geom_histogram(binwidth=1) + # scale_x_log10() + 
  labs(title=substitute(paste(bold("l"))),subtitle="Total no. of insect species") + theme_classic() + 
  xlab("Total no. of insect species per plant") + ylab("No. of plants")

h11 = ggplot(data=subset(d, Site=="JP"),aes(x=divHexp)) + geom_histogram() + ggtitle(substitute(paste(bold("b"),"  Otsu"))) +
  theme_classic() + xlab("Exponential Shannon diversity") + ylab("No. of plants")
h12 = ggplot(data=subset(d, Site=="JP"),aes(x=divD)) + geom_histogram() +
  theme_classic() + xlab("Simpson diversity") + ylab("No. of plants")

# Figure S2
h = (h5+h11) / (h6+h12)
ggsave(h, filename="../figs/TraitValues.pdf",width=6,height=4)


#############
# PVE barplot
res = c()
for(i in 1:2) {
  tmp = read.csv(paste0("./output/CHZoutS",i,".csv"),header=TRUE)
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
  tmp = read.csv(paste0("./output/JPNoutS",i,".csv"),header=TRUE)
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
  geom_text(data.frame(x=c(1.5,2.5),y=c(0.6,0.75)),mapping=aes(x=x,y=y),label="***",size=6) +
  labs(title=substitute(paste(bold("a"),"  Zurich")),subtitle="Herbivore damage (Leaf holes)")

bar2 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[10:12])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=c(1.5,2.5),y=c(0.25,0.25)),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Individual no. of external feeders")

bar3 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[7:9])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  labs(subtitle="Individual no. of internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveZH[13:15])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=c(1.5,2.5),y=c(0.25,0.25)),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Total no. of insect species")

barz = bar1 | bar2 | bar3 | bar4

# Japan
bar1 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[1:3])) + geom_col() +
  ylim(0,1) + ylab("PVE") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.30),mapping=aes(x=x,y=y),label="**",size=6) +
  geom_text(data.frame(x=2.5,y=0.30),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(title=substitute(paste(bold("b"),"  Otsu")),subtitle="Herbivore damage (Leaf area loss)")

bar2 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[10:12])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  labs(subtitle="Individual no. of external feeders")

bar3 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[7:9])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.2),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Individual no. of internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c(0,4,12)),y=pveJP[13:15])) + geom_col() +
  ylim(0,1) + ylab("") + xlab("No. of neighbors") + theme_classic() +
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=6) +
  labs(subtitle="Total no. of insect species")

barj = bar1 | bar2 | bar3 | bar4

# Figure S3
bar_pve = barz / barj
ggsave(bar_pve,filename="../figs/PVEsupp.pdf",width=12,height=8,bg="transparent")


# PVE barplot (J=4)
# Zurich
bar1 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[1:2])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("e"))), subtitle="Herbivore damage (Leaf holes)") +
  geom_text(data.frame(x=1.5,y=0.55),mapping=aes(x=x,y=y),label="***",size=8)

bar2 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[10:11])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("f"))), subtitle="Individual no. of external feeders") +
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=8) 

bar3 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[7:8])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("g"))), subtitle="Individual no. of internal feeders")

bar4 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveZH[13:14])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("h"))), subtitle="Total no. of insect species") + 
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=8) 

z = ((h1 | bar1) + plot_layout(widths=c(2.5,1))) / 
  ((h2 | bar2) + plot_layout(widths=c(2.5,1))) / 
  ((h3 | bar3) + plot_layout(widths=c(2.5,1))) / 
  ((h4 | bar4) + plot_layout(widths=c(2.5,1)))

# Japan
bar1 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[1:2])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("m"))), subtitle="Herbivore damage (Leaf area loss)") + 
  geom_text(data.frame(x=1.5,y=0.30),mapping=aes(x=x,y=y),label="**",size=8) 

bar2 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[10:11])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("n"))), subtitle="Individual no. of external feeders")

bar3 = ggplot(data=NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[7:8])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("o"))), subtitle="Individual no. of internal feeders") +
  geom_text(data.frame(x=1.5,y=0.2),mapping=aes(x=x,y=y),label="*",size=8)

bar4 = ggplot(NULL,aes(x=factor(c("focal","focal+neig."),levels=c("focal","focal+neig.")),y=pveJP[13:14])) + geom_col(width=0.7) +
  ylim(0,1) + ylab("PVE") + xlab("") + theme_classic() + labs(title=substitute(paste(bold("p"))), subtitle="Total no. of insect species") + 
  geom_text(data.frame(x=1.5,y=0.25),mapping=aes(x=x,y=y),label="*",size=8) 

j = ((h7 | bar1) + plot_layout(widths=c(2.5,1))) / 
  ((h8 | bar2) + plot_layout(widths=c(2.5,1))) / 
  ((h9 | bar3) + plot_layout(widths=c(2.5,1))) / 
  ((h10 | bar4) + plot_layout(widths=c(2.5,1)))

# main Figure 2
ggsave(z,filename="../figs/PVE_Zurich.pdf",width=9.5,height=9.5,bg="transparent")
ggsave(j,filename="../figs/PVE_Otsu.pdf",width=9.5,height=9.5,bg="transparent")
