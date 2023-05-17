#####################################
# Ordination analysis and bar plots #
#####################################

setwd("./insectData")

# load library
library(vegan)
library(ggplot2)
library(patchwork)

# load data on 2017 and 2018
d = read.csv("SurveyCombined4GWAS_max.csv")
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

pca_p = ggplot(data=df,aes(x=RDA1,y=RDA2,shape=Year,color=Site)) + geom_point(alpha=0.25) + 
  theme_classic() + geom_vline(xintercept=0,colour="grey",lty=2) + geom_hline(yintercept=0,colour="grey",lty=2) +
  xlab(paste0("RDA1 (",round(varp[1],2)*100,"%)")) + ylab(paste0("RDA2 (",round(varp[2],2)*100,"%)")) +
  scale_color_manual(values=c("red","blue"),labels=c("Otsu","Zurich")) + xlim(NA,0.08) + theme(legend.position=c(0.9,0.7))
  
# figure for supplementary materials
df = sort(apply(com2017[1601:3200,],2,sum), decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b1 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + geom_bar(stat="identity") + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Zurich 2017") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

df = sort(apply(com2017[1:1600,],2,sum), decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b2 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + geom_bar(stat="identity") + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Otsu 2017") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

df = sort(apply(com2018[1601:3200,],2,sum),decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b3 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + geom_bar(stat="identity") + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Zurich 2018") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

df = sort(apply(com2018[1:1600,],2,sum),decreasing=TRUE)[1:10]
df = data.frame(names=names(df),no=df)
b4 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + geom_bar(stat="identity") + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Otsu 2018") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

# load data on 2019
d2019 = read.csv("Survey20194GWAS_max.csv")

dz2019 = subset(d2019, Site=="ZH")
df = sort(apply(dz2019[,c("Pr","Mp","Bb","Le","Mummy","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Fi","Fo","Bh","Er")],2,sum),decreasing = T)[1:10]
df = data.frame(names=names(df),no=df)
b5 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + geom_bar(stat="identity") + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Zurich 2019") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

dj2019 = subset(d2019, Site=="JP")
df = sort(apply(dj2019[,c("Pr","Mp","Bb","Le","Mummy","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Fi","Fo","Bh","Er")],2,sum),decreasing = T)[1:10]
df = data.frame(names=names(df),no=df)
b6 = ggplot(data=df,aes(x=reorder(names,-no),y=no)) + geom_bar(stat="identity") + 
  theme_classic() + ylab("No. of individuals") + xlab("") + ggtitle("Otsu 2019") + 
  scale_y_log10() + theme(axis.text.x=element_text(angle=90))

p = (pca_p) / (b1+b2) / (b3+b4) / (b5+b6) + plot_annotation(tag_levels = "A") + plot_layout(heights=c(2,1,1,1))
ggsave(p,filename="InsectRDAandNo.pdf",width=8,height=10)

# trait values
h1 = ggplot(data=subset(d, Site=="ZH"),aes(x=Holes)) + geom_histogram() +
  theme_classic() + ggtitle("(A) Zurich") + xlab("No. of leaf holes") + ylab("No. of individuals")
h2 = ggplot(data=subset(d, Site=="ZH"),aes(x=chewer)) + geom_histogram() +
  theme_classic() + xlab("No. of external feeders") + ylab("No. of individuals")
h3 = ggplot(data=subset(d, Site=="ZH"),aes(x=sucker)) + geom_histogram() +
  theme_classic() + xlab("No. of internal feeders") + ylab("No. of individuals")
h4 = ggplot(data=subset(d, Site=="ZH"),aes(x=richness)) + geom_histogram() +
  theme_classic() + xlab("No. of species") + ylab("No. of individuals")
h5 = ggplot(data=subset(d, Site=="ZH"),aes(x=divHexp)) + geom_histogram() +
  theme_classic() + xlab("Exponential Shannon diversity") + ylab("No. of individuals")
h6 = ggplot(data=subset(d, Site=="ZH"),aes(x=divD)) + geom_histogram() +
  theme_classic() + xlab("Simpson diversity") + ylab("No. of individuals")

h7 = ggplot(data=subset(d, Site=="JP"),aes(x=Score)) + geom_histogram() +
  theme_classic() + ggtitle("(B) Otsu") + xlab("Leaf area loss") + ylab("No. of individuals")
h8 = ggplot(data=subset(d, Site=="JP"),aes(x=chewer)) + geom_histogram() +
  theme_classic() + xlab("No. of external feeders") + ylab("No. of individuals")
h9 = ggplot(data=subset(d, Site=="JP"),aes(x=sucker)) + geom_histogram() +
  theme_classic() + xlab("No. of internal feeders") + ylab("No. of individuals")
h10 = ggplot(data=subset(d, Site=="JP"),aes(x=richness)) + geom_histogram() +
  theme_classic() + xlab("No. of species") + ylab("No. of individuals")
h11 = ggplot(data=subset(d, Site=="JP"),aes(x=divHexp)) + geom_histogram() +
  theme_classic() + xlab("Exponential Shannon diversity") + ylab("No. of individuals")
h12 = ggplot(data=subset(d, Site=="JP"),aes(x=divD)) + geom_histogram() +
  theme_classic() + xlab("Simpson diversity") + ylab("No. of individuals")

h = (h1+h7) / (h2+h8) /(h3+h9) / (h4+h10) / (h5+h11) / (h6+h12)
ggsave(h, filename="TraitValues.pdf",width=6,height=10)

setwd("../")
