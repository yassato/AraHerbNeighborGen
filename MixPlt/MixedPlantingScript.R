###########################################
# GLM and figures for mixed planting data #
###########################################

# load library
library(emmeans)
library(tidyverse)
library(patchwork)
library(lmerTest)

setwd("./MixPlt")

# load data
d = read.csv("MixedPlantingData.csv",header=TRUE)
d = d[is.na(d$hole)==FALSE,] # excl. 4 dead individuals
d$Line = factor(d$Line,levels=c("Bg-2","Uod-1","Vastervik","Jm-0","Bla-1","Bro1-6"))

d$MonoPoly = as.factor(d$MonoPoly)
d$Year = as.factor(d$Year)
d$plotID = as.factor(d$plotID)

# for Table S5
library(lmerTest)
mainLMM = lmer(log(hole+1)~Year+Density+MonoPoly+Line+(1|plotID),offset=log(LeafLen),data=d)
fullLMM = lmer(log(hole+1)~Year+Density+MonoPoly+Line+MonoPoly*Density+Density*Line+MonoPoly*Line+(1|plotID),offset=log(LeafLen),data=d)
anova(mainLMM)
anova(fullLMM)

full_cont = lmer(log(hole+1)~Line*MonoPoly+MonoPoly*Density+Line*Density+Year+(1|plotID),offset=log(LeafLen),data=d)
summary(emmeans(full_cont,pairwise~MonoPoly|Line),adjust="sidak",mode="Satterthwaite")

# main figures
b1 = ggplot(d,aes(x=factor(Line:MonoPoly),y=hole/(LeafLen*0.1))) + 
  geom_boxplot(fill=rep(c("white","grey"),6),outlier.shape=NA) + 
  theme_classic() + scale_y_log10() + ylab("Herbivore damage") + xlab("") + theme(axis.text.x=element_text(angle=90)) +
  geom_jitter(shape=20, colour="black", alpha=0.25, position=position_jitter(0.2)) + geom_vline(xintercept=c(4.5,8.5),lty=1,col="grey") +
  scale_x_discrete(labels=c("Bg-2","Bg-2","Uod-1","Uod-1","Vastervik","Vastervik","Jm-0","Jm-0","Bla-1","Bla-1","Bro1-6","Bro1-6")) +
  geom_text(data.frame(x=c(1.5,5.5),y=rep(150,2)),mapping=aes(x=x,y=y),label="**",size=6) + 
  geom_text(data.frame(x=c(3.5,7.5),y=rep(150,2)),mapping=aes(x=x,y=y),label="*",size=6) +
  geom_label(data.frame(x=11,y=160),mapping=aes(x=x,y=y),label="Monoculture",size=3,fill="white",hjust=0) +
  geom_label(data.frame(x=11,y=100),mapping=aes(x=x,y=y),label="Mixture        ",size=3,fill=grey(0.5,0.5),hjust=0)

damhist = readRDS("../../figs/EffectSizeDamage.rds")
avgp = readRDS("../../figs/avgp.rds")
simeff = readRDS("../../figs/SimEffmain.rds")
damhist = damhist + ggtitle(substitute(paste(bold("a"))))
avgp = avgp + ggtitle(substitute(paste(bold("b"))))
simeff = simeff + ggtitle(substitute(paste(bold("c"))))

# composite main Figure 4
ab = ((damhist + avgp + simeff + plot_layout(widths=c(2,1,1))) / b1 + ggtitle(substitute(paste(bold("d"))))) + plot_layout(heights = c(1,1.5))
ggsave(ab,filename="monomixture.pdf",width=10.5,height=8)

# separated figures for supplementary materials
bs1 = ggplot(subset(d,Year=="2019"&Density=="Low"),aes(x=factor(Line:MonoPoly),y=hole/(LeafLen*0.1))) + 
  geom_boxplot(fill=rep(c("white","grey"),6),outlier.shape=NA) + geom_vline(xintercept=c(4.5,8.5),lty=1,col="grey") +
  theme_classic() + scale_y_log10(limits=c(0.45,200)) + ylab("Herbivore damage") + xlab("") + theme(axis.text.x=element_text(angle=90)) +
  geom_jitter(shape=20,colour="black",alpha=0.25,position=position_jitter(0.2)) + ggtitle(substitute(paste(bold("a"),"  Small patch in 2019"))) +
  scale_x_discrete(labels=c("Bg-2","Bg-2","Uod-1","Uod-1","Vastervik","Vastervik","Jm-0","Jm-0","Bla-1","Bla-1","Bro1-6","Bro1-6")) + 
  geom_text(data.frame(x=7.5,y=150),mapping=aes(x=x,y=y),label="*",size=6) +
  geom_text(data.frame(x=c(1:12),y=rep(1,12)),mapping=aes(x=x,y=y),label=paste0("(",aggregate(hole~MonoPoly+Line,subset(d,Year=="2019"&Density=="Low"),length)[,3],")"))

bs2 = ggplot(subset(d,Year=="2019"&Density=="High"),aes(x=factor(Line:MonoPoly),y=hole/(LeafLen*0.1))) + 
  geom_boxplot(fill=rep(c("white","grey"),6),outlier.shape=NA) + geom_vline(xintercept=c(4.5,8.5),lty=1,col="grey") +
  theme_classic() + scale_y_log10(limits=c(0.45,200)) + ylab("Herbivore damage") + xlab("") + theme(axis.text.x=element_text(angle=90)) +
  geom_jitter(shape=20,colour="black",alpha=0.25,position=position_jitter(0.2)) + ggtitle(substitute(paste(bold("a"),"  Large patch in 2019"))) + 
  scale_x_discrete(labels=c("Bg-2","Bg-2","Uod-1","Uod-1","Vastervik","Vastervik","Jm-0","Jm-0","Bla-1","Bla-1","Bro1-6","Bro1-6")) +
  geom_text(data.frame(x=1.5,y=150),mapping=aes(x=x,y=y),label="***",size=6) + geom_text(data.frame(x=c(3.5,7.5,11.5),y=rep(150,3)),mapping=aes(x=x,y=y),label="**",size=6) +
  geom_text(data.frame(x=c(1:12),y=rep(1,12)),mapping=aes(x=x,y=y),label=paste0("(",aggregate(hole~MonoPoly+Line,subset(d,Year=="2019"&Density=="High"),length)[,3],")"))

bs3 = ggplot(subset(d,Year=="2021"&Density=="Low"),aes(x=factor(Line:MonoPoly),y=hole/(LeafLen*0.1))) + 
  geom_boxplot(fill=rep(c("white","grey"),6),outlier.shape=NA) + geom_vline(xintercept=c(4.5,8.5),lty=1,col="grey") +
  theme_classic() + scale_y_log10(limits=c(0.45,200)) + ylab("Herbivore damage") + xlab("") + theme(axis.text.x=element_text(angle=90)) +
  geom_jitter(shape=20,colour="black",alpha=0.25,position=position_jitter(0.2)) + ggtitle(substitute(paste(bold("c"),"  Small patch in 2021"))) + 
  scale_x_discrete(labels=c("Bg-2","Bg-2","Uod-1","Uod-1","Vastervik","Vastervik","Jm-0","Jm-0","Bla-1","Bla-1","Bro1-6","Bro1-6")) +
  geom_text(data.frame(x=1.5,y=150),mapping=aes(x=x,y=y),label="*",size=6) + 
  geom_text(data.frame(x=c(1:12),y=rep(1,12)),mapping=aes(x=x,y=y),label=paste0("(",aggregate(hole~MonoPoly+Line,subset(d,Year=="2021"&Density=="Low"),length)[,3],")"))

bs4 = ggplot(subset(d,Year=="2021"&Density=="High"),aes(x=factor(Line:MonoPoly),y=hole/(LeafLen*0.1))) + 
  geom_boxplot(fill=rep(c("white","grey"),6),outlier.shape=NA) + geom_vline(xintercept=c(4.5,8.5),lty=1,col="grey") +
  theme_classic() + scale_y_log10(limits=c(0.45,200)) + ylab("Herbivore damage") + xlab("") + theme(axis.text.x=element_text(angle=90)) +
  geom_jitter(shape=20,colour="black",alpha=0.25,position=position_jitter(0.2)) + ggtitle(substitute(paste(bold("d"),"  Large patch in 2021"))) +
  scale_x_discrete(labels=c("Bg-2","Bg-2","Uod-1","Uod-1","Vastervik","Vastervik","Jm-0","Jm-0","Bla-1","Bla-1","Bro1-6","Bro1-6")) +
  geom_text(data.frame(x=c(1.5,5.5),y=rep(150,2)),mapping=aes(x=x,y=y),label="***",size=6) + geom_text(data.frame(x=3.5,y=150),mapping=aes(x=x,y=y),label="*",size=6) +
  geom_text(data.frame(x=c(1:12),y=rep(1,12)),mapping=aes(x=x,y=y),label=paste0("(",aggregate(hole~MonoPoly+Line,subset(d,Year=="2021"&Density=="High"),length)[,3],")"))

# Figure S12a-d
bs = bs1+bs2+bs3+bs4
ggsave(bs,filename="MixedPlantingSeparated.pdf",width=12,height=9)

##############
# for Table S6

# GLM for 2019 data
d2019 = subset(d, Year==2019)

# note: change * into + when testing main effects
full_low = lmer(log(hole+1)~MonoPoly*Line+(1|plotID),offset=log(LeafLen),data=subset(d2019,Density=="Low"))
anova(full_low)
emmeans(full_low,pairwise~MonoPoly|Line,adjust="sidak",mode="Satterthwaite")

full_high = lmer(log(hole+1)~MonoPoly*Line+(1|plotID),offset=log(LeafLen),data=subset(d2019,Density=="High"))
anova(full_high)
emmeans(full_high,pairwise~MonoPoly|Line,adjust="sidak",mode="Satterthwaite")


# GLM for 2021 data
d2021 = subset(d, Year==2021)

# note: change * into + when testing main effects
full_low = lmer(log(hole+1)~MonoPoly*Line+(1|plotID),offset=log(LeafLen),data=subset(d2021,Density=="Low"))
anova(full_low)
summary(emmeans(full_low,pairwise~MonoPoly|Line),adjust="sidak",mode="Satterthwaite")

full_high = lmer(log(hole+1)~MonoPoly*Line+(1|plotID),offset=log(LeafLen),data=subset(d2021,Density=="High"))
anova(full_high)
summary(emmeans(full_high,pairwise~MonoPoly|Line),adjust="sidak",mode="Satterthwaite")


setwd("../")
