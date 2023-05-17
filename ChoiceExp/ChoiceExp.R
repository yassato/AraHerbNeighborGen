##############################################
# GLM and figures for choice experiment data #
##############################################

library(ggplot2)
library(MASS)

setwd("./ChoiceExp")
d = read.csv("bioassay.csv",header=TRUE)
d$Line = factor(d$Line,levels=c("Bg-2","Uod-1","Vastervik","Jm-0","Bla-1","Bro1-6"))
d$Cup = factor(d$Cup)

d1 = subset(d,Set==1)
res = glm.nb(holes~factor(Cup)+Line,data=d1)
summary(res)
drop1(res,test="Chisq")

d2 = subset(d,Set==2)
res = glm.nb(holes~factor(Cup)+Line,data=d2)
summary(res)
drop1(res,test="Chisq")

d3 = subset(d,Set==3)
res = glm.nb(holes~factor(Cup)+Line,data=d3)
summary(res)
drop1(res,test="Chisq")

b = ggplot(d,aes(x=Line,y=holes+1)) + 
  geom_boxplot(outlier.shape=NA) + 
  theme_classic() + scale_y_log10() + ylab("No. of leaf holes + 1") + xlab("") + theme(axis.text.x=element_text(angle=90)) +
  geom_jitter(shape=20,colour="black",alpha=0.25,position=position_jitter(0.2)) + geom_vline(xintercept=c(2.5,4.5),lty=1,col="grey") + 
  geom_text(data.frame(x=1.5,y=30),mapping=aes(x=x,y=y),label="***",size=6) +
  geom_text(data.frame(x=3.5,y=30),mapping=aes(x=x,y=y),label="**",size=6)

saveRDS(b,file="../../figs/ChoiceExp_boxplot.rds")
ggsave(b,filename="choice.pdf",width=5,height=4)

setwd("../")
