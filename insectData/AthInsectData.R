###########################################
# Insect herbivory and abundance for GWAS #
###########################################

# load library
library(vegan)

setwd("./insectData")

#####################
# Merge 2017 and 2018

# load data
d17 = read.csv("AthInsectSurvey2017.csv")
d18 = read.csv("AthInsectSurvey2018.csv")

setdiff(colnames(d18), colnames(d17))
setdiff(colnames(d17), colnames(d18))
intersect(colnames(d18), colnames(d17))

# 2017 dataset
dJP17 = subset(d17, Site=="JP")
dJP17[,-c(1:12)][is.na(dJP17[,-c(1:12)])] = 0 # replace NA as zero
pheno_max = aggregate(cbind(Bolting,Score,Holes,Mines,Ps,Pa,Pp,Lc,Px,Pr_e,Pr_l,Ar,Tni,Syrphinae_l,Cv,Mp,Mp_w,Le,Le_w,Bb,Bb_w,Mummy,Wasps,Fi,Fo,Bh,Er)~IndivID, data=dJP17, max)
dJP17_max = cbind(dJP17[1:1600,1:12],pheno_max[,-1])

dZH17 = subset(d17, Site=="ZH")
dZH17[,-c(1:12)][is.na(dZH17[,-c(1:12)])] = 0 # replace NA as zero
pheno_max = aggregate(cbind(Bolting,Score,Holes,Mines,Ps,Pa,Pp,Lc,Px,Pr_e,Pr_l,Ar,Tni,Syrphinae_l,Cv,Mp,Mp_w,Le,Le_w,Bb,Bb_w,Mummy,Wasps,Fi,Fo,Bh,Er)~IndivID, data=dZH17, max)
dZH17_max = cbind(dZH17[1:1600,1:12],pheno_max[,-1])

# 2018 dataset
dJP18 = subset(d18, Site=="JP")
dJP18[,-c(1:12)][is.na(dJP18[,-c(1:12)])] = 0 # replace NA as zero
pheno_max = aggregate(cbind(Bolting,Score,Holes,Mines,Ps,Pa,Pp,Lc,Px,Pr_e,Pr_l,Ar,Tni,Syrphinae_l,Cv,Mp,Mp_w,Le,Le_w,Bb,Bb_w,Mummy,Wasps,Fi,Fo,Bh,Er)~IndivID, data=dJP18, max)
dJP18_max = cbind(dJP18[1:1600,1:12],pheno_max[,-1])

dZH18 = subset(d18, Site=="ZH")
dZH18[,-c(1:12)][is.na(dZH18[,-c(1:12)])] = 0 # replace NA as zero
pheno_max = aggregate(cbind(Bolting,Score,Holes,Mines,Ps,Pa,Pp,Lc,Px,Pr_e,Pr_l,Ar,Tni,Syrphinae_l,Cv,Mp,Mp_w,Le,Le_w,Bb,Bb_w,Mummy,Wasps,Fi,Fo,Bh,Er)~IndivID, data=dZH18, max)
dZH18_max = cbind(dZH18[1:1600,1:12],pheno_max[,-1])

d_max = rbind(dJP17_max,
              dZH17_max,
              dJP18_max,
              dZH18_max)

Pr = c()
for(i in 1:nrow(d_max)) {
  Pr = c(Pr, max(d_max$Pr_e[i],d_max$Pr_l[i]))
}

Mp = c()
for(i in 1:nrow(d_max)) {
  Mp = c(Mp, max(d_max$Mp[i],d_max$Mp_w[i]))
}

Bb = c()
for(i in 1:nrow(d_max)) {
  Bb = c(Bb, max(d_max$Bb[i],d_max$Bb_w[i]))
}

Le = c()
for(i in 1:nrow(d_max)) {
  Le = c(Le, max(d_max$Le[i],d_max$Le_w[i]))
}

d_max = data.frame(d_max,Pr)
d_max$Mp = Mp
d_max$Bb = Bb
d_max$Le = Le

com = d_max[,c("Pr","Mp","Bb","Le","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Mummy","Fi","Fo","Bh","Er")]
divH = diversity(com, index="shannon")
divHexp = diversity(exp(com), index="shannon")
divD = diversity(com, index="simpson")
richness = specnumber(com)
d_max = data.frame(d_max[,1:16], com, divH, divHexp, divD, richness)

PxPr = d_max$Px + d_max$Pr
PsPa = d_max$Ps + d_max$Pa
chewer = apply(d_max[,c("Ps","Pa","Pp","Lc","Px","Pr","Ar","Tni","Bh")],1,sum)
sucker = apply(d_max[,c("Mp","Bb","Le","Er","Fi","Fo","Mines")],1,sum)
d_max = data.frame(d_max, PxPr, PsPa, chewer, sucker)

write.csv(d_max,"SurveyCombined4GWAS_max.csv", row.names=FALSE)


######################
# formatting 2019 data
dJP19 = read.csv("AthInsectSurvey2019JP.csv")
dZH19 = read.csv("AthInsectSurvey2019ZH.csv")

setdiff(colnames(dJP19),colnames(dZH19))
setdiff(colnames(dZH19),colnames(dJP19))
intersect(colnames(dZH19),colnames(dJP19))
sum(colnames(dJP19)!=colnames(dZH19))

# 2019 dataset
dJP19[,-c(1:12)][is.na(dJP19[,-c(1:12)])] = 0 # replace NA as zero
pheno_max = aggregate(cbind(Bolting,Score,Holes,Mines,Ps,Pa,Pp,Lc,Px,Pr_e,Pr_l,Ar,Tni,Syrphinae_l,Cv,Mp,Mp_w,Le,Le_w,Bb,Bb_w,Mummy,Fi,Fo,Bh,Er)~IndivID, data=dJP19, max)
dJP19_max = cbind(dJP19[1:600,1:12],pheno_max[,-1])

dZH19[,-c(1:12)][is.na(dZH19[,-c(1:12)])] = 0 # replace NA as zero
pheno_max = aggregate(cbind(Bolting,Score,Holes,Mines,Ps,Pa,Pp,Lc,Px,Pr_e,Pr_l,Ar,Tni,Syrphinae_l,Cv,Mp,Mp_w,Le,Le_w,Bb,Bb_w,Mummy,Fi,Fo,Bh,Er)~IndivID, data=dZH19, max)
dZH19_max = cbind(dZH19[1:600,1:12],pheno_max[,-1])

d_max = rbind(dJP19_max, dZH19_max)

Pr = c()
for(i in 1:nrow(d_max)) {
  Pr = c(Pr, max(d_max$Pr_e[i],d_max$Pr_l[i]))
}

Mp = c()
for(i in 1:nrow(d_max)) {
  Mp = c(Mp, max(d_max$Mp[i],d_max$Mp_w[i]))
}

Bb = c()
for(i in 1:nrow(d_max)) {
  Bb = c(Bb, max(d_max$Bb[i],d_max$Bb_w[i]))
}

Le = c()
for(i in 1:nrow(d_max)) {
  Le = c(Le, max(d_max$Le[i],d_max$Le_w[i]))
}

d_max = data.frame(d_max,Pr)
d_max$Mp = Mp
d_max$Bb = Bb
d_max$Le = Le

com = d_max[,c("Pr","Mp","Bb","Le","Ps","Pa","Pp","Lc","Px","Ar","Tni","Syrphinae_l","Cv","Mummy","Fi","Fo","Bh","Er")]
divH = diversity(com, index="shannon")
divHexp = diversity(exp(com), index="shannon")
divD = diversity(com, index="simpson")
richness = specnumber(com)
d_max = data.frame(d_max[,1:16], com, divH, divHexp, divD, richness)

PxPr = d_max$Px + d_max$Pr
PsPa = d_max$Ps + d_max$Pa
chewer = apply(d_max[,c("Ps","Pa","Pp","Lc","Px","Pr","Ar","Tni","Bh")],1,sum)
sucker = apply(d_max[,c("Mp","Bb","Le","Er","Fi","Fo","Mines")],1,sum)
d_max = data.frame(d_max, PxPr, PsPa, chewer, sucker)

write.csv(d_max,"Survey20194GWAS_max.csv", row.names=FALSE)

setwd("../")
