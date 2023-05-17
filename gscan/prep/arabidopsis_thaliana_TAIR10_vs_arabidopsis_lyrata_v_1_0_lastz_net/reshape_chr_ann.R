##################################
# Reshaping csv into rehh inputs #
##################################

chr = read.csv("chr1.csv") # change chr when needed
# chr1 = read.csv("chr1.csv")
# chr2 = read.csv("chr2.csv")
# chr3 = read.csv("chr3.csv")
# chr4 = read.csv("chr4.csv")
# chr5 = read.csv("chr5.csv")

chr$Alyr = casefold(chr$Alyr,upper=TRUE)
chr$Atha = casefold(chr$Atha,upper=TRUE)

#omit non-ATGCs
chr = chr[(chr$Alyr!="N"),]
chr = chr[(chr$Alyr!="-"),]
chr = chr[(chr$Alyr!="Y"),]
chr = chr[(chr$Atha!="N"),]
chr = chr[(chr$Atha!="-"),]
chr = chr[(chr$Atha!="Y"),]


tab = table(chr$Pos)
dupl = as.numeric(names(tab[which(tab!=1)]))
sub1 = chr[(chr$Pos %in% as.numeric(names(tab[which(tab!=1)]))),]

maxTF = function(i) {
  sub = sub1[sub1$Pos==i,]
  maxTF = (sub$Score==max(sub$Score))
  return(as.numeric(sub[maxTF,]$X[1]))
}

library(parallel)
TFmat = mcmapply(maxTF, dupl, mc.cores=6L)

nondupl = chr[(chr$Pos %in% names(tab[which(tab==1)])),]
maxdupl = chr[(chr$X %in% TFmat),]
chr_new = rbind(nondupl,maxdupl)
chr_new = chr_new[order(chr_new$Pos),]
write.csv(chr_new,"chr1_ann.csv",row.names=FALSE) # change chr when needed

# output files, "chrX_ann.csv", should then be moved to data directory "../data".
